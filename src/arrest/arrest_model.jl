# Top-level arrest models: (arrest_state, progress, signals) -> arrest_level in [0,1].
# Each concrete model defines initial_arrest_state/advance_arrest/arrest_level.

abstract type AbstractArrestModel end

function initial_arrest_state end
function advance_arrest end
function arrest_level end

# recursively advances a model-defined NamedTuple state by rate*dt -- discrete
# step-loop hosts use this directly; an empty rate NamedTuple leaves the
# matching state branch unchanged.
step_state(state::NamedTuple, rate::NamedTuple{(),Tuple{}}, dt) = state
step_state(state::NamedTuple, rate::NamedTuple, dt) =
    NamedTuple{keys(state)}(map((s, r) -> step_state(s, r, dt), values(state), values(rate)))
step_state(state::Number, rate::Number, dt) = state + rate * dt
arrest_conditions(::AbstractArrestModel) = ()

# one induction controller + one breakage controller. level = induction*(1-breakage),
# fuzzy AND-NOT reducing to `in_diapause && !broken` at the {0,1} extremes.
struct ComposedArrest{I<:AbstractArrestController,B<:AbstractArrestController} <: AbstractArrestModel
    induction::I
    breakage::B
end
ComposedArrest(; induction, breakage=NeverController()) = ComposedArrest(induction, breakage)

initial_arrest_state(m::ComposedArrest) =
    (; induction=initial_controller_state(m.induction), breakage=initial_controller_state(m.breakage))

advance_arrest(m::ComposedArrest, arrest_state, progress, signals) =
    (; induction=controller_rate(m.induction, arrest_state.induction, progress, signals, m, arrest_state),
       breakage=controller_rate(m.breakage, arrest_state.breakage, progress, signals, m, arrest_state))

function arrest_level(m::ComposedArrest, arrest_state, progress, signals)
    ind = level(m.induction, arrest_state.induction, progress, signals, m, arrest_state)
    brk = level(m.breakage, arrest_state.breakage, progress, signals, m, arrest_state)
    ind * (1 - brk)
end

# built via pure Tuple recursion (_flatten_tuples), not Any[]/push! -- see
# the matching comment on composition.jl's trigger_conditions.
function arrest_conditions(m::ComposedArrest)
    _flatten_tuples(map(((:induction, m.induction), (:breakage, m.breakage))) do (k, ctrl)
        map(trigger_conditions(ctrl)) do cond
            (arrest_state, progress, signals) -> cond(getfield(arrest_state, k), progress, signals, m, arrest_state)
        end
    end)
end

# whole-model composition, same NamedTuple/max-min idiom as composition.jl --
# needed since Chortoicetes is "diapause (has breakage) OR quiescence (none)",
# not a single induction/breakage pair.
struct AnyArrestModel{T<:NamedTuple} <: AbstractArrestModel
    models::T
end
AnyArrestModel(; kw...) = AnyArrestModel(NamedTuple(kw))

struct AllArrestModel{T<:NamedTuple} <: AbstractArrestModel
    models::T
end
AllArrestModel(; kw...) = AllArrestModel(NamedTuple(kw))

const _CompositeArrestModel = Union{AnyArrestModel,AllArrestModel}

initial_arrest_state(m::_CompositeArrestModel) =
    NamedTuple{keys(m.models)}(map(initial_arrest_state, values(m.models)))

function advance_arrest(m::_CompositeArrestModel, arrest_state, progress, signals)
    ks = keys(m.models)
    NamedTuple{ks}(map(ks) do k
        advance_arrest(getfield(m.models, k), getfield(arrest_state, k), progress, signals)
    end)
end

function _sub_arrest_levels(m::_CompositeArrestModel, arrest_state, progress, signals)
    ks = keys(m.models)
    map(ks) do k
        arrest_level(getfield(m.models, k), getfield(arrest_state, k), progress, signals)
    end
end
# hard maximum/minimum -- a kink at the crossover regardless of child
# smoothing. A smooth version needs a SmoothingStrategy field here routed
# through HeatExchange.safe_max/safe_min.
arrest_level(m::AnyArrestModel, arrest_state, progress, signals) =
    maximum(_sub_arrest_levels(m, arrest_state, progress, signals))
arrest_level(m::AllArrestModel, arrest_state, progress, signals) =
    minimum(_sub_arrest_levels(m, arrest_state, progress, signals))

function arrest_conditions(m::_CompositeArrestModel)
    ks = keys(m.models)
    _flatten_tuples(map(ks) do k
        sub = getfield(m.models, k)
        map(arrest_conditions(sub)) do cond
            (arrest_state, progress, signals) -> cond(getfield(arrest_state, k), progress, signals)
        end
    end)
end
