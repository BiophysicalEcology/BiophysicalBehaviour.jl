# Named-NamedTuple composition (not Tuple like CombinedActivity/CombinedSurvival --
# sub-controllers carry state that must thread through ComponentArrays by name).
# Fuzzy-logic level combination: max for OR, min for AND (Zadeh), reducing to
# boolean any/all at the {0,1} extremes.

struct AnyController{T<:NamedTuple} <: AbstractArrestController
    controllers::T
end
AnyController(; kw...) = AnyController(NamedTuple(kw))

struct AllController{T<:NamedTuple} <: AbstractArrestController
    controllers::T
end
AllController(; kw...) = AllController(NamedTuple(kw))

const _CompositeController = Union{AnyController,AllController}

function initial_controller_state(c::_CompositeController)
    NamedTuple{keys(c.controllers)}(map(initial_controller_state, values(c.controllers)))
end

function controller_rate(c::_CompositeController, own_state, progress, signals, model, arrest_state)
    ks = keys(c.controllers)
    NamedTuple{ks}(map(ks) do k
        controller_rate(getfield(c.controllers, k), getfield(own_state, k), progress, signals, model, arrest_state)
    end)
end

function _sub_levels(c::_CompositeController, own_state, progress, signals, model, arrest_state)
    ks = keys(c.controllers)
    map(ks) do k
        controller_level(getfield(c.controllers, k), getfield(own_state, k), progress, signals, model, arrest_state)
    end
end
# hard maximum/minimum -- a kink at the crossover regardless of child
# smoothing. A smooth version needs a SmoothingStrategy field here routed
# through HeatExchange.safe_max/safe_min.
controller_level(c::AnyController, own_state, progress, signals, model, arrest_state) =
    maximum(_sub_levels(c, own_state, progress, signals, model, arrest_state))
controller_level(c::AllController, own_state, progress, signals, model, arrest_state) =
    minimum(_sub_levels(c, own_state, progress, signals, model, arrest_state))

register_callback(c::_CompositeController) = any(register_callback, values(c.controllers))

# built via pure Tuple recursion (_flatten_tuples), not Any[]/push! -- a host
# splicing this into its own fixed condition tuple (see e.g. egg_model's
# phases.jl) depends on the result being concretely typed for type-stable
# root-finding, not just "setup-time, so it doesn't matter" convenience.
function trigger_conditions(c::_CompositeController)
    ks = keys(c.controllers)
    _flatten_tuples(map(ks) do k
        sub = getfield(c.controllers, k)
        map(trigger_conditions(sub)) do cond
            (own_state, progress, signals, model, arrest_state) ->
                cond(getfield(own_state, k), progress, signals, model, arrest_state)
        end
    end)
end
