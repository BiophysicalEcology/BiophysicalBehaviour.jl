# =============================================================================
# Multi-part NLP formulation for `IPOPTControl` — the sole NLP strategy.
#
# `MultipartNLP` solves a genuine multi-part organism: one regulated core shared
# across all parts, per-part surface temperatures and effectors, and a single
# whole-organism respiration balance closed at the lung. A single `Body` is the
# one-part case, so there is no separate single-body formulation.
#
# The decision-variable layout is NOT hand-indexed. It is a nested `NamedTuple`
# *template* built from the organism's parts, with `Unitful` leaves, and
# `Flatten.flatten` / `Flatten.reconstruct` bridge between that named structure and
# the flat `Vector{Float64}` Ipopt optimises. Units live in the leaves' `Quantity`
# *types*: targeting `Real`, `flatten` rips out only each `Quantity`'s inner
# `Float64` magnitude (so the optimiser vector is homogeneous `Float64`), and
# `reconstruct` puts a `Float64` back into that typed slot — units round-trip
# through the type, never stripped or reattached by hand. Variables, bounds,
# scaling, and residuals are every one of them a template of the same shape
# flattened in the same order — so a bound can never drift out of alignment with
# the variable it bounds, and adding a part or a per-part effector is a change to
# the *structure*, never to a pile of `x[3 + 4i]` offsets. Enzyme differentiates
# straight through `reconstruct` (the reconstructed leaves are concretely typed),
# so the residual/objective read Unitful named fields off it.
#
# Variable template (`_variable_template`):
#   core_temperature            (K)          whole-organism
#   log_metabolic_heat_flow     (log W)      whole-organism (the log-transform invariant)
#   panting_rate                (—)          whole-organism
#   parts.<name>.skin_temperature        (K)
#   parts.<name>.insulation_temperature  (K)
#   parts.<name>.flesh_conductivity      (W/m/K)
#   parts.<name>.skin_wetness            (—)
#
# Residual template (`_residual_template`): per part a `surface_balance` (W) and a
# `skin_temperature` (K) residual, then the whole-organism respiration `balance`
# (W), then the `q10` inequality slack (last leaf → the one `≥ 0` constraint; all
# others are `= 0`).
#
# Geometry effectors (insulation_depth / axis_ratio) are deferred, exactly as the
# multi-part rule-based ladder defers Piloerect / Uncurl: no per-part body rebuild
# happens inside the residual, so the per-part `setups` are fixed across the solve.
#
# The multipart methods slot into the *same* direct-Ipopt + Enzyme callback
# machinery as the single-body strategies (ipopt.jl): `IPOPTSolverCache`,
# `_ipopt_solve!`, the five `Evaluate*` functors, and `_lagrangian` are all generic
# over the packed type and operate on the flat vector. This file only adds the
# `MultipartNLPPacked` methods those functors dispatch to. Physics is
# `HeatExchange.part_surface_residuals` per part, so zero physics is duplicated.
# =============================================================================

"""
    MultipartNLP <: HeatExchange.NLPStrategy

NLP strategy for a multi-part organism: shared regulated core, per-part surface
temperatures and effectors (flesh_conductivity, skin_wetness), whole-organism
respiration closed at the lung. Select via `IPOPTControl(; nlp_strategy = MultipartNLP())`.
Works for any `Organism` — a single `Body` is treated as the one-part case.
"""
struct MultipartNLP <: HeatExchange.NLPStrategy end

# --- Templates: one structural definition, flattened many ways ---------------
#
# Each helper builds the same nested shape with different leaf values. Because the
# field structure and order are identical, `Flatten.flatten(_, Real)` yields
# positionally-aligned vectors — the guarantee that bounds/scaling/initials line up
# with the variables they describe.

# Per-part leaf NamedTuple, four values in canonical order.
@inline _part_leaf(skin_temperature, insulation_temperature, flesh_conductivity, skin_wetness) =
    (; skin_temperature, insulation_temperature, flesh_conductivity, skin_wetness)

# Whole-organism scalars + a per-part-name NamedTuple of leaves.
@inline function _variable_structure(part_names::NTuple{N,Symbol}, core_temperature,
        log_metabolic_heat_flow, panting_rate, part_leaves::NTuple{N}) where {N}
    return (;
        core_temperature,
        log_metabolic_heat_flow,
        panting_rate,
        parts = NamedTuple{part_names}(part_leaves),
    )
end

# The unit-carrying prototype for `flatten`/`reconstruct`. Units live in the leaves'
# `Quantity` *types*: targeting `Real`, `Flatten.flatten` recurses into each
# `Quantity` and rips out only its inner `Float64` magnitude (so the flat optimiser
# vector is homogeneous `Float64` — the type-stable, Enzyme-safe form), and
# `Flatten.reconstruct` puts a `Float64` back into that typed slot, rebuilding the
# `Quantity` with its unit intact. Units are never stripped or reattached by hand;
# they round-trip through the type. Critically, the leaf units come from the *model's
# own quantities* (`temperature` and `conductivity` prototypes pulled from the
# organism), not hardcoded literals: whatever units the model uses become the
# optimiser's canonical units, and because the bounds/initial structures in
# `_inputs!` are built from those same model quantities, a bound can never be
# flattened in a different unit than the template reconstructs against — there is no
# canonical unit to drift out of, so nothing needs an `ustrip` to a fixed unit.
# `log_metabolic_heat_flow` / `panting_rate` / `skin_wetness` are genuinely
# dimensionless (`Float64`).
function _variable_template(part_names::NTuple{N,Symbol}, temperature, conductivity) where {N}
    t0 = zero(temperature)
    k0 = zero(conductivity)
    part_leaves = ntuple(_ -> _part_leaf(t0, t0, k0, 0.0), Val(N))
    return _variable_structure(part_names, t0, 0.0, 0.0, part_leaves)
end

# Residual structure: per-part surface + skin residuals, whole-organism balance,
# and the Q10 slack LAST (so it is the single `≥ 0` constraint; see `_constraint_bounds`).
# This defines the constraint order that `_heat_balance_residuals!` writes and that
# `_problem_size` counts — one place, so the two cannot disagree.
@inline function _residual_structure(part_names::NTuple{N,Symbol}, part_residuals::NTuple{N},
        whole_organism_balance, q10) where {N}
    return (;
        parts = NamedTuple{part_names}(part_residuals),
        whole_organism_balance,
        q10,
    )
end

# Zero-valued residual prototype, for counting constraints via Flatten.
@inline _residual_template(part_names::NTuple{N,Symbol}) where {N} =
    _residual_structure(part_names,
        ntuple(_ -> (; surface_balance = 0.0, skin_temperature = 0.0), Val(N)), 0.0, 0.0)

"""
    MultipartNLPPacked

Packed multi-part NLP problem: the fixed per-part `part_surface_residuals` setups,
the whole-organism respiration inputs, and the variable template (a nested
`NamedTuple` of `Float64` leaves) that maps Ipopt's flat vector to named per-part
fields via `Flatten.reconstruct`. Built once per organism + environment; the
decision variables arrive as the `x` vector at solve time.
"""
struct MultipartNLPPacked{S<:Tuple,Names,W,SM,T}
    setups::S               # per-part part_surface_residuals setup NamedTuples (fixed)
    part_names::Names        # NTuple{N,Symbol}, body order
    whole::W                 # whole-organism respiration inputs (resp_pars, lung_mass, atmos, …)
    smoothing::SM
    variable_template::T     # nested NamedTuple prototype (Unitful leaves)
end

_number_of_parts(p::MultipartNLPPacked) = length(p.setups)

# Reconstruct the named, Unitful decision variables from Ipopt's flat Float64
# vector. Targeting `Real`, `reconstruct` puts each `Float64` back into its
# `Quantity`-typed slot, so units come straight off the template's types — no
# manual reattachment. The reconstructed leaves are concretely typed (units are in
# the type), so this is type-stable and Enzyme-differentiable, and `effectors` (the
# raw `Vector{Float64}`) is indexed directly — nothing dynamically built in the AD
# path.
@inline _reconstruct_variables(packed::MultipartNLPPacked, effectors) =
    Flatten.reconstruct(packed.variable_template, effectors, Real)

# Build the packed problem from the organism + environment. Mirrors
# `solve_multipart_metabolic_rate`'s setup construction (same `part_surface_setups`),
# and keeps the whole-organism respiration inputs alongside for the NLP residual.
function HeatExchange.nlp_pack(::MultipartNLP, organism::Organism, environment,
        initial_skin_temperature, initial_insulation_temperature;
        smoothing::HeatExchange.SmoothingStrategy = HeatExchange.HardBound())
    body = HeatExchange.body(organism)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    core_temperature = HeatExchange.metabolism_pars(organism).core_temperature

    setups = part_surface_setups(
        organism, environment_pars, environment_vars,
        core_temperature, initial_skin_temperature, initial_insulation_temperature; smoothing,
    )
    metab_pars = HeatExchange.metabolism_pars(organism)
    flesh_conductivity = HeatExchange.conduction_pars_internal(organism).flesh_conductivity
    whole = (;
        resp_pars              = HeatExchange.respiration_pars(organism),
        lung_mass              = _lung_mass(body, lung_part(organism)),
        atmos                  = AtmosphericConditions(environment_vars),
        gas_fractions          = environment_pars.gas_fractions,
        air_temperature        = environment_vars.air_temperature,
        minimum_metabolic_heat = metab_pars.metabolic_heat_flow,
        respire                = HeatExchange.options(organism).respire,
        # The unit the log-metabolic variable reattaches to on the way out — pulled
        # from the model, so `exp(log_metabolic_heat_flow) * heat_flow_unit` inverts
        # the log-transform without hardcoding `u"W"`.
        heat_flow_unit         = oneunit(metab_pars.metabolic_heat_flow),
    )
    names = part_names(body)
    template = _variable_template(names, core_temperature, flesh_conductivity)
    return MultipartNLPPacked(setups, names, whole, smoothing, template)
end

# Variable / constraint counts derived from the templates, never hand-counted.
function _problem_size(p::MultipartNLPPacked)
    nvar = length(Flatten.flatten(p.variable_template, Real))
    ncon = length(Flatten.flatten(_residual_template(p.part_names), Real))
    return nvar, ncon
end

# =============================================================================
# Residuals — Float64 effectors in → named Unitful variables (via reconstruct) →
# physics → Float64 residuals out (flattened from the residual template).
#
# The single unit boundary: `_reconstruct_variables` attaches units at the top;
# each residual leaf is dimensionless-ised by its own `oneunit` at the bottom. Per-part surface balance and
# skin-temperature residuals come from `part_surface_residuals`; the whole-organism
# respiration balance closes metabolic − respiration = Σ per-part internal heat,
# exactly as `solve_coupled_metabolic_rate` does at its converged point.
# =============================================================================
function _heat_balance_residuals!(packed::MultipartNLPPacked, residuals, effectors, p)
    vars = _reconstruct_variables(packed, effectors)
    # Fields arrive Unitful straight off the template — no manual unit attachment.
    # (`metabolic_heat_flow` is the exception by construction: the optimiser carries
    # its *log*, a genuinely dimensionless quantity, so `exp` then a model-derived
    # `heat_flow_unit` is the log-transform's inverse, not a hardcoded reattachment.)
    core_temperature    = vars.core_temperature
    metabolic_heat_flow = exp(vars.log_metabolic_heat_flow) * packed.whole.heat_flow_unit
    panting_rate        = vars.panting_rate

    part_vars = values(vars.parts)
    N = length(part_vars)

    # Walk the parts with a plain `for` loop, accumulating Σ internal heat and Σ skin
    # temperature in loop-local scalars and writing each part's two residuals in place
    # via a running index. A plain loop (unlike a `map`/`ntuple` closure that
    # reassigns a captured accumulator) creates no `Core.Box`: the accumulators are
    # loop-carried SSA values, isbits, which Enzyme differentiates cleanly. Parts are
    # homogeneous here, so `part_vars[i]` / `packed.setups[i]` are type-stable.
    net_internal_total = zero(metabolic_heat_flow)
    skin_sum           = zero(core_temperature)
    k = 0
    @inbounds for i in 1:N
        v = part_vars[i]
        r = part_surface_residuals(
            packed.setups[i], core_temperature, v.skin_temperature,
            v.insulation_temperature, metabolic_heat_flow;
            k_flesh = v.flesh_conductivity, pant = panting_rate,
            skin_wetness = v.skin_wetness,
            resp_pars = packed.whole.resp_pars, smoothing = packed.smoothing,
        )
        # Each residual is dimensionless-ised by its *own* unit (`x / oneunit(x)`),
        # not stripped to a hardcoded one: the physics fixes each residual's unit
        # (W, K) deterministically across the solve, so this is the same number a
        # canonical `ustrip` would give but carries no unit literal.
        residuals[k += 1] = r.surface_balance / oneunit(r.surface_balance)
        residuals[k += 1] = r.residual_skin_temperature / oneunit(r.residual_skin_temperature)
        net_internal_total += r.net_metabolic_heat_internal
        skin_sum           += v.skin_temperature
    end

    # Whole-organism balance, closed at the lung. With respiration (the default) the
    # metabolic rate must close the respiration balance; without it the organism just
    # produces its internal heat directly — mirroring `solve_coupled_metabolic_rate`'s
    # `respire` branch so the NLP and the iterative solver constrain the same thing.
    skin_mean        = skin_sum / N
    lung_temperature = (core_temperature + skin_mean) / 2
    balance = if packed.whole.respire
        resp_pars_panted = setproperties(packed.whole.resp_pars, (; pant = panting_rate))
        respiration(
            MetabolicRates(; metabolic = metabolic_heat_flow, sum = net_internal_total,
                            minimum = packed.whole.minimum_metabolic_heat),
            resp_pars_panted, packed.whole.atmos, packed.whole.lung_mass,
            lung_temperature, packed.whole.air_temperature;
            gas_fractions = packed.whole.gas_fractions, O2conversion = Kleiber1961(),
            smoothing = packed.smoothing,
        ).balance
    else
        metabolic_heat_flow - net_internal_total
    end

    # Q10: metabolic_heat_flow − minimum_heat_flow · q10^((core − setpoint)/10 K) ≥ 0.
    # The Q10 exponent is an empirical correlation over a temperature difference: the
    # 10-degree reference interval is `10 * oneunit(ΔT)`, so the exponent is
    # dimensionless without stripping the difference to a hardcoded unit.
    ΔT_core = core_temperature - p.setpoint_temperature
    q10_slack = metabolic_heat_flow -
        p.minimum_heat_flow * p.q10 ^ (ΔT_core / (10 * oneunit(ΔT_core)))

    # The per-part residuals were written inside the fold above; the whole-organism
    # balance and the Q10 slack come last (the single `≥ 0` constraint — see
    # `_constraint_bounds`). This order matches `_residual_structure`, which
    # `_problem_size` counts, so the write and the count cannot disagree.
    @inbounds residuals[k += 1] = balance / oneunit(balance)
    @inbounds residuals[k += 1] = q10_slack / oneunit(q10_slack)
    return nothing
end

# =============================================================================
# Objective — a weighted sum of squared, normalised deviations (Tikhonov-style
# regularisation pulling each regulated quantity toward a reference). The objective
# is a sum of *terms*, one per *regulation target*: a target is the recipe (which
# quantity, toward what reference, at what scale and weight); the term is the number
# it contributes. The set of targets is *data*, not code: `_inputs!` authors a tuple
# of `RegulationTarget`s and the objective folds over it, naming no variable of its
# own. Each target co-locates its value-selector, reference, normalisation scale, and
# weight, so the three can never drift out of step with the quantity they regularise
# (the old failure mode: a flat `objective_parameters` bag whose fields had to be
# kept in lockstep with a separate hardcoded sum), and adding or dropping a regulated
# quantity is an edit to that one tuple.
# =============================================================================

"""
    RegulationTarget(select, reference, scale, weight)

One soft target of the NLP objective: hold `select(vars)` near `reference`,
measured in units of `scale`, weighted by `weight`. Its contribution — the
objective *term* — is `weight · ((select(vars) − reference) / scale)²`. `select`
maps the reconstructed decision variables to the regulated quantity — a top-level
function, hence a singleton that captures nothing, so a tuple of targets is
`isbits` and differentiates cleanly under Enzyme. `reference` and `scale` are
`Unitful` (or dimensionless) and share a dimension, so the ratio is dimensionless
and the term evaluates to a plain `Float64`.
"""
struct RegulationTarget{Select,Reference,Scale}
    select::Select
    reference::Reference
    scale::Scale
    weight::Float64
end

# The objective term this target contributes.
@inline _term(target::RegulationTarget, vars) =
    target.weight * ((target.select(vars) - target.reference) / target.scale)^2

# Unrolled, type-stable fold over the (heterogeneous) tuple of targets — a hand-
# written recursion rather than `sum(f, targets)` so the whole reduction stays on the
# plain-`+` Julia path Enzyme differentiates through: no `mapreduce` machinery, no
# dynamically built values. Bottoms out at the empty tuple.
@inline _sum_terms(::Tuple{}, vars) = 0.0
@inline _sum_terms(targets::Tuple, vars) =
    _term(first(targets), vars) + _sum_terms(Base.tail(targets), vars)

# Value-selectors: decode one regulated quantity from the reconstructed variables.
# Each is a top-level function (a singleton), so the targets tuple stays isbits.
# `metabolic_heat_flow` regulates in the log-variable's own dimensionless space:
# `exp(log_metabolic_heat_flow)` is the metabolic magnitude in the model's heat-flow
# unit, so its target's reference/scale are the matching dimensionless magnitudes and
# no unit is reattached here. The gradient and skin-wetness targets regulate the
# whole-organism mean over parts.
_regulated_core_temperature(vars) = vars.core_temperature
_regulated_metabolic_heat_flow(vars) = exp(vars.log_metabolic_heat_flow)
_regulated_panting_rate(vars) = vars.panting_rate
_regulated_core_skin_gradient(vars) = vars.core_temperature - _mean_skin_temperature(vars)
_regulated_skin_wetness(vars) = _mean_skin_wetness(vars)

# Whole-organism means over the per-part leaves. Plain loop with a loop-carried
# isbits accumulator (parts are homogeneous) — Enzyme-safe, nothing dynamically built.
@inline function _mean_skin_temperature(vars)
    parts = values(vars.parts)
    total = parts[1].skin_temperature
    @inbounds for i in 2:length(parts)
        total += parts[i].skin_temperature
    end
    return total / length(parts)
end
@inline function _mean_skin_wetness(vars)
    parts = values(vars.parts)
    total = parts[1].skin_wetness
    @inbounds for i in 2:length(parts)
        total += parts[i].skin_wetness
    end
    return total / length(parts)
end

function _objective_value(packed::MultipartNLPPacked, effectors, opt)
    vars = _reconstruct_variables(packed, effectors)
    return _sum_terms(opt.targets, vars)
end

# Per-variable Ipopt scaling, authored in the variable structure and flattened:
# temperatures ~300 K, log-metabolic/conductivity/pant ~O(1), skin_wetness ~O(0.01).
function _scaling(packed::MultipartNLPPacked)
    N = _number_of_parts(packed)
    part_leaves = ntuple(_ -> _part_leaf(1 / 300.0, 1 / 300.0, 1.0, 50.0), Val(N))
    scaling_structure = _variable_structure(packed.part_names, 1 / 300.0, 1.0, 1.0, part_leaves)
    return collect(Float64, Flatten.flatten(scaling_structure, Real))
end

# =============================================================================
# Per-call inputs (bounds, initial values, objective + NLP parameters).
#
# Bounds and initials are authored as templates of the variable structure and
# flattened — guaranteeing alignment with the variables. Whole-organism scalars come
# from `ThermoregulationLimits`; per-part bounds broadcast those whole-organism
# limits to every part (the uniform-tuning default). The objective is authored as a
# tuple of `RegulationTarget`s (references/scales Unitful — the objective folds over
# them doing dimensioned arithmetic, stripping only at the empirical boundary); the shared
# Ipopt callbacks pass `objective_parameters` / `nlp_parameters` through opaquely, so
# nothing downstream needs to change.
# =============================================================================
function _inputs!(lower_bounds, upper_bounds, initial_values, packed::MultipartNLPPacked,
                  organism, environment, init::NamedTuple)
    N       = _number_of_parts(packed)
    limits  = thermoregulation(organism)
    metab   = metabolism_pars(organism)
    internal = conduction_pars_internal(organism)
    evap     = evaporation_pars(organism)

    # Quantities stay Unitful; the bound/initial structures are built Unitful from the
    # model's own limits and `flatten(Real)` rips out the magnitudes. Because the
    # variable template's leaf units come from the same model quantities (see
    # `_variable_template`), these magnitudes reconstruct back into the identical
    # units — no unit is fixed by hand here.
    air_temperature      = environment.environment_vars.air_temperature
    setpoint_temperature = metab.core_temperature
    core_temperature_min = limits.core_temperature.reference
    core_temperature_max = limits.core_temperature.max
    skin_temperature_min = air_temperature - limits.skin_temperature_undershoot
    skin_temperature_max = core_temperature_max + limits.skin_temperature_core_overshoot
    heat_flow_min        = limits.minimum_heat_flow
    heat_flow_max        = heat_flow_min * limits.metabolic_heat_flow_max_multiplier
    flesh_conductivity_min = limits.flesh_conductivity.reference
    flesh_conductivity_max = limits.flesh_conductivity.max
    panting_rate_min     = Float64(limits.panting.pant.reference)
    panting_rate_max     = Float64(limits.panting.pant.max)
    skin_wetness_min     = Float64(limits.skin_wetness.reference)
    skin_wetness_max     = Float64(limits.skin_wetness.max)
    heat_flow_init       = clamp(init.metabolic_heat_flow, heat_flow_min, heat_flow_max)
    flesh_conductivity   = clamp(internal.flesh_conductivity, flesh_conductivity_min, flesh_conductivity_max)
    skin_wetness         = clamp(Float64(evap.skin_wetness), skin_wetness_min, skin_wetness_max)
    skin_temperature_init = clamp(init.skin_temperature, skin_temperature_min, skin_temperature_max)
    insulation_temp_init  = clamp(init.insulation_temperature, skin_temperature_min, skin_temperature_max)

    # The metabolic variable is stored as its log, so its bounds live in log space:
    # `log(heat_flow / oneunit(heat_flow))` is the log of the magnitude in the model's
    # own heat-flow unit — the same space `_regulated_metabolic_heat_flow` decodes into
    # and `packed.whole.heat_flow_unit` reattaches from, with no unit literal.
    heat_flow_unit    = oneunit(heat_flow_min)
    log_heat_flow_min = log(heat_flow_min / heat_flow_unit)
    log_heat_flow_max = log(heat_flow_max / heat_flow_unit)

    lower_leaves = ntuple(_ -> _part_leaf(skin_temperature_min, skin_temperature_min, flesh_conductivity_min, skin_wetness_min), Val(N))
    upper_leaves = ntuple(_ -> _part_leaf(skin_temperature_max, skin_temperature_max, flesh_conductivity_max, skin_wetness_max), Val(N))
    init_leaves  = ntuple(_ -> _part_leaf(skin_temperature_init, insulation_temp_init, flesh_conductivity, skin_wetness), Val(N))

    lower_structure = _variable_structure(packed.part_names, core_temperature_min, log_heat_flow_min, panting_rate_min, lower_leaves)
    upper_structure = _variable_structure(packed.part_names, core_temperature_max, log_heat_flow_max, panting_rate_max, upper_leaves)
    init_structure  = _variable_structure(packed.part_names,
        clamp(setpoint_temperature, core_temperature_min, core_temperature_max),
        clamp(log(heat_flow_init / heat_flow_unit), log_heat_flow_min, log_heat_flow_max),
        panting_rate_min, init_leaves)

    copyto!(lower_bounds,   Flatten.flatten(lower_structure, Real))
    copyto!(upper_bounds,   Flatten.flatten(upper_structure, Real))
    copyto!(initial_values, Flatten.flatten(init_structure, Real))

    # Objective normalisation scales. `minimum_normalisation_range` is a dimensionless
    # epsilon that `oneunit` gives each range's dimension, so a degenerate (zero-span)
    # limit never divides by zero.
    range_floor(span) = max(span, limits.minimum_normalisation_range * oneunit(span))
    core_temperature_range = range_floor(core_temperature_max - setpoint_temperature)
    heat_flow_range        = range_floor(heat_flow_max - heat_flow_min)
    panting_rate_range     = range_floor(panting_rate_max - panting_rate_min)
    skin_wetness_range     = range_floor(skin_wetness_max - skin_wetness_min)

    # The objective as data: one `RegulationTarget` per regulated quantity, each pairing
    # a value-selector with its reference, normalisation scale, and weight. `_objective_value`
    # folds over this tuple naming nothing — adding or dropping a regulated quantity is an
    # edit here alone. Selectors are top-level singletons, references/scales Unitful, so
    # the tuple is isbits and each target's term folds to a dimensionless Float64.
    targets = (
        RegulationTarget(_regulated_core_temperature,    setpoint_temperature, core_temperature_range,
                         Float64(limits.core_temperature_weight)),
        RegulationTarget(_regulated_metabolic_heat_flow, heat_flow_min / heat_flow_unit,
                         heat_flow_range / heat_flow_unit,
                         Float64(limits.metabolic_heat_weight)),
        RegulationTarget(_regulated_core_skin_gradient,  limits.target_core_skin_gradient, core_temperature_range,
                         Float64(limits.gradient_weight)),
        RegulationTarget(_regulated_panting_rate,        panting_rate_min,     panting_rate_range,
                         Float64(limits.panting_weight)),
        RegulationTarget(_regulated_skin_wetness,        skin_wetness_min,     skin_wetness_range,
                         Float64(limits.skin_wetness_weight)),
    )
    objective_parameters = (; targets)
    nlp_parameters = (;
        nlp_packed           = packed,
        minimum_heat_flow    = heat_flow_min,
        q10                  = Float64(metab.q10),
        setpoint_temperature = setpoint_temperature,
    )
    return objective_parameters, nlp_parameters
end

# =============================================================================
# Output assembly — reconstruct the converged `x` into named per-part fields plus
# the whole-organism quantities, and evaluate the physics once at the solution for
# the per-part heat flows. Shape mirrors `solve_multipart_metabolic_rate`'s output.
# =============================================================================
function _assemble(packed::MultipartNLPPacked, organism, environment, x_sol)
    vars = _reconstruct_variables(packed, x_sol)
    core_temperature    = vars.core_temperature
    metabolic_heat_flow = exp(vars.log_metabolic_heat_flow) * packed.whole.heat_flow_unit
    panting_rate        = vars.panting_rate

    part_vars = values(vars.parts)
    N = length(part_vars)
    part_results = ntuple(N) do i
        v = part_vars[i]
        skin_temperature       = v.skin_temperature
        insulation_temperature = v.insulation_temperature
        flesh_conductivity     = v.flesh_conductivity
        r = part_surface_residuals(
            packed.setups[i], core_temperature, skin_temperature, insulation_temperature,
            metabolic_heat_flow;
            k_flesh = flesh_conductivity, pant = panting_rate, skin_wetness = v.skin_wetness,
            resp_pars = packed.whole.resp_pars, smoothing = packed.smoothing,
        )
        (; skin_temperature, insulation_temperature, flesh_conductivity, v.skin_wetness,
           net_metabolic = r.net_metabolic_heat_internal, flows = r.balance)
    end
    parts = NamedTuple{packed.part_names}(part_results)

    skin_mean       = sum(p -> p.skin_temperature, part_results) / N
    insulation_mean = sum(p -> p.insulation_temperature, part_results) / N
    return (;
        metabolic_heat_flow,
        core_temperature,
        panting_rate,
        parts,
        net_metabolic_total = sum(p -> p.net_metabolic, part_results),
        skin_temperature = skin_mean,
        insulation_temperature = insulation_mean,
        lung_temperature = (core_temperature + skin_mean) / 2,
    )
end
