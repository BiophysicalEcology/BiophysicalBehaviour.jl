# =============================================================================
# Multi-part NLP formulation for `IPOPTControl` (Phase 7.5).
#
# The dorsal/ventral `WeightedMeanNLP` / `MultiSidedNLP` strategies (ipopt.jl)
# each hardcode a fixed variable layout for one lumped body. `MultipartNLP`
# generalises the NLP to a genuine multi-part organism: one regulated core shared
# across all parts, per-part surface temperatures and effectors, and a single
# whole-organism respiration balance closed at the lung.
#
# The decision-variable layout is NOT hand-indexed. It is a nested `NamedTuple`
# *template* built from the organism's parts, and `Flatten.flatten` /
# `Flatten.reconstruct` bridge between that named structure and the flat
# `Vector{Float64}` Ipopt optimises. Variables, bounds, scaling, and residuals are
# every one of them a template of the same shape flattened in the same order — so
# a bound can never drift out of alignment with the variable it bounds, and adding
# a part or a per-part effector is a change to the *structure*, never to a pile of
# `x[3 + 4i]` offsets. Enzyme differentiates straight through `reconstruct`, so the
# residual/objective read named fields off the reconstructed template.
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

# The structural prototype for `reconstruct`. All leaves are `Float64`: the flat
# optimiser vector Ipopt hands us is unitless (Float64 is the optimiser space), so
# the reconstructed template is unitless too and stays homogeneous — which keeps
# `reconstruct` type-stable (a heterogeneous `data[n]` would infer to a Union and
# break Enzyme). Units are a physics-boundary concern, attached where each named
# field enters the physics (`* u"K"`, `exp(_) * u"W"`, …), never carried in the
# optimiser vector.
function _variable_template(part_names::NTuple{N,Symbol}) where {N}
    part_leaves = ntuple(_ -> _part_leaf(0.0, 0.0, 0.0, 0.0), Val(N))
    return _variable_structure(part_names, 0.0, 0.0, 0.0, part_leaves)
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
    variable_template::T     # nested NamedTuple prototype (Float64 leaves)
end

_number_of_parts(p::MultipartNLPPacked) = length(p.setups)

# Reconstruct the named decision variables from Ipopt's flat Float64 vector. The
# template's Float64 leaves make `data[n]` concrete, so `reconstruct` is
# type-stable (and Enzyme-differentiable); `effectors` is indexed directly, so
# nothing is dynamically built in the AD path. Units are attached at the point each
# field enters the physics, not here.
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
    whole = (;
        resp_pars              = HeatExchange.respiration_pars(organism),
        lung_mass              = _lung_mass(body, lung_part(organism)),
        atmos                  = AtmosphericConditions(environment_vars),
        gas_fractions          = environment_pars.gas_fractions,
        air_temperature        = environment_vars.air_temperature,
        minimum_metabolic_heat = metab_pars.metabolic_heat_flow,
        respire                = HeatExchange.options(organism).respire,
    )
    names = part_names(body)
    return MultipartNLPPacked(setups, names, whole, smoothing, _variable_template(names))
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
# each residual leaf is `ustrip`ped at the bottom. Per-part surface balance and
# skin-temperature residuals come from `part_surface_residuals`; the whole-organism
# respiration balance closes metabolic − respiration = Σ per-part internal heat,
# exactly as `solve_coupled_metabolic_rate` does at its converged point.
# =============================================================================
function _heat_balance_residuals!(packed::MultipartNLPPacked, residuals, effectors, p)
    vars = _reconstruct_variables(packed, effectors)
    # Attach units as each named Float64 field enters the physics.
    core_temperature    = vars.core_temperature * u"K"
    metabolic_heat_flow = exp(vars.log_metabolic_heat_flow) * u"W"
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
            packed.setups[i], core_temperature, v.skin_temperature * u"K",
            v.insulation_temperature * u"K", metabolic_heat_flow;
            k_flesh = v.flesh_conductivity * u"W/m/K", pant = panting_rate,
            skin_wetness = v.skin_wetness,
            resp_pars = packed.whole.resp_pars, smoothing = packed.smoothing,
        )
        residuals[k += 1] = ustrip(u"W", r.surface_balance)
        residuals[k += 1] = ustrip(u"K", r.residual_skin_temperature)
        net_internal_total += r.net_metabolic_heat_internal
        skin_sum           += v.skin_temperature * u"K"
    end

    # Whole-organism respiration balance, closed at the lung (mean skin → lung temp).
    skin_mean        = skin_sum / N
    lung_temperature = (core_temperature + skin_mean) / 2
    resp_pars_panted = setproperties(packed.whole.resp_pars, (; pant = panting_rate))
    balance = respiration(
        MetabolicRates(; metabolic = metabolic_heat_flow, sum = net_internal_total,
                        minimum = packed.whole.minimum_metabolic_heat),
        resp_pars_panted, packed.whole.atmos, packed.whole.lung_mass,
        lung_temperature, packed.whole.air_temperature;
        gas_fractions = packed.whole.gas_fractions, O2conversion = Kleiber1961(),
        smoothing = packed.smoothing,
    ).balance

    # Q10: metabolic_heat_flow − minimum_heat_flow · q10^((core − setpoint)/10) ≥ 0.
    q10_slack = exp(vars.log_metabolic_heat_flow) -
        p.minimum_heat_flow * p.q10 ^ ((ustrip(u"K", core_temperature) - p.setpoint_temperature) / 10.0)

    # The per-part residuals were written inside the fold above; the whole-organism
    # balance and the Q10 slack come last (the single `≥ 0` constraint — see
    # `_constraint_bounds`). This order matches `_residual_structure`, which
    # `_problem_size` counts, so the write and the count cannot disagree.
    @inbounds residuals[k += 1] = ustrip(u"W", balance)
    @inbounds residuals[k += 1] = q10_slack
    return nothing
end

# =============================================================================
# Objective — the same five-term quadratic-penalty objective as the single-body
# strategies, over the whole-organism core / metabolic / panting decision variables
# and the *mean* per-part skin temperature and skin wetness. Reads named fields off
# the reconstructed template.
# =============================================================================
function _objective_value(packed::MultipartNLPPacked, effectors, opt)
    vars = _reconstruct_variables(packed, effectors)
    core_temperature    = vars.core_temperature          # already a Float64 in K's numeric space
    metabolic_heat_flow = exp(vars.log_metabolic_heat_flow)
    panting_rate        = vars.panting_rate

    part_vars = values(vars.parts)
    N = length(part_vars)
    skin_sum = zero(core_temperature)
    wet_sum  = zero(panting_rate)
    @inbounds for i in 1:N
        skin_sum += part_vars[i].skin_temperature
        wet_sum  += part_vars[i].skin_wetness
    end
    skin_mean         = skin_sum / N
    skin_wetness_mean = wet_sum / N

    return opt.core_temperature_penalty * ((core_temperature - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
           opt.metabolic_heat_penalty   * ((metabolic_heat_flow - opt.heat_flow_min_W) / opt.heat_flow_range)^2 +
           opt.gradient_penalty         * ((core_temperature - skin_mean - opt.target_gradient) / opt.gradient_range)^2 +
           opt.panting_penalty          * ((panting_rate - opt.panting_rate_min) / opt.panting_rate_range)^2 +
           opt.skin_wetness_penalty     * ((skin_wetness_mean - opt.skin_wetness_min) / opt.skin_wetness_range)^2
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
# limits to every part (the uniform-tuning default). The objective_parameters /
# nlp_parameters NamedTuples are structurally identical to the single-body path, so
# the shared Ipopt callbacks read them unchanged.
# =============================================================================
function _inputs!(lower_bounds, upper_bounds, initial_values, packed::MultipartNLPPacked,
                  organism, environment, init::NamedTuple)
    N       = _number_of_parts(packed)
    limits  = thermoregulation(organism)
    metab   = metabolism_pars(organism)
    internal = conduction_pars_internal(organism)
    evap     = evaporation_pars(organism)

    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metab.core_temperature)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - ustrip(u"K", limits.skin_temperature_undershoot)
    skin_temperature_max   = core_temperature_max + ustrip(u"K", limits.skin_temperature_core_overshoot)
    heat_flow_min_W        = ustrip(u"W", limits.minimum_heat_flow)
    heat_flow_max_W        = heat_flow_min_W * limits.metabolic_heat_flow_max_multiplier
    flesh_conductivity_min = ustrip(u"W/m/K", limits.flesh_conductivity.reference)
    flesh_conductivity_max = ustrip(u"W/m/K", limits.flesh_conductivity.max)
    panting_rate_min       = Float64(limits.panting.pant.reference)
    panting_rate_max       = Float64(limits.panting.pant.max)
    skin_wetness_min       = Float64(limits.skin_wetness.reference)
    skin_wetness_max       = Float64(limits.skin_wetness.max)
    heat_flow_init         = max(ustrip(u"W", init.metabolic_heat_flow), heat_flow_min_W)
    flesh_conductivity     = clamp(ustrip(u"W/m/K", internal.flesh_conductivity), flesh_conductivity_min, flesh_conductivity_max)
    skin_wetness           = clamp(Float64(evap.skin_wetness), skin_wetness_min, skin_wetness_max)
    skin_temperature_init  = clamp(ustrip(u"K", init.skin_temperature), skin_temperature_min, skin_temperature_max)
    insulation_temp_init   = clamp(ustrip(u"K", init.insulation_temperature), skin_temperature_min, skin_temperature_max)

    lower_leaves = ntuple(_ -> _part_leaf(skin_temperature_min, skin_temperature_min, flesh_conductivity_min, skin_wetness_min), Val(N))
    upper_leaves = ntuple(_ -> _part_leaf(skin_temperature_max, skin_temperature_max, flesh_conductivity_max, skin_wetness_max), Val(N))
    init_leaves  = ntuple(_ -> _part_leaf(skin_temperature_init, insulation_temp_init, flesh_conductivity, skin_wetness), Val(N))

    lower_structure = _variable_structure(packed.part_names, core_temperature_min, log(heat_flow_min_W), panting_rate_min, lower_leaves)
    upper_structure = _variable_structure(packed.part_names, core_temperature_max, log(heat_flow_max_W), panting_rate_max, upper_leaves)
    init_structure  = _variable_structure(packed.part_names,
        clamp(setpoint_temperature_K, core_temperature_min, core_temperature_max),
        clamp(log(heat_flow_init), log(heat_flow_min_W), log(heat_flow_max_W)),
        panting_rate_min, init_leaves)

    copyto!(lower_bounds,   Flatten.flatten(lower_structure, Real))
    copyto!(upper_bounds,   Flatten.flatten(upper_structure, Real))
    copyto!(initial_values, Flatten.flatten(init_structure, Real))

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, limits.minimum_normalisation_range)
    objective_parameters = (;
        setpoint_temperature_K   = Float64(setpoint_temperature_K),
        heat_flow_min_W          = Float64(heat_flow_min_W),
        core_temperature_penalty = Float64(limits.core_temperature_penalty),
        metabolic_heat_penalty   = Float64(limits.metabolic_heat_penalty),
        panting_penalty          = Float64(limits.panting_penalty),
        skin_wetness_penalty     = Float64(limits.skin_wetness_penalty),
        gradient_penalty         = Float64(limits.gradient_penalty),
        target_gradient          = Float64(limits.target_core_skin_gradient),
        core_temperature_range   = Float64(core_temperature_range),
        gradient_range           = Float64(core_temperature_range),
        panting_rate_min         = Float64(panting_rate_min),
        panting_rate_range       = max(panting_rate_max - panting_rate_min, 1e-6),
        skin_wetness_min         = Float64(skin_wetness_min),
        skin_wetness_range       = max(skin_wetness_max - skin_wetness_min, 1e-6),
        heat_flow_range          = max(heat_flow_max_W - heat_flow_min_W, 1.0),
    )
    nlp_parameters = (;
        nlp_packed           = packed,
        minimum_heat_flow    = Float64(heat_flow_min_W),
        q10                  = Float64(metab.q10),
        setpoint_temperature = Float64(setpoint_temperature_K),
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
    core_temperature    = vars.core_temperature * u"K"
    metabolic_heat_flow = exp(vars.log_metabolic_heat_flow) * u"W"
    panting_rate        = vars.panting_rate

    part_vars = values(vars.parts)
    N = length(part_vars)
    part_results = ntuple(N) do i
        v = part_vars[i]
        skin_temperature       = v.skin_temperature * u"K"
        insulation_temperature = v.insulation_temperature * u"K"
        flesh_conductivity     = v.flesh_conductivity * u"W/m/K"
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
