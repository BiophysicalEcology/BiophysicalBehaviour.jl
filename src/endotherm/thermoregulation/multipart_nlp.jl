# =============================================================================
# Multi-part NLP formulation for `IPOPTControl` (Phase 7.5).
#
# The dorsal/ventral `WeightedMeanNLP` / `MultiSidedNLP` strategies (ipopt.jl)
# each hardcode a fixed variable layout for one lumped body. `MultipartNLP`
# generalises the NLP to a genuine multi-part organism: one regulated core shared
# across all parts, per-part surface temperatures and effectors, and a single
# whole-organism respiration balance closed at the lung.
#
# It reuses the *same* direct-Ipopt + Enzyme callback machinery as the single-body
# strategies (ipopt.jl): `IPOPTSolverCache`, `_ipopt_solve!`, the five
# `Evaluate*` functors, and `_lagrangian` are all generic over the packed type and
# operate on a flat `Vector{Float64}`. This file only adds the `MultipartNLPPacked`
# methods those functors dispatch to — `_heat_balance_residuals!`, `_objective_value`,
# `_inputs!`, `_problem_size`, `_scaling`, `_assemble` — plus a `nlp_pack` method.
# The physics is `HeatExchange.part_surface_residuals` per part (the non-iterative
# twin of `solve_part_surface`), so zero physics is duplicated.
#
# Decision-variable layout (regular per-part stride, `N` = number of parts):
#   x[1]              core_temperature            (K)
#   x[2]              log_metabolic_heat_flow     (log W — the log-transform invariant)
#   x[3]              panting_rate                (dimensionless)
#   per part i, base = 3 + 4(i-1):
#     x[base+1]       skin_temperature            (K)
#     x[base+2]       insulation_temperature      (K)
#     x[base+3]       flesh_conductivity          (W/m/K)
#     x[base+4]       skin_wetness                (dimensionless)
# → 3 + 4N variables.
#
# Constraint layout:
#   per part i:  residuals[2i-1] surface balance (W), residuals[2i] skin temp (K)
#   residuals[2N+1]  whole-organism respiration balance (W)  [= 0]
#   residuals[2N+2]  Q10 inequality                          [>= 0]
# → 2N + 2 constraints (2N + 1 equality + 1 Q10 inequality).
#
# Geometry effectors (insulation_depth / axis_ratio) are deferred, exactly as the
# multi-part rule-based ladder defers Piloerect / Uncurl: no per-part body rebuild
# happens inside the residual, so the per-part `setups` are fixed across the solve.
# =============================================================================

"""
    MultipartNLP <: HeatExchange.NLPStrategy

NLP strategy for a multi-part organism: shared regulated core, per-part surface
temperatures and effectors (flesh_conductivity, skin_wetness), whole-organism
respiration closed at the lung. Select via `IPOPTControl(; nlp_strategy = MultipartNLP())`.
Works for any `Organism` — a single `Body` is treated as the one-part case.
"""
struct MultipartNLP <: HeatExchange.NLPStrategy end

"""
    MultipartNLPPacked

Packed multi-part NLP problem: the fixed per-part `solve_part_surface`/`part_surface_residuals`
setups plus the whole-organism respiration inputs. Decision variables (skin,
insulation, flesh_conductivity, skin_wetness per part; core, log-metabolic, pant
whole-organism) arrive as the Ipopt `x` vector at solve time, so nothing here
depends on them — the pack is built once per organism + environment.
"""
struct MultipartNLPPacked{S<:Tuple,Names,W,SM}
    setups::S               # per-part setup NamedTuples (fixed geometry/insulation/env)
    part_names::Names       # NTuple{N,Symbol}, body order
    whole::W                # whole-organism respiration inputs (resp_pars, lung_mass, atmos, …)
    smoothing::SM
end

_number_of_parts(p::MultipartNLPPacked) = length(p.setups)

# Build the packed problem from the organism + environment. Mirrors
# `solve_multipart_metabolic_rate`'s setup construction (same `part_surface_setups`),
# but keeps the whole-organism respiration inputs alongside for the NLP residual.
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
    return MultipartNLPPacked(setups, part_names(body), whole, smoothing)
end

# 3 whole-organism variables + 4 per part; 2 residuals per part + whole-organism
# balance + Q10 inequality.
_problem_size(p::MultipartNLPPacked) = (3 + 4 * _number_of_parts(p), 2 * _number_of_parts(p) + 2)

# =============================================================================
# Residuals — Float64 effectors in → Unitful physics → Float64 residuals out.
#
# The single boundary crossing: units are attached at the top (decode `x`),
# stripped at the bottom (`ustrip` each residual). Per-part surface balance and
# skin-temperature residuals come from `part_surface_residuals`; the whole-organism
# respiration balance closes metabolic − respiration = Σ per-part internal heat,
# exactly as `solve_coupled_metabolic_rate` does at its converged point.
# =============================================================================
function _heat_balance_residuals!(packed::MultipartNLPPacked, residuals, effectors, p)
    N = _number_of_parts(packed)
    core_temperature    = effectors[1] * u"K"
    metabolic_heat_flow = exp(effectors[2]) * u"W"
    panting_rate        = effectors[3]

    net_internal_total = zero(metabolic_heat_flow)
    skin_sum           = zero(core_temperature)
    @inbounds for i in 1:N
        base = 3 + 4 * (i - 1)
        skin_temperature       = effectors[base + 1] * u"K"
        insulation_temperature = effectors[base + 2] * u"K"
        flesh_conductivity     = effectors[base + 3] * u"W/m/K"
        skin_wetness           = effectors[base + 4]

        r = part_surface_residuals(
            packed.setups[i], core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow;
            k_flesh = flesh_conductivity, pant = panting_rate, skin_wetness,
            resp_pars = packed.whole.resp_pars, smoothing = packed.smoothing,
        )
        residuals[2 * (i - 1) + 1] = ustrip(u"W", r.surface_balance)
        residuals[2 * (i - 1) + 2] = ustrip(u"K", r.residual_skin_temperature)
        net_internal_total += r.net_metabolic_heat_internal
        skin_sum           += skin_temperature
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
    residuals[2 * N + 1] = ustrip(u"W", balance)
    # Q10: metabolic_heat_flow >= minimum_heat_flow · q10^((core − setpoint)/10)
    residuals[2 * N + 2] = exp(effectors[2]) -
        p.minimum_heat_flow * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
    return nothing
end

# =============================================================================
# Objective — the same five-term quadratic-penalty objective as the single-body
# strategies, over the whole-organism core / metabolic / panting decision variables
# and the *mean* per-part skin temperature and skin wetness.
# =============================================================================
function _objective_value(packed::MultipartNLPPacked, effectors, opt)
    N = _number_of_parts(packed)
    core_temperature    = effectors[1]
    metabolic_heat_flow = exp(effectors[2])
    panting_rate        = effectors[3]

    skin_sum = 0.0
    wet_sum  = 0.0
    @inbounds for i in 1:N
        base = 3 + 4 * (i - 1)
        skin_sum += effectors[base + 1]
        wet_sum  += effectors[base + 4]
    end
    skin_mean         = skin_sum / N
    skin_wetness_mean = wet_sum / N

    return opt.core_temperature_penalty * ((core_temperature - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
           opt.metabolic_heat_penalty   * ((metabolic_heat_flow - opt.heat_flow_min_W) / opt.heat_flow_range)^2 +
           opt.gradient_penalty         * ((core_temperature - skin_mean - opt.target_gradient) / opt.gradient_range)^2 +
           opt.panting_penalty          * ((panting_rate - opt.panting_rate_min) / opt.panting_rate_range)^2 +
           opt.skin_wetness_penalty     * ((skin_wetness_mean - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# Per-variable Ipopt scaling: temperatures ~300 K, log-metabolic/conductivity/pant
# ~O(1), skin_wetness ~O(0.01).
function _scaling(p::MultipartNLPPacked)
    N = _number_of_parts(p)
    x_scaling = ones(3 + 4 * N)
    x_scaling[1] = 1 / 300.0                 # core_temperature
    x_scaling[2] = 1.0                       # log_metabolic_heat_flow
    x_scaling[3] = 1.0                       # panting_rate
    @inbounds for i in 1:N
        base = 3 + 4 * (i - 1)
        x_scaling[base + 1] = 1 / 300.0      # skin_temperature
        x_scaling[base + 2] = 1 / 300.0      # insulation_temperature
        x_scaling[base + 3] = 1.0            # flesh_conductivity
        x_scaling[base + 4] = 50.0           # skin_wetness
    end
    return x_scaling
end

# =============================================================================
# Per-call inputs (bounds, initial values, objective + NLP parameters).
#
# Whole-organism scalars come from `ThermoregulationLimits`; per-part bounds
# broadcast those whole-organism limits to every part (the uniform-tuning default).
# The objective_parameters / nlp_parameters NamedTuples are structurally identical
# to the single-body path, so the shared Ipopt callbacks read them unchanged.
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
    flesh_conductivity     = ustrip(u"W/m/K", internal.flesh_conductivity)
    skin_wetness           = Float64(evap.skin_wetness)
    skin_temperature_init  = clamp(ustrip(u"K", init.skin_temperature), skin_temperature_min, skin_temperature_max)
    insulation_temp_init   = clamp(ustrip(u"K", init.insulation_temperature), skin_temperature_min, skin_temperature_max)

    # Whole-organism variables.
    lower_bounds[1] = core_temperature_min; upper_bounds[1] = core_temperature_max
    lower_bounds[2] = log(heat_flow_min_W); upper_bounds[2] = log(heat_flow_max_W)
    lower_bounds[3] = panting_rate_min;     upper_bounds[3] = panting_rate_max
    initial_values[1] = clamp(setpoint_temperature_K, core_temperature_min, core_temperature_max)
    initial_values[2] = clamp(log(heat_flow_init), lower_bounds[2], upper_bounds[2])
    initial_values[3] = clamp(panting_rate_min, panting_rate_min, panting_rate_max)

    # Per-part variables (whole-organism limits broadcast to every part).
    @inbounds for i in 1:N
        base = 3 + 4 * (i - 1)
        lower_bounds[base + 1] = skin_temperature_min;   upper_bounds[base + 1] = skin_temperature_max
        lower_bounds[base + 2] = skin_temperature_min;   upper_bounds[base + 2] = skin_temperature_max
        lower_bounds[base + 3] = flesh_conductivity_min; upper_bounds[base + 3] = flesh_conductivity_max
        lower_bounds[base + 4] = skin_wetness_min;       upper_bounds[base + 4] = skin_wetness_max
        initial_values[base + 1] = skin_temperature_init
        initial_values[base + 2] = insulation_temp_init
        initial_values[base + 3] = clamp(flesh_conductivity, flesh_conductivity_min, flesh_conductivity_max)
        initial_values[base + 4] = clamp(skin_wetness, skin_wetness_min, skin_wetness_max)
    end

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
# Output assembly — decode the converged `x` into a per-part NamedTuple plus the
# whole-organism quantities, and evaluate the physics once at the solution for the
# per-part heat flows. Shape mirrors `solve_multipart_metabolic_rate`'s output.
# =============================================================================
function _assemble(packed::MultipartNLPPacked, organism, environment, x_sol)
    N = _number_of_parts(packed)
    core_temperature    = x_sol[1] * u"K"
    metabolic_heat_flow = exp(x_sol[2]) * u"W"
    panting_rate        = x_sol[3]

    part_results = ntuple(N) do i
        base = 3 + 4 * (i - 1)
        skin_temperature       = x_sol[base + 1] * u"K"
        insulation_temperature = x_sol[base + 2] * u"K"
        flesh_conductivity     = x_sol[base + 3] * u"W/m/K"
        skin_wetness           = x_sol[base + 4]
        r = part_surface_residuals(
            packed.setups[i], core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow;
            k_flesh = flesh_conductivity, pant = panting_rate, skin_wetness,
            resp_pars = packed.whole.resp_pars, smoothing = packed.smoothing,
        )
        (; skin_temperature, insulation_temperature, flesh_conductivity, skin_wetness,
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
