# =============================================================================
# IPOPT-based endotherm thermoregulation control.
#
# Solver policy lives here: objective, variable bounds, initial effectors, and
# the Q10 metabolic-scaling inequality constraint.
# All physics (packing, per-iteration residuals, output assembly) is delegated
# to HeatExchange.nlp_pack / nlp_residuals / nlp_assemble_output.
#
# Two NLP formulations are supported via IPOPTControl.nlp_strategy:
#   WeightedMeanNLP  — 9 variables, 4 constraints (3 equality + 1 Q10)
#   MultiSidedNLP    — 11 variables, 6 constraints (5 equality + 1 Q10)
#
# Both execution paths talk to Ipopt.jl directly — no Optimization.jl /
# OptimizationIpopt wrappers. Going direct exposes `mult_g`, `mult_x_L`,
# `mult_x_U` on `IpoptProblem`, which is what makes dual warm-starting
# possible.
#
# Two execution paths:
#   1. `thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, init)`
#      — fresh per-call: allocates a single-use cache, runs one cold solve,
#      throws everything away. Use when each call has a structurally
#      different problem (different shape, different environment_vars layout).
#   2. `IPOPTSolverCache(control, organism, environment, init)`
#      then `thermoregulate(...; cache)` — warm-started: callback closures and
#      scratch buffers are reused, and primal+dual from the previous solve
#      are restored into the next `IpoptProblem` before `IpoptSolve` runs
#      with `warm_start_init_point = "yes"`. Use across a sweep where only
#      environment values change.
# =============================================================================

# =============================================================================
# Hard boundary between unitless optimizer world and Unitful scientific world.
#
#   `_objective_value` lives in optimizer world: pure Float64 arithmetic over
#       stripped reference values and dimensionless penalty weights.
#       Never sees a Unitful quantity.
#
#   `_heat_balance_residuals!` is the *only* function that crosses the
#       boundary: Float64 effectors in, Float64 residuals out, Unitful science
#       in between. Units are attached at the top, stripped at the bottom.
# =============================================================================

# =============================================================================
# WeightedMeanNLP residual / objective
#
# Variables (9):
#   x[1] = core_temperature (K)
#   x[2] = skin_temperature (K)
#   x[3] = insulation_temperature (K)
#   x[4] = log(metabolic_heat_flow) (log W)
#   x[5] = flesh_conductivity (W/m/K)
#   x[6] = panting_rate (dimensionless)
#   x[7] = skin_wetness (dimensionless)
#   x[8] = insulation_depth (m)
#   x[9] = axis_ratio_b (dimensionless)
#
# Constraints (4): residuals[1:3] == 0, residuals[4] >= 0 (Q10)
# =============================================================================

function _objective_value(::HeatExchange.WeightedMeanNLPPacked, effectors, opt)
    core_temperature    = effectors[1]
    skin_temperature    = effectors[2]
    metabolic_heat_flow = exp(effectors[4])
    panting_rate        = effectors[6]
    skin_wetness        = effectors[7]

    return opt.core_temperature_penalty * ((core_temperature - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
           opt.metabolic_heat_penalty   * ((metabolic_heat_flow - opt.heat_flow_min_W) / opt.heat_flow_range)^2 +
           opt.gradient_penalty         * ((core_temperature - skin_temperature - opt.target_gradient) / opt.gradient_range)^2 +
           opt.panting_penalty          * ((panting_rate - 1.0) / opt.panting_rate_range)^2 +
           opt.skin_wetness_penalty     * ((skin_wetness - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# Float64 effectors in → Unitful physics (via HeatExchange.nlp_residuals) → Float64 residuals out.
# residuals[1:3] == 0 (energy balance, internal conduction, skin temp).
# residuals[4]  >= 0 (Q10: metabolic_heat_flow >= minimum_heat_flow · q10^((T_core − T_setpoint)/10)).
function _heat_balance_residuals!(nlp_packed::HeatExchange.WeightedMeanNLPPacked, residuals, effectors, p)
    core_temperature       = effectors[1] * u"K"
    skin_temperature       = effectors[2] * u"K"
    insulation_temperature = effectors[3] * u"K"
    metabolic_heat_flow    = exp(effectors[4]) * u"W"
    flesh_conductivity     = effectors[5] * u"W/m/K"
    panting_rate           = effectors[6]
    skin_wetness           = effectors[7]
    insulation_depth       = effectors[8] * u"m"
    axis_ratio_b           = effectors[9]

    result = HeatExchange.nlp_residuals(nlp_packed, core_temperature, skin_temperature,
        insulation_temperature, metabolic_heat_flow, flesh_conductivity,
        panting_rate, skin_wetness, insulation_depth, axis_ratio_b)
    residuals[1] = ustrip(u"W", result.residuals[1])   # energy_balance
    residuals[2] = ustrip(u"W", result.residuals[2])   # internal_conduction
    residuals[3] = ustrip(u"K", result.residuals[3])   # skin_temperature
    residuals[4] = exp(effectors[4]) - p.minimum_heat_flow * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
    return nothing
end

# =============================================================================
# MultiSidedNLP residual / objective
#
# Variables (11):
#   x[1]  = core_temperature (K)
#   x[2]  = dorsal_skin_temperature (K)
#   x[3]  = dorsal_insulation_temperature (K)
#   x[4]  = ventral_skin_temperature (K)
#   x[5]  = ventral_insulation_temperature (K)
#   x[6]  = log(metabolic_heat_flow) (log W)
#   x[7]  = flesh_conductivity (W/m/K)
#   x[8]  = panting_rate (dimensionless)
#   x[9]  = skin_wetness (dimensionless)
#   x[10] = insulation_depth (m)
#   x[11] = axis_ratio_b (dimensionless)
#
# Constraints (6): residuals[1:5] == 0, residuals[6] >= 0 (Q10)
# =============================================================================

function _objective_value(::HeatExchange.MultiSidedNLPPacked, effectors, opt)
    core_temperature         = effectors[1]
    dorsal_skin_temperature  = effectors[2]
    ventral_skin_temperature = effectors[4]
    metabolic_heat_flow      = exp(effectors[6])
    panting_rate             = effectors[8]
    skin_wetness             = effectors[9]

    skin_temperature_mean = (dorsal_skin_temperature + ventral_skin_temperature) / 2

    return opt.core_temperature_penalty * ((core_temperature - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
           opt.metabolic_heat_penalty   * ((metabolic_heat_flow - opt.heat_flow_min_W) / opt.heat_flow_range)^2 +
           opt.gradient_penalty         * ((core_temperature - skin_temperature_mean - opt.target_gradient) / opt.gradient_range)^2 +
           opt.panting_penalty          * ((panting_rate - 1.0) / opt.panting_rate_range)^2 +
           opt.skin_wetness_penalty     * ((skin_wetness - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# residuals[1:5] == 0 (dorsal surface, dorsal skin, ventral surface, ventral skin, whole-organism).
# residuals[6]  >= 0 (Q10: metabolic_heat_flow >= minimum_heat_flow · q10^((T_core − T_setpoint)/10)).
function _heat_balance_residuals!(nlp_packed::HeatExchange.MultiSidedNLPPacked, residuals, effectors, p)
    # TODO automate this variable mapping
    core_temperature              = effectors[1]  * u"K"
    dorsal_skin_temperature       = effectors[2]  * u"K"
    dorsal_insulation_temperature = effectors[3]  * u"K"
    ventral_skin_temperature      = effectors[4]  * u"K"
    ventral_insulation_temperature = effectors[5] * u"K"
    metabolic_heat_flow           = exp(effectors[6]) * u"W"
    flesh_conductivity            = effectors[7]  * u"W/m/K"
    panting_rate                  = effectors[8]
    skin_wetness                  = effectors[9]
    insulation_depth              = effectors[10] * u"m"
    axis_ratio_b                  = effectors[11]

    result = HeatExchange.nlp_residuals(nlp_packed, core_temperature,
        dorsal_skin_temperature, dorsal_insulation_temperature,
        ventral_skin_temperature, ventral_insulation_temperature,
        metabolic_heat_flow, flesh_conductivity, panting_rate, skin_wetness,
        insulation_depth, axis_ratio_b)
    residuals[1] = ustrip(u"W", result.residuals[1])   # dorsal surface balance
    residuals[2] = ustrip(u"K", result.residuals[2])   # dorsal skin temperature
    residuals[3] = ustrip(u"W", result.residuals[3])   # ventral surface balance
    residuals[4] = ustrip(u"K", result.residuals[4])   # ventral skin temperature
    residuals[5] = ustrip(u"W", result.residuals[5])   # whole-organism balance
    residuals[6] = exp(effectors[6]) - p.minimum_heat_flow * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
    return nothing
end

# =============================================================================
# Lagrangian — `sigma * f(x) + dot(lambda, g(x))`. IPOPT's `evaluate_hessian` callback
# wants the Hessian of *this*, not of the objective alone. We define it as a
# single scalar function so a forward-over-reverse Enzyme pass (one per
# decision-variable direction) yields one Hessian column at a time.
#
# The residual buffer `residual_buffer` is passed in explicitly rather than allocated
# inside. Under the inner reverse pass we mark it as `Duplicated(residual_buffer, residual_adjoint)`
# so Enzyme writes adjoints into the supplied `residual_adjoint` rather than allocating
# a fresh shadow each call. Under the outer forward pass, both `residual_buffer` and
# `residual_adjoint` are themselves marked `Duplicated` so their forward tangents
# propagate. All four arrays live on the cache; the callbacks just thread
# them through.
#
# `p` is a NamedTuple with fields `sigma::Float64`, `lambda::Vector{Float64}`,
# `ipopt_parameters::IPOPTParameters`. The mutable reference lets the Lagrangian
# see the *current* `objective_parameters` / `nlp_parameters` set by `_inputs!`
# — no per-call bookkeeping inside the Hessian callback.
# =============================================================================

function _lagrangian(nlp_packed::HeatExchange.WeightedMeanNLPPacked, x, residual_buffer, p)
    _heat_balance_residuals!(nlp_packed, residual_buffer, x, p.ipopt_parameters.nlp_parameters)
    return p.sigma * _objective_value(nlp_packed, x, p.ipopt_parameters.objective_parameters) +
           p.lambda[1]*residual_buffer[1] + p.lambda[2]*residual_buffer[2] +
           p.lambda[3]*residual_buffer[3] + p.lambda[4]*residual_buffer[4]
end
function _lagrangian(nlp_packed::HeatExchange.MultiSidedNLPPacked, x, residual_buffer, p)
    _heat_balance_residuals!(nlp_packed, residual_buffer, x, p.ipopt_parameters.nlp_parameters)
    return p.sigma * _objective_value(nlp_packed, x, p.ipopt_parameters.objective_parameters) +
           p.lambda[1]*residual_buffer[1] + p.lambda[2]*residual_buffer[2] + p.lambda[3]*residual_buffer[3] +
           p.lambda[4]*residual_buffer[4] + p.lambda[5]*residual_buffer[5] + p.lambda[6]*residual_buffer[6]
end

# =============================================================================
# Per-call inputs (bounds, initial_values, objective_parameters, nlp_parameters)
#
# These helpers populate the mutable `IPOPTParameters` container and the
# bound/initial_values vectors that the cache reuses. They are also called from
# the fresh-per-call path with disposable buffers, so the bounds/initial_values
# computation lives in exactly one place. Each writes Float64 fields only — the
# NamedTuple field types are stable across calls, so assignment never reboxes.
# =============================================================================

function _inputs!(lower_bounds, upper_bounds, initial_values, nlp_packed,
                   organism, environment, init::NamedTuple)
    limits                         = thermoregulation(organism)
    metabolism_parameters          = metabolism_pars(organism)
    internal_conduction_parameters = conduction_pars_internal(organism)
    evaporation_parameters         = evaporation_pars(organism)

    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metabolism_parameters.core_temperature)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - 5.0
    skin_temperature_max   = core_temperature_max + 5.0
    heat_flow_min_W        = ustrip(u"W", limits.minimum_heat_flow)
    heat_flow_max_W        = heat_flow_min_W * 20.0
    flesh_conductivity_min = ustrip(u"W/m/K", limits.flesh_conductivity.reference)
    flesh_conductivity_max = ustrip(u"W/m/K", limits.flesh_conductivity.max)
    panting_rate_max       = Float64(limits.panting.pant.max)
    skin_wetness_min       = Float64(limits.skin_wetness.reference)
    skin_wetness_max       = Float64(limits.skin_wetness.max)
    insulation_depth_min   = ustrip(u"m", limits.insulation.dorsal.reference)
    insulation_depth_max   = ustrip(u"m", limits.insulation.dorsal.max)
    aspect_ratio_min       = Float64(limits.aspect_ratio_factor.reference)
    aspect_ratio_max       = organism.body.shape isa Sphere ? aspect_ratio_min : Float64(limits.aspect_ratio_factor.max)
    heat_flow_init         = max(ustrip(u"W", init.metabolic_heat_flow), heat_flow_min_W)
    flesh_conductivity     = ustrip(u"W/m/K", internal_conduction_parameters.flesh_conductivity)
    skin_wetness           = Float64(evaporation_parameters.skin_wetness)

    scalars = (;
        core_temperature_min, core_temperature_max,
        skin_temperature_min, skin_temperature_max,
        heat_flow_min_W, heat_flow_max_W,
        flesh_conductivity_min, flesh_conductivity_max,
        panting_rate_max,
        skin_wetness_min, skin_wetness_max,
        insulation_depth_min, insulation_depth_max,
        aspect_ratio_min, aspect_ratio_max,
        setpoint_temperature_K, heat_flow_init,
        flesh_conductivity, skin_wetness,
    )

    _write_layout!(lower_bounds, upper_bounds, initial_values, nlp_packed, init, scalars)

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
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
        panting_rate_range       = max(panting_rate_max - 1.0, 1e-6),
        skin_wetness_min         = Float64(skin_wetness_min),
        skin_wetness_range       = max(skin_wetness_max - skin_wetness_min, 1e-6),
        heat_flow_range          = max(heat_flow_max_W - heat_flow_min_W, 1.0),
    )
    nlp_parameters = (;
        nlp_packed,
        minimum_heat_flow    = Float64(heat_flow_min_W),
        q10                  = Float64(metabolism_parameters.q10),
        setpoint_temperature = Float64(setpoint_temperature_K),
    )
    return objective_parameters, nlp_parameters
end

# Strategy-specific variable-index layout. The two formulations differ only in
# which decision variables exist and at what positions; everything else lives
# in `_inputs!`. `s` is the NamedTuple of pre-computed scalars assembled there.
function _write_layout!(lower_bounds, upper_bounds, initial_values,
                         ::HeatExchange.WeightedMeanNLPPacked, init, s)
    lower_bounds[1] = s.core_temperature_min;   lower_bounds[2] = s.skin_temperature_min;   lower_bounds[3] = s.skin_temperature_min
    lower_bounds[4] = log(s.heat_flow_min_W);   lower_bounds[5] = s.flesh_conductivity_min; lower_bounds[6] = 1.0
    lower_bounds[7] = s.skin_wetness_min;       lower_bounds[8] = s.insulation_depth_min;   lower_bounds[9] = s.aspect_ratio_min
    upper_bounds[1] = s.core_temperature_max;   upper_bounds[2] = s.skin_temperature_max;   upper_bounds[3] = s.skin_temperature_max
    upper_bounds[4] = log(s.heat_flow_max_W);   upper_bounds[5] = s.flesh_conductivity_max; upper_bounds[6] = s.panting_rate_max
    upper_bounds[7] = s.skin_wetness_max;       upper_bounds[8] = s.insulation_depth_max;   upper_bounds[9] = s.aspect_ratio_max

    initial_values[1] = clamp(s.setpoint_temperature_K,                    lower_bounds[1], upper_bounds[1])
    initial_values[2] = clamp(ustrip(u"K", init.skin_temperature),         lower_bounds[2], upper_bounds[2])
    initial_values[3] = clamp(ustrip(u"K", init.insulation_temperature),   lower_bounds[3], upper_bounds[3])
    initial_values[4] = clamp(log(s.heat_flow_init),                       lower_bounds[4], upper_bounds[4])
    initial_values[5] = clamp(s.flesh_conductivity,                        lower_bounds[5], upper_bounds[5])
    initial_values[6] = clamp(1.0,                                         lower_bounds[6], upper_bounds[6])
    initial_values[7] = clamp(s.skin_wetness,                              lower_bounds[7], upper_bounds[7])
    initial_values[8] = clamp(s.insulation_depth_max,                      lower_bounds[8], upper_bounds[8])  # start with erected insulation
    initial_values[9] = clamp(s.aspect_ratio_min,                          lower_bounds[9], upper_bounds[9])  # start curled
    return nothing
end
function _write_layout!(lower_bounds, upper_bounds, initial_values,
                         nlp_packed::HeatExchange.MultiSidedNLPPacked, init, s)
    packed = nlp_packed.p

    lower_bounds[1]  = s.core_temperature_min
    lower_bounds[2]  = s.skin_temperature_min;  lower_bounds[3]  = s.skin_temperature_min    # dorsal skin, ins
    lower_bounds[4]  = s.skin_temperature_min;  lower_bounds[5]  = s.skin_temperature_min    # ventral skin, ins
    lower_bounds[6]  = log(s.heat_flow_min_W);  lower_bounds[7]  = s.flesh_conductivity_min
    lower_bounds[8]  = 1.0;                     lower_bounds[9]  = s.skin_wetness_min
    lower_bounds[10] = s.insulation_depth_min;  lower_bounds[11] = s.aspect_ratio_min
    upper_bounds[1]  = s.core_temperature_max
    upper_bounds[2]  = s.skin_temperature_max;  upper_bounds[3]  = s.skin_temperature_max
    upper_bounds[4]  = s.skin_temperature_max;  upper_bounds[5]  = s.skin_temperature_max
    upper_bounds[6]  = log(s.heat_flow_max_W);  upper_bounds[7]  = s.flesh_conductivity_max
    upper_bounds[8]  = s.panting_rate_max;      upper_bounds[9]  = s.skin_wetness_max
    upper_bounds[10] = s.insulation_depth_max;  upper_bounds[11] = s.aspect_ratio_max

    initial_values[1]  = clamp(s.setpoint_temperature_K,                                       lower_bounds[1],  upper_bounds[1])
    initial_values[2]  = clamp(ustrip(u"K", packed.initial_dorsal_skin_temperature),           lower_bounds[2],  upper_bounds[2])
    initial_values[3]  = clamp(ustrip(u"K", packed.initial_dorsal_insulation_temperature),     lower_bounds[3],  upper_bounds[3])
    initial_values[4]  = clamp(ustrip(u"K", packed.initial_ventral_skin_temperature),          lower_bounds[4],  upper_bounds[4])
    initial_values[5]  = clamp(ustrip(u"K", packed.initial_ventral_insulation_temperature),    lower_bounds[5],  upper_bounds[5])
    initial_values[6]  = clamp(log(s.heat_flow_init),                                          lower_bounds[6],  upper_bounds[6])
    initial_values[7]  = clamp(s.flesh_conductivity,                                           lower_bounds[7],  upper_bounds[7])
    initial_values[8]  = clamp(1.0,                                                            lower_bounds[8],  upper_bounds[8])
    initial_values[9]  = clamp(s.skin_wetness,                                                 lower_bounds[9],  upper_bounds[9])
    initial_values[10] = clamp(s.insulation_depth_max,                                         lower_bounds[10], upper_bounds[10])
    initial_values[11] = clamp(s.aspect_ratio_min,                                             lower_bounds[11], upper_bounds[11])
    return nothing
end

# =============================================================================
# Warm-start cache
#
# `IPOPTParameters` is a mutable struct that carries the two NamedTuples the
# Enzyme-AD callbacks consume (objective penalties + NLP physics parameters).
# The Ipopt callbacks hold a reference to this container and read its fields on
# every invocation. To set up a new solve we just write fresh NamedTuples into
# the two fields — same concrete types, no reboxing — and the callbacks pick
# them up automatically without rebuilding.
#
# Both fields are parametrically typed so accesses stay type-stable: as long
# as the per-call `objective_parameters` and `nlp_parameters` NamedTuples have the same key/type
# signature (guaranteed by the `_inputs!` helpers above), assignment is a
# bit-for-bit field overwrite. If a structurally different organism or
# environment is fed in (e.g. a different body shape), construct a new cache.
# =============================================================================

mutable struct IPOPTParameters{ObjectiveParametersType,NlpParametersType}
    objective_parameters::ObjectiveParametersType
    nlp_parameters::NlpParametersType
end

# Pre-allocated workspaces consumed by the Enzyme-AD calls inside the Ipopt
# callbacks. One instance per cache, sized at construction from (n_variables, n_constraints).
#
# Jacobian section (used by `evaluate_constraint_jacobian`):
#   jacobian_residual           length n_constraints, residuals primal
#   jacobian_residual_seed          length n_constraints, reverse-mode seed (unit vector for row j)
#   jacobian_x_derivative          length n_variables, gradient output → one Jacobian row
#
# Hessian section (used by `evaluate_hessian`, forward-over-reverse):
#   hessian_gradient         length n_variables, inner-Reverse gradient
#   hessian_gradient_tangent         length n_variables, outer-Forward dual → one Hessian column
#   hessian_seed       length n_variables, outer-Forward seed (unit vector for column j)
#   hessian_multipliers     length n_constraints, current IPOPT λ_g (copied per `evaluate_hessian` call)
#   hessian_residual          length n_constraints, Lagrangian residual buffer
#   hessian_residual_adjoint         length n_constraints, inner-Reverse adjoint of `hessian_residual`
#   hessian_residual_tangent      length n_constraints, outer-Forward tangent of `hessian_residual`
#   hessian_residual_adjoint_tangent     length n_constraints, outer-Forward tangent of `hessian_residual_adjoint`
struct IpoptCallbackBuffers
    jacobian_residual::Vector{Float64}
    jacobian_residual_seed::Vector{Float64}
    jacobian_x_derivative::Vector{Float64}
    hessian_gradient::Vector{Float64}
    hessian_gradient_tangent::Vector{Float64}
    hessian_seed::Vector{Float64}
    hessian_multipliers::Vector{Float64}
    hessian_residual::Vector{Float64}
    hessian_residual_adjoint::Vector{Float64}
    hessian_residual_tangent::Vector{Float64}
    hessian_residual_adjoint_tangent::Vector{Float64}
end

function IpoptCallbackBuffers(n_variables::Int, n_constraints::Int)
    IpoptCallbackBuffers(
        zeros(n_constraints), zeros(n_constraints), zeros(n_variables),               # jacobian buffers
        zeros(n_variables), zeros(n_variables), zeros(n_variables), zeros(n_constraints),     # hessian gradient / tangent / seed / multipliers
        zeros(n_constraints), zeros(n_constraints), zeros(n_constraints), zeros(n_constraints),     # hessian residual + adjoint + their forward tangents
    )
end

"""
    IPOPTSolverCache

    IPOPTSolverCache(control, organism, environment, init)

Reusable solver cache for the IPOPT thermoregulation path. Holds the
Ipopt callback functors, the per-call bound and `initial_values` buffers, the
reverse-mode Jacobian scratch vectors, the Julia-side primal+dual warm-start
vectors, and a mutable `IPOPTParameters` whose fields the functors read on
every solve.

Construct once per organism + smoothing + NLP-strategy combination:

    cache = IPOPTSolverCache(control, organism, environment, init)

where `init` is a NamedTuple of initial guesses (`metabolic_heat_flow`,
`skin_temperature`, `insulation_temperature`).

Then pass the cache as the trailing argument of `thermoregulate` to skip
the closure / scratch-buffer rebuilding and to enable primal+dual warm-
starting from the previous solve. The cache must not be shared between
threads.

Invalidate the cache (build a new one) if the organism's body shape, NLP
strategy, or smoothing strategy changes — the closures bake in those types.
Routine per-iteration changes (environment temperatures, q10, etc.) reuse
the cache without rebuilding.

## Warm-start scope

The cache warm-starts **both** primal and dual:

- *Primal*: the previous solve's effector vector (clamped into the new
  per-call bounds) is fed in as the next solve's initial `x`. The fresh
  path seeds only the three `init` fields from user-passed values; the
  cache seeds all 9 (11 for MultiSided), a strictly better starting point
  in smooth regimes.
- *Dual*: the Lagrange multipliers from the previous solve (`mult_g` on
  constraints, `mult_x_L` / `mult_x_U` on bounds) are restored into the
  next `IpoptProblem` before `IpoptSolve`, and Ipopt is told to use them
  via `warm_start_init_point=yes`. Bound multipliers for inactive bounds
  carry over harmlessly; for active bounds they accelerate convergence.

We go directly through `Ipopt.jl` rather than `OptimizationIpopt` because
the latter creates a fresh `IpoptProblem` on every `solve` and discards
the final multipliers — no Julia-side hook for them. Bounds still need to
be set at `CreateIpoptProblem` time, so we rebuild the C-side problem each
call with the tight per-call bounds; that allocation is small and Julia-
side state (`previous_primal`, `prev_mult_*`) carries the warm-start info across
rebuilds.
"""
struct IPOPTSolverCache{ParametersType<:IPOPTParameters, ObjectiveCallback, ConstraintsCallback, ObjectiveGradientCallback, ConstraintJacobianCallback, HessianCallback}
    ipopt_parameters::ParametersType
    # Cached callback functors. These read `ipopt_parameters.objective_parameters`
    # and `ipopt_parameters.nlp_parameters` at solve time, so updating those
    # NamedTuples retargets every callback without rebuilding the functors.
    # Stored on the cache so the same callable objects are passed to every
    # CreateIpoptProblem call — that way Julia compiles a single specialised
    # IpoptProblem type once and reuses it.
    evaluate_objective::ObjectiveCallback
    evaluate_constraints::ConstraintsCallback
    evaluate_objective_gradient::ObjectiveGradientCallback
    evaluate_constraint_jacobian::ConstraintJacobianCallback
    evaluate_hessian::HessianCallback
    # Per-call tight bounds — written by `_inputs!` each solve to match the
    # fresh `thermoregulate(..., organism, env, init)` path. Ipopt bakes
    # bounds into the C-side problem at CreateIpoptProblem time, so we rebuild
    # the C handle every solve (cheap — a single malloc + struct init) rather
    # than loosening the bounds to keep one handle alive.
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    initial_values::Vector{Float64}
    constraint_lower_bounds::Vector{Float64}
    constraint_upper_bounds::Vector{Float64}
    n_variables::Int                              # variables (9 weighted / 11 multisided)
    n_constraints::Int                              # constraints (4 weighted / 6 multisided)
    # All AD scratch buffers consumed by the Ipopt callbacks. See the
    # `IpoptCallbackBuffers` definition above.
    buffers::IpoptCallbackBuffers
    # Warm-start state, persisted Julia-side across IpoptProblem rebuilds.
    # On the first solve `has_previous_solve = false` so we feed `initial_values` as the primal and
    # let Ipopt initialise its own multipliers. After IpoptSolve we copy
    # `prob.x`, `prob.mult_g`, `prob.mult_x_L`, `prob.mult_x_U` into these
    # vectors. On subsequent solves we restore them before calling IpoptSolve
    # with `warm_start_init_point = "yes"`.
    previous_primal::Vector{Float64}
    previous_constraint_multipliers::Vector{Float64}
    previous_lower_bound_multipliers::Vector{Float64}
    previous_upper_bound_multipliers::Vector{Float64}
    has_previous_solve::Base.RefValue{Bool}
    # User-provided variable scaling (one entry per decision variable). See
    # the `_scaling` methods for the values. Passed to
    # Ipopt via `SetIpoptProblemScaling` together with `obj_scaling = 1.0`
    # and an all-ones `g_scaling`, after `CreateIpoptProblem` but before
    # `IpoptSolve`. Allocated once at cache construction and reused for
    # every solve.
    x_scaling::Vector{Float64}
    g_scaling::Vector{Float64}
end
function IPOPTSolverCache(control::IPOPTControl, organism::Organism, environment::NamedTuple, init::NamedTuple)
    nlp_packed = HeatExchange.nlp_pack(control.nlp_strategy, organism, environment,
                                       init.skin_temperature, init.insulation_temperature;
                                       smoothing = control.smoothing)
    return IPOPTSolverCache(nlp_packed, organism, environment, init)
end
function IPOPTSolverCache(nlp_packed, organism::Organism, environment::NamedTuple, init::NamedTuple)
    n_variables, n_constraints = _problem_size(nlp_packed)
    lower_bounds   = zeros(n_variables)
    upper_bounds   = zeros(n_variables)
    initial_values = zeros(n_variables)
    objective_parameters, nlp_parameters = _inputs!(lower_bounds, upper_bounds, initial_values,
        nlp_packed, organism, environment, init)
    ipopt_parameters = IPOPTParameters(objective_parameters, nlp_parameters)
    buffers = IpoptCallbackBuffers(n_variables, n_constraints)
    evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_constraint_jacobian, evaluate_hessian =
        _build_ipopt_callbacks(ipopt_parameters, n_variables, n_constraints, buffers)
    constraint_lower_bounds, constraint_upper_bounds = _constraint_bounds(n_constraints)
    return IPOPTSolverCache(ipopt_parameters,
        evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_constraint_jacobian, evaluate_hessian,
        lower_bounds, upper_bounds, initial_values,
        constraint_lower_bounds, constraint_upper_bounds,
        n_variables, n_constraints, buffers,
        zeros(n_variables), zeros(n_constraints), zeros(n_variables), zeros(n_variables), Ref(false),
        _scaling(nlp_packed), ones(n_constraints))
end

# Direct-Ipopt callbacks. Strategy-specific behaviour comes from method
# dispatch on `ipopt_parameters.nlp_parameters.nlp_packed` — passed as `Const`
# first arg to every Enzyme.autodiff call so Julia picks the correct
# `_objective_value`, `_heat_balance_residuals!`, `_lagrangian` method per
# `WeightedMeanNLPPacked` / `MultiSidedNLPPacked`. The rest of the wiring
# (Enzyme plumbing, sparsity, lower-triangle indexing) is shared.
#
# Functors read `ipopt_parameters.objective_parameters` /
# `ipopt_parameters.nlp_parameters` at solve time, so refreshing the parameters
# container retargets every callback without rebuilding them.
#
# Jacobian is treated as fully dense — for n_variables=9/11 sparse bookkeeping costs
# more than the saved entries. Sparsity convention: row-major dense, 1-based,
# entry k = (i-1)*n_variables + j ↔ (rows[k], cols[k]) = (i, j).
#
# Hessian uses forward-over-reverse Enzyme: the inner Reverse functor
# `LagrangianGradient` computes ∇ₓL into `buffers.hessian_gradient`, the outer
# Forward seeds e_j and pulls back column j of H into
# `buffers.hessian_gradient_tangent`. Sparsity is dense lower triangle
# (i ≥ j), entry k = i(i-1)/2 + j ↔ (i, j).
# --- Functor types for the five Ipopt callbacks ------------------------------
#
# Each callback is a `struct` with a call method. The struct fields stand in
# for what closures would capture (`ipopt_parameters`, `buffers`, problem dimensions).
# Two practical benefits over anonymous closures:
#   - the callable's type is named — easier to read in stack traces / `@code_*`.
#   - the captured fields are first-class struct fields with explicit types,
#     so type-stability is obvious by inspection rather than implied by capture.
# `LagrangianGradient` is the inner reverse-mode pass that `EvaluateHessian`
# differentiates forward over; it's its own functor for the same reasons.

struct EvaluateObjective{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
(f::EvaluateObjective)(x) =
    _objective_value(f.ipopt_parameters.nlp_parameters.nlp_packed, x, f.ipopt_parameters.objective_parameters)

struct EvaluateConstraints{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
(f::EvaluateConstraints)(x, g) =
    _heat_balance_residuals!(f.ipopt_parameters.nlp_parameters.nlp_packed, g, x, f.ipopt_parameters.nlp_parameters)

struct EvaluateObjectiveGradient{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
function (f::EvaluateObjectiveGradient)(x, grad_f)
    fill!(grad_f, 0)
    Enzyme.autodiff(Enzyme.Reverse, _objective_value,
                    Enzyme.Active,
                    Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
                    Enzyme.Duplicated(x, grad_f),
                    Enzyme.Const(f.ipopt_parameters.objective_parameters))
    return nothing
end

struct EvaluateConstraintJacobian{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
    buffers::IpoptCallbackBuffers
    n_variables::Int
    n_constraints::Int
end
function (f::EvaluateConstraintJacobian)(x, rows, cols, values)
    if values === nothing
        k = 1
        @inbounds for i in 1:f.n_constraints, j in 1:f.n_variables
            rows[k] = i; cols[k] = j; k += 1
        end
    else
        b = f.buffers
        @inbounds for i in 1:f.n_constraints
            fill!(b.jacobian_residual, 0); fill!(b.jacobian_residual_seed, 0); fill!(b.jacobian_x_derivative, 0); b.jacobian_residual_seed[i] = 1.0
            Enzyme.autodiff(Enzyme.Reverse, _heat_balance_residuals!,
                            Enzyme.Const,
                            Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
                            Enzyme.Duplicated(b.jacobian_residual, b.jacobian_residual_seed),
                            Enzyme.Duplicated(x, b.jacobian_x_derivative),
                            Enzyme.Const(f.ipopt_parameters.nlp_parameters))
            for j in 1:f.n_variables
                values[(i-1)*f.n_variables + j] = b.jacobian_x_derivative[j]
            end
        end
    end
    return nothing
end

# Inner reverse-mode pass: Enzyme.Reverse on `_lagrangian`. `x` and the
# residual buffer are both Duplicated so Enzyme writes adjoints into the
# supplied buffers rather than allocating shadows. Dispatched to the right
# `_lagrangian` method via the concrete packed type held on `ipopt_parameters`.
struct LagrangianGradient{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
function (f::LagrangianGradient)(g, x, residual_buffer, residual_adjoint, p)
    fill!(g, 0)
    fill!(residual_adjoint, 0)
    Enzyme.autodiff(Enzyme.Reverse, _lagrangian,
                    Enzyme.Active,
                    Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
                    Enzyme.Duplicated(x, g),
                    Enzyme.Duplicated(residual_buffer, residual_adjoint),
                    Enzyme.Const(p))
    return nothing
end

struct EvaluateHessian{P<:IPOPTParameters, LG<:LagrangianGradient} <: Function
    ipopt_parameters::P
    buffers::IpoptCallbackBuffers
    n_variables::Int
    lagrangian_gradient::LG
end
function (f::EvaluateHessian)(x, rows, cols, sigma, lambda, values)
    if values === nothing
        k = 1
        @inbounds for i in 1:f.n_variables, j in 1:i
            rows[k] = i; cols[k] = j; k += 1
        end
    else
        b = f.buffers
        # IPOPT passes the current objective factor (sigma) and constraint
        # multipliers; copy lambda into the cache so the Lagrangian functor
        # (via `p.sigma` / `p.lambda`) sees them under nested AD.
        copyto!(b.hessian_multipliers, lambda)
        p = (; sigma, lambda = b.hessian_multipliers, ipopt_parameters = f.ipopt_parameters)
        @inbounds for j in 1:f.n_variables
            fill!(b.hessian_seed, 0); b.hessian_seed[j] = 1.0
            fill!(b.hessian_gradient,    0); fill!(b.hessian_gradient_tangent,     0)
            fill!(b.hessian_residual,    0); fill!(b.hessian_residual_tangent,  0)
            fill!(b.hessian_residual_adjoint,   0); fill!(b.hessian_residual_adjoint_tangent, 0)
            # Outer Forward propagates tangents through everything
            # `lagrangian_gradient` reads/writes — including `residual_buffer`
            # and `residual_adjoint` — so every mutable buffer needs its
            # forward dual here.
            Enzyme.autodiff(Enzyme.Forward, f.lagrangian_gradient,
                            Enzyme.Const,
                            Enzyme.Duplicated(b.hessian_gradient, b.hessian_gradient_tangent),
                            Enzyme.Duplicated(x, b.hessian_seed),
                            Enzyme.Duplicated(b.hessian_residual,  b.hessian_residual_tangent),
                            Enzyme.Duplicated(b.hessian_residual_adjoint, b.hessian_residual_adjoint_tangent),
                            Enzyme.Const(p))
            for i in j:f.n_variables
                k = (i * (i - 1)) ÷ 2 + j
                values[k] = b.hessian_gradient_tangent[i]
            end
        end
    end
    return nothing
end

function _build_ipopt_callbacks(ipopt_parameters::IPOPTParameters, n_variables::Int, n_constraints::Int,
                                 buffers::IpoptCallbackBuffers)
    return (
        EvaluateObjective(ipopt_parameters),
        EvaluateConstraints(ipopt_parameters),
        EvaluateObjectiveGradient(ipopt_parameters),
        EvaluateConstraintJacobian(ipopt_parameters, buffers, n_variables, n_constraints),
        EvaluateHessian(ipopt_parameters, buffers, n_variables, LagrangianGradient(ipopt_parameters)),
    )
end

# Apply the project-wide Ipopt options to a fresh problem. Called once per
# `CreateIpoptProblem` — i.e. per cached solve — alongside the warm-start
# toggle. `print_level` is set per call because callers can flip `verbose`
# between solves; everything else is constant for the cache's lifetime.
function _apply_ipopt_options!(prob::Ipopt.IpoptProblem, verbose::Bool, warm_start::Bool)
    # `hessian_approximation` defaults to "exact" — we register a real
    # `evaluate_hessian` (Lagrangian Hessian via Enzyme forward-over-reverse), so no
    # override needed.
    Ipopt.AddIpoptNumOption(prob, "acceptable_tol",        1e-3)
    Ipopt.AddIpoptIntOption(prob, "acceptable_iter",       5)
    Ipopt.AddIpoptNumOption(prob, "tol",                   1e-4)
    Ipopt.AddIpoptIntOption(prob, "max_iter",              300)
    Ipopt.AddIpoptIntOption(prob, "print_level",           verbose ? 5 : 0)
    Ipopt.AddIpoptStrOption(prob, "warm_start_init_point", warm_start ? "yes" : "no")
    # Adaptive barrier-parameter strategy — Ipopt picks `mu` per-iteration
    # based on the optimality-error / complementarity ratio instead of the
    # default `monotone` schedule. On smooth NLPs (ours: ellipsoid heat
    # balance with continuous physical relationships) this typically
    # converges in fewer iterations because the barrier collapses faster
    # once the optimiser is in the right basin.
    Ipopt.AddIpoptStrOption(prob, "mu_strategy",           "adaptive")
    # User-provided variable scaling — set the actual factors via
    # `SetIpoptProblemScaling` after CreateIpoptProblem (see the `_scaling`
    # methods for values). Telling Ipopt
    # we'll supply our own scales skips its `gradient-based` heuristic pass
    # (an extra Jacobian eval at the initial point) and gives the solver a
    # well-conditioned internal problem from iteration zero.
    Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method",    "user-scaling")
    return prob
end

# ──────────────────────────────────────────────────────────────────────────
# TEMPORARY: variable scaling by raw index position.
#
# These methods hard-code Ipopt's per-variable scale factors keyed by integer
# index. The variable layout (`x[1] = core_T`, …) only exists as comments in
# the residual function headers — there is no shared schema between (a)
# `_inputs!` (which writes bounds and `initial_values`), (buffers) `_objective_value` /
# `_heat_balance_residuals!` (which decode `x[i]`), and (c) `_scaling`.
# Reordering / inserting / renaming a variable in any one of those three
# places without touching the other two silently mis-scales the problem.
#
# Proper fix: extend the ModelParameters.jl `Param` wrappers that already
# carry `bounds` / `units` / `val` (see `HeatExchange/src/traits.jl`) with a
# `scaling` field. Each decision variable's source `Param` then owns its
# scale alongside its bounds, and `_inputs!` reads scale and bounds from
# the same place — no positional duplication. Requires adding `scaling` to
# the relevant traits (`ThermoregulationLimits`, etc.) and threading it
# through `_inputs!` to fill `cache.x_scaling`. Worth doing the next time
# we touch the trait schema.
# ──────────────────────────────────────────────────────────────────────────
function _scaling(::HeatExchange.WeightedMeanNLPPacked)
    x_scaling = zeros(9)
    x_scaling[1] = 1/300.0     # core_T ~ 300 K
    x_scaling[2] = 1/300.0     # skin_temperature ~ 300 K
    x_scaling[3] = 1/300.0     # insulation_T ~ 300 K
    x_scaling[4] = 1.0         # log(metabolic_heat_flow) ~ 0 (W ~ 1)
    x_scaling[5] = 1.0         # flesh_conductivity ~ 1 W/m/K
    x_scaling[6] = 1.0         # panting ∈ [1, 15] → internal [1, 15]
    x_scaling[7] = 50.0        # skin_wetness ∈ [0, 0.05] → internal [0, 2.5]
    x_scaling[8] = 50.0        # insulation_depth ∈ [0, 0.05 n_constraints] → internal [0, 2.5]
    x_scaling[9] = 1.0         # axis_ratio ~ 1
    return x_scaling
end
function _scaling(::HeatExchange.MultiSidedNLPPacked)
    x_scaling = zeros(11)
    x_scaling[1] = 1/300.0     # core
    for i in 2:5
        x_scaling[i] = 1/300.0 # dorsal/ventral skin + insulation T
    end
    x_scaling[6]  = 1.0        # log_metabolic_heat_flow
    x_scaling[7]  = 1.0        # flesh_conductivity
    x_scaling[8]  = 1.0        # panting
    x_scaling[9]  = 50.0       # skin_wetness
    x_scaling[10] = 50.0       # insulation_depth
    x_scaling[11] = 1.0        # axis_ratio
    return x_scaling
end

# Decision-variable / constraint counts per NLP strategy. The only
# strategy-dependent dimensions; everything downstream (buffer sizes,
# sparsity, warm-start vectors) is allocated from these.
_problem_size(::HeatExchange.WeightedMeanNLPPacked) = (9, 4)
_problem_size(::HeatExchange.MultiSidedNLPPacked)   = (11, 6)

# Constraint bounds: residuals[1:end-1] == 0 (equalities), residuals[end] >= 0
# (Q10 inequality). Layout is identical for both NLP strategies.
function _constraint_bounds(n_constraints::Int)
    constraint_lower_bounds = zeros(n_constraints)
    constraint_upper_bounds = zeros(n_constraints)
    constraint_upper_bounds[end] = Inf
    return constraint_lower_bounds, constraint_upper_bounds
end

"""
    reset_warm_start!(cache::IPOPTSolverCache)

Discard the previous solution stored in the cache. The next solve will use
the cold initial guess derived from the `init` NamedTuple. Use this when the
trajectory takes a large discontinuous step (e.g. switching organism,
changing q10 regime) where seeding from the last solution would land the
optimiser in a bad basin.
"""
reset_warm_start!(cache::IPOPTSolverCache) = (cache.has_previous_solve[] = false; cache)

# =============================================================================
# Per-call input refresh.
# Overwrite the cache's tight-bound buffers and the two NamedTuples on
# `ipopt_parameters` in place. Dispatch on the packed type happens inside
# `_inputs!`, so the unused branch is dead-code-eliminated and no Union return
# appears on the hot path.
# =============================================================================

function _refresh_cache!(cache::IPOPTSolverCache,
                   nlp_packed, organism, environment, init::NamedTuple)
    objective_parameters, nlp_parameters = _inputs!(cache.lower_bounds, cache.upper_bounds, cache.initial_values, nlp_packed,
        organism, environment, init)
    cache.ipopt_parameters.objective_parameters = objective_parameters
    cache.ipopt_parameters.nlp_parameters = nlp_parameters
    return cache
end
_refresh_cache!(cache::Nothing, args...) = IPOPTSolverCache(args...)

# Direct-Ipopt solve with primal+dual warm-start.
#
# Each call rebuilds the C-side IpoptProblem so the new tight bounds in
# `cache.lower_bounds` / `cache.upper_bounds` take effect — same bounds the fresh
# `thermoregulate(..., organism, env, init)` path uses, so solutions
# stay bit-identical to that path's first iterate. CreateIpoptProblem is a
# malloc + a few struct inits; well under a percent of an IPOPT iteration's
# cost.
#
# Returns the IpoptProblem so the caller can pull the final `x` (and
# optionally inspect the multipliers via `cache.prev_*`).
function _ipopt_solve!(cache::IPOPTSolverCache, verbose::Bool)
    prob = Ipopt.CreateIpoptProblem(
        cache.n_variables, cache.lower_bounds, cache.upper_bounds,
        cache.n_constraints, cache.constraint_lower_bounds, cache.constraint_upper_bounds,
        cache.n_constraints * cache.n_variables,                       # nele_jac: dense
        (cache.n_variables * (cache.n_variables + 1)) ÷ 2,           # nele_hess: dense lower triangle
        cache.evaluate_objective, cache.evaluate_constraints, cache.evaluate_objective_gradient, cache.evaluate_constraint_jacobian,
        cache.evaluate_hessian,
    )
    _apply_ipopt_options!(prob, verbose, cache.has_previous_solve[])
    # User-scaling pairs with `nlp_scaling_method = "user-scaling"` in the
    # options above. `obj_scaling = 1.0` is fine — the objective is already
    # a sum of dimensionless squared penalties so it's intrinsically O(1).
    Ipopt.SetIpoptProblemScaling(prob, 1.0, cache.x_scaling, cache.g_scaling)

    if cache.has_previous_solve[]
        # Restore the previous primal+dual, with the primal clamped into the
        # new tight bounds. `skin_T_min = air_T − 5` shifts between solves;
        # an out-of-range initial x makes Ipopt spend its first iteration
        # snapping back to the feasible box. Dual multipliers can carry over
        # as-is: inactive bounds have ≈0 multipliers (no information lost),
        # active bounds get the previous shadow value as a good guess.
        @inbounds for i in 1:cache.n_variables
            prob.x[i] = clamp(cache.previous_primal[i], cache.lower_bounds[i], cache.upper_bounds[i])
        end
        copyto!(prob.mult_g,   cache.previous_constraint_multipliers)
        copyto!(prob.mult_x_L, cache.previous_lower_bound_multipliers)
        copyto!(prob.mult_x_U, cache.previous_upper_bound_multipliers)
    else
        # Cold start — seed primal from the inputs-helper defaults, let Ipopt
        # initialise its own multipliers (warm_start_init_point = "no" above).
        @inbounds for i in 1:cache.n_variables
            prob.x[i] = cache.initial_values[i]
        end
    end

    Ipopt.IpoptSolve(prob)

    # Persist primal+dual for the next solve.
    copyto!(cache.previous_primal,       prob.x)
    copyto!(cache.previous_constraint_multipliers,  prob.mult_g)
    copyto!(cache.previous_lower_bound_multipliers, prob.mult_x_L)
    copyto!(cache.previous_upper_bound_multipliers, prob.mult_x_U)
    cache.has_previous_solve[] = true
    return prob
end

# Unit reattachment + delegation to HeatExchange's output builder. Kept tiny
# and shape-agnostic so both fresh and cached paths share the assembly step.
# Dispatched on the packed type (which is 1:1 with the NLP strategy) so the
# variable layout per strategy stays adjacent to its decoding.
function _assemble(nlp_packed::HeatExchange.WeightedMeanNLPPacked, organism, environment, x_sol)
    return HeatExchange.nlp_assemble_output(nlp_packed, organism, environment,
        x_sol[1] * u"K", x_sol[2] * u"K", x_sol[3] * u"K", exp(x_sol[4]) * u"W",
        x_sol[5] * u"W/m/K", x_sol[6], x_sol[7], x_sol[8] * u"m", x_sol[9])
end
function _assemble(nlp_packed::HeatExchange.MultiSidedNLPPacked, organism, environment, x_sol)
    return HeatExchange.nlp_assemble_output(nlp_packed, organism, environment,
        x_sol[1]  * u"K",
        x_sol[2]  * u"K", x_sol[3]  * u"K",
        x_sol[4]  * u"K", x_sol[5]  * u"K",
        exp(x_sol[6]) * u"W",
        x_sol[7]  * u"W/m/K", x_sol[8], x_sol[9],
        x_sol[10] * u"m", x_sol[11])
end

# =============================================================================
# Main thermoregulate dispatch
# =============================================================================

"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, init;
                   cache=nothing, verbose=false)

Solve endotherm heat balance as a nonlinear program via IPOPT.

`init` is a NamedTuple of initial guesses with fields `metabolic_heat_flow`,
`skin_temperature`, and `insulation_temperature`.

The NLP formulation is selected by `IPOPTControl.nlp_strategy`:
- `WeightedMeanNLP()` (default): dorsal/ventral weighted-mean single-body formulation.
  Nine decision variables, four constraints (three equality + Q10 inequality).
- `MultiSidedNLP()`: explicit per-side heat balance. Eleven decision variables, six
  constraints (five equality + Q10 inequality). Per-side skin and insulation temperatures
  are independent decision variables, giving the solver more freedom to represent
  asymmetric environmental loading (sun above, ground below).

Both formulations penalise deviation from the setpoint core temperature and enforce
`metabolic_heat_flow >= heat_flow_min * Q10^((core_temperature − setpoint)/10)`.
Penalty weights in `ThermoregulationLimits`:
  - `core_temperature_penalty`  — lower → core_temperature rises sooner before effectors are exhausted
  - `metabolic_heat_penalty`   — small regularisation (default 0.1) preventing high-panting/high-metabolic_heat_flow
                                  degeneracy at cold temperatures; overridden at hot temperatures by Q10 constraint
  - `gradient_penalty`         — non-zero → penalise deviation from `target_core_skin_gradient`
  - `panting_penalty`          — lower → panting activates sooner
  - `skin_wetness_penalty`     — higher than `panting_penalty` → panting activates before sweating

# Arguments
- `organism`: must have `OrganismTraits` with `ThermoregulationLimits`
- `environment`: NamedTuple with `environment_pars` and `environment_vars`
- `init`: NamedTuple of initial guesses — `metabolic_heat_flow` (W),
  `skin_temperature` (K), `insulation_temperature` (K). Use
  `initial_physiological_state(organism, environment.environment_vars)` for defaults.
- `cache` (keyword, optional): an `IPOPTSolverCache` produced by
  `IPOPTSolverCache(control, …)`. Reuse across calls in a sweep to skip per-solve
  allocation and enable primal+dual warm-starting. Must match the organism's body-
  shape type, NLP strategy and smoothing strategy — invalidate the cache and rebuild
  if any of those change.
"""
function thermoregulate(
    ::Endotherm,
    control::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    init::NamedTuple;
    cache::Union{Nothing,IPOPTSolverCache} = nothing,
    verbose::Bool = false,
)
    nlp_packed = HeatExchange.nlp_pack(
        control.nlp_strategy, organism, environment, init.skin_temperature, init.insulation_temperature;
        smoothing = control.smoothing
    )
    # Either branch leaves the cache fully populated for `_ipopt_solve!`, 
    # with `_inputs!` called exactly once.
    cache = _refresh_cache!(cache, nlp_packed, organism, environment, init)
    prob = _ipopt_solve!(cache, verbose)
    return _assemble(nlp_packed, organism, environment, prob.x)
end
