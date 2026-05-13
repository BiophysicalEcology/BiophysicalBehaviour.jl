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
#   1. `thermoregulate(::Endotherm, ::IPOPTControl, organism, env, M, sT, iT)`
#      — fresh per-call: allocates a single-use cache, runs one cold solve,
#      throws everything away. Use when each call has a structurally
#      different problem (different shape, different env_vars layout).
#   2. `IPOPTSolverCache(control, organism, env, M, sT, iT)` then
#      `thermoregulate(..., cache)` — warm-started: callback closures and
#      scratch buffers are reused, and primal+dual from the previous solve
#      are restored into the next `IpoptProblem` before `IpoptSolve` runs
#      with `warm_start_init_point = "yes"`. Use across a sweep where only
#      environment values change.
# =============================================================================

# =============================================================================
# Hard boundary between unitless optimizer world and Unitful scientific world.
#
#   `_objective_value_*` lives in optimizer world: pure Float64 arithmetic over
#       stripped reference values and dimensionless penalty weights.
#       Never sees a Unitful quantity.
#
#   `_heat_balance_residuals_*!` is the *only* function that crosses the
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

function _objective_value_weighted(effectors, opt)
    opt.core_temperature_penalty * ((effectors[1] - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
    opt.metabolic_heat_penalty   * ((exp(effectors[4]) - opt.heat_flow_min_W)   / opt.heat_flow_range)^2        +
    opt.gradient_penalty         * ((effectors[1] - effectors[2] - opt.target_gradient) / opt.gradient_range)^2 +
    opt.panting_penalty          * ((effectors[6] - 1.0) / opt.panting_rate_range)^2                            +
    opt.skin_wetness_penalty     * ((effectors[7] - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# Float64 effectors in → Unitful physics (via HeatExchange.nlp_residuals) → Float64 residuals out.
# residuals[1:3] == 0 (energy balance, internal conduction, skin temp).
# residuals[4]  >= 0 (Q10: metabolic_heat_flow >= minimum_heat_flow · q10^((T_core − T_setpoint)/10)).
function _heat_balance_residuals_weighted!(residuals, effectors, p)
    core_temperature       = effectors[1] * u"K"
    skin_temperature       = effectors[2] * u"K"
    insulation_temperature = effectors[3] * u"K"
    metabolic_heat_flow    = exp(effectors[4]) * u"W"
    flesh_conductivity     = effectors[5] * u"W/m/K"
    panting_rate           = effectors[6]
    skin_wetness           = effectors[7]
    insulation_depth       = effectors[8] * u"m"
    axis_ratio_b           = effectors[9]

    r = HeatExchange.nlp_residuals(p.nlp_packed, core_temperature, skin_temperature,
        insulation_temperature, metabolic_heat_flow, flesh_conductivity,
        panting_rate, skin_wetness, insulation_depth, axis_ratio_b)
    residuals[1] = ustrip(u"W", r.residuals[1])   # energy_balance
    residuals[2] = ustrip(u"W", r.residuals[2])   # internal_conduction
    residuals[3] = ustrip(u"K", r.residuals[3])   # skin_temperature
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

function _objective_value_multisided(effectors, opt)
    skin_temp_mean = (effectors[2] + effectors[4]) / 2
    opt.core_temperature_penalty * ((effectors[1] - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
    opt.metabolic_heat_penalty   * ((exp(effectors[6]) - opt.heat_flow_min_W)   / opt.heat_flow_range)^2        +
    opt.gradient_penalty         * ((effectors[1] - skin_temp_mean - opt.target_gradient) / opt.gradient_range)^2 +
    opt.panting_penalty          * ((effectors[8] - 1.0) / opt.panting_rate_range)^2                            +
    opt.skin_wetness_penalty     * ((effectors[9] - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# residuals[1:5] == 0 (dorsal surface, dorsal skin, ventral surface, ventral skin, whole-organism).
# residuals[6]  >= 0 (Q10: metabolic_heat_flow >= minimum_heat_flow · q10^((T_core − T_setpoint)/10)).
function _heat_balance_residuals_multisided!(residuals, effectors, p)
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

    r = HeatExchange.nlp_residuals(p.nlp_packed, core_temperature,
        dorsal_skin_temperature, dorsal_insulation_temperature,
        ventral_skin_temperature, ventral_insulation_temperature,
        metabolic_heat_flow, flesh_conductivity, panting_rate, skin_wetness,
        insulation_depth, axis_ratio_b)
    residuals[1] = ustrip(u"W", r.residuals[1])   # dorsal surface balance
    residuals[2] = ustrip(u"K", r.residuals[2])   # dorsal skin temperature
    residuals[3] = ustrip(u"W", r.residuals[3])   # ventral surface balance
    residuals[4] = ustrip(u"K", r.residuals[4])   # ventral skin temperature
    residuals[5] = ustrip(u"W", r.residuals[5])   # whole-organism balance
    residuals[6] = exp(effectors[6]) - p.minimum_heat_flow * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
    return nothing
end

# =============================================================================
# Lagrangian — `sigma * f(x) + dot(lambda, g(x))`. IPOPT's `eval_h` callback
# wants the Hessian of *this*, not of the objective alone. We define it as a
# single scalar function so a forward-over-reverse Enzyme pass (one per
# decision-variable direction) yields one Hessian column at a time.
#
# The residual buffer `r_buf` is passed in explicitly rather than allocated
# inside. Under the inner reverse pass we mark it as `Duplicated(r_buf, r_adj)`
# so Enzyme writes adjoints into the supplied `r_adj` rather than allocating
# a fresh shadow each call. Under the outer forward pass, both `r_buf` and
# `r_adj` are themselves marked `Duplicated` so their forward tangents
# propagate. All four arrays live on the cache; the callbacks just thread
# them through.
#
# `p` is a NamedTuple with fields `sigma::Float64`, `lambda::Vector{Float64}`,
# `state::_IPOPTState`. The state reference lets the Lagrangian see the
# *current* `opt_pars` / `nlp_pars` set by `_*_inputs!` — no per-call
# bookkeeping inside the Hessian callback.
# =============================================================================

function _lagrangian_weighted(x, r_buf, p)
    _heat_balance_residuals_weighted!(r_buf, x, p.state.nlp_pars)
    return p.sigma * _objective_value_weighted(x, p.state.opt_pars) +
           p.lambda[1]*r_buf[1] + p.lambda[2]*r_buf[2] +
           p.lambda[3]*r_buf[3] + p.lambda[4]*r_buf[4]
end

function _lagrangian_multisided(x, r_buf, p)
    _heat_balance_residuals_multisided!(r_buf, x, p.state.nlp_pars)
    return p.sigma * _objective_value_multisided(x, p.state.opt_pars) +
           p.lambda[1]*r_buf[1] + p.lambda[2]*r_buf[2] + p.lambda[3]*r_buf[3] +
           p.lambda[4]*r_buf[4] + p.lambda[5]*r_buf[5] + p.lambda[6]*r_buf[6]
end

# =============================================================================
# Per-call inputs (bounds, u0, opt_pars, nlp_pars)
#
# These helpers populate the mutable state container and bound/u0 vectors that
# the cache reuses. They are also called from the fresh-per-call path with
# disposable buffers, so the bounds/u0 computation lives in exactly one place.
# Each writes Float64 fields only — the NamedTuple field types are stable
# across calls, so mutable-state field assignment never reboxes.
# =============================================================================

function _weighted_inputs!(lb, ub, u0, nlp_packed::HeatExchange.WeightedMeanNLPPacked,
                            organism, environment, limits, metab_pars, int_cond, evap_pars,
                            metabolic_heat_flow_init, skin_temperature_init,
                            insulation_temperature_init)
    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metab_pars.core_temperature)
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

    lb[1] = core_temperature_min;    lb[2] = skin_temperature_min;  lb[3] = skin_temperature_min
    lb[4] = log(heat_flow_min_W);    lb[5] = flesh_conductivity_min; lb[6] = 1.0
    lb[7] = skin_wetness_min;        lb[8] = insulation_depth_min;  lb[9] = aspect_ratio_min
    ub[1] = core_temperature_max;    ub[2] = skin_temperature_max;  ub[3] = skin_temperature_max
    ub[4] = log(heat_flow_max_W);    ub[5] = flesh_conductivity_max; ub[6] = panting_rate_max
    ub[7] = skin_wetness_max;        ub[8] = insulation_depth_max;  ub[9] = aspect_ratio_max

    heat_flow_init = max(ustrip(u"W", metabolic_heat_flow_init), heat_flow_min_W)
    u0[1] = clamp(setpoint_temperature_K,                    lb[1], ub[1])
    u0[2] = clamp(ustrip(u"K", skin_temperature_init),       lb[2], ub[2])
    u0[3] = clamp(ustrip(u"K", insulation_temperature_init), lb[3], ub[3])
    u0[4] = clamp(log(heat_flow_init),                       lb[4], ub[4])
    u0[5] = clamp(ustrip(u"W/m/K", int_cond.flesh_conductivity), lb[5], ub[5])
    u0[6] = clamp(1.0,                                       lb[6], ub[6])
    u0[7] = clamp(Float64(evap_pars.skin_wetness),           lb[7], ub[7])
    u0[8] = clamp(insulation_depth_max,                      lb[8], ub[8])   # start with erected insulation
    u0[9] = clamp(aspect_ratio_min,                          lb[9], ub[9])   # start curled

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
    opt_pars = (;
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
    nlp_pars = (;
        nlp_packed,
        minimum_heat_flow    = Float64(heat_flow_min_W),
        q10                  = Float64(metab_pars.q10),
        setpoint_temperature = Float64(setpoint_temperature_K),
    )
    return opt_pars, nlp_pars
end

function _multisided_inputs!(lb, ub, u0, nlp_packed::HeatExchange.MultiSidedNLPPacked,
                              organism, environment, limits, metab_pars, int_cond, evap_pars,
                              metabolic_heat_flow_init, _, _)
    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metab_pars.core_temperature)
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

    # 11-variable layout: see header table at top of MultiSidedNLP block.
    lb[1]  = core_temperature_min
    lb[2]  = skin_temperature_min; lb[3]  = skin_temperature_min   # dorsal skin, ins
    lb[4]  = skin_temperature_min; lb[5]  = skin_temperature_min   # ventral skin, ins
    lb[6]  = log(heat_flow_min_W); lb[7]  = flesh_conductivity_min
    lb[8]  = 1.0;                  lb[9]  = skin_wetness_min
    lb[10] = insulation_depth_min; lb[11] = aspect_ratio_min
    ub[1]  = core_temperature_max
    ub[2]  = skin_temperature_max; ub[3]  = skin_temperature_max
    ub[4]  = skin_temperature_max; ub[5]  = skin_temperature_max
    ub[6]  = log(heat_flow_max_W); ub[7]  = flesh_conductivity_max
    ub[8]  = panting_rate_max;     ub[9]  = skin_wetness_max
    ub[10] = insulation_depth_max; ub[11] = aspect_ratio_max

    pp = nlp_packed.p
    heat_flow_init = max(ustrip(u"W", metabolic_heat_flow_init), heat_flow_min_W)
    u0[1]  = clamp(setpoint_temperature_K,                                       lb[1],  ub[1])
    u0[2]  = clamp(ustrip(u"K", pp.initial_dorsal_skin_temperature),             lb[2],  ub[2])
    u0[3]  = clamp(ustrip(u"K", pp.initial_dorsal_insulation_temperature),       lb[3],  ub[3])
    u0[4]  = clamp(ustrip(u"K", pp.initial_ventral_skin_temperature),            lb[4],  ub[4])
    u0[5]  = clamp(ustrip(u"K", pp.initial_ventral_insulation_temperature),      lb[5],  ub[5])
    u0[6]  = clamp(log(heat_flow_init),                                          lb[6],  ub[6])
    u0[7]  = clamp(ustrip(u"W/m/K", int_cond.flesh_conductivity),                lb[7],  ub[7])
    u0[8]  = clamp(1.0,                                                          lb[8],  ub[8])
    u0[9]  = clamp(Float64(evap_pars.skin_wetness),                              lb[9],  ub[9])
    u0[10] = clamp(insulation_depth_max,                                         lb[10], ub[10])
    u0[11] = clamp(aspect_ratio_min,                                             lb[11], ub[11])

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
    opt_pars = (;
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
    nlp_pars = (;
        nlp_packed,
        minimum_heat_flow    = Float64(heat_flow_min_W),
        q10                  = Float64(metab_pars.q10),
        setpoint_temperature = Float64(setpoint_temperature_K),
    )
    return opt_pars, nlp_pars
end

# =============================================================================
# Warm-start cache
#
# `_IPOPTState` is a mutable struct that carries the two NamedTuples the
# Enzyme-AD callbacks consume (objective penalties + NLP physics parameters).
# The Ipopt callbacks hold a reference to this state and read the fields on
# every invocation. To set up a new solve we just write fresh NamedTuples
# into the two fields — same concrete types, no reboxing — and the callbacks
# pick them up automatically without rebuilding.
#
# Both fields are parametrically typed so accesses stay type-stable: as long
# as the per-call `opt_pars` and `nlp_pars` NamedTuples have the same key/type
# signature (guaranteed by the `_*_inputs!` helpers above), assignment is a
# bit-for-bit field overwrite. If a structurally different organism or
# environment is fed in (e.g. a different body shape), construct a new cache.
# =============================================================================

mutable struct _IPOPTState{OP,NP}
    opt_pars::OP
    nlp_pars::NP
end

# Pre-allocated workspaces consumed by the Enzyme-AD calls inside the Ipopt
# callbacks. One instance per cache, sized at construction from (n, m).
#
# Jacobian section (used by `eval_jac_g`):
#   jac_r           length m, residuals primal
#   jac_dr          length m, reverse-mode seed (unit vector for row j)
#   jac_dx          length n, gradient output → one Jacobian row
#
# Hessian section (used by `eval_h`, forward-over-reverse):
#   hess_g          length n, inner-Reverse gradient
#   hess_dg         length n, outer-Forward dual → one Hessian column
#   hess_seed       length n, outer-Forward seed (unit vector for column j)
#   hess_lambda     length m, current IPOPT λ_g (copied per `eval_h` call)
#   hess_r          length m, Lagrangian residual buffer
#   hess_dr         length m, inner-Reverse adjoint of `hess_r`
#   hess_r_fwd      length m, outer-Forward tangent of `hess_r`
#   hess_dr_fwd     length m, outer-Forward tangent of `hess_dr`
struct IpoptCallbackBuffers
    jac_r::Vector{Float64}
    jac_dr::Vector{Float64}
    jac_dx::Vector{Float64}
    hess_g::Vector{Float64}
    hess_dg::Vector{Float64}
    hess_seed::Vector{Float64}
    hess_lambda::Vector{Float64}
    hess_r::Vector{Float64}
    hess_dr::Vector{Float64}
    hess_r_fwd::Vector{Float64}
    hess_dr_fwd::Vector{Float64}
end

function IpoptCallbackBuffers(n::Int, m::Int)
    IpoptCallbackBuffers(
        zeros(m), zeros(m), zeros(n),               # jac_*
        zeros(n), zeros(n), zeros(n), zeros(m),     # hess_g/dg/seed/lambda
        zeros(m), zeros(m), zeros(m), zeros(m),     # hess_r/dr and forward duals
    )
end

"""
    IPOPTSolverCache

Reusable solver state for the IPOPT thermoregulation path. Holds the
Ipopt-callback closures, the per-call bound and `u0` buffers, the
reverse-mode Jacobian scratch vectors, the Julia-side primal+dual
warm-start vectors, and a mutable `_IPOPTState` whose fields the callbacks
read on every solve.

Construct once per organism + smoothing + NLP-strategy combination:

    cache = IPOPTSolverCache(control, organism, environment,
                              metabolic_heat_flow_init, skin_T_init, ins_T_init)

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
  path seeds only `M`, `skin_T`, and `ins_T` from user-passed values; the
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
side state (`prev_x`, `prev_mult_*`) carries the warm-start info across
rebuilds.
"""
mutable struct IPOPTSolverCache{S<:HeatExchange.NLPStrategy, ST<:_IPOPTState,
                                  EF, EG, EGF, EJG, EH}
    nlp_strategy::S
    state::ST
    # Cached callback closures. These read `state.opt_pars` / `state.nlp_pars`
    # at solve time, so updating the state retargets every callback without
    # rebuilding the closures. Stored on the cache so the same Function objects
    # are passed to every CreateIpoptProblem call — that way Julia compiles a
    # single specialised IpoptProblem type once and reuses it.
    eval_f::EF
    eval_g::EG
    eval_grad_f::EGF
    eval_jac_g::EJG
    eval_h::EH
    # Per-call tight bounds — written by `_*_inputs!` each solve to match the
    # fresh `thermoregulate(..., organism, env, M, sT, iT)` path. Ipopt bakes
    # bounds into the C-side problem at CreateIpoptProblem time, so we rebuild
    # the C handle every solve (cheap — a single malloc + struct init) rather
    # than loosening the bounds to keep one handle alive.
    lb::Vector{Float64}
    ub::Vector{Float64}
    u0::Vector{Float64}
    lcons::Vector{Float64}
    ucons::Vector{Float64}
    n::Int                              # variables (9 weighted / 11 multisided)
    m::Int                              # constraints (4 weighted / 6 multisided)
    # All AD scratch buffers consumed by the Ipopt callbacks. See the
    # `IpoptCallbackBuffers` definition above.
    buffers::IpoptCallbackBuffers
    # Warm-start state, persisted Julia-side across IpoptProblem rebuilds.
    # On the first solve `has_prev = false` so we feed `u0` as the primal and
    # let Ipopt initialise its own multipliers. After IpoptSolve we copy
    # `prob.x`, `prob.mult_g`, `prob.mult_x_L`, `prob.mult_x_U` into these
    # vectors. On subsequent solves we restore them before calling IpoptSolve
    # with `warm_start_init_point = "yes"`.
    prev_x::Vector{Float64}
    prev_mult_g::Vector{Float64}
    prev_mult_x_L::Vector{Float64}
    prev_mult_x_U::Vector{Float64}
    has_prev::Bool
    # User-provided variable scaling (one entry per decision variable). See
    # `_scaling_weighted` / `_scaling_multisided` for the values. Passed to
    # Ipopt via `SetIpoptProblemScaling` together with `obj_scaling = 1.0`
    # and an all-ones `g_scaling`, after `CreateIpoptProblem` but before
    # `IpoptSolve`. Allocated once at cache construction and reused for
    # every solve.
    x_scaling::Vector{Float64}
    g_scaling::Vector{Float64}
end

# Direct-Ipopt callbacks. Strategy-specific behaviour comes from the three
# top-level physics functions passed in (`obj_fn`, `res_fn!`, `lag_fn`) —
# the rest of the wiring (Enzyme calls, sparsity, lower-triangle indexing)
# is shared between WeightedMeanNLP and MultiSidedNLP.
#
# Closures read `state.opt_pars` / `state.nlp_pars` at solve time so
# refreshing the state retargets every callback without rebuilding them.
#
# Jacobian is treated as fully dense — for n=9/11 sparse bookkeeping costs
# more than the saved entries. Sparsity convention: row-major dense, 1-based,
# entry k = (i-1)*n + j ↔ (rows[k], cols[k]) = (i, j).
#
# Hessian uses forward-over-reverse Enzyme: the inner Reverse closure
# `lag_grad!` computes ∇ₓL into `b.hess_g`, the outer Forward seeds e_j and
# pulls back column j of H into `b.hess_dg`. Sparsity is dense lower
# triangle (i ≥ j), entry k = i(i-1)/2 + j ↔ (i, j).
function _build_ipopt_callbacks(obj_fn::F1, res_fn!::F2, lag_fn::F3,
                                 state::_IPOPTState, n::Int, m::Int,
                                 b::IpoptCallbackBuffers) where {F1, F2, F3}
    eval_f(x)    = obj_fn(x, state.opt_pars)
    eval_g(x, g) = res_fn!(g, x, state.nlp_pars)
    function eval_grad_f(x, grad_f)
        fill!(grad_f, 0)
        Enzyme.autodiff(Enzyme.Reverse, obj_fn,
                        Enzyme.Active,
                        Enzyme.Duplicated(x, grad_f),
                        Enzyme.Const(state.opt_pars))
        return nothing
    end
    function eval_jac_g(x, rows, cols, values)
        if values === nothing
            k = 1
            @inbounds for i in 1:m, j in 1:n
                rows[k] = i; cols[k] = j; k += 1
            end
        else
            @inbounds for i in 1:m
                fill!(b.jac_r, 0); fill!(b.jac_dr, 0); fill!(b.jac_dx, 0); b.jac_dr[i] = 1.0
                Enzyme.autodiff(Enzyme.Reverse, res_fn!,
                                Enzyme.Const,
                                Enzyme.Duplicated(b.jac_r, b.jac_dr),
                                Enzyme.Duplicated(x, b.jac_dx),
                                Enzyme.Const(state.nlp_pars))
                for j in 1:n
                    values[(i-1)*n + j] = b.jac_dx[j]
                end
            end
        end
        return nothing
    end
    # Inner reverse-gradient closure runs Enzyme.Reverse on `lag_fn(x, r, p)`,
    # marking both `x` and `r` as Duplicated so Enzyme writes adjoints into
    # preallocated buffers instead of allocating fresh shadows.
    function lag_grad!(g, x, r_buf, r_adj, p)
        fill!(g, 0)
        fill!(r_adj, 0)
        Enzyme.autodiff(Enzyme.Reverse, lag_fn,
                        Enzyme.Active,
                        Enzyme.Duplicated(x, g),
                        Enzyme.Duplicated(r_buf, r_adj),
                        Enzyme.Const(p))
        return nothing
    end
    function eval_h(x, rows, cols, obj_factor, lambda, values)
        if values === nothing
            k = 1
            @inbounds for i in 1:n, j in 1:i
                rows[k] = i; cols[k] = j; k += 1
            end
        else
            # IPOPT passes the current obj_factor (`sigma`) and constraint
            # multipliers; copy into the cache so the Lagrangian closure (via
            # `p.sigma` / `p.lambda`) sees them under nested AD.
            copyto!(b.hess_lambda, lambda)
            p = (; sigma = obj_factor, lambda = b.hess_lambda, state)
            @inbounds for j in 1:n
                fill!(b.hess_seed, 0); b.hess_seed[j] = 1.0
                fill!(b.hess_g,    0); fill!(b.hess_dg,     0)
                fill!(b.hess_r,    0); fill!(b.hess_r_fwd,  0)
                fill!(b.hess_dr,   0); fill!(b.hess_dr_fwd, 0)
                # Outer Forward propagates tangents through everything
                # `lag_grad!` reads/writes — including `r` and `r_adj` —
                # so every mutable buffer needs its forward dual here.
                Enzyme.autodiff(Enzyme.Forward, lag_grad!,
                                Enzyme.Const,
                                Enzyme.Duplicated(b.hess_g, b.hess_dg),
                                Enzyme.Duplicated(x, b.hess_seed),
                                Enzyme.Duplicated(b.hess_r,  b.hess_r_fwd),
                                Enzyme.Duplicated(b.hess_dr, b.hess_dr_fwd),
                                Enzyme.Const(p))
                for i in j:n
                    k = (i * (i - 1)) ÷ 2 + j
                    values[k] = b.hess_dg[i]
                end
            end
        end
        return nothing
    end
    return eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h
end

# Apply the project-wide Ipopt options to a fresh problem. Called once per
# `CreateIpoptProblem` — i.e. per cached solve — alongside the warm-start
# toggle. `print_level` is set per call because callers can flip `verbose`
# between solves; everything else is constant for the cache's lifetime.
function _apply_ipopt_options!(prob::Ipopt.IpoptProblem, verbose::Bool, warm_start::Bool)
    # `hessian_approximation` defaults to "exact" — we register a real
    # `eval_h` (Lagrangian Hessian via Enzyme forward-over-reverse), so no
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
    # `SetIpoptProblemScaling` after CreateIpoptProblem (see
    # `_scaling_weighted` / `_scaling_multisided` for values). Telling Ipopt
    # we'll supply our own scales skips its `gradient-based` heuristic pass
    # (an extra Jacobian eval at the initial point) and gives the solver a
    # well-conditioned internal problem from iteration zero.
    Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method",    "user-scaling")
    return prob
end

# ──────────────────────────────────────────────────────────────────────────
# TEMPORARY: variable scaling by raw index position.
#
# These two functions hard-code Ipopt's per-variable scale factors keyed by
# integer index. The variable layout (`x[1] = core_T`, …) only exists as
# comments in the residual function headers — there is no shared schema
# between (a) `_*_inputs!` (which writes bounds and `u0`), (b)
# `_objective_value_*` / `_heat_balance_residuals_*` (which decode `x[i]`),
# and (c) these scaling functions. Reordering / inserting / renaming a
# variable in any one of those three places without touching the other two
# silently mis-scales the problem.
#
# Proper fix: extend the ModelParameters.jl `Param` wrappers that already
# carry `bounds` / `units` / `val` (see `HeatExchange/src/traits.jl`) with a
# `scaling` field. Each decision variable's source `Param` then owns its
# scale alongside its bounds, and `_*_inputs!` reads scale and bounds from
# the same place — no positional duplication. Requires adding `scaling` to
# the relevant traits (`ThermoregulationLimits`, etc.) and threading it
# through `_*_inputs!` to fill `cache.x_scaling`. Worth doing the next time
# we touch the trait schema.
# ──────────────────────────────────────────────────────────────────────────
function _scaling_weighted()
    x_scaling = zeros(9)
    x_scaling[1] = 1/300.0     # core_T ~ 300 K
    x_scaling[2] = 1/300.0     # skin_T ~ 300 K
    x_scaling[3] = 1/300.0     # insulation_T ~ 300 K
    x_scaling[4] = 1.0         # log(metabolic_heat_flow) ~ 0 (W ~ 1)
    x_scaling[5] = 1.0         # flesh_conductivity ~ 1 W/m/K
    x_scaling[6] = 1.0         # panting ∈ [1, 15] → internal [1, 15]
    x_scaling[7] = 50.0        # skin_wetness ∈ [0, 0.05] → internal [0, 2.5]
    x_scaling[8] = 50.0        # insulation_depth ∈ [0, 0.05 m] → internal [0, 2.5]
    x_scaling[9] = 1.0         # axis_ratio ~ 1
    return x_scaling
end

function _scaling_multisided()
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

# Construct the cache for a given strategy. The first solve still primes type
# inference, but every subsequent solve reuses every long-lived allocation.
function IPOPTSolverCache(control::IPOPTControl, organism::Organism, environment::NamedTuple,
                          metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init)
    metab_pars = metabolism_pars(organism)
    limits     = thermoregulation(organism)
    int_cond   = conduction_pars_internal(organism)
    evap_pars  = evaporation_pars(organism)

    nlp_packed = HeatExchange.nlp_pack(control.nlp_strategy, organism, environment,
                                       skin_temperature_init, insulation_temperature_init;
                                       smoothing = control.smoothing)

    return _build_cache(control.nlp_strategy, nlp_packed, organism, environment,
                        limits, metab_pars, int_cond, evap_pars,
                        metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init)
end

function _build_cache(strategy::HeatExchange.WeightedMeanNLP,
                      nlp_packed::HeatExchange.WeightedMeanNLPPacked,
                      organism, environment, limits, metab_pars, int_cond, evap_pars,
                      M_init, sT_init, iT_init)
    n, m = 9, 4
    lb = zeros(n); ub = zeros(n); u0 = zeros(n)
    opt_pars, nlp_pars = _weighted_inputs!(lb, ub, u0, nlp_packed,
        organism, environment, limits, metab_pars, int_cond, evap_pars,
        M_init, sT_init, iT_init)
    state    = _IPOPTState(opt_pars, nlp_pars)
    buffers  = IpoptCallbackBuffers(n, m)
    eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h =
        _build_ipopt_callbacks(_objective_value_weighted,
                                _heat_balance_residuals_weighted!,
                                _lagrangian_weighted,
                                state, n, m, buffers)
    return IPOPTSolverCache(strategy, state,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h,
        lb, ub, u0, [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, Inf],
        n, m, buffers,
        zeros(n), zeros(m), zeros(n), zeros(n), false,
        _scaling_weighted(), ones(m))
end

function _build_cache(strategy::HeatExchange.MultiSidedNLP,
                      nlp_packed::HeatExchange.MultiSidedNLPPacked,
                      organism, environment, limits, metab_pars, int_cond, evap_pars,
                      M_init, sT_init, iT_init)
    n, m = 11, 6
    lb = zeros(n); ub = zeros(n); u0 = zeros(n)
    opt_pars, nlp_pars = _multisided_inputs!(lb, ub, u0, nlp_packed,
        organism, environment, limits, metab_pars, int_cond, evap_pars,
        M_init, sT_init, iT_init)
    state    = _IPOPTState(opt_pars, nlp_pars)
    buffers  = IpoptCallbackBuffers(n, m)
    eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h =
        _build_ipopt_callbacks(_objective_value_multisided,
                                _heat_balance_residuals_multisided!,
                                _lagrangian_multisided,
                                state, n, m, buffers)
    return IPOPTSolverCache(strategy, state,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h,
        lb, ub, u0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, Inf],
        n, m, buffers,
        zeros(n), zeros(m), zeros(n), zeros(n), false,
        _scaling_multisided(), ones(m))
end

"""
    reset_warm_start!(cache::IPOPTSolverCache)

Discard the previous solution stored in the cache. The next solve will use
the cold initial guess derived from `metabolic_heat_flow_init`,
`skin_temperature_init` and `insulation_temperature_init`. Use this when the
trajectory takes a large discontinuous step (e.g. switching organism,
changing q10 regime) where seeding from the last solution would land the
optimiser in a bad basin.
"""
reset_warm_start!(cache::IPOPTSolverCache) = (cache.has_prev = false; cache)

# =============================================================================
# Per-call solve paths.
# Both paths delegate to `_ipopt_solve!`. The cached path refreshes the
# mutable `_IPOPTState` and tight-bound buffers in place, then calls into
# `_ipopt_solve!` which restores the previous primal+dual warm-start state
# before `IpoptSolve`. The fresh path builds a single-use cache, runs one
# cold `_ipopt_solve!`, and lets the GC collect the cache.
# =============================================================================

# Refresh the cache's state and bounds for a new solve, then run IPOPT.
# Splits inputs computation by strategy via dispatch on cache.nlp_strategy so
# the unused branch is dead-code-eliminated and no Union return appears on the
# hot path.
function _solve_cached!(cache::IPOPTSolverCache{<:HeatExchange.WeightedMeanNLP},
                        nlp_packed, organism, environment, limits, metab_pars,
                        int_cond, evap_pars, M_init, sT_init, iT_init; verbose)
    opt_pars, nlp_pars = _weighted_inputs!(cache.lb, cache.ub, cache.u0, nlp_packed,
        organism, environment, limits, metab_pars, int_cond, evap_pars,
        M_init, sT_init, iT_init)
    cache.state.opt_pars = opt_pars
    cache.state.nlp_pars = nlp_pars
    return _ipopt_solve!(cache, verbose)
end

function _solve_cached!(cache::IPOPTSolverCache{<:HeatExchange.MultiSidedNLP},
                        nlp_packed, organism, environment, limits, metab_pars,
                        int_cond, evap_pars, M_init, sT_init, iT_init; verbose)
    opt_pars, nlp_pars = _multisided_inputs!(cache.lb, cache.ub, cache.u0, nlp_packed,
        organism, environment, limits, metab_pars, int_cond, evap_pars,
        M_init, sT_init, iT_init)
    cache.state.opt_pars = opt_pars
    cache.state.nlp_pars = nlp_pars
    return _ipopt_solve!(cache, verbose)
end

# Direct-Ipopt solve with primal+dual warm-start.
#
# Each call rebuilds the C-side IpoptProblem so the new tight bounds in
# `cache.lb` / `cache.ub` take effect — same bounds the fresh
# `thermoregulate(..., organism, env, M, sT, iT)` path uses, so solutions
# stay bit-identical to that path's first iterate. CreateIpoptProblem is a
# malloc + a few struct inits; well under a percent of an IPOPT iteration's
# cost.
#
# Returns the IpoptProblem so the caller can pull the final `x` (and
# optionally inspect the multipliers via `cache.prev_*`).
function _ipopt_solve!(cache::IPOPTSolverCache, verbose::Bool)
    prob = Ipopt.CreateIpoptProblem(
        cache.n, cache.lb, cache.ub,
        cache.m, cache.lcons, cache.ucons,
        cache.m * cache.n,                       # nele_jac: dense
        (cache.n * (cache.n + 1)) ÷ 2,           # nele_hess: dense lower triangle
        cache.eval_f, cache.eval_g, cache.eval_grad_f, cache.eval_jac_g,
        cache.eval_h,
    )
    _apply_ipopt_options!(prob, verbose, cache.has_prev)
    # User-scaling pairs with `nlp_scaling_method = "user-scaling"` in the
    # options above. `obj_scaling = 1.0` is fine — the objective is already
    # a sum of dimensionless squared penalties so it's intrinsically O(1).
    Ipopt.SetIpoptProblemScaling(prob, 1.0, cache.x_scaling, cache.g_scaling)

    if cache.has_prev
        # Restore the previous primal+dual, with the primal clamped into the
        # new tight bounds. `skin_T_min = air_T − 5` shifts between solves;
        # an out-of-range initial x makes Ipopt spend its first iteration
        # snapping back to the feasible box. Dual multipliers can carry over
        # as-is: inactive bounds have ≈0 multipliers (no information lost),
        # active bounds get the previous shadow value as a good guess.
        @inbounds for i in 1:cache.n
            prob.x[i] = clamp(cache.prev_x[i], cache.lb[i], cache.ub[i])
        end
        copyto!(prob.mult_g,   cache.prev_mult_g)
        copyto!(prob.mult_x_L, cache.prev_mult_x_L)
        copyto!(prob.mult_x_U, cache.prev_mult_x_U)
    else
        # Cold start — seed primal from the inputs-helper defaults, let Ipopt
        # initialise its own multipliers (warm_start_init_point = "no" above).
        @inbounds for i in 1:cache.n
            prob.x[i] = cache.u0[i]
        end
    end

    Ipopt.IpoptSolve(prob)

    # Persist primal+dual for the next solve.
    copyto!(cache.prev_x,       prob.x)
    copyto!(cache.prev_mult_g,  prob.mult_g)
    copyto!(cache.prev_mult_x_L, prob.mult_x_L)
    copyto!(cache.prev_mult_x_U, prob.mult_x_U)
    cache.has_prev = true
    return prob
end

# Fresh-per-call entrypoint — uses the same direct-Ipopt machinery as the
# cache, but allocates a single-use cache, runs one solve, and lets the GC
# clean up. `has_prev` starts at `false` so `_ipopt_solve!` treats it as a
# cold start (no warm-start multipliers, default Ipopt heuristic for initial
# duals). Bounds are the tight per-call values from `_*_inputs!`, matching
# the cached path's first iterate bit-for-bit.
function _run_ipopt(
    nlp_packed::HeatExchange.WeightedMeanNLPPacked,
    organism, environment, limits, metab_pars, int_cond, evap_pars,
    metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init;
    verbose,
)
    cache = _build_cache(HeatExchange.WeightedMeanNLP(), nlp_packed, organism, environment,
        limits, metab_pars, int_cond, evap_pars,
        metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init)
    prob = _ipopt_solve!(cache, verbose)
    return _assemble_weighted(nlp_packed, organism, environment, prob.x)
end

function _run_ipopt(
    nlp_packed::HeatExchange.MultiSidedNLPPacked,
    organism, environment, limits, metab_pars, int_cond, evap_pars,
    metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init;
    verbose,
)
    cache = _build_cache(HeatExchange.MultiSidedNLP(), nlp_packed, organism, environment,
        limits, metab_pars, int_cond, evap_pars,
        metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init)
    prob = _ipopt_solve!(cache, verbose)
    return _assemble_multisided(nlp_packed, organism, environment, prob.x)
end

# Unit reattachment + delegation to HeatExchange's output builder. Kept tiny
# and shape-agnostic so both fresh and cached paths share the assembly step.
function _assemble_weighted(nlp_packed, organism, environment, x_sol)
    return HeatExchange.nlp_assemble_output(nlp_packed, organism, environment,
        x_sol[1] * u"K", x_sol[2] * u"K", x_sol[3] * u"K", exp(x_sol[4]) * u"W",
        x_sol[5] * u"W/m/K", x_sol[6], x_sol[7], x_sol[8] * u"m", x_sol[9])
end

function _assemble_multisided(nlp_packed, organism, environment, x_sol)
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
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment,
                   metabolic_heat_flow_init, skin_temperature_init,
                   insulation_temperature_init; verbose=false)
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment,
                   metabolic_heat_flow_init, skin_temperature_init,
                   insulation_temperature_init, cache::IPOPTSolverCache;
                   verbose=false)

Solve endotherm heat balance as a nonlinear program via IPOPT.

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
- `metabolic_heat_flow_init`: initial guess for metabolic heat generation (W)
- `skin_temperature_init`: initial guess for skin temperature (K)
- `insulation_temperature_init`: initial guess for insulation surface temperature (K)
- `cache` (optional): an `IPOPTSolverCache` produced by `IPOPTSolverCache(control, …)`.
  Reuse across calls in a sweep to skip ~190 KiB of per-solve setup. Must match the
  organism's body-shape type, NLP strategy and smoothing strategy — invalidate the
  cache and rebuild if any of those change.
"""
function thermoregulate(
    ::Endotherm,
    control::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    metabolic_heat_flow_init,
    skin_temperature_init,
    insulation_temperature_init;
    verbose = false,
)
    metab_pars = metabolism_pars(organism)
    limits     = thermoregulation(organism)
    int_cond   = conduction_pars_internal(organism)
    evap_pars  = evaporation_pars(organism)

    nlp_packed = HeatExchange.nlp_pack(control.nlp_strategy, organism, environment,
                                       skin_temperature_init, insulation_temperature_init;
                                       smoothing = control.smoothing)

    _run_ipopt(nlp_packed, organism, environment, limits, metab_pars, int_cond, evap_pars,
               metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init;
               verbose)
end

function thermoregulate(
    ::Endotherm,
    control::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    metabolic_heat_flow_init,
    skin_temperature_init,
    insulation_temperature_init,
    cache::IPOPTSolverCache;
    verbose = false,
)
    metab_pars = metabolism_pars(organism)
    limits     = thermoregulation(organism)
    int_cond   = conduction_pars_internal(organism)
    evap_pars  = evaporation_pars(organism)

    nlp_packed = HeatExchange.nlp_pack(control.nlp_strategy, organism, environment,
                                       skin_temperature_init, insulation_temperature_init;
                                       smoothing = control.smoothing)

    prob = _solve_cached!(cache, nlp_packed, organism, environment, limits, metab_pars,
        int_cond, evap_pars, metabolic_heat_flow_init, skin_temperature_init,
        insulation_temperature_init; verbose)

    return cache.nlp_strategy isa HeatExchange.WeightedMeanNLP ?
        _assemble_weighted(nlp_packed, organism, environment, prob.x) :
        _assemble_multisided(nlp_packed, organism, environment, prob.x)
end
