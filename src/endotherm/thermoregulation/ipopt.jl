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
# WeightedMeanNLP helpers
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

function _run_ipopt(
    nlp_packed::HeatExchange.WeightedMeanNLPPacked,
    organism, environment, limits, metab_pars, int_cond, evap_pars,
    metabolic_heat_flow_init, skin_temperature_init, insulation_temperature_init;
    verbose,
)
    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metab_pars.core_temperature)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - 5.0
    skin_temperature_max   = core_temperature_max + 5.0
    heat_flow_min_W        = ustrip(u"W", limits.minimum_heat_flow)
    heat_flow_max_W        = heat_flow_min_W * 20.0
    log_heat_flow_min      = log(heat_flow_min_W)
    log_heat_flow_max      = log(heat_flow_max_W)
    flesh_conductivity_min = ustrip(u"W/m/K", limits.flesh_conductivity.reference)
    flesh_conductivity_max = ustrip(u"W/m/K", limits.flesh_conductivity.max)
    panting_rate_min       = 1.0
    panting_rate_max       = Float64(limits.panting.pant.max)
    skin_wetness_min       = Float64(limits.skin_wetness.reference)
    skin_wetness_max       = Float64(limits.skin_wetness.max)
    insulation_depth_min   = ustrip(u"m", limits.insulation.dorsal.reference)
    insulation_depth_max   = ustrip(u"m", limits.insulation.dorsal.max)
    axis_ratio_min         = Float64(limits.axis_ratio_factor.reference)
    axis_ratio_max         = organism.body.shape isa Sphere ? axis_ratio_min : Float64(limits.axis_ratio_factor.max)

    lower_bounds = [core_temperature_min, skin_temperature_min, skin_temperature_min, log_heat_flow_min,
                    flesh_conductivity_min, panting_rate_min, skin_wetness_min,
                    insulation_depth_min, axis_ratio_min]
    upper_bounds = [core_temperature_max, skin_temperature_max, skin_temperature_max, log_heat_flow_max,
                    flesh_conductivity_max, panting_rate_max, skin_wetness_max,
                    insulation_depth_max, axis_ratio_max]

    heat_flow_init = max(ustrip(u"W", metabolic_heat_flow_init), heat_flow_min_W)
    initial_effectors = clamp.(
        [setpoint_temperature_K,
         ustrip(u"K", skin_temperature_init),
         ustrip(u"K", insulation_temperature_init),
         log(heat_flow_init),
         ustrip(u"W/m/K", int_cond.flesh_conductivity),
         1.0,
         Float64(evap_pars.skin_wetness),
         insulation_depth_max,   # start with erected insulation
         axis_ratio_min],      # start curled
        lower_bounds, upper_bounds,
    )

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
    opt_pars = (;
        setpoint_temperature_K,
        heat_flow_min_W,
        core_temperature_penalty = limits.core_temperature_penalty,
        metabolic_heat_penalty   = limits.metabolic_heat_penalty,
        panting_penalty          = limits.panting_penalty,
        skin_wetness_penalty     = limits.skin_wetness_penalty,
        gradient_penalty         = limits.gradient_penalty,
        target_gradient          = limits.target_core_skin_gradient,
        core_temperature_range,
        gradient_range           = core_temperature_range,
        panting_rate_range       = max(panting_rate_max - 1.0, 1e-6),
        skin_wetness_min,
        skin_wetness_range       = max(skin_wetness_max - skin_wetness_min, 1e-6),
        heat_flow_range          = max(heat_flow_max_W - heat_flow_min_W, 1.0),
    )

    nlp_pars = (;
        nlp_packed,
        minimum_heat_flow    = heat_flow_min_W,
        q10                  = Float64(metab_pars.q10),
        setpoint_temperature = setpoint_temperature_K,
    )

    obj_fn(x, _)     = _objective_value_weighted(x, opt_pars)
    res_fn!(r, x, _) = _heat_balance_residuals_weighted!(r, x, nlp_pars)

    # Enzyme reverse-mode AD. Bypass DifferentiationInterface — the params
    # NamedTuples are already captured in the closures, and DI's vector packing
    # of `p` doesn't compose with Unitful values.
    # hess! and cons_h! are registered (required by IpoptOptimizer) but not
    # called at runtime when hessian_approximation="limited-memory" is set.
    function grad_fn!(g, x, _)
        fill!(g, 0)
        Enzyme.autodiff(Enzyme.Reverse, _objective_value_weighted,
                        Enzyme.Active,
                        Enzyme.Duplicated(x, g),
                        Enzyme.Const(opt_pars))
        return nothing
    end
    function cons_j_fn!(J, x, _)
        # Forward mode: 9 passes (one per input) instead of reverse's 4.
        # Reverse mode currently produces silent NaN for cols 2 (skin_T)
        # and 3 (ins_T) in this call chain — even though every called
        # function is type-stable per JET. Forward mode gives correct
        # finite values matching finite differences. The reverse-mode
        # NaN is an Enzyme limitation unrelated to type stability; revisit
        # when the Enzyme issue is resolved upstream.
        m = size(J, 1)
        n = length(x)
        r  = zeros(m)
        dr = zeros(m)
        dx = zeros(n)
        for j in 1:n
            fill!(r, 0); fill!(dr, 0); fill!(dx, 0); dx[j] = 1.0
            Enzyme.autodiff(Enzyme.Forward, _heat_balance_residuals_weighted!,
                            Enzyme.Const,
                            Enzyme.Duplicated(r, dr),
                            Enzyme.Duplicated(x, dx),
                            Enzyme.Const(nlp_pars))
            @views J[:, j] .= dr
        end
        return nothing
    end
    hess_fn!(H, x, _)        = (fill!(H, 0); nothing)
    cons_h_fn!(res, x, _)    = (foreach(r -> fill!(r, 0), res); nothing)

    optimization_func = OptimizationFunction(obj_fn, SciMLBase.NoAD();
        cons     = res_fn!,
        grad     = grad_fn!,
        hess     = hess_fn!,
        cons_j   = cons_j_fn!,
        cons_h   = cons_h_fn!,
    )
    optimization_prob = OptimizationProblem(optimization_func, initial_effectors, nothing;
        lb     = lower_bounds,
        ub     = upper_bounds,
        lcons  = [0.0, 0.0, 0.0, 0.0],
        ucons  = [0.0, 0.0, 0.0, Inf],
    )
    ipopt_sol = solve(optimization_prob,
        IpoptOptimizer(;
            hessian_approximation = "limited-memory",
            acceptable_tol        = 1e-3,
            acceptable_iter       = 5,
        );
        verbose,
        reltol   = 1e-4,
        maxiters = 300,
    )

    x_sol = ipopt_sol.u
    core_temperature       = x_sol[1] * u"K"
    skin_temperature       = x_sol[2] * u"K"
    insulation_temperature = x_sol[3] * u"K"
    metabolic_heat_flow    = exp(x_sol[4]) * u"W"
    flesh_conductivity     = x_sol[5] * u"W/m/K"
    panting_rate           = x_sol[6]
    skin_wetness           = x_sol[7]
    insulation_depth       = x_sol[8] * u"m"
    axis_ratio_b           = x_sol[9]

    return HeatExchange.nlp_assemble_output(nlp_pars.nlp_packed, organism, environment,
        core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow,
        flesh_conductivity, panting_rate, skin_wetness, insulation_depth, axis_ratio_b)
end

# =============================================================================
# MultiSidedNLP helpers
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

function _run_ipopt(
    nlp_packed::HeatExchange.MultiSidedNLPPacked,
    organism, environment, limits, metab_pars, int_cond, evap_pars,
    metabolic_heat_flow_init, ::Any, ::Any;
    verbose,
)
    air_temperature_K      = ustrip(u"K", environment.environment_vars.air_temperature)
    setpoint_temperature_K = ustrip(u"K", metab_pars.core_temperature)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - 5.0
    skin_temperature_max   = core_temperature_max + 5.0
    heat_flow_min_W        = ustrip(u"W", limits.minimum_heat_flow)
    heat_flow_max_W        = heat_flow_min_W * 20.0
    log_heat_flow_min      = log(heat_flow_min_W)
    log_heat_flow_max      = log(heat_flow_max_W)
    flesh_conductivity_min = ustrip(u"W/m/K", limits.flesh_conductivity.reference)
    flesh_conductivity_max = ustrip(u"W/m/K", limits.flesh_conductivity.max)
    panting_rate_min       = 1.0
    panting_rate_max       = Float64(limits.panting.pant.max)
    skin_wetness_min       = Float64(limits.skin_wetness.reference)
    skin_wetness_max       = Float64(limits.skin_wetness.max)
    insulation_depth_min   = ustrip(u"m", limits.insulation.dorsal.reference)
    insulation_depth_max   = ustrip(u"m", limits.insulation.dorsal.max)
    axis_ratio_min         = Float64(limits.axis_ratio_factor.reference)
    axis_ratio_max         = organism.body.shape isa Sphere ? axis_ratio_min : Float64(limits.axis_ratio_factor.max)

    lower_bounds = [core_temperature_min,
                    skin_temperature_min, skin_temperature_min,   # dorsal skin, ins
                    skin_temperature_min, skin_temperature_min,   # ventral skin, ins
                    log_heat_flow_min,
                    flesh_conductivity_min, panting_rate_min, skin_wetness_min,
                    insulation_depth_min, axis_ratio_min]
    upper_bounds = [core_temperature_max,
                    skin_temperature_max, skin_temperature_max,
                    skin_temperature_max, skin_temperature_max,
                    log_heat_flow_max,
                    flesh_conductivity_max, panting_rate_max, skin_wetness_max,
                    insulation_depth_max, axis_ratio_max]

    pp = nlp_packed.p
    heat_flow_init = max(ustrip(u"W", metabolic_heat_flow_init), heat_flow_min_W)
    initial_effectors = clamp.(
        [setpoint_temperature_K,
         ustrip(u"K", pp.initial_dorsal_skin_temperature),
         ustrip(u"K", pp.initial_dorsal_insulation_temperature),
         ustrip(u"K", pp.initial_ventral_skin_temperature),
         ustrip(u"K", pp.initial_ventral_insulation_temperature),
         log(heat_flow_init),
         ustrip(u"W/m/K", int_cond.flesh_conductivity),
         1.0,
         Float64(evap_pars.skin_wetness),
         insulation_depth_max,   # start with erected insulation
         axis_ratio_min],      # start curled
        lower_bounds, upper_bounds,
    )

    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
    opt_pars = (;
        setpoint_temperature_K,
        heat_flow_min_W,
        core_temperature_penalty = limits.core_temperature_penalty,
        metabolic_heat_penalty   = limits.metabolic_heat_penalty,
        panting_penalty          = limits.panting_penalty,
        skin_wetness_penalty     = limits.skin_wetness_penalty,
        gradient_penalty         = limits.gradient_penalty,
        target_gradient          = limits.target_core_skin_gradient,
        core_temperature_range,
        gradient_range           = core_temperature_range,
        panting_rate_range       = max(panting_rate_max - 1.0, 1e-6),
        skin_wetness_min,
        skin_wetness_range       = max(skin_wetness_max - skin_wetness_min, 1e-6),
        heat_flow_range          = max(heat_flow_max_W - heat_flow_min_W, 1.0),
    )

    nlp_pars = (;
        nlp_packed,
        minimum_heat_flow    = heat_flow_min_W,
        q10                  = Float64(metab_pars.q10),
        setpoint_temperature = setpoint_temperature_K,
    )

    obj_fn(x, _)     = _objective_value_multisided(x, opt_pars)
    res_fn!(r, x, _) = _heat_balance_residuals_multisided!(r, x, nlp_pars)

    function grad_fn!(g, x, _)
        fill!(g, 0)
        Enzyme.autodiff(Enzyme.Reverse, _objective_value_multisided,
                        Enzyme.Active,
                        Enzyme.Duplicated(x, g),
                        Enzyme.Const(opt_pars))
        return nothing
    end
    function cons_j_fn!(J, x, _)
        m = size(J, 1)
        n = length(x)
        r  = zeros(m)
        dr = zeros(m)
        dx = zeros(n)
        for j in 1:n
            fill!(r, 0); fill!(dr, 0); fill!(dx, 0); dx[j] = 1.0
            Enzyme.autodiff(Enzyme.Forward, _heat_balance_residuals_multisided!,
                            Enzyme.Const,
                            Enzyme.Duplicated(r, dr),
                            Enzyme.Duplicated(x, dx),
                            Enzyme.Const(nlp_pars))
            @views J[:, j] .= dr
        end
        return nothing
    end
    hess_fn!(H, x, _)     = (fill!(H, 0); nothing)
    cons_h_fn!(res, x, _) = (foreach(r -> fill!(r, 0), res); nothing)

    optimization_func = OptimizationFunction(obj_fn, SciMLBase.NoAD();
        cons     = res_fn!,
        grad     = grad_fn!,
        hess     = hess_fn!,
        cons_j   = cons_j_fn!,
        cons_h   = cons_h_fn!,
    )
    optimization_prob = OptimizationProblem(optimization_func, initial_effectors, nothing;
        lb     = lower_bounds,
        ub     = upper_bounds,
        lcons  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ucons  = [0.0, 0.0, 0.0, 0.0, 0.0, Inf],
    )
    ipopt_sol = solve(optimization_prob,
        IpoptOptimizer(;
            hessian_approximation = "limited-memory",
            acceptable_tol        = 1e-3,
            acceptable_iter       = 5,
        );
        verbose,
        reltol   = 1e-4,
        maxiters = 300,
    )

    x_sol = ipopt_sol.u
    core_temperature              = x_sol[1]  * u"K"
    dorsal_skin_temperature       = x_sol[2]  * u"K"
    dorsal_insulation_temperature = x_sol[3]  * u"K"
    ventral_skin_temperature      = x_sol[4]  * u"K"
    ventral_insulation_temperature = x_sol[5] * u"K"
    metabolic_heat_flow           = exp(x_sol[6]) * u"W"
    flesh_conductivity            = x_sol[7]  * u"W/m/K"
    panting_rate                  = x_sol[8]
    skin_wetness                  = x_sol[9]
    insulation_depth              = x_sol[10] * u"m"
    axis_ratio_b                  = x_sol[11]

    return HeatExchange.nlp_assemble_output(nlp_pars.nlp_packed, organism, environment,
        core_temperature,
        dorsal_skin_temperature, dorsal_insulation_temperature,
        ventral_skin_temperature, ventral_insulation_temperature,
        metabolic_heat_flow, flesh_conductivity, panting_rate, skin_wetness,
        insulation_depth, axis_ratio_b)
end

# =============================================================================
# Main thermoregulate dispatch
# =============================================================================

"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment,
                   metabolic_heat_flow_init, skin_temperature_init,
                   insulation_temperature_init; verbose=false)

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
