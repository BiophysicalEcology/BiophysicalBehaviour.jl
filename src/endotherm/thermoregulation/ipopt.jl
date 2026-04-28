# =============================================================================
# IPOPT-based endotherm thermoregulation control.
#
# Solver policy lives here: objective, variable bounds, initial effectors, and
# the Q10 metabolic-scaling inequality constraint.
# All physics (packing, per-iteration residuals, output assembly) is delegated
# to HeatExchange.nlp_pack / nlp_residuals / nlp_assemble_output.
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

# Pure Float64 objective. `opt` carries only stripped scalars and dimensionless
# weights/ranges. metabolic_heat_penalty is a regularisation term that breaks
# degeneracy in cold conditions; the Q10 inequality constraint overrides it in
# hot conditions. gradient_penalty (default 0) penalises deviation from the
# target core–skin difference. skin_wetness_penalty > panting_penalty makes
# panting activate before sweating.
function _objective_value(effectors, opt)
    opt.core_temperature_penalty * ((effectors[1] - opt.setpoint_temperature_K) / opt.core_temperature_range)^2 +
    opt.metabolic_heat_penalty   * ((exp(effectors[4]) - opt.heat_flow_min_W)   / opt.heat_flow_range)^2        +
    opt.gradient_penalty         * ((effectors[1] - effectors[2] - opt.target_gradient) / opt.gradient_range)^2 +
    opt.panting_penalty          * ((effectors[6] - 1.0) / opt.panting_rate_range)^2                            +
    opt.skin_wetness_penalty     * ((effectors[7] - opt.skin_wetness_min) / opt.skin_wetness_range)^2
end

# Float64 effectors in → Unitful physics (via HeatExchange.nlp_residuals) → Float64 residuals out.
# residuals[1:3] = 0 (energy balance, internal conduction, skin temp).
# residuals[4]  >= 0 (Q10: metabolic_heat_flow >= Q_gen_min · q10^((T_core − T_setpoint)/10)).
function _heat_balance_residuals!(residuals, effectors, p)
    core_temperature       = effectors[1] * u"K"
    skin_temperature       = effectors[2] * u"K"
    insulation_temperature = effectors[3] * u"K"
    metabolic_heat_flow    = exp(effectors[4]) * u"W"   # effectors[4] = log(metabolic_heat_flow)
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
    residuals[4] = exp(effectors[4]) - p.Q_gen_min * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
    return nothing
end

"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment,
                   metabolic_heat_flow_init, skin_temperature_init,
                   insulation_temperature_init; verbose=false)

Solve endotherm heat balance as a nonlinear program via IPOPT.

Uses nine NLP decision variables: six physiological effectors (metabolic heat
generation, flesh conductivity, panting rate, skin wetness, insulation depth,
body shape) plus three temperature state variables (core, skin, insulation surface).
Three heat-balance equality constraints from `HeatExchange.nlp_residuals` enforce
physical consistency so that the temperatures are outcomes of the effectors, not
independent choices. A fourth inequality constraint enforces Q10 metabolic scaling:
`metabolic_heat_flow >= heat_flow_min * Q10^((core_temperature − setpoint)/10)`.

`metabolic_heat_flow` is not penalised in the objective — it is a state variable driven by the
energy balance. The `core_temperature` term implicitly sets `metabolic_heat_flow`: in the cold it
rises freely to close the deficit; in the heat it stays near `heat_flow_min` because extra heat
would push `core_temperature` above setpoint. All penalty terms are normalised to [0,1].
Penalty weights in `ThermoregulationLimits`:
  - `core_temperature_penalty`  — lower → core_temperature rises sooner before effectors are exhausted
  - `metabolic_heat_penalty`   — small regularisation (default 0.1) preventing high-panting/high-metabolic_heat_flow
                                  degeneracy at cold temperatures; overridden at hot temperatures by
                                  the Q10 inequality constraint which forces metabolic_heat_flow to rise
  - `gradient_penalty`         — non-zero → penalise deviation from `target_core_skin_gradient`;
                                  activates vasodilation/evaporation before core_temperature moves (default 0)
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
    ::IPOPTControl,
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
    T_air      = environment.environment_vars.air_temperature
    T_setpoint = metab_pars.core_temperature

    # ---- Pre-solve physics packing ------------------------------------------
    nlp_packed = HeatExchange.nlp_pack(WeightedMeanNLP(), organism, environment,
                                       skin_temperature_init, insulation_temperature_init)

    # ---- Boundary: strip into pure Float64 optimizer-world parameters --------
    # x = [core_temperature, skin_temperature, insulation_temperature,    ← temperature states
    #       log(metabolic_heat_flow), flesh_conductivity, panting_rate,   ← effectors
    #       skin_wetness, insulation_depth, axis_ratio_b]                  ← effectors
    # The three equality constraints make temperatures outcomes of the six effectors.
    air_temperature_K      = ustrip(u"K", T_air)
    setpoint_temperature_K = ustrip(u"K", T_setpoint)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - 5.0
    skin_temperature_max   = core_temperature_max + 5.0
    heat_flow_min_W        = ustrip(u"W", limits.Q_minimum_ref)
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
    aspect_ratio_min       = Float64(limits.aspect_ratio_factor.reference)
    aspect_ratio_max       = organism.body.shape isa Sphere ? aspect_ratio_min : Float64(limits.aspect_ratio_factor.max)

    lower_bounds = [core_temperature_min, skin_temperature_min, skin_temperature_min, log_heat_flow_min,
                    flesh_conductivity_min, panting_rate_min, skin_wetness_min,
                    insulation_depth_min, aspect_ratio_min]
    upper_bounds = [core_temperature_max, skin_temperature_max, skin_temperature_max, log_heat_flow_max,
                    flesh_conductivity_max, panting_rate_max, skin_wetness_max,
                    insulation_depth_max, aspect_ratio_max]

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
         aspect_ratio_min],      # start curled
        lower_bounds, upper_bounds,
    )

    # gradient_range reuses core_temperature_range: max plausible deviation of
    # (core_temperature − skin_temperature) from target equals the max core_temperature
    # excursion above setpoint.
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

    # nlp_packed carries physics; Q_gen_min/q10/setpoint_temperature are Q10 solver policy.
    nlp_pars = (;
        nlp_packed,
        Q_gen_min            = heat_flow_min_W,
        q10                  = Float64(metab_pars.q10),
        setpoint_temperature = setpoint_temperature_K,
    )

    # ---- Closures bridging IPOPT to the boundary helpers --------------------
    # Each closure captures only the world it needs: `_objective_value` sees
    # `opt_pars`, `_heat_balance_residuals!` sees `nlp_pars`. `nothing` is
    # passed as Optimization's `p` since the closures already carry the params.
    obj_fn(x, _) = _objective_value(x, opt_pars)
    res_fn!(r, x, _) = _heat_balance_residuals!(r, x, nlp_pars)

    # Bypass DifferentiationInterface (incompatible with Unitful params NamedTuple).
    # hess! and cons_h! are registered (required by IpoptOptimizer) but not called
    # at runtime when hessian_approximation="limited-memory" is set — L-BFGS builds
    # the Hessian approximation from gradient differences, eliminating O(n²) evaluations.
    grad_fn!(g, x, _) = FiniteDiff.finite_difference_gradient!(g, e -> _objective_value(e, opt_pars), x)
    hess_fn!(H, x, _) = (H .= FiniteDiff.finite_difference_hessian(e -> _objective_value(e, opt_pars), x))
    cons_j_fn!(J, x, _) = (J .= FiniteDiff.finite_difference_jacobian(
        e -> (r = zeros(eltype(e), 4); _heat_balance_residuals!(r, e, nlp_pars); r),
        x,
    ))
    function cons_h_fn!(res, x, _)
        for i in eachindex(res)
            res[i] .= FiniteDiff.finite_difference_hessian(
                e -> (r = zeros(4); _heat_balance_residuals!(r, e, nlp_pars); r[i]),
                x,
            )
        end
    end

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
        ucons  = [0.0, 0.0, 0.0, Inf],  # residuals[4] >= 0: metabolic_heat_flow >= Q10-scaled minimum
    )
    # hessian_approximation and tolerance options go to the IpoptOptimizer struct (not solve kwargs).
    # reltol and maxiters are common interface args forwarded via solve.
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

    # ---- Boundary: attach units to solution, delegate output assembly --------
    x_sol = ipopt_sol.u
    core_temperature       = x_sol[1] * u"K"
    skin_temperature       = x_sol[2] * u"K"
    insulation_temperature = x_sol[3] * u"K"
    metabolic_heat_flow    = exp(x_sol[4]) * u"W"   # x_sol[4] = log(metabolic_heat_flow)
    flesh_conductivity     = x_sol[5] * u"W/m/K"
    panting_rate           = x_sol[6]
    skin_wetness           = x_sol[7]
    insulation_depth       = x_sol[8] * u"m"
    axis_ratio_b           = x_sol[9]

    return HeatExchange.nlp_assemble_output(nlp_pars.nlp_packed, organism, environment,
        core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow,
        flesh_conductivity, panting_rate, skin_wetness, insulation_depth, axis_ratio_b)
end
