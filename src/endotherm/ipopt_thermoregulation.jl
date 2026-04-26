"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, generated_heat_flow_init, skin_temperature_init, insulation_temperature_init)

Solve endotherm heat balance as a nonlinear program via IPOPT.

Uses nine NLP decision variables: six physiological effectors (metabolic heat
generation, flesh conductivity, panting rate, skin wetness, insulation depth,
body shape) plus three temperature state variables (core, skin, insulation surface).
Three heat-balance equality constraints from `HeatExchange.heat_balance` enforce
physical consistency so that the temperatures are outcomes of the effectors, not
independent choices. A fourth inequality constraint enforces Q10 metabolic scaling:
`generated_heat_flow >= heat_flow_min * Q10^((core_temperature − setpoint)/10)`.

`generated_heat_flow` is not penalised in the objective — it is a state variable driven by the
energy balance. The `core_temperature` term implicitly sets `generated_heat_flow`: in the cold it
rises freely to close the deficit; in the heat it stays near `heat_flow_min` because extra heat
would push `core_temperature` above setpoint. All penalty terms are normalised to [0,1].
Penalty weights in `ThermoregulationLimits`:
  - `core_temperature_penalty`  — lower → core_temperature rises sooner before effectors are exhausted
  - `metabolic_heat_penalty`   — small regularisation (default 0.1) preventing high-panting/high-generated_heat_flow
                                  degeneracy at cold temperatures; overridden at hot temperatures by
                                  the Q10 inequality constraint which forces generated_heat_flow to rise
  - `gradient_penalty`         — non-zero → penalise deviation from `target_core_skin_gradient`;
                                  activates vasodilation/evaporation before core_temperature moves (default 0)
  - `panting_penalty`          — lower → panting activates sooner
  - `skin_wetness_penalty`     — higher than `panting_penalty` → panting activates before sweating

Uses a dorsal/ventral mean-weighted body where dorsal_weight = sky_view_factor +
vegetation_view_factor. This matches the weighting used by `solve_metabolic_rate`.

# Arguments
- `organism`: must have `OrganismTraits` with `ThermoregulationLimits`
- `environment`: NamedTuple with `environment_pars` and `environment_vars`
- `generated_heat_flow_init`: initial guess for metabolic heat generation (W)
- `skin_temperature_init`: initial guess for skin temperature (K)
- `insulation_temperature_init`: initial guess for insulation surface temperature (K)
"""
function thermoregulate(
    ::Endotherm,
    ::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow_init,
    skin_temperature_init,
    insulation_temperature_init;
    verbose = false,
)
    env_pars  = stripparams(environment.environment_pars)
    env_vars  = environment.environment_vars

    ins_pars   = insulation_pars(organism)
    ext_cond   = conduction_pars_external(organism)
    int_cond   = conduction_pars_internal(organism)
    rad_pars   = radiation_pars(organism)
    evap_pars  = evaporation_pars(organism)
    resp_pars  = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)
    limits     = thermoregulation(organism)

    T_air       = env_vars.air_temperature
    T_setpoint  = metab_pars.core_temperature
    T_vegetation = env_vars.reference_air_temperature

    # ---- View factor geometry (mirrors solve_metabolic_rate) -----------------
    vegetation_view_factor = rad_pars.sky_view_factor * env_vars.shade
    sky_view_factor        = rad_pars.sky_view_factor - vegetation_view_factor
    ground_view_factor     = 1 - sky_view_factor - vegetation_view_factor
    bush_view_factor       = rad_pars.bush_view_factor
    dorsal_weight          = sky_view_factor + vegetation_view_factor
    ventral_weight         = 1 - dorsal_weight
    ventral_fraction       = rad_pars.ventral_fraction

    # ---- Mean-weighted insulation geometry -----------------------------------
    fat = Fat(int_cond.fat_fraction, int_cond.fat_density)

    mean_insulation_depth       = ins_pars.dorsal.depth       * dorsal_weight + ins_pars.ventral.depth       * ventral_weight
    mean_fibre_diameter         = ins_pars.dorsal.diameter    * dorsal_weight + ins_pars.ventral.diameter    * ventral_weight
    mean_fibre_density          = ins_pars.dorsal.density     * dorsal_weight + ins_pars.ventral.density     * ventral_weight
    mean_fibre_length           = ins_pars.dorsal.length      * dorsal_weight + ins_pars.ventral.length      * ventral_weight
    mean_fibre_reflectance      = ins_pars.dorsal.reflectance * dorsal_weight + ins_pars.ventral.reflectance * ventral_weight
    mean_fibre_conductivity     = ins_pars.dorsal.conductivity * dorsal_weight + ins_pars.ventral.conductivity * ventral_weight

    mean_fibre_props = FibreProperties(;
        diameter     = mean_fibre_diameter,
        length       = mean_fibre_length,
        density      = mean_fibre_density,
        depth        = mean_insulation_depth,
        reflectance  = mean_fibre_reflectance,
        conductivity = mean_fibre_conductivity,
    )
    mean_ins_pars = InsulationParameters(;
        dorsal                  = mean_fibre_props,
        ventral                 = mean_fibre_props,
        depth_compressed        = ins_pars.depth_compressed,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )
    mean_fur  = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
    mean_body = Body(organism.body.shape, CompositeInsulation(mean_fur, fat))

    # ---- Solar heat flows ----------------------------------------------------
    reference_total_area      = BiophysicalGeometry.total_area(organism.body)
    reference_conduction_area = reference_total_area * ext_cond.conduction_fraction
    absorptivities            = Absorptivities(rad_pars, env_pars)
    solar_view_factors        = ViewFactors(sky_view_factor, ground_view_factor, 0.0, 0.0)
    solar_conds               = SolarConditions(env_vars)
    reference_silhouette_area = silhouette_area(organism.body, rad_pars.solar_orientation)
    solar_result              = solar(organism.body, absorptivities, solar_view_factors,
                                      solar_conds, reference_silhouette_area, reference_conduction_area)

    if solar_result.solar_flow > 0.0u"W"
        dorsal_solar_flow  = 2.0 * solar_result.direct_flow + solar_result.solar_sky_flow * 2.0
        ventral_solar_flow = (solar_result.solar_substrate_flow /
                              (1.0 - sky_view_factor - vegetation_view_factor)) *
                             (1.0 - 2.0 * ext_cond.conduction_fraction)
    else
        dorsal_solar_flow  = 0.0u"W"
        ventral_solar_flow = 0.0u"W"
    end
    mean_solar_flow = dorsal_solar_flow * dorsal_weight + ventral_solar_flow * ventral_weight

    # ---- Mean view factors: dorsal sees sky+veg, ventral sees ground+bush ----
    mean_view_factors = ViewFactors(
        sky_view_factor        * 2.0 * dorsal_weight,
        ground_view_factor     * 2.0 * ventral_weight,
        bush_view_factor       * 2.0 * ventral_weight,
        vegetation_view_factor * 2.0,
    )

    # ---- Mean substrate conductance (ventral side only) ----------------------
    mean_body_total_area      = BiophysicalGeometry.total_area(mean_body)
    ventral_conduction_area   = mean_body_total_area * ext_cond.conduction_fraction * 2
    ventral_conduction_coeff  = (ventral_conduction_area * env_vars.substrate_conductivity) /
                                env_pars.conduction_depth
    mean_conduction_coeff     = ventral_conduction_coeff * ventral_weight

    # ---- Pack environment and traits for heat_balance ------------------------
    mean_insulation_temperature_init = insulation_temperature_init * 0.7 + skin_temperature_init * 0.3
    insulation_props_init = insulation_properties(mean_ins_pars, mean_insulation_temperature_init, ventral_fraction)

    env_temperatures = EnvironmentTemperatures(
        T_air,
        env_vars.sky_temperature,
        env_vars.ground_temperature,
        T_vegetation,
        env_vars.bush_temperature,
        env_vars.substrate_temperature,
    )
    atmos = AtmosphericConditions(env_vars)
    heat_balance_env = (;
        temperature            = env_temperatures,
        view_factors           = mean_view_factors,
        atmos,
        fluid                  = env_pars.fluid,
        solar_flow             = mean_solar_flow,
        gas_fractions          = env_pars.gas_fractions,
        convection_enhancement = env_pars.convection_enhancement,
    )
    mean_body_emissivity = rad_pars.body_emissivity_dorsal * dorsal_weight +
                           rad_pars.body_emissivity_ventral * ventral_weight
    heat_balance_traits = (;
        fat_conductivity   = int_cond.fat_conductivity,
        flesh_conductivity = int_cond.flesh_conductivity,
        ϵ_body             = mean_body_emissivity,
        skin_wetness       = evap_pars.skin_wetness,
        insulation_wetness = evap_pars.insulation_wetness,
        bare_skin_fraction = evap_pars.bare_skin_fraction,
        eye_fraction       = evap_pars.eye_fraction,
    )
    geometry_vars = GeometryVariables(;
        side                    = :dorsal,
        conductance_coefficient = mean_conduction_coeff,
        ventral_fraction,
        conduction_fraction     = ext_cond.conduction_fraction,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )

    # ---- Decision variable bounds (SI, unitless Float64) --------------------
    # x = [core_temperature, skin_temperature, insulation_temperature,    ← temperature states
    #       log(generated_heat_flow), flesh_conductivity, panting_rate,   ← effectors
    #       skin_wetness, insulation_depth, aspect_ratio_factor]           ← effectors
    # The three equality constraints make temperatures outcomes of the six effectors.
    air_temperature_K      = ustrip(u"K", T_air)
    setpoint_temperature_K = ustrip(u"K", T_setpoint)
    core_temperature_min   = ustrip(u"K", limits.core_temperature.reference)
    core_temperature_max   = ustrip(u"K", limits.core_temperature.max)
    skin_temperature_min   = air_temperature_K - 5.0
    skin_temperature_max   = core_temperature_max + 5.0
    heat_flow_min          = ustrip(u"W", limits.Q_minimum_ref)
    heat_flow_max          = ustrip(u"W", limits.Q_minimum_ref) * 20.0
    log_heat_flow_min      = log(heat_flow_min)
    log_heat_flow_max      = log(heat_flow_max)
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

    heat_flow_init = max(ustrip(u"W", generated_heat_flow_init), heat_flow_min)
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

    # ---- Normalisation ranges for objective penalties -----------------------
    heat_flow_range        = max(heat_flow_max - heat_flow_min, 1.0)
    core_temperature_range = max(core_temperature_max - setpoint_temperature_K, 1e-6)
    panting_rate_range     = max(panting_rate_max - 1.0, 1e-6)
    skin_wetness_range     = max(skin_wetness_max - skin_wetness_min, 1e-6)
    # gradient_range reuses core_temperature_range: max plausible deviation of
    # (core_temperature − skin_temperature) from target equals the max core_temperature
    # excursion above setpoint.

    nlp_pars = (;
        setpoint_temperature     = setpoint_temperature_K,
        core_temperature_penalty = limits.core_temperature_penalty,
        metabolic_heat_penalty   = limits.metabolic_heat_penalty,
        panting_penalty          = limits.panting_penalty,
        skin_wetness_penalty     = limits.skin_wetness_penalty,
        gradient_penalty         = limits.gradient_penalty,
        target_gradient          = limits.target_core_skin_gradient,
        core_temperature_range,
        gradient_range           = core_temperature_range,
        panting_rate_range,
        skin_wetness_min,
        skin_wetness_range,
        heat_flow_min,             # used by Q10 constraint and metabolic regularisation
        heat_flow_range,
        q10                      = metab_pars.q10,
        mean_body,
        mean_ins_pars,
        mean_fibre_props,
        fat,
        body_shape     = organism.body.shape,
        body_is_sphere = organism.body.shape isa Sphere,
        insulation_props_init,
        geometry_vars,
        heat_balance_env,
        heat_balance_traits,
        resp_pars,
        ventral_fraction,
        substrate_conductivity = env_vars.substrate_conductivity,
        conduction_depth       = env_pars.conduction_depth,
        ventral_weight,
        ext_cond,
    )

    # ---- Objective -----------------------------------------------------------
    # metabolic_heat_penalty is a regularisation term (default 0.1) that breaks
    # degeneracy in cold conditions: without it, the optimizer can freely combine high
    # generated_heat_flow + high panting and satisfy the energy balance equally well.
    # The Q10 inequality constraint (4th residual) overrides it in hot conditions.
    # gradient_penalty (default 0) adds an optional term penalising deviation from the
    # target core–skin temperature difference, activating vasodilation/evaporation
    # before absolute core_temperature deviation becomes the primary signal.
    # skin_wetness_penalty > panting_penalty → panting activates before sweating.
    function objective(effectors, p)
        p.core_temperature_penalty * ((effectors[1] - p.setpoint_temperature) / p.core_temperature_range)^2 +
        p.metabolic_heat_penalty   * ((exp(effectors[4]) - p.heat_flow_min) / p.heat_flow_range)^2         +
        p.gradient_penalty         * ((effectors[1] - effectors[2] - p.target_gradient) / p.gradient_range)^2 +
        p.panting_penalty          * ((effectors[6] - 1.0) / p.panting_rate_range)^2                       +
        p.skin_wetness_penalty     * ((effectors[7] - p.skin_wetness_min) / p.skin_wetness_range)^2
    end

    # ---- Constraints: three equality + one Q10 inequality residual ----------
    # residuals[1:3] = 0  (equality: energy balance, internal conduction, skin temp)
    # residuals[4]  >= 0  (inequality: generated_heat_flow >= heat_flow_min * Q10^((core_temperature − setpoint)/10))
    function heat_balance_residuals!(residuals, effectors, p)
        core_temperature       = effectors[1] * u"K"
        skin_temperature       = effectors[2] * u"K"
        insulation_temperature = effectors[3] * u"K"
        generated_heat_flow    = exp(effectors[4]) * u"W"   # effectors[4] = log(generated_heat_flow)
        flesh_conductivity     = effectors[5] * u"W/m/K"
        panting_rate           = effectors[6]
        skin_wetness           = effectors[7]
        insulation_depth       = effectors[8] * u"m"
        aspect_ratio_factor    = effectors[9]

        # Rebuild insulation and body from piloerection and aspect_ratio_factor effectors
        trial_fibre_props  = setproperties(p.mean_fibre_props; depth = insulation_depth)
        trial_ins_pars     = setproperties(p.mean_ins_pars; dorsal = trial_fibre_props, ventral = trial_fibre_props)
        trial_fur          = Fur(insulation_depth, p.mean_fibre_props.diameter, p.mean_fibre_props.density)
        trial_shape        = p.body_is_sphere ? p.body_shape : setproperties(p.body_shape; aspect_ratio_b = aspect_ratio_factor)
        trial_body         = Body(trial_shape, CompositeInsulation(trial_fur, p.fat))

        # Recompute temperature-dependent insulation properties
        mean_insulation_temperature = insulation_temperature * 0.7 + skin_temperature * 0.3
        trial_insulation_props      = insulation_properties(trial_ins_pars, mean_insulation_temperature, p.ventral_fraction)

        # Update conductance for new body size (conduction area changes with shape)
        trial_conduction_area  = BiophysicalGeometry.total_area(trial_body) * p.ext_cond.conduction_fraction * 2
        trial_conduction_coeff = (trial_conduction_area * p.substrate_conductivity) / p.conduction_depth * p.ventral_weight
        trial_geometry_vars    = setproperties(p.geometry_vars; conductance_coefficient = trial_conduction_coeff)

        balance = HeatExchange.heat_balance(
            core_temperature, skin_temperature, insulation_temperature, generated_heat_flow;
            body             = trial_body,
            insulation_pars  = trial_ins_pars,
            insulation       = trial_insulation_props,
            geometry_vars    = trial_geometry_vars,
            environment_vars = p.heat_balance_env,
            traits           = p.heat_balance_traits,
            resp_pars        = p.resp_pars,
            k_flesh          = flesh_conductivity,
            pant             = panting_rate,
            skin_wetness,
        )
        residuals[1] = ustrip(u"W", balance.residual_energy_balance)
        residuals[2] = ustrip(u"W", balance.residual_internal_conduction)
        residuals[3] = ustrip(u"K", balance.residual_skin_temperature)
        q10_minimum  = p.heat_flow_min * p.q10 ^ ((effectors[1] - p.setpoint_temperature) / 10.0)
        residuals[4] = exp(effectors[4]) - q10_minimum   # >= 0 enforced via ucons[4] = Inf
        return nothing
    end

    # Bypass DifferentiationInterface (incompatible with Unitful params NamedTuple).
    # hess! and cons_h! are registered (required by IpoptOptimizer) but not called
    # at runtime when hessian_approximation="limited-memory" is set — L-BFGS builds
    # the Hessian approximation from gradient differences, eliminating O(n²) evaluations.
    function objective_gradient!(grad, effectors, p)
        FiniteDiff.finite_difference_gradient!(grad, e -> objective(e, p), effectors)
    end
    function objective_hessian!(H, effectors, p)
        H .= FiniteDiff.finite_difference_hessian(e -> objective(e, p), effectors)
    end
    function constraint_jacobian!(J, effectors, p)
        J .= FiniteDiff.finite_difference_jacobian(
            e -> (res = zeros(eltype(e), 4); heat_balance_residuals!(res, e, p); res),
            effectors,
        )
    end
    function constraint_hessians!(res, effectors, p)
        for i in eachindex(res)
            res[i] .= FiniteDiff.finite_difference_hessian(
                e -> begin
                    r = zeros(4)
                    heat_balance_residuals!(r, e, p)
                    r[i]
                end,
                effectors,
            )
        end
    end

    optimization_func = OptimizationFunction(objective, SciMLBase.NoAD();
        cons     = heat_balance_residuals!,
        grad     = objective_gradient!,
        hess     = objective_hessian!,
        cons_j   = constraint_jacobian!,
        cons_h   = constraint_hessians!,
    )
    optimization_prob = OptimizationProblem(optimization_func, initial_effectors, nlp_pars;
        lb     = lower_bounds,
        ub     = upper_bounds,
        lcons  = [0.0, 0.0, 0.0, 0.0],
        ucons  = [0.0, 0.0, 0.0, Inf],  # residuals[4] >= 0: generated_heat_flow >= Q10-scaled minimum
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

    # ---- Reconstruct full solution from solved effector values ---------------
    x_sol = ipopt_sol.u
    core_temperature       = x_sol[1] * u"K"
    skin_temperature       = x_sol[2] * u"K"
    insulation_temperature = x_sol[3] * u"K"
    generated_heat_flow    = exp(x_sol[4]) * u"W"   # x_sol[4] = log(generated_heat_flow)
    flesh_conductivity     = x_sol[5] * u"W/m/K"
    panting_rate           = x_sol[6]
    skin_wetness           = x_sol[7]
    insulation_depth       = x_sol[8] * u"m"
    aspect_ratio_factor    = x_sol[9]

    sol_fibre_props    = setproperties(mean_fibre_props; depth = insulation_depth)
    sol_ins_pars       = setproperties(mean_ins_pars; dorsal = sol_fibre_props, ventral = sol_fibre_props)
    sol_fur            = Fur(insulation_depth, mean_fibre_props.diameter, mean_fibre_props.density)
    sol_shape          = organism.body.shape isa Sphere ?
                         organism.body.shape :
                         setproperties(organism.body.shape; aspect_ratio_b = aspect_ratio_factor)
    sol_body           = Body(sol_shape, CompositeInsulation(sol_fur, fat))
    sol_mean_insulation_temperature = insulation_temperature * 0.7 + skin_temperature * 0.3
    sol_insulation_props = insulation_properties(sol_ins_pars, sol_mean_insulation_temperature, ventral_fraction)
    sol_conduction_area  = BiophysicalGeometry.total_area(sol_body) * ext_cond.conduction_fraction * 2
    sol_conduction_coeff = (sol_conduction_area * env_vars.substrate_conductivity) /
                           env_pars.conduction_depth * ventral_weight
    sol_geometry_vars    = setproperties(geometry_vars; conductance_coefficient = sol_conduction_coeff)

    heat_balance_result = HeatExchange.heat_balance(
        core_temperature, skin_temperature, insulation_temperature, generated_heat_flow;
        body             = sol_body,
        insulation_pars  = sol_ins_pars,
        insulation       = sol_insulation_props,
        geometry_vars    = sol_geometry_vars,
        environment_vars = heat_balance_env,
        traits           = heat_balance_traits,
        resp_pars,
        k_flesh          = flesh_conductivity,
        pant             = panting_rate,
        skin_wetness,
    )

    lung_temperature = (core_temperature + skin_temperature) / 2

    # ---- Longwave flows (mirrors solve_metabolic_rate simplified approach) ---
    σ                        = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    sol_total_area           = BiophysicalGeometry.total_area(sol_body)
    sol_ventral_radiant_area = sol_total_area * (1 - ext_cond.conduction_fraction)
    dorsal_longwave_emission  = 2 * sky_view_factor * σ * rad_pars.body_emissivity_dorsal *
                                sol_total_area * insulation_temperature^4
    ventral_longwave_emission = 2 * ground_view_factor * σ * rad_pars.body_emissivity_ventral *
                                sol_ventral_radiant_area * insulation_temperature^4
    longwave_flow_out = dorsal_longwave_emission * dorsal_weight + ventral_longwave_emission * ventral_weight
    longwave_flow_in  = longwave_flow_out - heat_balance_result.radiation_heat_flow

    # ---- Respiration mass flows at solution ----------------------------------
    panting_resp_pars    = setproperties(resp_pars; pant = panting_rate)
    respiration_result   = HeatExchange.respiration(
        MetabolicRates(; metabolic = generated_heat_flow, sum = generated_heat_flow, minimum = 0.0u"W"),
        panting_resp_pars,
        heat_balance_env.atmos,
        organism.body.shape.mass,
        lung_temperature,
        T_air;
        gas_fractions = env_pars.gas_fractions,
        O2conversion  = Kleiber1961(),
    )
    latent_heat_vaporisation = FluidProperties.enthalpy_of_vaporisation(T_air)
    m_sweat = u"g/hr"(heat_balance_result.skin_evaporation_heat_flow / latent_heat_vaporisation)
    m_evap  = u"g/hr"(respiration_result.respiration_mass + m_sweat)

    # ---- Assemble outputs matching solve_metabolic_rate structure ------------
    thermoregulation_out = (;
        core_temperature,
        skin_temperature,
        insulation_temperature,
        lung_temperature,
        skin_temperature_dorsal        = skin_temperature,
        skin_temperature_ventral       = skin_temperature,
        insulation_temperature_dorsal  = insulation_temperature,
        insulation_temperature_ventral = insulation_temperature,
        aspect_ratio_factor,
        pant                           = panting_rate,
        skin_wetness,
        flesh_conductivity,
        insulation_conductivity_effective  = sol_insulation_props.conductivities.average,
        insulation_conductivity_dorsal     = sol_insulation_props.conductivities.dorsal,
        insulation_conductivity_ventral    = sol_insulation_props.conductivities.ventral,
        insulation_conductivity_compressed = sol_insulation_props.conductivity_compressed,
        insulation_depth_dorsal  = sol_ins_pars.dorsal.depth,
        insulation_depth_ventral = sol_ins_pars.ventral.depth,
        metab_pars.q10,
    )

    sol_total_area_out   = BiophysicalGeometry.total_area(sol_body)
    area_skin_out        = skin_area(sol_body)
    area_evaporation_out = evaporation_area(sol_body)
    area_convection_out  = sol_total_area_out * (1 - ext_cond.conduction_fraction)
    area_silhouette_out  = silhouette_area(sol_body, rad_pars.solar_orientation)

    morphology = (;
        total_area               = sol_total_area_out,
        area_skin                = area_skin_out,
        area_evaporation         = area_evaporation_out,
        area_convection          = area_convection_out,
        area_conduction          = reference_conduction_area / 2,
        area_silhouette          = area_silhouette_out,
        sky_view_factor,
        ground_view_factor,
        volume                   = sol_body.geometry.volume,
        volume_flesh             = flesh_volume(sol_body),
        characteristic_dimension = HeatExchange.characteristic_dimension(HeatExchange.VolumeCubeRoot(), sol_body),
        fat_mass                 = sol_body.shape.mass * fat.fraction,
        sol_body.geometry.length...,
    )

    energy_flows = (;
        solar_flow            = heat_balance_result.solar_heat_flow,
        longwave_flow_in,
        generated_heat_flow,
        evaporation_heat_flow = heat_balance_result.skin_evaporation_heat_flow +
                                heat_balance_result.insulation_evaporation_heat_flow +
                                heat_balance_result.respiration_heat_flow,
        longwave_flow_out,
        convection_heat_flow  = heat_balance_result.convection_heat_flow,
        conduction_flow       = heat_balance_result.conduction_heat_flow,
        balance               = heat_balance_result.residual_energy_balance,
        ntry                  = 1,
        success               = true,
    )

    mass_flows = (;
        air_flow             = respiration_result.air_flow,
        oxygen_flow_standard = respiration_result.oxygen_flow_standard,
        m_evap,
        respiration_mass     = respiration_result.respiration_mass,
        m_sweat,
        molar_fluxes_in      = respiration_result.molar_fluxes_in,
        molar_fluxes_out     = respiration_result.molar_fluxes_out,
    )

    return (; thermoregulation = thermoregulation_out, morphology, energy_flows, mass_flows)
end
