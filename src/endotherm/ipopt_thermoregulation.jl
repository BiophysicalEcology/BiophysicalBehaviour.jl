"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, Q_gen, T_skin, T_insulation)

Solve endotherm heat balance as a nonlinear program via IPOPT.

Treats T_core, T_skin, T_ins, Q_gen, k_flesh, pant, skin_wetness, ins_depth,
and shape_b as decision variables, covering all physiological effectors in the
rule-based solver. Three heat-balance equality constraints (from
`HeatExchange.heat_balance`) are imposed; the objective minimises deviation
from the setpoint core temperature, metabolic heat above the minimum, and
evaporative effector use.

Uses a single mean-weighted body (dorsal×dmult + ventral×vmult) where
dmult = sky_view_factor + vegetation_view_factor. This matches the weighting
used by `solve_metabolic_rate`.

# Arguments
- `organism`: must have `OrganismTraits` with `BehavioralTraits{ThermoregulationLimits{IPOPTControl}}`
- `environment`: NamedTuple with `environment_pars` and `environment_vars`
- `Q_gen`: initial guess for metabolic heat generation (W)
- `T_skin`: initial guess for skin temperature (K)
- `T_insulation`: initial guess for insulation surface temperature (K)
"""
function thermoregulate(
    ::Endotherm,
    ::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    Q_gen_init,
    T_skin_init,
    T_ins_init;
    verbose = false,
)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars

    ins_pars  = insulation_pars(organism)
    ext_cond  = conduction_pars_external(organism)
    int_cond  = conduction_pars_internal(organism)
    rad_pars  = radiation_pars(organism)
    evap_pars = evaporation_pars(organism)
    resp_pars = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)
    limits    = thermoregulation(organism)

    T_air = environment_vars.air_temperature
    T_setpoint = metab_pars.core_temperature
    vegetation_temperature = environment_vars.reference_air_temperature

    # View factor geometry — mirrors solve_metabolic_rate
    vegetation_factor = rad_pars.sky_view_factor * environment_vars.shade
    sky_factor_ref    = rad_pars.sky_view_factor - vegetation_factor
    ground_factor_ref = 1 - sky_factor_ref - vegetation_factor
    bush_factor_ref   = rad_pars.bush_view_factor
    dmult = sky_factor_ref + vegetation_factor
    vmult = 1 - dmult

    # Mean-weighted insulation geometry
    fat  = Fat(int_cond.fat_fraction, int_cond.fat_density)
    pven = rad_pars.ventral_fraction

    mean_depth       = ins_pars.dorsal.depth       * dmult + ins_pars.ventral.depth       * vmult
    mean_diam        = ins_pars.dorsal.diameter    * dmult + ins_pars.ventral.diameter    * vmult
    mean_density     = ins_pars.dorsal.density     * dmult + ins_pars.ventral.density     * vmult
    mean_length      = ins_pars.dorsal.length      * dmult + ins_pars.ventral.length      * vmult
    mean_reflectance = ins_pars.dorsal.reflectance * dmult + ins_pars.ventral.reflectance * vmult
    mean_k_fibre     = ins_pars.dorsal.conductivity * dmult + ins_pars.ventral.conductivity * vmult

    mean_fibres = FibreProperties(;
        diameter     = mean_diam,
        length       = mean_length,
        density      = mean_density,
        depth        = mean_depth,
        reflectance  = mean_reflectance,
        conductivity = mean_k_fibre,
    )
    mean_ins_pars = InsulationParameters(;
        dorsal                  = mean_fibres,
        ventral                 = mean_fibres,
        depth_compressed        = ins_pars.depth_compressed,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )
    mean_fur = Fur(mean_depth, mean_diam, mean_density)
    body = Body(organism.body.shape, CompositeInsulation(mean_fur, fat))

    # Solar flows
    total_area_orig = BiophysicalGeometry.total_area(organism.body)
    area_cond_orig  = total_area_orig * ext_cond.conduction_fraction
    absorptivities  = Absorptivities(rad_pars, environment_pars)
    vf_solar        = ViewFactors(sky_factor_ref, ground_factor_ref, 0.0, 0.0)
    solar_conds     = SolarConditions(environment_vars)
    area_sil_orig   = silhouette_area(organism.body, rad_pars.solar_orientation)
    solar_out       = solar(organism.body, absorptivities, vf_solar, solar_conds,
                            area_sil_orig, area_cond_orig)

    if solar_out.solar_flow > 0.0u"W"
        dorsal_solar  = 2.0 * solar_out.direct_flow + solar_out.solar_sky_flow * 2.0
        ventral_solar = (solar_out.solar_substrate_flow /
                         (1.0 - sky_factor_ref - vegetation_factor)) *
                        (1.0 - 2.0 * ext_cond.conduction_fraction)
    else
        dorsal_solar  = 0.0u"W"
        ventral_solar = 0.0u"W"
    end
    solar_mean = dorsal_solar * dmult + ventral_solar * vmult

    # Mean view factors: dorsal sees sky+veg, ventral sees ground+bush
    vf_mean = ViewFactors(
        sky_factor_ref    * 2.0 * dmult,
        ground_factor_ref * 2.0 * vmult,
        bush_factor_ref   * 2.0 * vmult,
        vegetation_factor * 2.0,
    )

    # Mean conductance: only ventral side conducts to substrate
    body_total_area   = BiophysicalGeometry.total_area(body)
    area_cond_ventral = body_total_area * ext_cond.conduction_fraction * 2
    ventral_conductance = (area_cond_ventral * environment_vars.substrate_conductivity) /
                          environment_pars.conduction_depth
    mean_conductance = ventral_conductance * vmult

    # Insulation properties at initial temperature guess
    T_ins_avg  = T_ins_init * 0.7 + T_skin_init * 0.3
    insulation = insulation_properties(mean_ins_pars, T_ins_avg, pven)

    # Pack environment and traits for heat_balance
    temperature = EnvironmentTemperatures(
        T_air,
        environment_vars.sky_temperature,
        environment_vars.ground_temperature,
        vegetation_temperature,
        environment_vars.bush_temperature,
        environment_vars.substrate_temperature,
    )
    atmos = AtmosphericConditions(environment_vars)
    env_packed = (;
        temperature,
        view_factors           = vf_mean,
        atmos,
        fluid                  = environment_pars.fluid,
        solar_flow             = solar_mean,
        gas_fractions          = environment_pars.gas_fractions,
        convection_enhancement = environment_pars.convection_enhancement,
    )
    ϵ_body = rad_pars.body_emissivity_dorsal * dmult + rad_pars.body_emissivity_ventral * vmult
    traits_packed = (;
        fat_conductivity   = int_cond.fat_conductivity,
        flesh_conductivity = int_cond.flesh_conductivity,
        ϵ_body,
        skin_wetness       = evap_pars.skin_wetness,
        insulation_wetness = evap_pars.insulation_wetness,
        bare_skin_fraction = evap_pars.bare_skin_fraction,
        eye_fraction       = evap_pars.eye_fraction,
    )
    geometry_vars = GeometryVariables(;
        side                    = :dorsal,
        conductance_coefficient = mean_conductance,
        ventral_fraction        = pven,
        conduction_fraction     = ext_cond.conduction_fraction,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )

    # ---- Decision variable bounds (all SI, unitless Float64) ---------------
    # x = [T_core, T_skin, T_ins, Q_gen, k_flesh, pant, skin_wetness, ins_depth, shape_b]
    T_air_K   = ustrip(u"K", T_air)
    T_set_K   = ustrip(u"K", T_setpoint)
    T_core_lo = ustrip(u"K", limits.T_core.reference)
    T_core_hi = ustrip(u"K", limits.T_core.max)
    T_skin_lo = T_air_K - 5.0
    T_skin_hi = T_core_hi + 5.0
    Q_lo  = ustrip(u"W", limits.Q_minimum_ref)
    Q_hi  = ustrip(u"W", limits.Q_minimum_ref) * 20.0
    k_lo  = ustrip(u"W/m/K", limits.k_flesh.reference)
    k_hi  = ustrip(u"W/m/K", limits.k_flesh.max)
    p_lo  = 1.0
    p_hi  = Float64(limits.panting.pant.max)
    sw_lo = Float64(limits.skin_wetness.reference)
    sw_hi = Float64(limits.skin_wetness.max)
    ins_lo = ustrip(u"m", limits.insulation.dorsal.reference)
    ins_hi = ustrip(u"m", limits.insulation.dorsal.max)
    sh_lo  = Float64(limits.shape_b.reference)
    sh_hi  = organism.body.shape isa Sphere ? sh_lo : Float64(limits.shape_b.max)

    lb = [T_core_lo, T_skin_lo, T_skin_lo, Q_lo, k_lo, p_lo, sw_lo, ins_lo, sh_lo]
    ub = [T_core_hi, T_skin_hi, T_skin_hi, Q_hi, k_hi, p_hi, sw_hi, ins_hi, sh_hi]

    x0 = clamp.(
        [T_set_K,
         ustrip(u"K", T_skin_init),
         ustrip(u"K", T_ins_init),
         ustrip(u"W", Q_gen_init),
         ustrip(u"W/m/K", int_cond.flesh_conductivity),
         1.0,
         Float64(evap_pars.skin_wetness),
         ustrip(u"m", limits.insulation.dorsal.max),  # start with erected insulation
         sh_lo],                                        # start curled
        lb, ub,
    )

    Q_range = max(Q_hi - Q_lo, 1.0)
    params = (;
        T_setpoint = T_set_K,
        w_pant     = 1.0,
        w_sweat    = 1.0,
        w_qgen     = 0.01,
        Q_lo,
        Q_range,
        body,
        ins_pars   = mean_ins_pars,
        mean_fibres,
        fat,
        shape_pars = organism.body.shape,
        is_sphere  = organism.body.shape isa Sphere,
        insulation,
        geometry_vars,
        env        = env_packed,
        traits     = traits_packed,
        resp_pars,
        pven,
        substrate_conductivity = environment_vars.substrate_conductivity,
        conduction_depth       = environment_pars.conduction_depth,
        vmult,
        ext_cond,
    )

    # ---- Objective: stay at setpoint, minimise Q_gen above minimum, minimise evaporative cost
    function f(x, p)
        (x[1] - p.T_setpoint)^2 +
        p.w_qgen  * ((x[4] - p.Q_lo) / p.Q_range)^2 +
        p.w_pant  * (x[6] - 1.0)^2 +
        p.w_sweat * x[7]^2
    end

    # ---- Constraints: three heat-balance residuals = 0 --------------------
    function cons!(res, x, p)
        T_core    = x[1] * u"K"
        T_skin    = x[2] * u"K"
        T_ins     = x[3] * u"K"
        Q_gen     = x[4] * u"W"
        k_flesh   = x[5] * u"W/m/K"
        pant_v    = x[6]
        skin_w    = x[7]
        ins_d     = x[8] * u"m"
        shape_b_v = x[9]

        # Rebuild insulation and body from piloerection and shape_b decision variables
        new_fibres   = setproperties(p.mean_fibres; depth = ins_d)
        new_ins_pars = setproperties(p.ins_pars; dorsal = new_fibres, ventral = new_fibres)
        new_fur      = Fur(ins_d, p.mean_fibres.diameter, p.mean_fibres.density)
        new_shape    = p.is_sphere ? p.shape_pars : setproperties(p.shape_pars; b = shape_b_v)
        new_body     = Body(new_shape, CompositeInsulation(new_fur, p.fat))

        # Recompute temperature-dependent insulation properties
        T_ins_avg      = T_ins * 0.7 + T_skin * 0.3
        new_insulation = insulation_properties(new_ins_pars, T_ins_avg, p.pven)

        # Update conductance coefficient for new body size (substrate conduction area changes with shape)
        new_area_cond   = BiophysicalGeometry.total_area(new_body) * p.ext_cond.conduction_fraction * 2
        new_conductance = (new_area_cond * p.substrate_conductivity) / p.conduction_depth * p.vmult
        new_geometry_vars = setproperties(p.geometry_vars; conductance_coefficient = new_conductance)

        out = HeatExchange.heat_balance(
            T_core, T_skin, T_ins, Q_gen;
            body             = new_body,
            insulation_pars  = new_ins_pars,
            insulation       = new_insulation,
            geometry_vars    = new_geometry_vars,
            environment_vars = p.env,
            traits           = p.traits,
            resp_pars        = p.resp_pars,
            k_flesh          = k_flesh,
            pant             = pant_v,
            skin_wetness     = skin_w,
        )
        res[1] = ustrip(u"W", out.residual_energy_balance)
        res[2] = ustrip(u"W", out.residual_internal_conduction)
        res[3] = ustrip(u"K", out.residual_skin_temperature)
        return nothing
    end

    # Bypass DifferentiationInterface (incompatible with Unitful params NamedTuple).
    # Use out-of-place FiniteDiff forms to avoid JacobianCache argument requirement.
    function grad!(g, x, p)
        FiniteDiff.finite_difference_gradient!(g, x_ -> f(x_, p), x)
    end
    function hess!(H, x, p)
        H .= FiniteDiff.finite_difference_hessian(x_ -> f(x_, p), x)
    end
    function cons_j!(J, x, p)
        J .= FiniteDiff.finite_difference_jacobian(
            x_ -> (res = zeros(eltype(x_), 3); cons!(res, x_, p); res),
            x,
        )
    end
    function cons_h!(res, x, p)
        for i in eachindex(res)
            res[i] .= FiniteDiff.finite_difference_hessian(
                x_ -> begin
                    r = zeros(3)
                    cons!(r, x_, p)
                    r[i]
                end,
                x,
            )
        end
    end

    optf = OptimizationFunction(f, SciMLBase.NoAD();
        cons = cons!, grad = grad!, hess = hess!, cons_j = cons_j!, cons_h = cons_h!,
    )
    prob = OptimizationProblem(optf, x0, params;
        lb = lb, ub = ub,
        lcons = zeros(3), ucons = zeros(3),
    )
    sol = solve(prob, IpoptOptimizer(); verbose)

    # ---- Evaluate heat_balance at solution for full output -----------------
    xsol = sol.u
    T_core_sol    = xsol[1] * u"K"
    T_skin_sol    = xsol[2] * u"K"
    T_ins_sol     = xsol[3] * u"K"
    Q_gen_sol     = xsol[4] * u"W"
    k_flesh_sol   = xsol[5] * u"W/m/K"
    pant_sol      = xsol[6]
    skin_w_sol    = xsol[7]
    ins_depth_sol = xsol[8] * u"m"
    shape_b_sol   = xsol[9]

    # Rebuild body and insulation at solution values
    sol_fibres    = setproperties(mean_fibres; depth = ins_depth_sol)
    sol_ins_pars  = setproperties(mean_ins_pars; dorsal = sol_fibres, ventral = sol_fibres)
    sol_fur       = Fur(ins_depth_sol, mean_fibres.diameter, mean_fibres.density)
    sol_shape     = organism.body.shape isa Sphere ?
                    organism.body.shape :
                    setproperties(organism.body.shape; b = shape_b_sol)
    body_sol      = Body(sol_shape, CompositeInsulation(sol_fur, fat))
    T_ins_avg_sol = T_ins_sol * 0.7 + T_skin_sol * 0.3
    insulation_sol = insulation_properties(sol_ins_pars, T_ins_avg_sol, pven)
    area_cond_sol  = BiophysicalGeometry.total_area(body_sol) * ext_cond.conduction_fraction * 2
    conductance_sol = (area_cond_sol * environment_vars.substrate_conductivity) /
                      environment_pars.conduction_depth * vmult
    geometry_vars_sol = setproperties(geometry_vars; conductance_coefficient = conductance_sol)

    hb = HeatExchange.heat_balance(
        T_core_sol, T_skin_sol, T_ins_sol, Q_gen_sol;
        body             = body_sol,
        insulation_pars  = sol_ins_pars,
        insulation       = insulation_sol,
        geometry_vars    = geometry_vars_sol,
        environment_vars = env_packed,
        traits           = traits_packed,
        resp_pars        = resp_pars,
        k_flesh          = k_flesh_sol,
        pant             = pant_sol,
        skin_wetness     = skin_w_sol,
    )

    lung_temperature = (T_core_sol + T_skin_sol) / 2

    # Longwave flows (mirrors solve_metabolic_rate simplified approach)
    σ_const        = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    area_total_lw  = BiophysicalGeometry.total_area(body_sol)
    area_radiant_v = area_total_lw * (1 - ext_cond.conduction_fraction)
    dorsal_lw_out  = 2 * sky_factor_ref * σ_const * rad_pars.body_emissivity_dorsal *
                     area_total_lw * T_ins_sol^4
    ventral_lw_out = 2 * ground_factor_ref * σ_const * rad_pars.body_emissivity_ventral *
                     area_radiant_v * T_ins_sol^4
    longwave_flow_out = dorsal_lw_out * dmult + ventral_lw_out * vmult
    longwave_flow_in  = longwave_flow_out - hb.radiation_heat_flow

    # Respiration mass flows at solution point
    resp_pars_eff = setproperties(resp_pars; pant = pant_sol)
    resp_out = HeatExchange.respiration(
        MetabolicRates(; metabolic = Q_gen_sol, sum = Q_gen_sol, minimum = 0.0u"W"),
        resp_pars_eff,
        env_packed.atmos,
        organism.body.shape.mass,
        lung_temperature,
        T_air;
        gas_fractions = environment_pars.gas_fractions,
        O2conversion  = Kleiber1961(),
    )
    latent_heat_vap = FluidProperties.enthalpy_of_vaporisation(T_air)
    m_sweat = u"g/hr"(hb.skin_evaporation_heat_flow / latent_heat_vap)
    m_evap  = u"g/hr"(resp_out.respiration_mass + m_sweat)

    thermoregulation_out = (;
        core_temperature               = T_core_sol,
        skin_temperature               = T_skin_sol,
        insulation_temperature         = T_ins_sol,
        lung_temperature,
        skin_temperature_dorsal        = T_skin_sol,
        skin_temperature_ventral       = T_skin_sol,
        insulation_temperature_dorsal  = T_ins_sol,
        insulation_temperature_ventral = T_ins_sol,
        shape_b                        = shape_b_sol,
        pant                           = pant_sol,
        skin_wetness                   = skin_w_sol,
        flesh_conductivity             = k_flesh_sol,
        insulation_conductivity_effective  = insulation_sol.conductivities.average,
        insulation_conductivity_dorsal     = insulation_sol.conductivities.dorsal,
        insulation_conductivity_ventral    = insulation_sol.conductivities.ventral,
        insulation_conductivity_compressed = insulation_sol.conductivity_compressed,
        insulation_depth_dorsal  = sol_ins_pars.dorsal.depth,
        insulation_depth_ventral = sol_ins_pars.ventral.depth,
        metab_pars.q10,
    )

    total_area_out = BiophysicalGeometry.total_area(body_sol)
    area_skin_out  = skin_area(body_sol)
    area_evap_out  = evaporation_area(body_sol)
    area_conv_out  = total_area_out * (1 - ext_cond.conduction_fraction)
    area_sil_out   = silhouette_area(body_sol, rad_pars.solar_orientation)
    volume_out     = body_sol.geometry.volume

    morphology = (;
        total_area         = total_area_out,
        area_skin          = area_skin_out,
        area_evaporation   = area_evap_out,
        area_convection    = area_conv_out,
        area_conduction    = area_cond_orig / 2,
        area_silhouette    = area_sil_out,
        sky_view_factor    = sky_factor_ref,
        ground_view_factor = ground_factor_ref,
        volume             = volume_out,
        volume_flesh       = flesh_volume(body_sol),
        characteristic_dimension = HeatExchange.characteristic_dimension(HeatExchange.VolumeCubeRoot(), body_sol),
        fat_mass           = body_sol.shape.mass * fat.fraction,
        body_sol.geometry.length...,
    )

    energy_flows = (;
        solar_flow            = hb.solar_heat_flow,
        longwave_flow_in,
        generated_heat_flow   = Q_gen_sol,
        evaporation_heat_flow = hb.skin_evaporation_heat_flow +
                                hb.insulation_evaporation_heat_flow +
                                hb.respiration_heat_flow,
        longwave_flow_out,
        convection_heat_flow  = hb.convection_heat_flow,
        conduction_flow       = hb.conduction_heat_flow,
        balance               = hb.residual_energy_balance,
        ntry                  = 1,
        success               = true,
    )

    mass_flows = (;
        air_flow             = resp_out.air_flow,
        oxygen_flow_standard = resp_out.oxygen_flow_standard,
        m_evap,
        respiration_mass     = resp_out.respiration_mass,
        m_sweat,
        molar_fluxes_in      = resp_out.molar_fluxes_in,
        molar_fluxes_out     = resp_out.molar_fluxes_out,
    )

    return (; thermoregulation = thermoregulation_out, morphology, energy_flows, mass_flows)
end
