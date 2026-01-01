function example_environment_vars(;
    T_air=u"K"((20.0)u"°C"),
    rh=0.05,
    wind_speed=0.1u"m/s",
    P_atmos=101325.0u"Pa",
    zenith_angle=20.0u"°",
    k_substrate=2.79u"W/m/K",
    global_radiation=0.0u"W/m^2",
    diffuse_fraction=0.0,
    shade=0,
)
    EnvironmentalVars(;
        T_air,
        T_air_reference=T_air,
        T_sky=T_air,
        T_ground=T_air,
        T_substrate=T_air,
        T_bush=T_air,
        T_vegetation=T_air,
        rh,
        wind_speed,
        P_atmos,
        zenith_angle,
        k_substrate,
        global_radiation,
        diffuse_fraction,
        shade,
    )
end

function example_environment_pars(;
    α_ground=0.8,
    ϵ_ground=1.0,
    ϵ_sky=1.0,
    elevation=0.0u"m",
    fluid=0,
    fN2=0.7902,
    fO2=0.2095,
    fCO2=0.000412,
    convection_enhancement=1.0,
)
    EnvironmentalPars(;
        α_ground,
        ϵ_ground,
        ϵ_sky,
        elevation,
        fluid,
        fN2,
        fO2,
        fCO2,
        convection_enhancement,
    )
end

function example_shape_pars(;
    mass=65.0u"kg",
    ρ_flesh=1000.0u"kg/m^3",
    shape_b=1.1,
    shape_c=1.1,
)
    Ellipsoid(mass, ρ_flesh, shape_b, shape_c)
end

function example_conduction_pars_external(;
    conduction_fraction=0.0,
)
    ExternalConductionParameters(;
        conduction_fraction,
    )
end

function example_conduction_pars_internal(;
    fat_fraction=0.0,
    k_flesh=0.9u"W/m/K",
    k_fat=0.23u"W/m/K",
    ρ_fat=901.0u"kg/m^3",
)
    InternalConductionParameters(;
        fat_fraction,
        k_flesh,
        k_fat,
        ρ_fat,
    )
end

function example_radiation_pars(;
    α_body_dorsal=0.8,
    α_body_ventral=0.8,
    ϵ_body_dorsal=0.99,
    ϵ_body_ventral=0.99,
    F_sky=0.5,
    F_ground=0.5,
    F_bush=0.0,
    ventral_fraction=0.5,
    solar_orientation=Intermediate(),
)
    RadiationParameters(;
        α_body_dorsal,
        α_body_ventral,
        ϵ_body_dorsal,
        ϵ_body_ventral,
        F_sky,
        F_ground,
        F_bush,
        ventral_fraction,
        solar_orientation,
    )
end

function example_evaporation_pars(;
    skin_wetness=0.005,
    insulation_wetness=0.0,
    eye_fraction=0.0,
    bare_skin_fraction=0.0,
    insulation_fraction=1.0,
)
    EvaporationParameters(;
        skin_wetness,
        insulation_wetness,
        eye_fraction,
        bare_skin_fraction,
        insulation_fraction,
    )
end

function example_hydraulic_pars(;
    water_potential=0.0u"J/kg",
    hydraulic_conductance=0.0u"kg / (m^2 * s * (J/kg))",
    specific_hydration=0.000304u"m^3 / (m^3 * (J/kg))",
)
    HydraulicParameters(;
        water_potential,
        hydraulic_conductance,
        specific_hydration,
    )
end

function example_respiration_pars(;
    fO2_extract=0.2,
    pant=1.0,
    rq=0.8,
    Δ_breath=0.0u"K",
    rh_exit=1.0,
)
    RespirationParameters(;
        fO2_extract,
        pant,
        rq,
        Δ_breath,
        rh_exit,
    )
end

function example_metabolism_pars(;
    T_core=u"K"((37.0)u"°C"),
    Q_metabolism=77.61842u"W",
    q10=2.0,
    model=Kleiber(),
)
    MetabolismParameters(;
        T_core,
        Q_metabolism,
        q10,
        model,
    )
end

function example_insulation_pars(;
    insulation_conductivity_dorsal=nothing,
    insulation_conductivity_ventral=nothing,
    fibre_diameter_dorsal=30e-06u"m",
    fibre_diameter_ventral=30e-06u"m",
    fibre_length_dorsal=23.9e-03u"m",
    fibre_length_ventral=23.9e-03u"m",
    insulation_depth_dorsal=2e-03u"m",
    insulation_depth_ventral=2e-03u"m",
    fibre_density_dorsal=3000e+04u"1/m^2",
    fibre_density_ventral=3000e+04u"1/m^2",
    insulation_reflectance_dorsal=0.2,
    insulation_reflectance_ventral=0.2,
    insulation_depth_compressed=2e-03u"m",
    fibre_conductivity=0.209u"W/m/K",
    longwave_depth_fraction=1.0,
)
    InsulationParameters(;
        insulation_conductivity_dorsal,
        insulation_conductivity_ventral,
        fibre_diameter_dorsal,
        fibre_diameter_ventral,
        fibre_length_dorsal,
        fibre_length_ventral,
        insulation_depth_dorsal,
        insulation_depth_ventral,
        fibre_density_dorsal,
        fibre_density_ventral,
        insulation_reflectance_dorsal,
        insulation_reflectance_ventral,
        insulation_depth_compressed,
        fibre_conductivity,
        longwave_depth_fraction,
    )
end

function example_model_pars(;
    respire=true,
    simulsol_tolerance=1e-3u"K",
    resp_tolerance=1e-5,
)
    EndoModelPars(;
        respire,
        simulsol_tolerance,
        resp_tolerance,
    )
end

function example_endotherm_thermoregulation_pars(;
    thermoregulation_mode = 1,
    tolerance = 0.005,
    max_iterations = 1000,

    Q_minimum = 77.61842u"W",
    Q_minimum_ref = 77.61842u"W",

    insulation_depth_dorsal = 2e-03u"m",
    insulation_depth_ventral = 2e-03u"m",
    insulation_depth_dorsal_max = 2e-03u"m",
    insulation_depth_ventral_max = 2e-03u"m",
    insulation_depth_dorsal_ref = 2e-03u"m",
    insulation_depth_ventral_ref = 2e-03u"m",
    insulation_step = 0.0,

    shape_b = 1.1,
    shape_b_step = 0.1,
    shape_b_max = 5.0,

    k_flesh = 0.9u"W/m/K",
    k_flesh_step = 0.1u"W/m/K",
    k_flesh_max = 2.8u"W/m/K",

    T_core = (37.0 + 273.15)u"K",
    T_core_step = 0.1u"K",
    T_core_max = (39.0 + 273.15)u"K",
    T_core_ref = (37.0 + 273.15)u"K",

    pant = 1.0,
    pant_step = 0.1,
    pant_max = 10.0,
    pant_cost = 0.0u"W",
    pant_multiplier = 1.05,

    skin_wetness = 0.005,
    skin_wetness_step = 0.001,
    skin_wetness_max = 1.0,
)
    EndothermThermoregulationParameters(;
        thermoregulation_mode,
        tolerance,
        max_iterations,
        
        Q_minimum,
        Q_minimum_ref,

        insulation_depth_dorsal,
        insulation_depth_ventral,
        insulation_depth_dorsal_max,
        insulation_depth_ventral_max,
        insulation_depth_dorsal_ref,
        insulation_depth_ventral_ref,
        insulation_step,

        shape_b,
        shape_b_step,
        shape_b_max,

        k_flesh,
        k_flesh_step,
        k_flesh_max,

        T_core,
        T_core_step,
        T_core_max,
        T_core_ref,

        pant,
        pant_step,
        pant_max,
        pant_cost,
        pant_multiplier,

        skin_wetness,
        skin_wetness_step,
        skin_wetness_max,
    )
end
