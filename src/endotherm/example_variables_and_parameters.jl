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
    gasfrac=FluidProperties.GasFractions(),
    convection_enhancement=1.0,
)
    EnvironmentalPars(;
        α_ground,
        ϵ_ground,
        ϵ_sky,
        elevation,
        fluid,
        gasfrac,
        convection_enhancement,
    )
end

function example_ellipsoid_shape_pars(;
    mass=65.0u"kg",
    ρ_flesh=1000.0u"kg/m^3",
    shape_b=1.1,
    shape_c=1.1,
)
    Ellipsoid(mass, ρ_flesh, shape_b, shape_c)
end

# Alias for convenience
example_shape_pars(; kwargs...) = example_ellipsoid_shape_pars(; kwargs...)

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

function example_metabolic_rate_options(;
    respire=true,
    simulsol_tolerance=1e-3u"K",
    resp_tolerance=1e-5,
)
    SolveMetabolicRateOptions(;
        respire,
        simulsol_tolerance,
        resp_tolerance,
    )
end

"""
    example_thermoregulation_limits(; kwargs...)

Create example `ThermoregulationLimits` with sensible defaults.
"""
function example_thermoregulation_limits(;
    # Control
    thermoregulation_mode=CoreFirst(),
    tolerance=0.005,
    max_iterations=1000,
    # Metabolic reference
    Q_minimum_ref=77.61842u"W",
    # Insulation (piloerection)
    insulation_depth_dorsal=2e-03u"m",
    insulation_depth_ventral=2e-03u"m",
    insulation_depth_dorsal_max=2e-03u"m",
    insulation_depth_ventral_max=2e-03u"m",
    insulation_depth_dorsal_ref=2e-03u"m",
    insulation_depth_ventral_ref=2e-03u"m",
    insulation_step=0.0,
    # Shape (uncurl)
    shape_b=1.1,
    shape_b_step=0.1,
    shape_b_max=5.0,
    # Tissue conductivity (vasodilation)
    k_flesh=0.9u"W/m/K",
    k_flesh_step=0.1u"W/m/K",
    k_flesh_max=2.8u"W/m/K",
    # Core temperature (hyperthermia)
    T_core=(37.0 + 273.15)u"K",
    T_core_step=0.1u"K",
    T_core_max=(39.0 + 273.15)u"K",
    T_core_ref=(37.0 + 273.15)u"K",
    # Panting
    pant_current=1.0,
    pant_step=0.1,
    pant_max=10.0,
    pant_cost=0.0u"W",
    pant_multiplier=1.05,
    # Sweating
    skin_wetness=0.005,
    skin_wetness_step=0.001,
    skin_wetness_max=1.0,
)
    control = RuleBasedSequentialControl(;
        mode=thermoregulation_mode,
        tolerance,
        max_iterations,
    )

    insulation = InsulationLimits(;
        dorsal=SteppedParameter(;
            current=insulation_depth_dorsal,
            reference=insulation_depth_dorsal_ref,
            max=insulation_depth_dorsal_max,
            step=insulation_step,
        ),
        ventral=SteppedParameter(;
            current=insulation_depth_ventral,
            reference=insulation_depth_ventral_ref,
            max=insulation_depth_ventral_max,
            step=insulation_step,
        ),
    )

    shape_b_param = SteppedParameter(;
        current=shape_b,
        max=shape_b_max,
        step=shape_b_step,
    )

    k_flesh_param = SteppedParameter(;
        current=k_flesh,
        max=k_flesh_max,
        step=k_flesh_step,
    )

    T_core_param = SteppedParameter(;
        current=T_core,
        reference=T_core_ref,
        max=T_core_max,
        step=T_core_step,
    )

    panting = PantingLimits(;
        pant=SteppedParameter(;
            current=pant_current,
            max=pant_max,
            step=pant_step,
        ),
        cost=pant_cost,
        multiplier=pant_multiplier,
        T_core_ref=T_core_ref,
    )

    skin_wetness_param = SteppedParameter(;
        current=skin_wetness,
        max=skin_wetness_max,
        step=skin_wetness_step,
    )

    ThermoregulationLimits(;
        control,
        Q_minimum_ref,
        insulation,
        shape_b=shape_b_param,
        k_flesh=k_flesh_param,
        T_core=T_core_param,
        panting,
        skin_wetness=skin_wetness_param,
    )
end

"""
    example_behavioral_traits(; kwargs...)

Create example `BehavioralTraits` with sensible defaults.
"""
function example_behavioral_traits(;
    thermoregulation=example_thermoregulation_limits(),
    activity=Diurnal(),
)
    BehavioralTraits(; thermoregulation, activity)
end

"""
    example_organism_traits(; thermal_strategy=Endotherm(), physiology=nothing, behavior=nothing, kwargs...)

Create example `OrganismTraits` combining thermal strategy, physiology and behavior.

If `physiology` is not provided, creates a default `HeatExchangeTraits`.
If `behavior` is not provided, creates a default `BehavioralTraits`.

Additional keyword arguments are passed to `example_heat_exchange_traits`.
"""
function example_organism_traits(;
    thermal_strategy=Endotherm(),
    physiology=nothing,
    behavior=nothing,
    kwargs...,
)
    phys = isnothing(physiology) ? example_heat_exchange_traits(; kwargs...) : physiology
    behav = isnothing(behavior) ? example_behavioral_traits() : behavior
    OrganismTraits(thermal_strategy, phys, behav)
end

"""
    example_heat_exchange_traits(; kwargs...)

Create example `HeatExchangeTraits` with sensible defaults.
"""
function example_heat_exchange_traits(;
    shape_pars=example_shape_pars(),
    insulation_pars=example_insulation_pars(),
    conduction_pars_external=example_conduction_pars_external(),
    conduction_pars_internal=example_conduction_pars_internal(),
    convection_pars=ConvectionParameters(),
    radiation_pars=example_radiation_pars(),
    evaporation_pars=example_evaporation_pars(),
    hydraulic_pars=example_hydraulic_pars(),
    respiration_pars=example_respiration_pars(),
    metabolism_pars=example_metabolism_pars(),
    options=example_metabolic_rate_options(),
)
    HeatExchange.HeatExchangeTraits(;
        shape_pars,
        insulation_pars,
        conduction_pars_external,
        conduction_pars_internal,
        convection_pars,
        radiation_pars,
        evaporation_pars,
        hydraulic_pars,
        respiration_pars,
        metabolism_pars,
        options,
    )
end
