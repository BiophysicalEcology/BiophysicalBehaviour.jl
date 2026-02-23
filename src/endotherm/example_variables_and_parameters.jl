function example_environment_vars(;
    air_temperature=u"K"((20.0)u"°C"),
    relative_humidity=0.05,
    wind_speed=0.1u"m/s",
    atmospheric_pressure=101325.0u"Pa",
    zenith_angle=20.0u"°",
    k_substrate=2.79u"W/m/K",
    global_radiation=0.0u"W/m^2",
    diffuse_fraction=0.0,
    shade=0,
)
    EnvironmentalVars(;
        air_temperature,
        reference_air_temperature=air_temperature,
        sky_temperature=air_temperature,
        ground_temperature=air_temperature,
        substrate_temperature=air_temperature,
        bush_temperature=air_temperature,
        vegetation_temperature=air_temperature,
        relative_humidity,
        wind_speed,
        atmospheric_pressure,
        zenith_angle,
        k_substrate,
        global_radiation,
        diffuse_fraction,
        shade,
    )
end

function example_environment_pars(;
    ground_albedo=0.8,
    ground_emissivity=1.0,
    sky_emissivity=1.0,
    elevation=0.0u"m",
    fluid=0,
    gas_fractions=FluidProperties.GasFractions(),
    convection_enhancement=1.0,
)
    EnvironmentalPars(;
        ground_albedo,
        ground_emissivity,
        sky_emissivity,
        elevation,
        fluid,
        gas_fractions,
        convection_enhancement,
    )
end

function example_ellipsoid_shape_pars(;
    mass=65.0u"kg",
    ρ_flesh=1000.0u"kg/m^3",
    shape_coefficient_b=1.1,
    shape_coefficient_c=1.1,
)
    Ellipsoid(mass, ρ_flesh, shape_coefficient_b, shape_coefficient_c)
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
    flesh_conductivity=0.9u"W/m/K",
    fat_conductivity=0.23u"W/m/K",
    fat_density=901.0u"kg/m^3",
)
    InternalConductionParameters(;
        fat_fraction,
        flesh_conductivity,
        fat_conductivity,
        fat_density,
    )
end

function example_radiation_pars(;
    body_absorptivity_dorsal=0.8,
    body_absorptivity_ventral=0.8,
    body_emissivity_dorsal=0.99,
    body_emissivity_ventral=0.99,
    sky_view_factor=0.5,
    ground_view_factor=0.5,
    bush_view_factor=0.0,
    ventral_fraction=0.5,
    solar_orientation=Intermediate(),
)
    RadiationParameters(;
        body_absorptivity_dorsal,
        body_absorptivity_ventral,
        body_emissivity_dorsal,
        body_emissivity_ventral,
        sky_view_factor,
        ground_view_factor,
        bush_view_factor,
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
    oxygen_extraction_efficiency=0.2,
    panting_rate=1.0,
    respiratory_quotient=0.8,
    exhaled_temperature_offset=0.0u"K",
    exhaled_relative_humidity=1.0,
)
    RespirationParameters(;
        oxygen_extraction_efficiency,
        pant=panting_rate,
        respiratory_quotient,
        exhaled_temperature_offset,
        exhaled_relative_humidity,
    )
end

function example_metabolism_pars(;
    core_temperature=u"K"((37.0)u"°C"),
    metabolic_flux=77.61842u"W",
    q10=2.0,
    model=Kleiber(),
)
    MetabolismParameters(;
        core_temperature,
        metabolic_flux,
        q10,
        model,
    )
end

function example_insulation_pars(;
    # Dorsal fibre properties
    fibre_diameter_dorsal=30e-06u"m",
    fibre_length_dorsal=23.9e-03u"m",
    fibre_density_dorsal=3000e+04u"1/m^2",
    insulation_depth_dorsal=2e-03u"m",
    insulation_reflectance_dorsal=0.2,
    fibre_conductivity_dorsal=0.209u"W/m/K",
    # Ventral fibre properties
    fibre_diameter_ventral=30e-06u"m",
    fibre_length_ventral=23.9e-03u"m",
    fibre_density_ventral=3000e+04u"1/m^2",
    insulation_depth_ventral=2e-03u"m",
    insulation_reflectance_ventral=0.2,
    fibre_conductivity_ventral=0.209u"W/m/K",
    # Global insulation parameters
    insulation_depth_compressed=2e-03u"m",
    longwave_depth_fraction=1.0,
)
    dorsal = HeatExchange.FibreProperties(;
        diameter=fibre_diameter_dorsal,
        length=fibre_length_dorsal,
        density=fibre_density_dorsal,
        depth=insulation_depth_dorsal,
        reflectance=insulation_reflectance_dorsal,
        conductivity=fibre_conductivity_dorsal,
    )
    ventral = HeatExchange.FibreProperties(;
        diameter=fibre_diameter_ventral,
        length=fibre_length_ventral,
        density=fibre_density_ventral,
        depth=insulation_depth_ventral,
        reflectance=insulation_reflectance_ventral,
        conductivity=fibre_conductivity_ventral,
    )
    InsulationParameters(;
        dorsal,
        ventral,
        depth_compressed=insulation_depth_compressed,
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
    minimum_metabolic_flux_ref=77.61842u"W",
    # Insulation (piloerection)
    insulation_depth_dorsal=2e-03u"m",
    insulation_depth_ventral=2e-03u"m",
    insulation_depth_dorsal_max=2e-03u"m",
    insulation_depth_ventral_max=2e-03u"m",
    insulation_depth_dorsal_ref=2e-03u"m",
    insulation_depth_ventral_ref=2e-03u"m",
    insulation_step=0.0,
    # Shape (uncurl)
    shape_coefficient_b=1.1,
    shape_coefficient_b_step=0.1,
    shape_coefficient_b_max=5.0,
    # Tissue conductivity (vasodilation)
    flesh_conductivity=0.9u"W/m/K",
    flesh_conductivity_step=0.1u"W/m/K",
    flesh_conductivity_max=2.8u"W/m/K",
    # Core temperature (hyperthermia)
    core_temperature=(37.0 + 273.15)u"K",
    core_temperature_step=0.1u"K",
    core_temperature_max=(39.0 + 273.15)u"K",
    core_temperature_ref=(37.0 + 273.15)u"K",
    # Panting
    panting_rate_current=1.0,
    panting_rate_step=0.1,
    panting_rate_max=10.0,
    panting_cost=0.0u"W",
    panting_multiplier=1.05,
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

    shape_coefficient_b_param = SteppedParameter(;
        current=shape_coefficient_b,
        max=shape_coefficient_b_max,
        step=shape_coefficient_b_step,
    )

    flesh_conductivity_param = SteppedParameter(;
        current=flesh_conductivity,
        max=flesh_conductivity_max,
        step=flesh_conductivity_step,
    )

    core_temperature_param = SteppedParameter(;
        current=core_temperature,
        reference=core_temperature_ref,
        max=core_temperature_max,
        step=core_temperature_step,
    )

    panting = PantingLimits(;
        panting_rate=SteppedParameter(;
            current=panting_rate_current,
            max=panting_rate_max,
            step=panting_rate_step,
        ),
        cost=panting_cost,
        multiplier=panting_multiplier,
        core_temperature_ref=core_temperature_ref,
    )

    skin_wetness_param = SteppedParameter(;
        current=skin_wetness,
        max=skin_wetness_max,
        step=skin_wetness_step,
    )

    ThermoregulationLimits(;
        control,
        minimum_metabolic_flux_ref,
        insulation,
        shape_coefficient_b=shape_coefficient_b_param,
        flesh_conductivity=flesh_conductivity_param,
        core_temperature=core_temperature_param,
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
