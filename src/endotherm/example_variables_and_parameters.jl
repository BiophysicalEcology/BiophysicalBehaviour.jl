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
    minimum_heat_flow=77.61842u"W",
    # Insulation (piloerection)
    insulation_depth_dorsal=2e-03u"m",
    insulation_depth_ventral=2e-03u"m",
    insulation_depth_dorsal_max=2e-03u"m",
    insulation_depth_ventral_max=2e-03u"m",
    insulation_depth_dorsal_ref=2e-03u"m",
    insulation_depth_ventral_ref=2e-03u"m",
    insulation_step=0.0,
    # Shape (uncurl)
    aspect_ratio_factor=1.1,
    aspect_ratio_step=0.1,
    aspect_ratio_max=5.0,
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

    aspect_ratio_param = SteppedParameter(;
        current=aspect_ratio_factor,
        max=aspect_ratio_max,
        step=aspect_ratio_step,
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
        pant=SteppedParameter(;
            current=pant_current,
            max=pant_max,
            step=pant_step,
        ),
        cost=pant_cost,
        multiplier=pant_multiplier,
        core_temperature_ref=core_temperature_ref,
    )

    skin_wetness_param = SteppedParameter(;
        current=skin_wetness,
        max=skin_wetness_max,
        step=skin_wetness_step,
    )

    ThermoregulationLimits(;
        control,
        minimum_heat_flow,
        insulation,
        aspect_ratio_factor=aspect_ratio_param,
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
    activity_period=Diurnal(),
)
    BehavioralTraits(; thermoregulation, activity_period)
end

"""
    example_organism_traits(; thermal_strategy=Endotherm(), heat_exchange=nothing, behavior=nothing, kwargs...)

Create example `OrganismTraits` combining thermal strategy, heat exchange traits and behavior.

If `heat_exchange` is not provided, creates a default `HeatExchangeTraits`.
If `behavior` is not provided, creates a default `BehavioralTraits`.

Additional keyword arguments are passed to `example_heat_exchange_traits`.
"""
function example_organism_traits(;
    thermal_strategy=Endotherm(),
    heat_exchange=nothing,
    behavior=nothing,
    kwargs...,
)
    he = isnothing(heat_exchange) ? example_heat_exchange_traits(; kwargs...) : heat_exchange
    behav = isnothing(behavior) ? example_behavioral_traits() : behavior
    OrganismTraits(thermal_strategy, he, behav)
end
