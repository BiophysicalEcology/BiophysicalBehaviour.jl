# =============================================================================
# Mode dispatch helpers
# =============================================================================

"""
    simultaneous_pant(mode::AbstractThermoregulationMode) -> Bool

Return true if the mode allows panting as an effector.
"""
simultaneous_pant(::CoreFirst) = false
simultaneous_pant(::CoreAndPantingFirst) = true
simultaneous_pant(::CorePantingSweatingFirst) = true

"""
    simultaneous_sweat(mode::AbstractThermoregulationMode) -> Bool

Return true if the mode allows sweating as an effector.
"""
simultaneous_sweat(::CoreFirst) = false
simultaneous_sweat(::CoreAndPantingFirst) = false
simultaneous_sweat(::CorePantingSweatingFirst) = true

# =============================================================================
# Thermoregulation entry points
# =============================================================================

"""
    thermoregulate(organism, environment, generated_heat_flow, skin_temperature, insulation_temperature)

Run the thermoregulation loop to find heat balance.

Dispatches on the organism's thermal strategy (`Endotherm`, `Ectotherm`, `Heterotherm`).

The organism must have `OrganismTraits` containing `BehavioralTraits` with
thermoregulation limits. Metabolic rate options are extracted from the organism's
HeatExchangeTraits.

Returns the result from `solve_metabolic_rate`.
"""
function thermoregulate(
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    thermoregulate(
        thermal_strategy(organism),
        organism,
        environment,
        generated_heat_flow,
        skin_temperature,
        insulation_temperature,
    )
end

"""
    thermoregulate(::Endotherm, organism, environment, generated_heat_flow, skin_temperature, insulation_temperature)

Run the endotherm thermoregulation loop to find heat balance.

Dispatches on the control strategy from the organism's thermoregulation limits.
"""
function thermoregulate(
    ::Endotherm,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    thermoregulate(
        Endotherm(),
        control_strategy(organism),
        organism,
        environment,
        generated_heat_flow,
        skin_temperature,
        insulation_temperature,
    )
end

"""
    thermoregulate(::Endotherm, ::RuleBasedSequentialControl, organism, environment, generated_heat_flow, skin_temperature, insulation_temperature)

Run the endotherm thermoregulation loop using rule-based sequential control.

Applies thermoregulation behaviors in order:
1. Reduce insulation (piloerection)
2. Uncurl (increase surface area)
3. Vasodilate (increase tissue conductivity)
4. Hyperthermia (allow core temperature to rise)
5. Pant (evaporative cooling via respiration)
6. Sweat (evaporative cooling via skin)

Returns the final `endotherm_out` result from `solve_metabolic_rate`.
"""
function thermoregulate(
    ::Endotherm,
    ::RuleBasedSequentialControl,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    # Extract thermoregulation limits from organism's behavioral traits
    limits = thermoregulation(organism)
    endotherm_out = nothing

    # Extract control parameters
    (; mode, tolerance, max_iterations) = limits.control

    # Extract current limits (will be updated during loop)
    insulation_limits = limits.insulation
    shape_coefficient_b_limits = limits.shape_coefficient_b
    flesh_conductivity_limits = limits.flesh_conductivity
    core_temperature_limits = limits.core_temperature
    panting_limits = limits.panting
    skin_wetness_limits = limits.skin_wetness

    # Check if starting in piloerect state
    if insulation_limits.dorsal.step > 0.0 &&
       (insulation_limits.dorsal.current + insulation_limits.ventral.current) > 0u"mm"
        # Start with erect insulation (set to max)
        insulation_limits = @set insulation_limits.dorsal.current = insulation_limits.dorsal.max
        insulation_limits = @set insulation_limits.ventral.current = insulation_limits.ventral.max
        # Apply to organism with step=0 (just set to max, don't decrement)
        zero_step_limits = @set insulation_limits.dorsal.step = 0.0
        zero_step_limits = @set zero_step_limits.ventral.step = 0.0
        insulation_limits, organism = piloerect(organism, zero_step_limits)
        # Restore original step values
        insulation_limits = @set insulation_limits.dorsal.step = limits.insulation.dorsal.step
        insulation_limits = @set insulation_limits.ventral.step = limits.insulation.ventral.step
    end

    endotherm_out = solve_metabolic_rate(organism, environment, skin_temperature, insulation_temperature)
    skin_temperature = endotherm_out.thermoregulation.skin_temperature
    insulation_temperature = endotherm_out.thermoregulation.insulation_temperature
    generated_heat_flow = endotherm_out.energy_flows.generated_heat_flow

    # Current minimum_metabolic_heat_flow (may be modified by panting/hyperthermia)
    minimum_metabolic_heat_flow = limits.minimum_metabolic_heat_flow_ref

    iteration = 0

    # Start of thermoregulation loop
    while generated_heat_flow < minimum_metabolic_heat_flow * (1 - tolerance)
        iteration += 1
        if iteration > max_iterations
            @warn "max_iterations exceeded"
            return endotherm_out
        end

        # -------------------------------------------------------------------------
        # 1. Reduce insulation (piloerection)
        # -------------------------------------------------------------------------
        if (insulation_limits.dorsal.current > insulation_limits.dorsal.reference) &&
           (insulation_limits.ventral.current > insulation_limits.ventral.reference)

            insulation_limits, organism = piloerect(organism, insulation_limits)

        # -------------------------------------------------------------------------
        # 2. Uncurl (increase surface area)
        # -------------------------------------------------------------------------
        elseif shape_coefficient_b_limits.current < shape_coefficient_b_limits.max
            shape_coefficient_b_limits, organism = uncurl(organism, shape_coefficient_b_limits)

        # -------------------------------------------------------------------------
        # 3. Vasodilate (increase flesh_conductivity)
        # -------------------------------------------------------------------------
        elseif flesh_conductivity_limits.current < flesh_conductivity_limits.max
            flesh_conductivity_limits, organism = vasodilate(organism, flesh_conductivity_limits)

        # -------------------------------------------------------------------------
        # 4. Allow core temperature to rise (and possibly pant and sweat in parallel)
        # -------------------------------------------------------------------------
        elseif core_temperature_limits.current < core_temperature_limits.max
            core_temperature_limits, minimum_metabolic_heat_flow, organism = hyperthermia(
                organism, core_temperature_limits, panting_limits.cost
            )
            if simultaneous_pant(mode) && panting_limits.panting_rate.current < panting_limits.panting_rate.max
                # Pant in parallel to allowing core temperature to rise
                panting_limits, minimum_metabolic_heat_flow, organism = pant(organism, panting_limits)
            end
            if simultaneous_sweat(mode)
                # Sweat in parallel to allowing core temperature to rise and panting
                if (skin_wetness_limits.current > skin_wetness_limits.max) ||
                   (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return endotherm_out
                end
                skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
            end

        # -------------------------------------------------------------------------
        # 5. Pant to dump heat evaporatively (and possibly sweat in parallel)
        # -------------------------------------------------------------------------
        elseif panting_limits.panting_rate.current < panting_limits.panting_rate.max
            panting_limits, minimum_metabolic_heat_flow, organism = pant(organism, panting_limits)
            if simultaneous_sweat(mode)
                if (skin_wetness_limits.current > skin_wetness_limits.max) ||
                   (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return endotherm_out
                end
                # Sweat in parallel to panting
                skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
            end

        # -------------------------------------------------------------------------
        # 6. Sweat to dump heat evaporatively
        # -------------------------------------------------------------------------
        else
            if (skin_wetness_limits.current > skin_wetness_limits.max) ||
               (skin_wetness_limits.step <= 0)
                return endotherm_out
            end
            skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
        end

        endotherm_out = solve_metabolic_rate(organism, environment, skin_temperature, insulation_temperature)
        skin_temperature = endotherm_out.thermoregulation.skin_temperature
        insulation_temperature = endotherm_out.thermoregulation.insulation_temperature
        generated_heat_flow = endotherm_out.energy_flows.generated_heat_flow
    end

    return endotherm_out
end

"""
    thermoregulate(::Endotherm, ::PDEControl, organism, environment, generated_heat_flow, skin_temperature, insulation_temperature)

Run the endotherm thermoregulation loop using PDE-based control.

!!! warning
    This control strategy is not yet implemented.
"""
function thermoregulate(
    ::Endotherm,
    ::PDEControl,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    error("PDEControl thermoregulation not yet implemented")
end

# Placeholder for future implementations
function thermoregulate(
    ::Ectotherm,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    error("Ectotherm thermoregulation not yet implemented")
end

function thermoregulate(
    ::Heterotherm,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    error("Heterotherm thermoregulation not yet implemented")
end
