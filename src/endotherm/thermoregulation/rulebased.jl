# =============================================================================
# Rule-based sequential thermoregulation control.
#
# Six effector functions adjust insulation depth, body shape, tissue
# conductivity, core temperature setpoint, panting rate, and skin wetness.
# The `thermoregulate(::Endotherm, ::RuleBasedSequentialControl, ...)` method
# applies them in priority order until the heat balance is satisfied.
# =============================================================================

# ----------------------------------------------------------------------------
# Effectors
# ----------------------------------------------------------------------------

"""
    piloerect(organism::Organism, insulation_limits::InsulationLimits)

Reduce insulation depth toward reference values (flatten fur/feathers).

Returns updated `InsulationLimits` and `organism`.
"""
function piloerect(organism::Organism, insulation_limits::InsulationLimits)
    insulation_pars = HeatExchange.insulation_pars(organism)
    shape_pars = HeatExchange.shape_pars(organism)
    fibre_length_dorsal = insulation_pars.dorsal.length
    fibre_length_ventral = insulation_pars.ventral.length

    # Decrement each insulation depth towards its reference
    insulation_depth_dorsal = max(
        insulation_limits.dorsal.reference,
        insulation_limits.dorsal.current - insulation_limits.dorsal.step * fibre_length_dorsal
    )
    insulation_depth_ventral = max(
        insulation_limits.ventral.reference,
        insulation_limits.ventral.current - insulation_limits.ventral.step * fibre_length_ventral
    )

    insulation_limits = @set insulation_limits.dorsal.current = insulation_depth_dorsal
    insulation_limits = @set insulation_limits.ventral.current = insulation_depth_ventral

    # Compute mean insulation properties (ventral_fraction-weighted)
    ventral_frac = HeatExchange.radiation_pars(organism).ventral_fraction
    mean_insulation_depth = insulation_depth_dorsal * (1 - ventral_frac) +
                            insulation_depth_ventral * ventral_frac
    mean_fibre_diameter = insulation_pars.dorsal.diameter * (1 - ventral_frac) +
                          insulation_pars.ventral.diameter * ventral_frac
    mean_fibre_density = insulation_pars.dorsal.density * (1 - ventral_frac) +
                         insulation_pars.ventral.density * ventral_frac

    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    geometry = rebuild_body(shape_pars, mean_insulation_depth, mean_fibre_diameter,
                            mean_fibre_density, fat)

    organism = @set organism.traits.heat_exchange.insulation_pars.dorsal.depth = insulation_depth_dorsal
    organism = @set organism.traits.heat_exchange.insulation_pars.ventral.depth = insulation_depth_ventral
    organism = @set organism.body = geometry

    return insulation_limits, organism
end

"""
    uncurl(organism::Organism, aspect_ratio_limits::SteppedParameter)

Increase body aspect ratio (uncurl from ball to elongated).

Returns updated `SteppedParameter` and `organism`.
"""
function uncurl(organism::Organism, aspect_ratio_limits::SteppedParameter)
    shape_pars = HeatExchange.shape_pars(organism)

    # No meaning to uncurl a sphere
    if shape_pars isa Sphere
        aspect_ratio_limits = @set aspect_ratio_limits.current = aspect_ratio_limits.max
        return aspect_ratio_limits, organism
    end

    aspect_ratio_factor = min(aspect_ratio_limits.current + aspect_ratio_limits.step, aspect_ratio_limits.max)
    aspect_ratio_limits = @set aspect_ratio_limits.current = aspect_ratio_factor

    new_shape_pars = @set shape_pars.axis_ratio_b = aspect_ratio_factor
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = rebuild_body(new_shape_pars, fur, fat)
    organism = @set organism.body = geometry

    return aspect_ratio_limits, organism
end

"""
    vasodilate(organism::Organism, flesh_conductivity_limits::SteppedParameter)

Increase tissue thermal conductivity (vasodilation).

Returns updated `SteppedParameter` and `organism`.
"""
function vasodilate(organism::Organism, flesh_conductivity_limits::SteppedParameter)
    flesh_conductivity = min(flesh_conductivity_limits.current + flesh_conductivity_limits.step, flesh_conductivity_limits.max)
    flesh_conductivity_limits = @set flesh_conductivity_limits.current = flesh_conductivity

    organism = @set organism.traits.heat_exchange.conduction_pars_internal.flesh_conductivity = flesh_conductivity

    return flesh_conductivity_limits, organism
end

"""
    hyperthermia(organism::Organism, core_temperature_limits::SteppedParameter, pant_cost)

Allow core temperature to rise (hyperthermia).

Returns updated `SteppedParameter`, new minimum_heat_flow, and `organism`.
"""
function hyperthermia(organism::Organism, core_temperature_limits::SteppedParameter, pant_cost)
    minimum_heat_flow = thermoregulation(organism).minimum_heat_flow
    core_temp = min(core_temperature_limits.current + core_temperature_limits.step, core_temperature_limits.max)
    core_temperature_limits = @set core_temperature_limits.current = core_temp

    metabolism = HeatExchange.metabolism_pars(organism)
    minimum_heat_flow = (minimum_heat_flow + pant_cost) *
                q10_scale(metabolism.q10, core_temp, core_temperature_limits.reference)

    organism = @set organism.traits.heat_exchange.metabolism_pars.core_temperature = core_temp
    organism = @set organism.traits.heat_exchange.metabolism_pars.metabolic_heat_flow = minimum_heat_flow

    return core_temperature_limits, minimum_heat_flow, organism
end

"""
    pant(organism::Organism, panting_limits::PantingLimits)

Increase panting rate for evaporative cooling.

Returns updated `PantingLimits`, new minimum_heat_flow, and `organism`.
"""
function pant(organism::Organism, panting_limits::PantingLimits)
    minimum_heat_flow = thermoregulation(organism).minimum_heat_flow
    pant_rate_limits = panting_limits.pant
    pant_rate = min(pant_rate_limits.current + pant_rate_limits.step, pant_rate_limits.max)

    pant_cost = ((pant_rate - 1) / (pant_rate_limits.max + 1e-6 - 1)) *
                (panting_limits.multiplier - 1) * minimum_heat_flow

    panting_limits = @set panting_limits.pant.current = pant_rate
    panting_limits = @set panting_limits.cost = pant_cost

    metabolism = HeatExchange.metabolism_pars(organism)
    minimum_heat_flow = (minimum_heat_flow + pant_cost) *
                q10_scale(metabolism.q10, metabolism.core_temperature, panting_limits.core_temperature_ref)

    organism = @set organism.traits.heat_exchange.metabolism_pars.metabolic_heat_flow = minimum_heat_flow
    organism = @set organism.traits.heat_exchange.respiration_pars.pant = pant_rate

    return panting_limits, minimum_heat_flow, organism
end

"""
    sweat(organism::Organism, skin_wetness_limits::SteppedParameter)

Increase skin wetness for evaporative cooling (sweating).

Returns updated `SteppedParameter` and `organism`.
"""
function sweat(organism::Organism, skin_wetness_limits::SteppedParameter)
    skin_wetness = min(skin_wetness_limits.current + skin_wetness_limits.step, skin_wetness_limits.max)
    skin_wetness_limits = @set skin_wetness_limits.current = skin_wetness

    organism = @set organism.traits.heat_exchange.evaporation_pars.skin_wetness = skin_wetness

    return skin_wetness_limits, organism
end

# ----------------------------------------------------------------------------
# Sequential control loop
# ----------------------------------------------------------------------------

"""
    thermoregulate(::Endotherm, ::RuleBasedSequentialControl, organism, environment, init)

Run the endotherm thermoregulation loop using rule-based sequential control.

`init` is a NamedTuple of initial guesses (`metabolic_heat_flow`,
`skin_temperature`, `insulation_temperature`). Only `skin_temperature` and
`insulation_temperature` seed the first `solve_metabolic_rate` call;
`metabolic_heat_flow` is computed from scratch on every iteration.

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
    init::NamedTuple,
)
    (; skin_temperature, insulation_temperature) = init

    limits = thermoregulation(organism)
    endotherm_out = nothing

    (; mode, tolerance, max_iterations) = limits.control

    insulation_limits        = limits.insulation
    aspect_ratio_limits      = limits.aspect_ratio_factor
    flesh_conductivity_limits = limits.flesh_conductivity
    core_temperature_limits  = limits.core_temperature
    panting_limits           = limits.panting
    skin_wetness_limits      = limits.skin_wetness

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

    endotherm_out    = solve_metabolic_rate(organism, environment, skin_temperature, insulation_temperature)
    skin_temperature = endotherm_out.thermoregulation.skin_temperature
    insulation_temperature = endotherm_out.thermoregulation.insulation_temperature
    metabolic_heat_flow = endotherm_out.energy_flows.metabolic_heat_flow

    minimum_heat_flow = limits.minimum_heat_flow

    iteration = 0

    while metabolic_heat_flow < minimum_heat_flow * (1 - tolerance)
        iteration += 1
        if iteration > max_iterations
            @warn "max_iterations exceeded"
            return endotherm_out
        end

        if (insulation_limits.dorsal.current > insulation_limits.dorsal.reference) &&
           (insulation_limits.ventral.current > insulation_limits.ventral.reference)

            insulation_limits, organism = piloerect(organism, insulation_limits)

        elseif aspect_ratio_limits.current < aspect_ratio_limits.max
            aspect_ratio_limits, organism = uncurl(organism, aspect_ratio_limits)

        elseif flesh_conductivity_limits.current < flesh_conductivity_limits.max
            flesh_conductivity_limits, organism = vasodilate(organism, flesh_conductivity_limits)

        elseif core_temperature_limits.current < core_temperature_limits.max
            core_temperature_limits, minimum_heat_flow, organism = hyperthermia(
                organism, core_temperature_limits, panting_limits.cost
            )
            if simultaneous_pant(mode) && panting_limits.pant.current < panting_limits.pant.max
                panting_limits, minimum_heat_flow, organism = pant(organism, panting_limits)
            end
            if simultaneous_sweat(mode)
                if (skin_wetness_limits.current > skin_wetness_limits.max) ||
                   (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return endotherm_out
                end
                skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
            end

        elseif panting_limits.pant.current < panting_limits.pant.max
            panting_limits, minimum_heat_flow, organism = pant(organism, panting_limits)
            if simultaneous_sweat(mode)
                if (skin_wetness_limits.current > skin_wetness_limits.max) ||
                   (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return endotherm_out
                end
                skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
            end

        else
            if (skin_wetness_limits.current > skin_wetness_limits.max) ||
               (skin_wetness_limits.step <= 0)
                return endotherm_out
            end
            skin_wetness_limits, organism = sweat(organism, skin_wetness_limits)
        end

        endotherm_out    = solve_metabolic_rate(organism, environment, skin_temperature, insulation_temperature)
        skin_temperature = endotherm_out.thermoregulation.skin_temperature
        insulation_temperature = endotherm_out.thermoregulation.insulation_temperature
        metabolic_heat_flow = endotherm_out.energy_flows.metabolic_heat_flow
    end

    return endotherm_out
end
