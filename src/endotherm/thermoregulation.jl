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

    # Update limits
    insulation_limits = @set insulation_limits.dorsal.current = insulation_depth_dorsal
    insulation_limits = @set insulation_limits.ventral.current = insulation_depth_ventral

    # Compute mean insulation properties
    ventral_frac = HeatExchange.radiation_pars(organism).ventral_fraction
    mean_insulation_depth = insulation_depth_dorsal * (1 - ventral_frac) +
                            insulation_depth_ventral * ventral_frac
    mean_fibre_diameter = insulation_pars.dorsal.diameter * (1 - ventral_frac) +
                          insulation_pars.ventral.diameter * ventral_frac
    mean_fibre_density = insulation_pars.dorsal.density * (1 - ventral_frac) +
                         insulation_pars.ventral.density * ventral_frac

    # Rebuild fur and geometry
    fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    geometry = Body(shape_pars, CompositeInsulation(fur, fat))

    # Update organism
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

    new_shape_pars = @set shape_pars.aspect_ratio_b = aspect_ratio_factor
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = Body(new_shape_pars, CompositeInsulation(fur, fat))
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

Returns updated `SteppedParameter`, new Q_minimum, and `organism`.
"""
function hyperthermia(organism::Organism, core_temperature_limits::SteppedParameter, pant_cost)
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    core_temp = min(core_temperature_limits.current + core_temperature_limits.step, core_temperature_limits.max)
    core_temperature_limits = @set core_temperature_limits.current = core_temp

    metabolism = HeatExchange.metabolism_pars(organism)
    q10mult = metabolism.q10^((ustrip(u"K", core_temp - core_temperature_limits.reference)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    organism = @set organism.traits.heat_exchange.metabolism_pars.core_temperature = core_temp
    organism = @set organism.traits.heat_exchange.metabolism_pars.metabolic_heat_flow = Q_minimum

    return core_temperature_limits, Q_minimum, organism
end

"""
    pant(organism::Organism, panting_limits::PantingLimits)

Increase panting rate for evaporative cooling.

Returns updated `PantingLimits`, new Q_minimum, and `organism`.
"""
function pant(organism::Organism, panting_limits::PantingLimits)
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    pant_rate_limits = panting_limits.pant
    pant_rate = min(pant_rate_limits.current + pant_rate_limits.step, pant_rate_limits.max)

    pant_cost = ((pant_rate - 1) / (pant_rate_limits.max + 1e-6 - 1)) *
                (panting_limits.multiplier - 1) * Q_minimum_ref

    panting_limits = @set panting_limits.pant.current = pant_rate
    panting_limits = @set panting_limits.cost = pant_cost

    metabolism = HeatExchange.metabolism_pars(organism)
    q10mult = metabolism.q10^((ustrip(u"K", metabolism.core_temperature - panting_limits.core_temperature_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    organism = @set organism.traits.heat_exchange.metabolism_pars.metabolic_heat_flow = Q_minimum
    organism = @set organism.traits.heat_exchange.respiration_pars.pant = pant_rate

    return panting_limits, Q_minimum, organism
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
