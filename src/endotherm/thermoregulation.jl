"""
    piloerect(organism::Organism, insulation_limits::InsulationLimits)

Reduce insulation depth toward reference values (flatten fur/feathers).

Returns updated `InsulationLimits` and `organism`.
"""
function piloerect(organism::Organism, insulation_limits::InsulationLimits)
    phys = organism.traits.physiology
    insulation_pars = HeatExchange.insulationpars(phys)
    shape_pars = HeatExchange.shapepars(phys)
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
    dorsal_limits′ = setproperties(insulation_limits.dorsal; current=insulation_depth_dorsal)
    ventral_limits′ = setproperties(insulation_limits.ventral; current=insulation_depth_ventral)
    insulation_limits′ = InsulationLimits(; dorsal=dorsal_limits′, ventral=ventral_limits′)

    # Update organism insulation parameters (nested FibreProperties)
    dorsal′ = setproperties(insulation_pars.dorsal; depth=insulation_depth_dorsal)
    ventral′ = setproperties(insulation_pars.ventral; depth=insulation_depth_ventral)
    insulation_pars′ = setproperties(insulation_pars; dorsal=dorsal′, ventral=ventral′)

    # Compute mean insulation properties
    ventral_frac = HeatExchange.radiationpars(phys).ventral_fraction
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

    # Update organism immutably
    phys′ = setproperties(phys; insulation_pars=insulation_pars′)
    traits′ = setproperties(organism.traits; physiology=phys′)
    organism′ = setproperties(organism; traits=traits′, body=geometry)

    return insulation_limits′, organism′
end

"""
    uncurl(organism::Organism, shape_b_limits::SteppedParameter)

Increase body shape parameter (uncurl from ball to elongated).

Returns updated `SteppedParameter` and `organism`.
"""
function uncurl(organism::Organism, shape_b_limits::SteppedParameter)
    phys = organism.traits.physiology
    shape_pars = HeatExchange.shapepars(phys)

    # No meaning to uncurl a sphere
    if shape_pars isa Sphere
        shape_b_limits′ = setproperties(shape_b_limits; current=shape_b_limits.max)
        return shape_b_limits′, organism
    end

    shape_b = min(shape_b_limits.current + shape_b_limits.step, shape_b_limits.max)
    shape_b_limits′ = setproperties(shape_b_limits; current=shape_b)

    new_shape_pars = setproperties(shape_pars; b=shape_b)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = Body(new_shape_pars, CompositeInsulation(fur, fat))
    organism′ = setproperties(organism; body=geometry)

    return shape_b_limits′, organism′
end

"""
    vasodilate(organism::Organism, k_flesh_limits::SteppedParameter)

Increase tissue thermal conductivity (vasodilation).

Returns updated `SteppedParameter` and `organism`.
"""
function vasodilate(organism::Organism, k_flesh_limits::SteppedParameter)
    phys = organism.traits.physiology
    k_flesh = min(k_flesh_limits.current + k_flesh_limits.step, k_flesh_limits.max)
    k_flesh_limits′ = setproperties(k_flesh_limits; current=k_flesh)

    conduction_pars_internal = setproperties(
        HeatExchange.conductionpars_internal(phys); k_flesh
    )
    phys′ = setproperties(phys; conduction_pars_internal)
    traits′ = setproperties(organism.traits; physiology=phys′)
    organism′ = setproperties(organism; traits=traits′)

    return k_flesh_limits′, organism′
end

"""
    hyperthermia(organism::Organism, T_core_limits::SteppedParameter, pant_cost)

Allow core temperature to rise (hyperthermia).

Returns updated `SteppedParameter`, new Q_minimum, and `organism`.
"""
function hyperthermia(organism::Organism, T_core_limits::SteppedParameter, pant_cost)
    phys = organism.traits.physiology
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    T_core = min(T_core_limits.current + T_core_limits.step, T_core_limits.max)
    T_core_limits′ = setproperties(T_core_limits; current=T_core)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", T_core - T_core_limits.reference)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = setproperties(
        metabolism; T_core=T_core, Q_metabolism=Q_minimum
    )
    phys′ = setproperties(phys; metabolism_pars=new_metabolism)
    traits′ = setproperties(organism.traits; physiology=phys′)
    organism′ = setproperties(organism; traits=traits′)

    return T_core_limits′, Q_minimum, organism′
end

"""
    pant(organism::Organism, panting_limits::PantingLimits)

Increase panting rate for evaporative cooling.

Returns updated `PantingLimits`, new Q_minimum, and `organism`.
"""
function pant(organism::Organism, panting_limits::PantingLimits)
    phys = organism.traits.physiology
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    pant_rate_limits = panting_limits.pant
    pant_rate = min(pant_rate_limits.current + pant_rate_limits.step, pant_rate_limits.max)

    pant_cost = ((pant_rate - 1) / (pant_rate_limits.max + 1e-6 - 1)) *
                (panting_limits.multiplier - 1) * Q_minimum_ref

    pant_rate_limits′ = setproperties(pant_rate_limits; current=pant_rate)
    panting_limits′ = setproperties(panting_limits; pant=pant_rate_limits′, cost=pant_cost)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", metabolism.T_core - panting_limits.T_core_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = setproperties(metabolism; Q_metabolism=Q_minimum)
    new_respiration = setproperties(
        HeatExchange.respirationpars(phys); pant=pant_rate
    )
    phys′ = setproperties(phys;
        metabolism_pars=new_metabolism,
        respiration_pars=new_respiration
    )
    traits′ = setproperties(organism.traits; physiology=phys′)
    organism′ = setproperties(organism; traits=traits′)

    return panting_limits′, Q_minimum, organism′
end

"""
    sweat(organism::Organism, skin_wetness_limits::SteppedParameter)

Increase skin wetness for evaporative cooling (sweating).

Returns updated `SteppedParameter` and `organism`.
"""
function sweat(organism::Organism, skin_wetness_limits::SteppedParameter)
    phys = organism.traits.physiology
    skin_wetness = min(skin_wetness_limits.current + skin_wetness_limits.step, skin_wetness_limits.max)
    skin_wetness_limits′ = setproperties(skin_wetness_limits; current=skin_wetness)

    new_evaporation = setproperties(HeatExchange.evaporationpars(phys); skin_wetness)
    phys′ = setproperties(phys; evaporation_pars=new_evaporation)
    traits′ = setproperties(organism.traits; physiology=phys′)
    organism′ = setproperties(organism; traits=traits′)

    return skin_wetness_limits′, organism′
end
