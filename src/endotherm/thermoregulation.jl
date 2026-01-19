# Helper to get physiology traits from organism (works with both HeatExchangeTraits and OrganismTraits)
_physiology(t::OrganismTraits) = t.physiology
_physiology(t::HeatExchange.HeatExchangeTraits) = t

# Helper to update physiology within traits
function _update_physiology(t::OrganismTraits, new_phys)
    ConstructionBase.setproperties(t; physiology=new_phys)
end
function _update_physiology(::HeatExchange.HeatExchangeTraits, new_phys)
    new_phys
end

"""
    piloerect(organism::Organism, insulation_limits::InsulationLimits)

Reduce insulation depth toward reference values (flatten fur/feathers).

Returns updated `InsulationLimits` and `organism`.
"""
function piloerect(organism::Organism, insulation_limits::InsulationLimits)
    phys = _physiology(organism.traits)
    insulation_pars = HeatExchange.insulationpars(phys)
    shape_pars = HeatExchange.shapepars(phys)
    fibre_length_dorsal = insulation_pars.fibre_length_dorsal
    fibre_length_ventral = insulation_pars.fibre_length_ventral

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
    dorsal_limits′ = ConstructionBase.setproperties(insulation_limits.dorsal; current=insulation_depth_dorsal)
    ventral_limits′ = ConstructionBase.setproperties(insulation_limits.ventral; current=insulation_depth_ventral)
    insulation_limits′ = InsulationLimits(; dorsal=dorsal_limits′, ventral=ventral_limits′)

    # Update organism insulation parameters
    insulation_pars′ = ConstructionBase.setproperties(
        insulation_pars;
        insulation_depth_dorsal=insulation_depth_dorsal,
        insulation_depth_ventral=insulation_depth_ventral
    )

    # Compute mean insulation properties
    ventral_frac = HeatExchange.radiationpars(phys).ventral_fraction
    mean_insulation_depth = insulation_depth_dorsal * (1 - ventral_frac) +
                            insulation_depth_ventral * ventral_frac
    mean_fibre_diameter = insulation_pars.fibre_diameter_dorsal * (1 - ventral_frac) +
                          insulation_pars.fibre_diameter_ventral * ventral_frac
    mean_fibre_density = insulation_pars.fibre_density_dorsal * (1 - ventral_frac) +
                         insulation_pars.fibre_density_ventral * ventral_frac

    # Rebuild fur and geometry
    fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    geometry = Body(shape_pars, CompositeInsulation(fur, fat))

    # Update organism immutably
    phys′ = ConstructionBase.setproperties(phys; insulation_pars=insulation_pars′)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′, body=geometry)

    return insulation_limits′, organism′
end

"""
    uncurl(organism::Organism, shape_b_limits::SteppedParameter)

Increase body shape parameter (uncurl from ball to elongated).

Returns updated `SteppedParameter` and `organism`.
"""
function uncurl(organism::Organism, shape_b_limits::SteppedParameter)
    phys = _physiology(organism.traits)
    shape_pars = HeatExchange.shapepars(phys)

    # No meaning to uncurl a sphere
    if shape_pars isa Sphere
        shape_b_limits′ = ConstructionBase.setproperties(shape_b_limits; current=shape_b_limits.max)
        return shape_b_limits′, organism
    end

    shape_b = min(shape_b_limits.current + shape_b_limits.step, shape_b_limits.max)
    shape_b_limits′ = ConstructionBase.setproperties(shape_b_limits; current=shape_b)

    new_shape_pars = ConstructionBase.setproperties(shape_pars; b=shape_b)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = Body(new_shape_pars, CompositeInsulation(fur, fat))
    organism′ = ConstructionBase.setproperties(organism; body=geometry)

    return shape_b_limits′, organism′
end

"""
    vasodilate(organism::Organism, k_flesh_limits::SteppedParameter)

Increase tissue thermal conductivity (vasodilation).

Returns updated `SteppedParameter` and `organism`.
"""
function vasodilate(organism::Organism, k_flesh_limits::SteppedParameter)
    phys = _physiology(organism.traits)
    k_flesh = min(k_flesh_limits.current + k_flesh_limits.step, k_flesh_limits.max)
    k_flesh_limits′ = ConstructionBase.setproperties(k_flesh_limits; current=k_flesh)

    new_conduction = ConstructionBase.setproperties(
        HeatExchange.conductionpars_internal(phys); k_flesh=k_flesh
    )
    phys′ = ConstructionBase.setproperties(phys; conduction_pars_internal=new_conduction)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return k_flesh_limits′, organism′
end

"""
    hyperthermia(organism::Organism, T_core_limits::SteppedParameter, pant_cost)

Allow core temperature to rise (hyperthermia).

Returns updated `SteppedParameter`, new Q_minimum, and `organism`.
"""
function hyperthermia(organism::Organism, T_core_limits::SteppedParameter, pant_cost)
    phys = _physiology(organism.traits)
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    T_core = min(T_core_limits.current + T_core_limits.step, T_core_limits.max)
    T_core_limits′ = ConstructionBase.setproperties(T_core_limits; current=T_core)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", T_core - T_core_limits.reference)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism; T_core=T_core, Q_metabolism=Q_minimum
    )
    phys′ = ConstructionBase.setproperties(phys; metabolism_pars=new_metabolism)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return T_core_limits′, Q_minimum, organism′
end

"""
    pant(organism::Organism, panting_limits::PantingLimits)

Increase panting rate for evaporative cooling.

Returns updated `PantingLimits`, new Q_minimum, and `organism`.
"""
function pant(organism::Organism, panting_limits::PantingLimits)
    phys = _physiology(organism.traits)
    Q_minimum_ref = thermoregulation(organism).Q_minimum_ref
    pant_rate_limits = panting_limits.pant
    pant_rate = min(pant_rate_limits.current + pant_rate_limits.step, pant_rate_limits.max)

    pant_cost = ((pant_rate - 1) / (pant_rate_limits.max + 1e-6 - 1)) *
                (panting_limits.multiplier - 1) * Q_minimum_ref

    pant_rate_limits′ = ConstructionBase.setproperties(pant_rate_limits; current=pant_rate)
    panting_limits′ = ConstructionBase.setproperties(panting_limits; pant=pant_rate_limits′, cost=pant_cost)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", metabolism.T_core - panting_limits.T_core_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism; Q_metabolism=Q_minimum
    )
    new_respiration = ConstructionBase.setproperties(
        HeatExchange.respirationpars(phys); pant=pant_rate
    )
    phys′ = ConstructionBase.setproperties(
        phys;
        metabolism_pars=new_metabolism,
        respiration_pars=new_respiration
    )
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return panting_limits′, Q_minimum, organism′
end

"""
    sweat(organism::Organism, skin_wetness_limits::SteppedParameter)

Increase skin wetness for evaporative cooling (sweating).

Returns updated `SteppedParameter` and `organism`.
"""
function sweat(organism::Organism, skin_wetness_limits::SteppedParameter)
    phys = _physiology(organism.traits)
    skin_wetness = min(skin_wetness_limits.current + skin_wetness_limits.step, skin_wetness_limits.max)
    skin_wetness_limits′ = ConstructionBase.setproperties(skin_wetness_limits; current=skin_wetness)

    new_evaporation = ConstructionBase.setproperties(
        HeatExchange.evaporationpars(phys); skin_wetness=skin_wetness
    )
    phys′ = ConstructionBase.setproperties(phys; evaporation_pars=new_evaporation)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return skin_wetness_limits′, organism′
end
