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
    piloerect(limits::InsulationLimits, organism::Organism)

Reduce insulation depth toward reference values (flatten fur/feathers).

Returns updated `InsulationLimits` and `organism`.
"""
function piloerect(limits::InsulationLimits, organism::Organism)
    phys = _physiology(organism.traits)
    insulation_pars = HeatExchange.insulationpars(phys)
    shape_pars = HeatExchange.shapepars(phys)
    fibre_length_dorsal = insulation_pars.fibre_length_dorsal
    fibre_length_ventral = insulation_pars.fibre_length_ventral

    # Decrement each insulation depth towards its reference
    new_dorsal_current = max(
        limits.dorsal.reference,
        limits.dorsal.current - limits.dorsal.step * fibre_length_dorsal
    )
    new_ventral_current = max(
        limits.ventral.reference,
        limits.ventral.current - limits.ventral.step * fibre_length_ventral
    )

    # Update limits
    new_dorsal = ConstructionBase.setproperties(limits.dorsal; current=new_dorsal_current)
    new_ventral = ConstructionBase.setproperties(limits.ventral; current=new_ventral_current)
    limits′ = InsulationLimits(; dorsal=new_dorsal, ventral=new_ventral)

    # Update organism insulation parameters
    insulation_pars′ = ConstructionBase.setproperties(
        insulation_pars;
        insulation_depth_dorsal=new_dorsal_current,
        insulation_depth_ventral=new_ventral_current
    )

    # Compute mean insulation properties
    ventral_frac = HeatExchange.radiationpars(phys).ventral_fraction
    mean_insulation_depth = new_dorsal_current * (1 - ventral_frac) +
                            new_ventral_current * ventral_frac
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

    return limits′, organism′
end

"""
    uncurl(limits::SteppedParameter, organism::Organism)

Increase body shape parameter (uncurl from ball to elongated).

Returns updated `SteppedParameter` and `organism`.
"""
function uncurl(limits::SteppedParameter, organism::Organism)
    phys = _physiology(organism.traits)
    shape_pars = HeatExchange.shapepars(phys)

    # No meaning to uncurl a sphere
    if shape_pars isa Sphere
        limits′ = ConstructionBase.setproperties(limits; current=limits.max)
        return limits′, organism
    end

    new_current = min(limits.current + limits.step, limits.max)
    limits′ = ConstructionBase.setproperties(limits; current=new_current)

    new_shape_pars = ConstructionBase.setproperties(shape_pars; b=new_current)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = Body(new_shape_pars, CompositeInsulation(fur, fat))
    organism′ = ConstructionBase.setproperties(organism; body=geometry)

    return limits′, organism′
end

"""
    vasodilate(limits::SteppedParameter, organism::Organism)

Increase tissue thermal conductivity (vasodilation).

Returns updated `SteppedParameter` and `organism`.
"""
function vasodilate(limits::SteppedParameter, organism::Organism)
    phys = _physiology(organism.traits)
    new_current = min(limits.current + limits.step, limits.max)
    limits′ = ConstructionBase.setproperties(limits; current=new_current)

    new_conduction = ConstructionBase.setproperties(
        HeatExchange.conductionpars_internal(phys); k_flesh=new_current
    )
    phys′ = ConstructionBase.setproperties(phys; conduction_pars_internal=new_conduction)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return limits′, organism′
end

"""
    hyperthermia(limits::SteppedParameter, Q_minimum_ref, pant_cost, organism::Organism)

Allow core temperature to rise (hyperthermia).

Returns updated `SteppedParameter`, new Q_minimum, and `organism`.
"""
function hyperthermia(limits::SteppedParameter, Q_minimum_ref, pant_cost, organism::Organism)
    phys = _physiology(organism.traits)
    new_current = min(limits.current + limits.step, limits.max)
    limits′ = ConstructionBase.setproperties(limits; current=new_current)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", new_current - limits.reference)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism; T_core=new_current, Q_metabolism=Q_minimum
    )
    phys′ = ConstructionBase.setproperties(phys; metabolism_pars=new_metabolism)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return limits′, Q_minimum, organism′
end

"""
    pant(limits::PantingLimits, Q_minimum_ref, organism::Organism)

Increase panting rate for evaporative cooling.

Returns updated `PantingLimits`, new Q_minimum, and `organism`.
"""
function pant(limits::PantingLimits, Q_minimum_ref, organism::Organism)
    phys = _physiology(organism.traits)
    pant_param = limits.pant
    new_pant = min(pant_param.current + pant_param.step, pant_param.max)

    pant_cost = ((new_pant - 1) / (pant_param.max + 1e-6 - 1)) *
                (limits.multiplier - 1) * Q_minimum_ref

    pant_param′ = ConstructionBase.setproperties(pant_param; current=new_pant)
    limits′ = ConstructionBase.setproperties(limits; pant=pant_param′, cost=pant_cost)

    metabolism = HeatExchange.metabolismpars(phys)
    q10mult = metabolism.q10^((ustrip(u"K", metabolism.T_core - limits.T_core_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism; Q_metabolism=Q_minimum
    )
    new_respiration = ConstructionBase.setproperties(
        HeatExchange.respirationpars(phys); pant=new_pant
    )
    phys′ = ConstructionBase.setproperties(
        phys;
        metabolism_pars=new_metabolism,
        respiration_pars=new_respiration
    )
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return limits′, Q_minimum, organism′
end

"""
    sweat(limits::SteppedParameter, organism::Organism)

Increase skin wetness for evaporative cooling (sweating).

Returns updated `SteppedParameter` and `organism`.
"""
function sweat(limits::SteppedParameter, organism::Organism)
    phys = _physiology(organism.traits)
    new_current = min(limits.current + limits.step, limits.max)
    limits′ = ConstructionBase.setproperties(limits; current=new_current)

    new_evaporation = ConstructionBase.setproperties(
        HeatExchange.evaporationpars(phys); skin_wetness=new_current
    )
    phys′ = ConstructionBase.setproperties(phys; evaporation_pars=new_evaporation)
    traits′ = _update_physiology(organism.traits, phys′)
    organism′ = ConstructionBase.setproperties(organism; traits=traits′)

    return limits′, organism′
end
