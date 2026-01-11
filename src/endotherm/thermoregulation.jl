function piloerect(
    insulation_depth_dorsal,
    insulation_depth_ventral,
    insulation_depth_dorsal_ref,
    insulation_depth_ventral_ref,
    insulation_step,
    organism
)
    insulation_pars = organism.traits.insulation_pars
    shape_pars = organism.traits.shape_pars
    fibre_length_dorsal = insulation_pars.fibre_length_dorsal
    fibre_length_ventral = insulation_pars.fibre_length_ventral

    # Decrement each insulation depth towards its reference
    insulation_depth_dorsal = max(
        insulation_depth_dorsal_ref,
        insulation_depth_dorsal - insulation_step * fibre_length_dorsal
    )
    insulation_depth_ventral = max(
        insulation_depth_ventral_ref,
        insulation_depth_ventral - insulation_step * fibre_length_ventral
    )

    insulation_pars = ConstructionBase.setproperties(
        insulation_pars,
        insulation_depth_dorsal = insulation_depth_dorsal,
        insulation_depth_ventral = insulation_depth_ventral
    )

    # Compute mean insulation properties
    ventral_frac = organism.traits.radiation_pars.ventral_fraction
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
    traits′ = ConstructionBase.setproperties(organism.traits, insulation_pars = insulation_pars)
    organism′ = ConstructionBase.setproperties(organism, traits = traits′, body = geometry)

    return insulation_depth_dorsal, insulation_depth_ventral, organism′
end

function uncurl(shape_b, shape_b_step, shape_b_max, organism)
    shape_pars = organism.traits.shape_pars

    # No meaning to uncurl a sphere
    if shape_pars isa Sphere
        return shape_b_max, organism
    end

    shape_b = min(shape_b + shape_b_step, shape_b_max)

    new_shape_pars = ConstructionBase.setproperties(shape_pars, b=shape_b)
    fat = BiophysicalGeometry.inner_insulation(organism.body.insulation)
    fur = BiophysicalGeometry.outer_insulation(organism.body.insulation)
    geometry = Body(new_shape_pars, CompositeInsulation(fur, fat))
    organism′ = ConstructionBase.setproperties(organism, body=geometry)

    return shape_b, organism′
end

function vasodilate(k_flesh, k_flesh_step, k_flesh_max, organism)
    k_flesh = min(k_flesh + k_flesh_step, k_flesh_max)

    new_conduction = ConstructionBase.setproperties(
        organism.traits.conduction_pars_internal, k_flesh=k_flesh
    )
    new_traits = ConstructionBase.setproperties(
        organism.traits, conduction_pars_internal=new_conduction
    )
    organism′ = ConstructionBase.setproperties(organism, traits=new_traits)

    return k_flesh, organism′
end

function hyperthermia(
    T_core, T_core_step, T_core_max, T_core_ref, Q_minimum_ref, pant_cost, organism
)
    # Update T_core
    T_core = min(T_core + T_core_step, T_core_max)

    metabolism = organism.traits.metabolism_pars
    q10mult = metabolism.q10^((ustrip(u"K", T_core - T_core_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism, T_core=T_core, Q_metabolism=Q_minimum
    )
    new_traits = ConstructionBase.setproperties(
        organism.traits, metabolism_pars=new_metabolism
    )
    organism′ = ConstructionBase.setproperties(organism, traits=new_traits)

    return T_core, Q_minimum, organism′
end

function pant(pant, pant_step, pant_max, T_core_ref, Q_minimum_ref, pant_multiplier, organism)
    pant = min(pant + pant_step, pant_max)

    pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                (pant_multiplier - 1) * Q_minimum_ref

    metabolism = organism.traits.metabolism_pars
    q10mult = metabolism.q10^((ustrip(u"K", metabolism.T_core - T_core_ref)) / 10)
    Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

    new_metabolism = ConstructionBase.setproperties(
        metabolism, Q_metabolism=Q_minimum
    )
    new_respiration = ConstructionBase.setproperties(
        organism.traits.respiration_pars, pant=pant
    )
    new_traits = ConstructionBase.setproperties(
        organism.traits,
        metabolism_pars=new_metabolism,
        respiration_pars=new_respiration
    )
    organism′ = ConstructionBase.setproperties(organism, traits=new_traits)

    return pant, pant_cost, Q_minimum, organism′
end

function sweat(skin_wetness, skin_wetness_step, skin_wetness_max, organism)
    skin_wetness = min(skin_wetness + skin_wetness_step, skin_wetness_max)

    new_evaporation = ConstructionBase.setproperties(
        organism.traits.evaporation_pars, skin_wetness=skin_wetness
    )
    new_traits = ConstructionBase.setproperties(
        organism.traits, evaporation_pars=new_evaporation
    )
    organism′ = ConstructionBase.setproperties(organism, traits=new_traits)

    return skin_wetness, organism′
end
