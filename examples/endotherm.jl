using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using FluidProperties
using ConstructionBase
using ModelParameters
using Unitful, UnitfulMoles

# simulation with example parameters and variables
# which are NicheMapR's endoR/endoR_devel defaults
shape_pars = example_shape_pars()
insulation_pars = example_insulation_pars()
radiation_pars = example_radiation_pars()
metabolism_pars = example_metabolism_pars()
evaporation_pars = example_evaporation_pars()
respiration_pars = example_respiration_pars()

# set up geometry
conduction_pars_internal = example_conduction_pars_internal()
fat = Fat(conduction_pars_internal.fat_fraction, conduction_pars_internal.ρ_fat)
mean_insulation_depth = insulation_pars.insulation_depth_dorsal * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.insulation_depth_ventral * radiation_pars.ventral_fraction
mean_fibre_diameter = insulation_pars.fibre_diameter_dorsal * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.fibre_diameter_ventral * radiation_pars.ventral_fraction
mean_fibre_density = insulation_pars.fibre_density_dorsal * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.fibre_density_ventral * radiation_pars.ventral_fraction
fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

# Create physiology traits (HeatExchangeTraits)
physiology_traits = HeatExchangeTraits(
    shape_pars,
    insulation_pars,
    example_conduction_pars_external(),
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    example_hydraulic_pars(),
    respiration_pars,
    metabolism_pars,
    example_metabolic_rate_options(),
)

# Create thermoregulation limits
core_temperature_ref = metabolism_pars.core_temperature
thermoregulation_limits = ThermoregulationLimits(;
    control=RuleBasedSequentialControl(;
        mode=CoreFirst(),
        tolerance=0.005,
        max_iterations=200,
    ),
    Q_minimum_ref=metabolism_pars.Q_metabolism,
    insulation=InsulationLimits(;
        dorsal=SteppedParameter(;
            current=insulation_pars.insulation_depth_dorsal,
            reference=insulation_pars.insulation_depth_dorsal,
            max=insulation_pars.insulation_depth_dorsal,
            step=0.0,
        ),
        ventral=SteppedParameter(;
            current=insulation_pars.insulation_depth_ventral,
            reference=insulation_pars.insulation_depth_ventral,
            max=insulation_pars.insulation_depth_ventral,
            step=0.0,
        ),
    ),
    shape_b=SteppedParameter(;
        current=shape_pars.b,
        max=5.0,
        step=0.1,
    ),
    k_flesh=SteppedParameter(;
        current=conduction_pars_internal.k_flesh,
        max=2.8u"W/m/K",
        step=0.1u"W/m/K",
    ),
    core_temperature=SteppedParameter(;
        current=core_temperature_ref,
        reference=core_temperature_ref,
        max=core_temperature_ref + 5.0u"K",
        step=0.1u"K",
    ),
    panting=PantingLimits(;
        pant=SteppedParameter(;
            current=respiration_pars.pant,
            max=15.0,
            step=0.01,
        ),
        cost=0.0u"W",
        multiplier=1.0,
        core_temperature_ref=core_temperature_ref,
    ),
    skin_wetness=SteppedParameter(;
        current=evaporation_pars.skin_wetness,
        max=0.5,
        step=0.01,
    ),
)

# Combine physiology and behavior into OrganismTraits
behavioral_traits = BehavioralTraits(;
    thermoregulation=thermoregulation_limits,
    activity=Diurnal(),
)
organism_traits = OrganismTraits(Endotherm(), physiology_traits, behavioral_traits)

organism = Organism(geometry, organism_traits)

environment_vars = example_environment_vars()
environment_pars = example_environment_pars()
environment = (; environment_pars, environment_vars)

# initial conditions
skin_temperature = metabolism_pars.core_temperature - 3.0u"K"
insulation_temperature = environment_vars.air_temperature
Q_gen = 0.0u"W"

endotherm_out = thermoregulate(
    organism,
    environment,
    Q_gen,
    skin_temperature,
    insulation_temperature,
)
thermoregulation = endotherm_out.thermoregulation
morphology = endotherm_out.morphology
energy_fluxes = endotherm_out.energy_fluxes
mass_fluxes = endotherm_out.mass_fluxes
