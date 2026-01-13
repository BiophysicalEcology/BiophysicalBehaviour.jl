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

traits = HeatExchangeTraits(
    shape_pars,
    insulation_pars,
    example_conduction_pars_external(),
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    example_evaporation_pars(),
    example_hydraulic_pars(),
    example_respiration_pars(),
    metabolism_pars,
)

organism = Organism(geometry, traits)

environment_vars = example_environment_vars()
environment_pars = example_environment_pars()
environment = (; environment_pars, environment_vars)

# initial conditions
T_skin = metabolism_pars.T_core - 3.0u"K"
T_insulation = environment_vars.T_air
Q_gen = 0.0u"W"

thermoregulation_limits = example_thermoregulation_limits()

model_pars = example_model_pars()

endotherm_out = endotherm_thermoregulation_original(
    Q_gen,
    T_skin,
    T_insulation,
    organism,
    thermoregulation_limits,
    environment,
    model_pars
);
thermoregulation = endotherm_out.thermoregulation
morphology = endotherm_out.morphology
energy_fluxes = endotherm_out.energy_fluxes
mass_fluxes = endotherm_out.mass_fluxes
