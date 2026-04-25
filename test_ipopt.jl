using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

traits    = example_organism_traits()
ins_pars  = example_insulation_pars()
int_cond  = example_conduction_pars_internal()
rad_pars  = example_radiation_pars()
shape_pars = example_shape_pars()

pven = rad_pars.ventral_fraction
fat  = Fat(int_cond.fat_fraction, int_cond.fat_density)
mean_depth   = ins_pars.dorsal.depth    * (1-pven) + ins_pars.ventral.depth    * pven
mean_diam    = ins_pars.dorsal.diameter * (1-pven) + ins_pars.ventral.diameter * pven
mean_density = ins_pars.dorsal.density  * (1-pven) + ins_pars.ventral.density  * pven
fur  = Fur(mean_depth, mean_diam, mean_density)
body = Body(shape_pars, CompositeInsulation(fur, fat))

organism = Organism(body, traits)

env_pars = example_environment_pars()
env_vars = example_environment_vars()
environment = (; environment_pars=env_pars, environment_vars=env_vars)

# Reference solution from rule-based solver
ref = thermoregulate(organism, environment, 0.0u"W", 293.15u"K", 293.15u"K")
println("RuleBasedSequentialControl:")
println("  T_core  = ", ref.thermoregulation.core_temperature)
println("  T_skin  = ", ref.thermoregulation.skin_temperature)
println("  T_ins   = ", ref.thermoregulation.insulation_temperature)
println("  Q_gen   = ", ref.energy_flows.generated_heat_flow)

Q_init    = ref.energy_flows.generated_heat_flow
T_sk_init = ref.thermoregulation.skin_temperature
T_in_init = ref.thermoregulation.insulation_temperature

println("\nRunning IPOPTControl...")
# IPOPTControl dispatch reads bounds from organism's thermoregulation limits;
# the control field itself doesn't matter for this dispatch path.
sol = thermoregulate(Endotherm(), IPOPTControl(), organism, environment, Q_init, T_sk_init, T_in_init)
println("IPOPTControl:")
println("  T_core  = ", sol.thermoregulation.core_temperature)
println("  T_skin  = ", sol.thermoregulation.skin_temperature)
println("  T_ins   = ", sol.thermoregulation.insulation_temperature)
println("  Q_gen   = ", sol.energy_flows.generated_heat_flow)
println("  balance = ", sol.energy_flows.balance)
