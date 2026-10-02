# A bird: rules and optimisation

The budgerigar, *Melopsittacus undulatus*, is a 34 g parrot of the Australian arid zone. Weathers and
Schoenbaechler (1976) measured its metabolic rate, body temperature and evaporative water loss in a chamber from
0 to 45 °C, and Kearney et al. (2021) used those measurements to test the NicheMapR endotherm model. This
tutorial runs both controllers of this package on the same bird, against the same data.

```@setup budgerigar
using Main.FigureHelpers
using CairoMakie
```

## The bird

An ellipsoid of 33.7 g, with feathers 23 mm long lying about 6 mm deep at rest, defending 38 °C. Its basal
metabolic rate comes from the allometry of McKechnie and Wolf (2004).

```@example budgerigar
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, FluidProperties, Unitful

shape = example_ellipsoid_shape_pars(; mass = 33.7u"g", axis_ratio_b = 1.1, axis_ratio_c = 1.1)
plumage = example_insulation_pars(;
    fibre_diameter_dorsal = 30.0u"μm", fibre_diameter_ventral = 30.0u"μm",
    fibre_length_dorsal = 23.1u"mm", fibre_length_ventral = 22.7u"mm",
    insulation_depth_dorsal = 5.9u"mm", insulation_depth_ventral = 5.7u"mm",
    insulation_depth_compressed = 5.7u"mm",
    fibre_density_dorsal = 5000e4u"1/m^2", fibre_density_ventral = 5000e4u"1/m^2",
    insulation_reflectance_dorsal = 0.248, insulation_reflectance_ventral = 0.351,
)
basal = metabolic_rate(McKechnieWolf(), shape.mass)
core, core_max = u"K"(38.0u"°C"), u"K"(43.0u"°C")
physiology_traits = example_heat_exchange_traits(;
    shape_pars = shape,
    insulation_pars = plumage,
    respiration_pars = example_respiration_pars(; oxygen_extraction_efficiency = 0.25, exhaled_temperature_offset = 5.0u"K"),
    evaporation_pars = example_evaporation_pars(; skin_wetness = 0.005),
    metabolism_pars = example_metabolism_pars(; core_temperature = core, q10 = 2.0, metabolic_heat_flow = basal),
)
feathers = FibrousLayer((plumage.dorsal.depth + plumage.ventral.depth) / 2, plumage.dorsal.diameter, plumage.dorsal.density)
bird_body = Body(shape, CompositeInsulation(feathers, FatLayer(0.0, 901.0u"kg/m^3")))
basal
```

What it can do. Its feathers can be raised to 70 % of their length. It can stretch out, vasodilate, let its core
warm by 5 °C and pant up to 15 times its resting ventilation. A bird has no sweat glands, but water does
evaporate through its skin, and up to 5 % of the skin is allowed to be wet.

```@example budgerigar
function budgerigar(control; fluffing = true, weights...)
    raised(side) = fluffing ? 0.7 * side.length : side.depth
    limits = ThermoregulationLimits(;
        control,
        minimum_heat_flow = basal,
        insulation = InsulationLimits(;
            dorsal = SteppedParameter(; current = raised(plumage.dorsal), reference = plumage.dorsal.depth,
                                      max = raised(plumage.dorsal), step = fluffing ? 0.1 : 0.0),
            ventral = SteppedParameter(; current = raised(plumage.ventral), reference = plumage.ventral.depth,
                                       max = raised(plumage.ventral), step = fluffing ? 0.1 : 0.0),
        ),
        axis_ratio_factor = SteppedParameter(; current = 1.1, max = fluffing ? 5.0 : 1.1, step = 0.1),
        flesh_conductivity = SteppedParameter(; current = 0.9u"W/m/K", max = 2.8u"W/m/K", step = 0.1u"W/m/K"),
        core_temperature = SteppedParameter(; current = core, reference = core, max = core_max, step = 0.1u"K"),
        panting = PantingLimits(; pant = SteppedParameter(; current = 1.0, max = 15.0, step = 0.01),
                                multiplier = 1.0, core_temperature_ref = core),
        skin_wetness = SteppedParameter(; current = 0.005, max = 0.05, step = 0.0025),
        weights...,
    )
    Organism(bird_body, OrganismTraits(Endotherm(), physiology_traits, BehavioralTraits(; thermoregulation = limits)))
end
nothing # hide
```

`fluffing = false` gives a bird that cannot change its plumage or its posture. It is there because the
optimiser cannot change them either, see [Thermoregulation by optimisation](../manual/optimisation.md), and a
fair comparison of the two controllers needs a rule-based bird with the same abilities.

## The chamber

Still air, in the dark. The humidity follows the experiment: 15 % relative humidity below 30 °C, and above it
the vapour density of air at 40 °C and 30 %.

```@example budgerigar
vapour_density = wet_air_properties(40.0u"°C", 0.3, 101325.0u"Pa").vapour_density
humidity(T) = T < u"K"(30.0u"°C") ? 0.15 : min(1.0, vapour_density / wet_air_properties(T, 1.0, 101325.0u"Pa").vapour_density)
chamber(air_temperature) = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(; air_temperature, relative_humidity = humidity(air_temperature), wind_speed = 0.1u"m/s"),
)
air = [u"K"(T * u"°C") for T in 0.0:1.0:46.0]
nothing # hide
```

## By rules

The mode is `CorePantingSweatingFirst`: the bird pants, and wets its skin, as its core temperature rises.

```@example budgerigar
function sweep_rules(bird)
    map(air) do air_temperature
        init = (; metabolic_heat_flow = 0.0u"W", skin_temperature = core - 3.0u"K", insulation_temperature = air_temperature)
        thermoregulate(bird, chamber(air_temperature), init)
    end
end
rule_control = RuleBasedSequentialControl(; mode = CorePantingSweatingFirst())
rules = sweep_rules(budgerigar(rule_control))
rules_fixed = sweep_rules(budgerigar(rule_control; fluffing = false))
nothing # hide
```

## By optimisation

The weights say what this bird minds. Its core temperature is allowed to drift, as it does in the data, so the
weight on it is low. Panting is made costly relative to cutaneous evaporation, and a rise in metabolic rate
more costly still.

```@example budgerigar
weights = (; core_temperature_weight = 0.1, panting_weight = 5.0, skin_wetness_weight = 0.1,
           flesh_conductivity_weight = 0.05, metabolic_heat_weight = 10.0)
optimising_bird = budgerigar(IPOPTControl(); weights...)

function sweep_optimiser(bird)
    control = control_strategy(bird)
    init = (; metabolic_heat_flow = basal, skin_temperature = core - 3.0u"K", insulation_temperature = first(air))
    cache = IPOPTSolverCache(control, bird, chamber(first(air)), init)
    map(air) do air_temperature
        thermoregulate(Endotherm(), control, bird, chamber(air_temperature), init; cache)
    end
end
optimised = sweep_optimiser(optimising_bird)

is_feasible(result) = all(result.parts) do part
    abs(part.flows.residual_energy_balance - part.flows.residual_internal_conduction) < 0.001u"W" &&
        abs(part.flows.residual_skin_temperature) < 0.01u"K"
end
all(is_feasible, optimised)
```

## Against the measurements

```@example budgerigar
using CSV, DataFrames
data_directory = joinpath(pkgdir(BiophysicalBehaviour), "test", "data", "budgerigar")
observed(file) = CSV.read(joinpath(data_directory, file), DataFrame)
metabolism_observed = observed("Weathers1976Fig1.csv")
temperature_observed = observed("Weathers1976Fig2.csv")
water_observed = observed("Weathers1976Fig3.csv")

mass_g = ustrip(u"g", shape.mass)
observed_watts = [ustrip(u"W", HeatExchange.O2_to_Joules(Typical(), (v * mass_g)u"ml/hr", 0.8)) for v in metabolism_observed.mlO2gh]
evaporation_of(result) = sum(part.flows.skin_evaporation_heat_flow for part in result.parts) + result.parts.body.flows.respiration_heat_flow
x = ustrip.(u"°C", air) # hide
fig = Figure(size = (760, 620)) # hide
grey, blue, orange = :grey35, :steelblue, :darkorange # hide
ax = Axis(fig[1, 1]; ylabel = "Metabolic rate (W)") # hide
scatter!(ax, metabolism_observed.Tair, observed_watts; color = (:black, 0.35), markersize = 6, label = "observed") # hide
lines!(ax, x, [watts(r.energy_flows.metabolic_heat_flow) for r in rules]; color = grey, linewidth = 2, label = "rules") # hide
lines!(ax, x, [watts(r.energy_flows.metabolic_heat_flow) for r in rules_fixed]; color = blue, linewidth = 2, label = "rules, plumage fixed") # hide
lines!(ax, x, [watts(r.metabolic_heat_flow) for r in optimised]; color = orange, linewidth = 2, linestyle = :dash, label = "optimisation") # hide
axislegend(ax; position = :rt, labelsize = 10) # hide
ax = Axis(fig[1, 2]; ylabel = "Core temperature (°C)") # hide
scatter!(ax, temperature_observed.Tair, temperature_observed.Tb; color = (:black, 0.35), markersize = 6) # hide
lines!(ax, x, [ustrip(u"°C", r.thermoregulation.core_temperature) for r in rules]; color = grey, linewidth = 2) # hide
lines!(ax, x, [ustrip(u"°C", r.thermoregulation.core_temperature) for r in rules_fixed]; color = blue, linewidth = 2) # hide
lines!(ax, x, [ustrip(u"°C", r.core_temperature) for r in optimised]; color = orange, linewidth = 2, linestyle = :dash) # hide
ax = Axis(fig[2, 1]; xlabel = "Air temperature (°C)", ylabel = "Evaporative water loss (g/h)") # hide
scatter!(ax, water_observed.Tair, water_observed.mgH2Ogh .* mass_g ./ 1000; color = (:black, 0.35), markersize = 6) # hide
lines!(ax, x, [ustrip(u"g/hr", r.mass_flows.m_evap) for r in rules]; color = grey, linewidth = 2) # hide
lines!(ax, x, [ustrip(u"g/hr", r.mass_flows.m_evap) for r in rules_fixed]; color = blue, linewidth = 2) # hide
latent = HeatExchange.enthalpy_of_vaporisation.(air) # hide
lines!(ax, x, [ustrip(u"g/hr", evaporation_of(r) / L) for (r, L) in zip(optimised, latent)]; color = orange, linewidth = 2, linestyle = :dash) # hide
ax = Axis(fig[2, 2]; xlabel = "Air temperature (°C)", ylabel = "Panting multiplier") # hide
lines!(ax, x, [r.thermoregulation.pant for r in rules]; color = grey, linewidth = 2) # hide
lines!(ax, x, [r.thermoregulation.pant for r in rules_fixed]; color = blue, linewidth = 2) # hide
lines!(ax, x, [r.panting_rate for r in optimised]; color = orange, linewidth = 2, linestyle = :dash) # hide
fig # hide
```

For the optimiser, the water loss plotted is the latent heat lost from the skin and in breathing divided by the
latent heat of vaporisation. Its heat lost in breathing includes a small sensible part, so the line is a slight
overestimate.

Three things can be read from the figure.

**In the cold, plumage is most of the answer.** The rule-based bird with its feathers raised to 16 mm needs
about a quarter less heat than the bird whose feathers stay at 6 mm, and lies below the observations where the
fixed plumage lies above them. Real budgerigars are between the two. The optimiser follows the fixed-plumage
line, because it has the same plumage.

**In the thermoneutral zone the controllers agree**, since there the minimum metabolic rate is the answer
whatever is done to reach it.

**In the heat they differ in how, not in how much.** The rules, in this mode, raise core temperature, panting
and skin wetness in step, after going to full vasodilation. The optimiser lets the core rise sooner, to 39.6 °C in
air at 30 °C where the rule-based bird is still near 38 °C, because the weight on core temperature is low. It
vasodilates only part of the way and pants somewhat more. Both reproduce the observed rise in body temperature
and in water loss above 30 °C, and the observed upturn of metabolic rate with the ``Q_{10}`` effect of a warmer
core. Below 25 °C all three predict several times the water loss that was measured. One candidate cause, not
tested here, is that the exhaled air of a small bird in the cold is cooler and drier than is assumed.

**The optimiser changes its mind at 40 °C.** Between 39 and 40 °C its solution jumps: the core drops by 3 °C,
vasodilation is abandoned and panting doubles. Nothing in the bird changes there. The problem has more than one
local optimum, and the solver has moved from one to another, see
[Thermoregulation by optimisation](../manual/optimisation.md#What-to-check). The rule-based solution has no such
jumps, because its order is fixed.

## How the responses are used

```@example budgerigar
fig = Figure(size = (760, 320)) # hide
ax = Axis(fig[1, 1]; xlabel = "Air temperature (°C)", ylabel = "Flesh conductivity (W/m/K)") # hide
lines!(ax, x, [ustrip(u"W/m/K", r.thermoregulation.flesh_conductivity) for r in rules_fixed]; color = blue, linewidth = 2, label = "rules, plumage fixed") # hide
lines!(ax, x, [ustrip(u"W/m/K", r.parts.body.flesh_conductivity) for r in optimised]; color = orange, linewidth = 2, linestyle = :dash, label = "optimisation") # hide
axislegend(ax; position = :lt, labelsize = 10) # hide
ax = Axis(fig[1, 2]; xlabel = "Air temperature (°C)", ylabel = "Skin wetness (%)") # hide
lines!(ax, x, [100 * r.thermoregulation.skin_wetness for r in rules_fixed]; color = blue, linewidth = 2) # hide
lines!(ax, x, [100 * r.parts.body.skin_wetness for r in optimised]; color = orange, linewidth = 2, linestyle = :dash) # hide
ax = Axis(fig[1, 3]; xlabel = "Air temperature (°C)", ylabel = "Feather depth (mm)") # hide
lines!(ax, x, [ustrip(u"mm", r.thermoregulation.insulation_depth) for r in rules]; color = grey, linewidth = 2, label = "rules") # hide
lines!(ax, x, [ustrip(u"mm", r.thermoregulation.insulation_depth) for r in rules_fixed]; color = blue, linewidth = 2) # hide
axislegend(ax; position = :rt, labelsize = 10) # hide
fig # hide
```

The rule-based bird goes to full vasodilation before anything else, because vasodilation is early in the
order. The optimiser uses a little, because the weights give it a small cost and the alternatives are cheap.
Which is nearer the truth is a question about budgerigars. The value of having both is that the question can
be put.

## The range of the comparison

The sweep stops at 46 °C. In tests above that the optimiser did not return a usable solution for this bird: the
solutions it found had a metabolic rate well above the rule-based one, and by 50 °C the heat budget was not
satisfied at the point returned. The rule-based controller continues, to a panting multiplier above 9 at
50 °C. See the limits listed in [Thermoregulation by optimisation](../manual/optimisation.md#Present-limits).
