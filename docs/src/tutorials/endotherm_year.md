# A desert mammal through the year

[A mammal across air temperatures](mammal.md) puts an endotherm in a metabolic chamber. This tutorial puts one
outdoors: a 1 kg mammal at Palm Springs, California, in the microclimates of [A lizard's day](lizard.md), hour
by hour through the middle day of each month. What does the desert cost it in energy and water, first sitting
in the open and then choosing where to be?

```@setup endotherm_year
using Main.FigureHelpers
using Main.ExampleEnvironments
using CairoMakie
```

## The animal

A 1 kg ellipsoid three times as long as it is wide, with 2 mm of fur, defending 38 °C. Its core can rise by
5 °C, it can pant to 15 times its resting ventilation, and all its skin can be wet. It pants and sweats as its
core warms, the mode `CorePantingSweatingFirst`.

```@example endotherm_year
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

shape = example_shape_pars(; mass = 1.0u"kg", axis_ratio_b = 3.0, axis_ratio_c = 3.0)
basal = metabolic_rate(McKechnieWolf(), shape.mass)
core = u"K"(38.0u"°C")
physiology_traits = example_heat_exchange_traits(;
    shape_pars = shape,
    metabolism_pars = example_metabolism_pars(; core_temperature = core, q10 = 2.0, metabolic_heat_flow = basal),
)
limits = example_thermoregulation_limits(;
    thermoregulation_mode = CorePantingSweatingFirst(),
    max_iterations = 200,
    minimum_heat_flow = basal,
    axis_ratio_factor = 3.0,
    core_temperature = core, core_temperature_ref = core, core_temperature_max = core + 5.0u"K",
    pant_max = 15.0, pant_step = 0.5, pant_multiplier = 1.0,
    skin_wetness_max = 1.0, skin_wetness_step = 0.05,
)
fibres = insulation_pars(physiology_traits).dorsal
coat = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3"))
animal = Organism(Body(shape, coat),
    OrganismTraits(Endotherm(), physiology_traits, BehavioralTraits(; thermoregulation = limits)))
basal
```

## In the open

The microclimates are those built in [A lizard's day](lizard.md). [`interpolate_environment`](@ref) gives the
conditions at a position, here the default position of a set of positional limits: on the surface, at the
lowest height, in no shade.

```@example endotherm_year
environments = palm_springs_environments() # as built in A lizard's day
ground = example_environment_pars(; elevation = 130.0u"m", ground_albedo = 0.33)
in_the_open = example_ectotherm_behavioral_limits()

exposed = map(1:288) do step
    environment_vars = interpolate_environment(environments, step, in_the_open, ground)
    init = BiophysicalBehaviour.initial_physiological_state(animal, environment_vars)
    thermoregulate(animal, (; environment_pars = ground, environment_vars), init)
end
nothing # hide
```

```@example endotherm_year
by_hour_and_month(f, results) = [f(results[(m - 1) * 24 + h]) for h in 1:24, m in 1:12]
months = (1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]) # hide
function year_maps(panels; size = (760, 560)) # hide
    fig = Figure(; size) # hide
    for (i, (title, values, colormap)) in enumerate(panels) # hide
        row, col = fldmod1(i, 2) # hide
        ax = Axis(fig[row, 2col - 1]; xlabel = "Month", ylabel = "Hour", title, xticks = months, titlesize = 12) # hide
        hm = heatmap!(ax, 1:12, 0:23, permutedims(values); colormap) # hide
        Colorbar(fig[row, 2col], hm) # hide
    end # hide
    fig # hide
end # hide
year_maps(( # hide
    ("Metabolic rate (multiple of basal)", by_hour_and_month(r -> r.energy_flows.metabolic_heat_flow / basal, exposed), :viridis), # hide
    ("Evaporative water loss (g/h)", by_hour_and_month(r -> ustrip(u"g/hr", r.mass_flows.m_evap), exposed), :dense), # hide
    ("Core temperature (°C)", by_hour_and_month(r -> ustrip(u"°C", r.thermoregulation.core_temperature), exposed), :thermal), # hide
    ("Panting multiplier", by_hour_and_month(r -> r.thermoregulation.pant, exposed), :amp), # hide
)) # hide
```

On winter nights the animal needs several times its basal rate. On summer days it is at the minimum, with its
core above the setpoint, panting, and losing water at a rate it could not sustain: a 1 kg animal holds about
650 g of water.

## Choosing where to be

An endotherm has the positions an ectotherm has, and uses them for a different reason: each degree avoided by
moving is water not evaporated, and each degree gained by shelter on a cold night is food not burned. Given
available environments and a set of positional limits, [`thermoregulate`](@ref) for an endotherm first selects
a position and then thermoregulates physiologically there.

The selection is on *operative temperature*, the temperature of the animal as a passive object with no
metabolism, against two thresholds. Above the upper one it tries shade, then height, then the burrow. Below the
lower one it leaves shade, and retreats if colder than its critical minimum. Outside its activity period it is
in its burrow. See [An endotherm that chooses where to be](../manual/ectotherm.md#An-endotherm-that-chooses-where-to-be).

Here the animal is diurnal, comfortable at operative temperatures from 10 to 38 °C, with a burrow from 30 cm
down:

```@example endotherm_year
positions = example_ectotherm_behavioral_limits(;
    active_temperature_min = u"K"(10.0u"°C"),
    active_temperature_max = u"K"(38.0u"°C"),
    critical_temperature_min = u"K"(-10.0u"°C"),
    critical_temperature_max = u"K"(45.0u"°C"),
    can_solar_orient = false,
    can_press_to_ground = false,
    depth_min_underground = 13,
    burrow_shade_mode = MinShadeOnly(),
)
environments.depths[13]
```

```@example endotherm_year
function choose(animal, positions)
    depth = 1
    map(1:288) do step
        out = thermoregulate(Endotherm(), RuleBasedSequentialControl(), animal, environments, positions, ground, step, depth)
        depth = out.depth_node
        out
    end
end
chosen = choose(animal, positions)
keys(first(chosen))
```

Each hour's output has the position chosen, the operative temperature `Te` there, and under `endotherm_out`
the physiological result.

```@example endotherm_year
year_maps(( # hide
    ("Shade used on the surface (%)", by_hour_and_month(r -> r.depth_node == 1 ? 100 * r.shade : NaN, chosen), :Greens), # hide
    ("Depth (cm)", by_hour_and_month(r -> ustrip(u"cm", environments.depths[r.depth_node]), chosen), :turbid), # hide
    ("Metabolic rate (multiple of basal)", by_hour_and_month(r -> r.endotherm_out.energy_flows.metabolic_heat_flow / basal, chosen), :viridis), # hide
    ("Evaporative water loss (g/h)", by_hour_and_month(r -> ustrip(u"g/hr", r.endotherm_out.mass_flows.m_evap), chosen), :dense), # hide
)) # hide
```

Blank cells in the first panel are hours underground.

## What the choice is worth

Summed over the middle day of each month:

```@example endotherm_year
daily(f, results) = [sum(f(results[(m - 1) * 24 + h]) for h in 1:24) for m in 1:12]
energy(results, f) = ustrip.(u"kJ", daily(r -> f(r).energy_flows.metabolic_heat_flow * 1u"hr", results))
water(results, f) = ustrip.(u"g", daily(r -> f(r).mass_flows.m_evap * 1u"hr", results))
fig = Figure(size = (760, 340)) # hide
ax = Axis(fig[1, 1]; xlabel = "Month", ylabel = "Energy (kJ/day)", xticks = months) # hide
scatterlines!(ax, 1:12, energy(exposed, identity); linewidth = 2, color = :grey35, label = "in the open") # hide
scatterlines!(ax, 1:12, energy(chosen, r -> r.endotherm_out); linewidth = 2, color = :darkorange, label = "choosing") # hide
hlines!(ax, [ustrip(u"kJ", basal * 24u"hr")]; color = :grey, linestyle = :dash) # hide
axislegend(ax; position = :ct, labelsize = 10) # hide
ax = Axis(fig[1, 2]; xlabel = "Month", ylabel = "Evaporative water loss (g/day)", xticks = months) # hide
scatterlines!(ax, 1:12, water(exposed, identity); linewidth = 2, color = :grey35) # hide
scatterlines!(ax, 1:12, water(chosen, r -> r.endotherm_out); linewidth = 2, color = :darkorange) # hide
fig # hide
```

```@example endotherm_year
markdown_table(["", "In the open", "Choosing where to be"], [ # hide
    ("energy, mean of the twelve days (kJ/day)", sum(energy(exposed, identity)) / 12, sum(energy(chosen, r -> r.endotherm_out)) / 12), # hide
    ("as a multiple of basal", sum(energy(exposed, identity)) / 12 / ustrip(u"kJ", basal * 24u"hr"), sum(energy(chosen, r -> r.endotherm_out)) / 12 / ustrip(u"kJ", basal * 24u"hr")), # hide
    ("water in July (g/day)", water(exposed, identity)[7], water(chosen, r -> r.endotherm_out)[7]), # hide
    ("highest core temperature (°C)", maximum(ustrip(u"°C", r.thermoregulation.core_temperature) for r in exposed), maximum(ustrip(u"°C", r.endotherm_out.thermoregulation.core_temperature) for r in chosen)), # hide
]) # hide
```

The dashed line is the basal rate through a day. The burrow saves energy through the whole year, by keeping the
animal out of the cold night air and from under the cold night sky. Shade and the burrow together save water
through the summer.

The animal is simple, and the result is a demonstration of the calculation. Its burrow is ventilated soil air
with no nest, it does not huddle or raise its fur, and its activity costs nothing. Each is a change to the
organism or its limits.

## From here to a map

The same loop runs for any place with a microclimate.
[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) solves Microclimate.jl from
gridded climate and terrain data, for a point or every cell of a raster, see
[From gridded climate data](../manual/environments.md#From-gridded-climate-data). Running an animal on each cell
turns the figures above into maps of the energy and water a place costs it.

The model is at steady state in every hour. A 1 kg animal takes a good part of an hour to warm or cool, and a
larger one much longer, so the peaks above are those of an animal with no thermal inertia. Behaviour during
heat budgets through time, alternating bouts of activity and rest as body temperature rises and falls, is in
development.
