# Get started

BiophysicalBehaviour.jl computes what an organism does to regulate its body temperature, given the heat budget
of [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl).

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl")
```

```@setup get_started
using Main.FigureHelpers
using CairoMakie
```

## An endotherm

An organism here is a body and an [`OrganismTraits`](@ref): a thermal strategy, the heat-exchange traits of
HeatExchange.jl, and the behavioural traits that say what it can change. The example constructors give a 65 kg
ellipsoid with 2 mm of fur, defending a core temperature of 37 °C with a minimum metabolic rate of 77.6 W:

```@example get_started
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

traits = example_organism_traits()
physiology_traits = heat_exchange(traits)
fibres = insulation_pars(physiology_traits).dorsal
fur = FibrousLayer(fibres.depth, fibres.diameter, fibres.density)
fat = FatLayer(0.0, 901.0u"kg/m^3")
mammal = Organism(Body(shape_pars(physiology_traits), CompositeInsulation(fur, fat)), traits)
nothing # hide
```

The environment is that of HeatExchange.jl: the conditions of the hour and the properties of the site.

```@example get_started
environment(air_temperature) = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature)),
)
nothing # hide
```

[`thermoregulate`](@ref) takes the organism, the environment and a first guess of the temperatures of the skin
and the fur surface. In air at 5 °C:

```@example get_started
cold = environment(5.0u"°C")
init = BiophysicalBehaviour.initial_physiological_state(mammal, cold.environment_vars)
out = thermoregulate(mammal, cold, init)
out.energy_flows.metabolic_heat_flow
```

The heat budget balances at a metabolic rate above the minimum, so there is nothing to regulate: the animal
makes more heat, and the answer is that of `solve_metabolic_rate` in HeatExchange.jl.

At 35 °C the heat budget would balance only at a metabolic rate *below* the minimum. The animal makes more heat
than it can lose, and must act:

```@example get_started
hot = environment(35.0u"°C")
out = thermoregulate(mammal, hot, BiophysicalBehaviour.initial_physiological_state(mammal, hot.environment_vars))
state = out.thermoregulation
markdown_table(["Response", "At rest", "At 35 °C"], [ # hide
    ("axis ratio of the body", 1.1, state.axis_ratio_b), # hide
    ("flesh conductivity", 0.9u"W/m/K", state.flesh_conductivity), # hide
    ("core temperature", 37.0u"°C", celsius(state.core_temperature)), # hide
    ("panting multiplier", 1.0, state.pant), # hide
    ("skin wetness", 0.005, state.skin_wetness), # hide
    ("metabolic rate", 77.6u"W", out.energy_flows.metabolic_heat_flow), # hide
    ("evaporative water loss", "", out.mass_flows.m_evap), # hide
]) # hide
```

It has stretched out (the axis ratio is the length of the body relative to its width), sent blood to its skin,
let its core temperature rise by 2 °C and is breathing at several times the resting rate. Each was tried in
that order, one step at a time, until the heat budget balanced. See [Endotherm thermoregulation by rules](manual/endotherm_rules.md), and
[A mammal across air temperatures](tutorials/mammal.md) for the whole curve.

## An ectotherm

An ectotherm regulates mostly by where it is, so it needs somewhere to go.
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) computes the conditions above and
below the ground through a day. Two runs of it, one in the open and one in deep shade, bound the environments
within reach:

```@example get_started
using Microclimate

model = MicroModel(;
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
)
function microclimate(shade)
    inputs = MicroInputs(;
        site = example_site(),
        soil_profile = example_soil_profile(),
        environment_minmax = example_monthly_weather(),
        environment_daily = example_daily_environment(; shade = fill(shade, 12)),
        environment_hourly = example_hourly_environment(),
        initial_soil_moisture = fill(0.1, length(model.depths)),
    )
    solve(MicroProblem(model, inputs))
end
environments = AvailableEnvironments(microclimate(0.0), microclimate(0.9), 0.0, 0.9, model.depths, model.heights)
nothing # hide
```

These are the example inputs of Microclimate.jl: the middle day of each month at Madison, Wisconsin. The
lizard is the 40 g desert iguana of the HeatExchange.jl documentation, with the default behavioural traits: it
is diurnal, forages between 24 and 34 °C, prefers 30 °C, and can use shade and a burrow.

```@example get_started
lizard_traits = example_ectotherm_organism_traits()
lizard = Organism(Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()), lizard_traits)
limits = thermoregulation(lizard)
site = example_environment_pars(; elevation = 226.0u"m")
nothing # hide
```

For an ectotherm [`thermoregulate`](@ref) is called once for each hour. Early afternoon in July is hour 14 of
day 7:

```@example get_started
afternoon = (7 - 1) * 24 + 14
hour = thermoregulate(lizard, environments, limits, site, afternoon, 1)
(celsius(hour.core_temperature), hour.shade, hour.state)
```

Through the whole day, passing the depth chosen in one hour on to the next:

```@example get_started
function simulate(organism, environments, site, steps)
    limits = thermoregulation(organism)
    depth = limits.depth.reference
    commenced = false
    map(steps) do step
        (step - 1) % 24 == 0 && (commenced = false)
        out = thermoregulate(organism, environments, limits, site, step, depth; activity_commenced = commenced)
        depth = out.depth_node
        commenced = commenced || !(out.state isa Resting)
        out
    end
end

july = (7 - 1) * 24 .+ (1:24)
day = simulate(lizard, environments, site, july)
fig = Figure(size = (720, 560)) # hide
hours = 0:23 # hide
ax = Axis(fig[1, 1]; ylabel = "Temperature (°C)") # hide
open_air = environments.min_shade_result # hide
state_bands!(ax, hours, [s.state isa Active ? 2 : s.state isa Basking ? 1 : 0 for s in day]; y = 0.0, height = 2.0) # hide
hspan!(ax, 24, 34; color = (:grey, 0.15)) # hide
lines!(ax, hours, ustrip.(u"°C", open_air.profile.air_temperature[july, 1]); color = :steelblue, linestyle = :dash, label = "air at 1 cm, open") # hide
lines!(ax, hours, ustrip.(u"°C", open_air.soil_temperature[july, 1]); color = :sienna, linestyle = :dash, label = "ground surface, open") # hide
lines!(ax, hours, [ustrip(u"°C", s.core_temperature) for s in day]; color = :black, linewidth = 3, label = "body") # hide
axislegend(ax; position = :lt, labelsize = 11) # hide
ax2 = Axis(fig[2, 1]; ylabel = "Shade (%)") # hide
barplot!(ax2, hours, [100 * s.shade * (s.depth_node == 1) for s in day]; color = :seagreen) # hide
ax3 = Axis(fig[3, 1]; xlabel = "Hour", ylabel = "Depth (cm)", yreversed = true) # hide
barplot!(ax3, hours, [ustrip(u"cm", environments.depths[s.depth_node]) for s in day]; color = :sienna) # hide
hidexdecorations!(ax; grid = false) # hide
hidexdecorations!(ax2; grid = false) # hide
linkxaxes!(ax, ax2, ax3) # hide
rowsize!(fig.layout, 1, Relative(0.6)) # hide
fig # hide
```

The grey band is the range of body temperatures for activity. The coloured strip along the bottom of the top
panel is the state of the animal: blue at rest, orange basking, red active. See
[Ectotherm thermoregulation](manual/ectotherm.md) for the sequence of decisions, and
[A lizard's day](tutorials/lizard.md) for a year in the desert.

## Where next

- [Behaviour as control](manual/control.md) for the idea behind the package, and
  [Gradients](manual/gradients.md) for what a response changes.
- [Get started](https://biophysicalecology.github.io/HeatExchange.jl/dev/get_started) in the documentation of
  HeatExchange.jl for the heat budget underneath.
- [Thermoregulation by optimisation](manual/optimisation.md) for the alternative to rules.
- [Bodies of many parts](manual/multipart.md) for animals with a trunk, a head and limbs.
