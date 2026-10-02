# A lizard's day

The desert iguana, *Dipsosaurus dorsalis*, is active at body temperatures that would kill most lizards, in a
place where the ground reaches 60 °C. Porter et al. (1973) used it to show that the times and places of an
animal's activity could be computed from its heat budget and the microclimates open to it. This tutorial repeats
that calculation: a 40 g lizard at Palm Springs, California, on the middle day of each month.

It needs two more packages than the rest of this documentation,
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) for the environments and the Statistics
standard library.

```@setup lizard
using Main.FigureHelpers
using CairoMakie
```

## The microclimates

The site is at 33.8° N and 130 m, on pale sand. The weather is monthly climate: mean daily minima and maxima of
air temperature, wind speed and humidity, under clear skies.

```@example lizard
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Microclimate, Unitful
using Statistics: mean

elevation = 130.0u"m"
site = example_site(; latitude = 33.8u"°", longitude = -116.5u"°", elevation,
    albedo = 0.33, roughness_height = 0.0005u"m")

temperature_min = [3.7, 6.0, 8.2, 11.2, 15.0, 19.1, 22.6, 22.2, 18.8, 13.5, 7.5, 3.5]u"°C"
temperature_max = [19.3, 21.8, 23.6, 27.0, 31.0, 35.9, 39.2, 38.4, 35.4, 30.3, 23.9, 19.7]u"°C"
weather = example_monthly_weather(;
    reference_temperature_min = temperature_min,
    reference_temperature_max = temperature_max,
    reference_wind_speed_min = [0.24, 0.26, 0.29, 0.31, 0.31, 0.31, 0.31, 0.29, 0.28, 0.25, 0.24, 0.23]u"m/s",
    reference_wind_speed_max = [2.40, 2.62, 2.91, 3.06, 3.06, 3.06, 3.06, 2.91, 2.76, 2.47, 2.40, 2.33]u"m/s",
    reference_humidity_min = [34.0, 33.7, 33.4, 31.1, 31.4, 31.1, 33.4, 35.2, 34.2, 33.0, 33.1, 34.0] ./ 100,
    reference_humidity_max = [95.6, 94.3, 89.4, 83.3, 82.7, 83.2, 86.2, 89.1, 90.8, 92.2, 94.8, 99.3] ./ 100,
    cloud_cover_min = fill(0.0, 12),
    cloud_cover_max = fill(0.0, 12),
)
nothing # hide
```

The lizard can climb into creosote bushes, so the air is described at five heights up to 2 m, the last being
the height of the weather data. The soil is described at the default depths of Microclimate.jl, down to 2 m.

```@example lizard
heights = [1.0, 50.0, 100.0, 150.0, 200.0]u"cm"
model = MicroModel(; heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
)
depths = model.depths
deep_soil_temperature = mean([mean(u"K".(temperature_min)), mean(u"K".(temperature_max))])

function microclimate(shade)
    inputs = MicroInputs(;
        site,
        soil_profile = example_soil_profile(depths; bulk_density = 1.5u"Mg/m^3"),
        environment_minmax = weather,
        environment_daily = example_daily_environment(;
            shade = fill(shade, 12),
            rainfall = fill(0.0u"kg/m^2", 12),
            deep_soil_temperature = fill(deep_soil_temperature, 12),
        ),
        environment_hourly = example_hourly_environment(; elevation),
        initial_soil_moisture = fill(0.2, length(depths)),
    )
    solve(MicroProblem(model, inputs))
end

open_ground = microclimate(0.0)
shaded_ground = microclimate(0.9)
environments = AvailableEnvironments(open_ground, shaded_ground, 0.0, 0.9, depths, heights)
nothing # hide
```

What the lizard has to work with, in January and July:

```@example lizard
month_hours(month) = (month - 1) * 24 .+ (1:24)
fig = Figure(size = (760, 340)) # hide
for (i, (month, name)) in enumerate(((1, "January"), (7, "July"))) # hide
    ax = Axis(fig[1, i]; xlabel = "Hour", ylabel = "Temperature (°C)", title = name, limits = (nothing, (-5, 70))) # hide
    t = month_hours(month) # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.soil_temperature[t, 1]); color = :sienna, linewidth = 2, label = "ground surface") # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.profile.air_temperature[t, 1]); color = :steelblue, linewidth = 2, label = "air, 1 cm") # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.profile.air_temperature[t, 5]); color = :steelblue, linestyle = :dash, label = "air, 2 m") # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.soil_temperature[t, 7]); color = :black, linestyle = :dot, label = "soil, 10 cm") # hide
    hspan!(ax, 38, 43; color = (:red, 0.12)) # hide
    i == 1 && axislegend(ax; position = :lt, labelsize = 10) # hide
end # hide
fig # hide
```

The red band is the range of body temperatures at which the lizard is active.

## The lizard

The traits are those of Porter et al. (1973). The lizard is active between 38 and 43 °C and prefers 38.5 °C.
It changes colour, from a solar absorptivity of 0.8 when cool to 0.6 when hot. It turns broadside to the sun to
bask. It climbs. It does not seek shade on the ground, but it has a burrow, no shallower than 2.5 cm, under open
ground.

```@example lizard
desert_iguana(; changes...) = Organism(Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()), example_ectotherm_organism_traits(;
    activity_period = CombinedActivity(Diurnal(), Crepuscular()),
    target_temperature = u"K"(38.5u"°C"),
    active_temperature_min = u"K"(38.0u"°C"),
    active_temperature_max = u"K"(43.0u"°C"),
    basking_temperature_min = u"K"(34.0u"°C"),
    emerge_temperature_min = u"K"(15.0u"°C"),
    critical_temperature_min = u"K"(3.0u"°C"),
    critical_temperature_max = u"K"(44.0u"°C"),
    can_climb = true,
    can_seek_shade = false,
    can_retreat_underground = true,
    depth_min_underground = 3,
    burrow_shade_mode = MinShadeOnly(),
    can_solar_orient = true,
    can_press_to_ground = false,
    can_change_absorptivity = true,
    absorptivity_min = 0.6,
    absorptivity_max = 0.8,
    absorptivity_step = 0.003,
    heat_exchange = example_ectotherm_heat_exchange_traits(;
        conduction_pars_external = example_ectotherm_conduction_pars_external(; conduction_fraction = 0.0),
        evaporation_pars = example_ectotherm_evaporation_pars(; eye_fraction = 0.0003, skin_wetness = 0.001),
        radiation_pars = example_ectotherm_radiation_pars(;
            body_absorptivity_dorsal = 0.8, body_absorptivity_ventral = 0.8,
            body_emissivity_dorsal = 0.95, body_emissivity_ventral = 0.95,
            solar_orientation = Intermediate()),
        respiration_pars = example_ectotherm_respiration_pars(; mouth_fraction = 0.0),
    ),
    changes...,
))
iguana = desert_iguana()
ground = example_environment_pars(; elevation, ground_albedo = 0.33)
nothing # hide
```

## The year

[`thermoregulate`](@ref) is called for each of the 288 hours. The depth of each hour is passed to the next,
and a flag records whether the lizard has yet been out on the day:

```@example lizard
function simulate(organism, environments, ground, steps)
    limits = thermoregulation(organism)
    depth = limits.depth.reference
    commenced = false
    map(steps) do step
        (step - 1) % 24 == 0 && (commenced = false)
        out = thermoregulate(organism, environments, limits, ground, step, depth; activity_commenced = commenced)
        depth = out.depth_node
        commenced = commenced || !(out.state isa Resting)
        out
    end
end

year = simulate(iguana, environments, ground, 1:288)
nothing # hide
```

### A day in spring and a day in summer

```@example lizard
state_code(s) = s.state isa Active ? 2 : s.state isa Basking ? 1 : 0
position_cm(s) = s.depth_node > 1 ? -ustrip(u"cm", depths[s.depth_node]) :
                 s.height > heights[1] ? ustrip(u"cm", s.height) : 0.0
fig = Figure(size = (760, 560)) # hide
for (i, (month, name)) in enumerate(((4, "April"), (7, "July"))) # hide
    t = month_hours(month) # hide
    day = year[t] # hide
    ax = Axis(fig[1, i]; ylabel = "Temperature (°C)", title = name, limits = (nothing, (0, 62))) # hide
    state_bands!(ax, 0:23, state_code.(day); y = 0.0, height = 2.5) # hide
    hspan!(ax, 38, 43; color = (:red, 0.12)) # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.soil_temperature[t, 1]); color = :sienna, linestyle = :dash, label = "ground surface") # hide
    lines!(ax, 0:23, ustrip.(u"°C", open_ground.profile.air_temperature[t, 1]); color = :steelblue, linestyle = :dash, label = "air, 1 cm") # hide
    lines!(ax, 0:23, [ustrip(u"°C", s.core_temperature) for s in day]; color = :black, linewidth = 3, label = "body") # hide
    i == 1 && axislegend(ax; position = :lt, labelsize = 10) # hide
    hidexdecorations!(ax; grid = false) # hide
    ax2 = Axis(fig[2, i]; ylabel = "Height or depth (cm)", limits = (nothing, (-12, 110))) # hide
    barplot!(ax2, 0:23, position_cm.(day); color = [p < 0 ? :sienna : :seagreen for p in position_cm.(day)]) # hide
    hlines!(ax2, [0.0]; color = :black) # hide
    hidexdecorations!(ax2; grid = false) # hide
    ax3 = Axis(fig[3, i]; xlabel = "Hour", ylabel = "Absorptivity", limits = (nothing, (0.55, 0.85))) # hide
    scatterlines!(ax3, 0:23, [s.absorptivity for s in day]; color = :grey30) # hide
end # hide
rowsize!(fig.layout, 1, Relative(0.5)) # hide
fig # hide
```

The strip at the foot of the top panels is the state of the lizard: blue at rest, orange basking, red active.

In April the lizard comes out in mid-morning, basks, and is active through the middle of the day, paling as it
warms. It never needs to leave the ground. In July it is out soon after sunrise and too hot on the ground by
mid-morning. It climbs, where the air is cooler and moves faster, and when that fails it goes underground
through the middle of the day. It comes out again in the late afternoon, up in the bushes first and then on the
ground. This is the pattern that Porter et al. (1973) computed and observed: one period of activity in spring,
two in summer.

### The year at a glance

```@example lizard
body = [ustrip(u"°C", year[(m - 1) * 24 + h].core_temperature) for h in 1:24, m in 1:12]
state = [state_code(year[(m - 1) * 24 + h]) for h in 1:24, m in 1:12]
position = [position_cm(year[(m - 1) * 24 + h]) for h in 1:24, m in 1:12]
months = (1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]) # hide
fig = Figure(size = (760, 330)) # hide
ax = Axis(fig[1, 1]; xlabel = "Month", ylabel = "Hour", title = "Body temperature (°C)", xticks = months) # hide
hm = heatmap!(ax, 1:12, 0:23, permutedims(body); colormap = :thermal) # hide
Colorbar(fig[2, 1], hm; vertical = false) # hide
ax = Axis(fig[1, 2]; xlabel = "Month", title = "State", xticks = months) # hide
heatmap!(ax, 1:12, 0:23, permutedims(state); colormap = [RGBf(0.45, 0.55, 0.75), RGBf(0.98, 0.75, 0.35), RGBf(0.80, 0.30, 0.25)], colorrange = (0, 2)) # hide
Label(fig[2, 2], "blue: resting, orange: basking, red: active"; fontsize = 11, tellwidth = false) # hide
ax = Axis(fig[1, 3]; xlabel = "Month", title = "Height (+) or depth (−), cm", xticks = months) # hide
hm = heatmap!(ax, 1:12, 0:23, permutedims(position); colormap = :broc, colorrange = (-100, 100)) # hide
Colorbar(fig[2, 3], hm; vertical = false) # hide
fig # hide
```

```@example lizard
hours_in(code) = [count(==(code), state[:, m]) for m in 1:12]
fig, ax = figure_axis("Month", "Hours on the middle day"; size = (700, 320), xticks = months) # hide
barplot!(ax, repeat(1:12, 2), vcat(hours_in(2), hours_in(1)); stack = repeat(1:2; inner = 12), # hide
    color = repeat([RGBf(0.80, 0.30, 0.25), RGBf(0.98, 0.75, 0.35)]; inner = 12)) # hide
fig # hide
```

Red is active and orange basking. The lizard is not active at all from November to February: nowhere within
its reach is warm enough. That seasonal limit, and the two-peaked days of summer, both follow from two
threshold traits, 38 and 43 °C, and the physics.

## What the behaviours are for

Each behaviour can be taken away, see [Ectotherm thermoregulation](../manual/ectotherm.md):

```@example lizard
variants = (
    "all behaviours" => iguana,
    "no climbing" => desert_iguana(; can_climb = false),
    "no colour change" => desert_iguana(; can_change_absorptivity = false),
    "no burrow" => desert_iguana(; can_retreat_underground = false),
)
rows = map(variants) do (label, animal)
    result = simulate(animal, environments, ground, 1:288)
    (label, count(s -> s.state isa Active, result), celsius(maximum(s.core_temperature for s in result)))
end
markdown_table(["Behaviours", "Hours active, of 288", "Highest body temperature"], rows) # hide
```

Climbing buys hours of activity in summer. The burrow buys survival: without it the body temperature on the
hottest afternoon is above the critical maximum of 44 °C.

## Where next

The body temperatures and activity states computed here are the inputs to models of what the animal can do
with its time: [ThermalPhysiology.jl](https://github.com/BiophysicalEcology/ThermalPhysiology.jl) for
performance and survival as functions of body temperature, and energy and water budgets driven by hours of
activity. For a site described by real weather in place of monthly means, see the documentation of
Microclimate.jl.
