# Activity and available environments

An organism that regulates its body temperature by moving needs a description of where it can go. This page
describes that description, [`AvailableEnvironments`](@ref), and the rules that say when the organism is out in
it at all.

```@setup environments
using Main.FigureHelpers
using Main.ExampleEnvironments
using CairoMakie
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
```

## Two microclimates

At any hour the places within reach of a small animal differ in three ways: how much of the sky is hidden by
vegetation, how far they are above the ground, and how far below it. A run of
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) gives the conditions at one level of
shade, at a set of heights in the air and a set of depths in the soil, for every hour. Two runs, at the least
and the most shade available, bound the rest:

```julia
environments = AvailableEnvironments(open_result, shaded_result, 0.0, 0.9, depths, heights)
```

```@example environments
environments = madison_environments() # the microclimates of Get started
nothing # hide
```

This is the arrangement of NicheMapR, in which the microclimate model is run for `minshade` and `maxshade` and
the ectotherm model moves between its `metout`/`soil` and `shadmet`/`shadsoil` tables (Kearney and Porter
2020).

| Field | Content |
|:--|:--|
| `min_shade_result`, `max_shade_result` | the two results of `Microclimate.solve` |
| `min_shade_fraction`, `max_shade_fraction` | the shade of each, from 0 to 1 |
| `depths` | the depths of the soil nodes. Node 1 is the surface |
| `heights` | the heights of the air nodes. Node 1 is the height of the animal on the ground |

Position is three numbers: a shade fraction, a height node and a depth node. They are held as the `current`
values of three [`SteppedParameter`](@ref)s in [`EctothermBehavioralLimits`](@ref).

Any object with the fields that are read will do in place of a `MicroResult`. The tests of the package build
one from the output tables of NicheMapR.

## From a position to an environment

[`interpolate_environment`](@ref) turns a position and an hour into the `EnvironmentalVars` that the heat
budget needs. Here is July on the surface in the open, on the surface in half shade, and 20 cm down:

```@example environments
limits = example_ectotherm_behavioral_limits()
site = example_environment_pars(; elevation = 226.0u"m")
step = (7 - 1) * 24 + 14

using Setfield: @set
half_shade = @set limits.shade.current = 0.45
underground = @set limits.depth.current = 11

rows = map((; open = limits, half_shade, underground)) do position
    interpolate_environment(environments, step, position, site)
end
fields = (:air_temperature, :ground_temperature, :sky_temperature, :relative_humidity, :wind_speed, :global_radiation, :shade) # hide
show_value(x) = x isa Unitful.Temperature ? celsius(x) : x # hide
markdown_table(["Variable", "Open", "Half shade", "20 cm down"], # hide
    [("`$f`", (show_value(getfield(r, f)) for r in rows)...) for f in fields]) # hide
```

Above the ground the air temperature, wind speed and the temperatures of the sky and the ground are
interpolated linearly between the two runs, according to where the current shade lies between the two shade
fractions. Solar radiation is passed on as it is in the open, with the shade fraction beside it, and
HeatExchange.jl applies the shade. Relative humidity is recomputed so that the vapour pressure of the air is
the same in the shade as in the open.

Below the ground there is no sun and almost no wind, and the air, sky and ground are all at the temperature of
the soil at that depth. Shade underground is not interpolated: a [`BurrowShadeMode`](@ref) says whether the
burrow is in the open ([`MinShadeOnly`](@ref)), in the shade ([`MaxShadeOnly`](@ref)) or wherever is tolerable
([`AdaptiveBurrowShade`](@ref)).

## From gridded climate data

[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) runs Microclimate.jl from
gridded climate and terrain datasets, for a list of points or a whole raster. Its output is a stack of rasters
with dimensions for point, depth and height. A small adapter turns the output for one point into the object
that [`AvailableEnvironments`](@ref) reads:

```julia
using MicroclimateMapper, Microclimate, RasterDataSources, Dates

depths = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 0.5, 1.0]u"m"
output_layers = (
    LayerSpec(:soil_temperature, :soil), LayerSpec(:soil_humidity, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:air_temperature, :profile), LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:global_radiation, :scalar), LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:diffuse_fraction, :scalar), LayerSpec(:reference_temperature, :scalar),
    LayerSpec(:pressure, :scalar), LayerSpec(:zenith_angle, :solar),
)
model = MicroMapModel(;
    micro_model = MicroModel(; depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model()),
    dem_source = CRUCL2, weather_source = CRUCL2,   # or TerraClimate, SILO, BARRA, ...
    output_layers,
)
problem = (; model, points = [geocode("Palm Springs, CA")], dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
           soil_profile = example_soil_profile(depths), init = (; soil_moisture = fill(0.2, length(depths))))
open_output = solve(MicroVectorProblem(; problem...))
shaded_output = solve(MicroVectorProblem(; problem..., data = (; shade = 0.9)))

point_environment(output) = (;
    pressure = collect(output.pressure[point = 1]),
    reference_temperature = collect(output.reference_temperature[point = 1]),
    global_radiation = collect(output.global_radiation[point = 1]),
    diffuse_fraction = collect(output.diffuse_fraction[point = 1]),
    sky_temperature = collect(output.sky_temperature[point = 1]),
    soil_temperature = collect(output.soil_temperature[point = 1]),
    soil_humidity = collect(output.soil_humidity[point = 1]),
    soil_thermal_conductivity = collect(output.soil_thermal_conductivity[point = 1]),
    profile = (;
        air_temperature = collect(output.air_temperature[point = 1]),
        relative_humidity = collect(output.relative_humidity[point = 1]),
        wind_speed = collect(output.wind_speed[point = 1]),
    ),
    solar_radiation = (; zenith_angle = collect(output.zenith_angle[point = 1])),
)

environments = AvailableEnvironments(point_environment(open_output), point_environment(shaded_output),
                                     0.0, 0.9, depths, heights)
```

This block is not run here, since it downloads climate data. Changing `weather_source` runs the same animal at
the same place against another dataset, and a raster in place of the points gives an environment, and so an
animal, for every cell.

## The thermal landscape

Solving the heat budget of the same lizard at every position gives the body temperatures available to it. This
is the map that the controller searches:

```@example environments
lizard = Organism(Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()), example_ectotherm_organism_traits())
july = (7 - 1) * 24 .+ (1:24)
shades = 0.0:0.1:0.9
body_temperature(position, step) =
    solve_body_temperature(lizard, interpolate_environment(environments, step, position, site), site)
above = [ustrip(u"°C", body_temperature((@set limits.shade.current = s), t)) for t in july, s in shades]
nodes = 2:14
below = [ustrip(u"°C", environments.min_shade_result.soil_temperature[t, n]) for t in july, n in nodes]
fig = Figure(size = (720, 520)) # hide
ax = Axis(fig[1, 1]; ylabel = "Shade (%)", title = "Body temperature on the surface (°C)") # hide
hm = heatmap!(ax, 0:23, 100 .* shades, above; colormap = :thermal, colorrange = (5, 60)) # hide
contour!(ax, 0:23, 100 .* shades, above; levels = [24.0, 34.0], color = :white, linewidth = 2) # hide
ax2 = Axis(fig[2, 1]; xlabel = "Hour", ylabel = "Depth (cm)", yreversed = true, title = "Soil temperature in the open (°C)") # hide
depths_cm = ustrip.(u"cm", environments.depths[nodes]) # hide
heatmap!(ax2, 0:23, depths_cm, below; colormap = :thermal, colorrange = (5, 60)) # hide
contour!(ax2, 0:23, depths_cm, below; levels = [24.0, 34.0], color = :white, linewidth = 2) # hide
Colorbar(fig[1:2, 2], hm) # hide
hidexdecorations!(ax; grid = false) # hide
fig # hide
```

The white contours are the limits of the range for activity, 24 and 34 °C. Between them, on the surface, the
animal can forage. Outside them it must be somewhere else, and the lower panel shows what the soil offers.

## When to be active

Whether an organism is abroad at a given hour is decided first by its [`ActivityPeriod`](@ref), from the
zenith angle of the sun and the solar radiation, and not by its state. This is open-loop control, see
[Behaviour as control](control.md#Open-and-closed-loops).

| Type | Active when |
|:--|:--|
| [`Diurnal`](@ref) | the sun is above the horizon and there is sunlight |
| [`Nocturnal`](@ref) | the sun is below the horizon, or there is no sunlight |
| [`Crepuscular`](@ref) | the sun is within 5° of the horizon |
| [`CombinedActivity`](@ref) | any of the periods it holds, for example `CombinedActivity(Diurnal(), Crepuscular())` |
| [`ResponsiveActivity`](@ref) | a function of the zenith angle and solar radiation returns `true` |

```@example environments
open_air = environments.min_shade_result
zenith = open_air.solar_radiation.zenith_angle[july]
sunlight = open_air.global_radiation[july]
periods = (Diurnal(), Nocturnal(), Crepuscular())
fig = Figure(size = (720, 200)) # hide
ax = Axis(fig[1, 1]; xlabel = "Hour", yticks = (1:3, ["Diurnal", "Nocturnal", "Crepuscular"])) # hide
for (i, period) in enumerate(periods) # hide
    active = [is_active(period, z, s) for (z, s) in zip(zenith, sunlight)] # hide
    for (h, a) in zip(0:23, active) # hide
        a && poly!(ax, Rect2f(h - 0.5, i - 0.35, 1.0, 0.7); color = :grey35) # hide
    end # hide
end # hide
xlims!(ax, -0.5, 23.5) # hide
ylims!(ax, 0.4, 3.6) # hide
fig # hide
```

Outside its activity period the animal is at rest in its retreat. Within it, the animal may still be at rest,
if no position within reach brings its body temperature into range.

## The state of the organism

The outcome for each hour is an [`OrganismState`](@ref), decided by where the body temperature lies between the
[threshold traits](states_traits.md):

| State | Condition | NicheMapR `ACT` |
|:--|:--|:--|
| [`Resting`](@ref) | outside the activity period, underground, or above ground with a body temperature outside the basking and activity ranges | 0 |
| [`Basking`](@ref) | body temperature from `basking_temperature_min` up to `active_temperature_min` | 1 |
| [`Active`](@ref) | body temperature from `active_temperature_min` to `active_temperature_max` | 2 |

Hours in the `Active` state are the time available for feeding and for everything else that needs the animal
to be out and moving. They are the link from this package to models of energy and water budgets and of
population growth.
