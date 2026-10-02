# A mammal across air temperatures

The tutorial
[An endotherm: metabolic rate](https://biophysicalecology.github.io/HeatExchange.jl/dev/tutorials/endotherm)
in the documentation of HeatExchange.jl follows a 65 kg animal from 0 °C upward, and stops where the metabolic
rate needed to hold its core temperature falls to the least it can produce. This tutorial carries on from
there, through the thermoneutral zone and out the other side, and compares the result with `endoR_devel` of
NicheMapR. It reproduces the calculation behind Fig. 2 of Kearney et al. (2021).

```@setup mammal
using Main.FigureHelpers
using CairoMakie
```

## The animal

The defaults of `endoR`: a 65 kg ellipsoid a little longer than it is wide, with 2 mm of fur, defending 37 °C
with a basal metabolic rate of 77.6 W. That is roughly a resting, lightly clothed human.

```@example mammal
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

physiology_traits = example_heat_exchange_traits()
limits = example_thermoregulation_limits()
fibres = insulation_pars(physiology_traits).dorsal
insulation = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density),
                                 FatLayer(0.0, 901.0u"kg/m^3"))
mammal = Organism(
    Body(shape_pars(physiology_traits), insulation),
    OrganismTraits(Endotherm(), physiology_traits, BehavioralTraits(; thermoregulation = limits)),
)
nothing # hide
```

What it can do about heat, from [`ThermoregulationLimits`](@ref):

```@example mammal
markdown_table(["Response", "From", "To", "Step"], [ # hide
    ("uncurl: ratio of long to short axis", limits.axis_ratio_factor.current, limits.axis_ratio_factor.max, limits.axis_ratio_factor.step), # hide
    ("vasodilate: flesh conductivity", limits.flesh_conductivity.current, limits.flesh_conductivity.max, limits.flesh_conductivity.step), # hide
    ("raise core temperature", celsius(limits.core_temperature.current), celsius(limits.core_temperature.max), limits.core_temperature.step), # hide
    ("pant: multiple of resting ventilation", limits.panting.pant.current, limits.panting.pant.max, limits.panting.pant.step), # hide
    ("sweat: wet fraction of skin", limits.skin_wetness.current, limits.skin_wetness.max, limits.skin_wetness.step), # hide
]) # hide
```

Its fur is not raised or flattened: the depth has no range. The controller is the default,
`RuleBasedSequentialControl` with `CoreFirst`, which is `TREGMODE = 1`.

## The sweep

The environment is a metabolic chamber: still air at 0.1 m/s and 5 % relative humidity, no sun, and walls at
air temperature.

```@example mammal
function respond(animal, air_temperature)
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature))
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    thermoregulate(animal, environment, BiophysicalBehaviour.initial_physiological_state(animal, environment_vars))
end

air = 0.0:1.0:50.0
sweep = [respond(mammal, T * u"°C") for T in air]
nothing # hide
```

```@example mammal
using CSV, DataFrames
reference = CSV.read(joinpath(@__DIR__, "..", "data", "endoR_thermoreg_reference.csv"), DataFrame)
fig = Figure(size = (760, 760)) # hide
state(f) = [f(out.thermoregulation) for out in sweep] # hide
panels = ( # hide
    ("Metabolic rate (W)", [watts(out.energy_flows.metabolic_heat_flow) for out in sweep], reference.metabolic_W), # hide
    ("Evaporative water loss (g/h)", [ustrip(u"g/hr", out.mass_flows.m_evap) for out in sweep], reference.respiratory_water_g_h .+ reference.cutaneous_water_g_h), # hide
    ("Axis ratio", state(s -> s.axis_ratio_b), reference.shape_b), # hide
    ("Flesh conductivity (W/m/K)", state(s -> ustrip(u"W/m/K", s.flesh_conductivity)), reference.flesh_conductivity), # hide
    ("Core temperature (°C)", state(s -> ustrip(u"°C", s.core_temperature)), reference.core_C), # hide
    ("Panting multiplier", state(s -> s.pant), reference.pant), # hide
    ("Skin temperature (°C)", state(s -> ustrip(u"°C", s.skin_temperature)), reference.skin_C), # hide
    ("Skin wetness (%)", state(s -> 100 * s.skin_wetness), reference.skin_wetness_pct), # hide
) # hide
for (i, (label, here, there)) in enumerate(panels) # hide
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel = i > 6 ? "Air temperature (°C)" : "", ylabel = label) # hide
    lines!(ax, air, here; linewidth = 2, color = :black, label = "BiophysicalBehaviour.jl") # hide
    scatter!(ax, reference.air_temperature_C, there; color = :firebrick, markersize = 7, label = "NicheMapR endoR_devel") # hide
    i == 1 && axislegend(ax; position = :rt, labelsize = 10) # hide
end # hide
fig # hide
```

The curve has the form that every textbook of thermal physiology draws, and here nothing about it was drawn.
From the left:

**Below the lower critical temperature** the animal is in its most heat-conserving state and the heat budget
gives the metabolic rate directly. It falls in a straight line as the air warms, with a slope set by the
insulation and the surface area.

**The thermoneutral zone** begins where that line meets the basal rate. From there the controller holds the
metabolic rate at the minimum by spending one response at a time. First posture: the animal uncurls, increasing
its surface area. When it is fully stretched, vasodilation: flesh conductivity rises and the skin warms towards
the core.

**Above the upper critical temperature** those are exhausted. The core temperature is let rise by 2 °C, and
with it the metabolic rate, by the ``Q_{10}`` effect. Then panting begins, and evaporative water loss climbs
steeply. Sweating, last in the order, has only begun by 50 °C.

## The critical temperatures are outputs

The limits of the thermoneutral zone can be read from the sweep. They are properties of the animal *and* the
chamber, see [States, thresholds and traits](../manual/states_traits.md):

```@example mammal
function critical_temperatures(animal; kw...)
    function respond_in(air_temperature)
        environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature), kw...)
        environment = (; environment_pars = example_environment_pars(), environment_vars)
        thermoregulate(animal, environment, BiophysicalBehaviour.initial_physiological_state(animal, environment_vars))
    end
    outputs = [respond_in(T * u"°C") for T in air]
    resting = thermoregulation(animal)
    basal = resting.minimum_heat_flow
    lower = findfirst(out -> out.energy_flows.metabolic_heat_flow < 1.01 * basal, outputs)
    upper = findfirst(out -> out.thermoregulation.core_temperature > resting.core_temperature.reference, outputs)
    (air[lower] * u"°C", air[upper] * u"°C")
end

markdown_table(["Conditions", "Lower critical temperature", "Upper critical temperature"], [ # hide
    ("still air, 0.1 m/s", critical_temperatures(mammal)...), # hide
    ("wind, 2 m/s", critical_temperatures(mammal; wind_speed = 2.0u"m/s")...), # hide
    ("humid, 80 %", critical_temperatures(mammal; relative_humidity = 0.8)...), # hide
]) # hide
```

Wind moves both limits upward by several degrees, and humid air moves them down. Neither is a trait of the
animal.

## Comparison with NicheMapR

The points in the figure above are `endoR_devel(THERMOREG = 1, TREGMODE = 1)` with the same limits, written by
`docs/src/data/nichemapr_reference.R`. At the air temperatures that both were run for:

```@example mammal
at(T) = sweep[findfirst(==(T), air)]
rows = map(eachrow(reference)[1:3:end]) do r
    out = at(r.air_temperature_C)
    s = out.thermoregulation
    (r.air_temperature_C,
     round(watts(out.energy_flows.metabolic_heat_flow); digits = 1), round(r.metabolic_W; digits = 1),
     round(s.axis_ratio_b; digits = 1), r.shape_b,
     round(ustrip(u"W/m/K", s.flesh_conductivity); digits = 1), r.flesh_conductivity,
     round(s.pant; digits = 1), round(r.pant; digits = 1))
end
markdown_table(["Air (°C)", "Metabolic rate (W)", "NicheMapR", "Axis ratio", "NicheMapR", "Flesh conductivity", "NicheMapR", "Panting", "NicheMapR"], rows) # hide
```

```@example mammal
difference = [100 * (watts(at(r.air_temperature_C).energy_flows.metabolic_heat_flow) / r.metabolic_W - 1) for r in eachrow(reference)]
(mean_absolute = sum(abs, difference) / length(difference), largest = maximum(abs, difference))
```

The differences in metabolic rate, in per cent, have two sources. One is the heat budget itself, about 1 % in
the cold, discussed in the documentation of HeatExchange.jl. The other is the loop: both stop at the first
state in which the required rate exceeds the minimum less the tolerance, and since each step of a response
changes the required rate by a few watts, two implementations that differ by a fraction of a watt before a step
can stop one step apart.

## Changing the animal

The responses are ranges of parameters, and each can be changed. A thicker coat, and one that can be raised:

```@example mammal
function furred(depth, raised)
    depth, raised = u"m"(depth), u"m"(raised)
    traits = example_heat_exchange_traits(; insulation_pars = example_insulation_pars(;
        insulation_depth_dorsal = depth, insulation_depth_ventral = depth, insulation_depth_compressed = depth))
    regulation = example_thermoregulation_limits(;
        insulation_depth_dorsal = raised, insulation_depth_ventral = raised,
        insulation_depth_dorsal_ref = depth, insulation_depth_ventral_ref = depth,
        insulation_depth_dorsal_max = raised, insulation_depth_ventral_max = raised,
        insulation_step = depth == raised ? 0.0 : 0.05)
    fibres = insulation_pars(traits).dorsal
    coat = CompositeInsulation(FibrousLayer(raised, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3"))
    Organism(Body(shape_pars(traits), coat), OrganismTraits(Endotherm(), traits, BehavioralTraits(; thermoregulation = regulation)))
end

coats = ("2 mm, fixed" => mammal, "10 mm, fixed" => furred(10.0u"mm", 10.0u"mm"), "10 mm, raised to 20 mm" => furred(10.0u"mm", 20.0u"mm"))
fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate (W)"; size = (700, 400)) # hide
ax2 = Axis(fig[1, 2]; xlabel = "Air temperature (°C)", ylabel = "Fur depth (mm)") # hide
for (label, animal) in coats # hide
    outputs = [respond(animal, T * u"°C") for T in air] # hide
    lines!(ax, air, [watts(out.energy_flows.metabolic_heat_flow) for out in outputs]; linewidth = 2, label) # hide
    lines!(ax2, air, [ustrip(u"mm", out.thermoregulation.insulation_depth) for out in outputs]; linewidth = 2) # hide
end # hide
axislegend(ax; position = :rt, labelsize = 10) # hide
fig # hide
```

A coat that can be raised extends the thermoneutral zone downward: the animal with 10 mm of fur that it can
double has the lower critical temperature of the thicker coat and the upper critical temperature of the thinner
one.
