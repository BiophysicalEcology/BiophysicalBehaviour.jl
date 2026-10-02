# Endotherm thermoregulation by rules

An endotherm holds its core temperature. In the cold it does so by making heat, and the heat budget of
HeatExchange.jl gives the amount. In the heat it cannot make less than its minimum, and must lose more. This
page describes the sequence of responses by which it does that, which is the thermoregulatory loop of the
NicheMapR endotherm model (Kearney et al. 2021).

```@setup endotherm_rules
using Main.FigureHelpers
using CairoMakie
```

## The aim

`solve_metabolic_rate` returns the metabolic heat production ``Q_{\text{gen}}`` at which the heat budget
balances for a given core temperature. The animal has a least rate at which it can produce heat,
``Q_{\text{min}}``, its basal or resting rate. The controller acts while

```math
Q_{\text{gen}} < Q_{\text{min}}\,(1 - \text{tolerance})
```

that is, while holding the core temperature would need the animal to produce less heat than it can. The state
is ``Q_{\text{gen}}``, the reference ``Q_{\text{min}}``, and each action increases the heat that the animal
can lose, so that the ``Q_{\text{gen}}`` required rises towards the reference.

```julia
thermoregulate(organism, environment, init)
```

| Argument | Meaning |
|:--|:--|
| `organism` | an `Organism` whose traits have the [`Endotherm`](@ref) strategy and a [`ThermoregulationLimits`](@ref) |
| `environment` | `(; environment_pars, environment_vars)`, as for HeatExchange.jl |
| `init` | first guesses `(; metabolic_heat_flow, skin_temperature, insulation_temperature)` for the solver |

## The sequence

```@example endotherm_rules
ladder_diagram(["flatten fur", "uncurl", "vasodilate", "raise core\ntemperature", "pant", "sweat"]; # hide
    title = "Required heat production below the minimum", colour = RGBf(0.95, 0.80, 0.78)) # hide
```

| | Response | Function | Effector type | What it changes |
|:--|:--|:--|:--|:--|
| 1 | flatten the fur or feathers | [`piloerect`](@ref) | [`Piloerect`](@ref) | insulation depth, down to its `reference` |
| 2 | stretch out | [`uncurl`](@ref) | [`Uncurl`](@ref) | the ratio of the long axis of the body to the short, up to its `max` |
| 3 | send blood to the skin | [`vasodilate`](@ref) | [`Vasodilate`](@ref) | flesh conductivity, up to its `max` |
| 4 | let the core warm | [`hyperthermia`](@ref) | [`Hyperthermia`](@ref) | core temperature, up to its `max` |
| 5 | pant | [`pant`](@ref) | [`Pant`](@ref) | the multiplier on ventilation |
| 6 | sweat | [`sweat`](@ref) | [`Sweat`](@ref) | the fraction of the skin that is wet |

The loop starts with the fur fully raised, the condition that conserves most heat. If the heat budget balances
there at or above ``Q_{\text{min}}``, that is the answer, and the animal is at or below its lower critical
temperature. Otherwise the first response that is not exhausted is moved one step, the heat budget is solved
again, and the test is repeated.

The order is that of cost. The first three cost nothing but the loss of their use for anything else. A rise in
core temperature costs performance and a higher metabolic rate. Panting and sweating cost water.

### The costs that feed back

Two responses raise the reference as well as the state, see [Behaviour as control](control.md). A warmer core
has a higher metabolic rate, and panting is work. After each step of either,

```math
Q_{\text{min}} = (Q_{\text{basal}} + Q_{\text{pant}})\, Q_{10}^{(T_c - T_{c,\text{ref}})/10}
```

``Q_{\text{pant}}`` rises linearly with the panting rate to `(multiplier - 1)` times the basal rate at
full panting, see [`PantingLimits`](@ref). ``Q_{10}`` is that of the `MetabolismParameters` of HeatExchange.jl.

## Modes

Many animals do not wait for their core temperature to reach its limit before panting, or for panting to reach
its limit before sweating. The [`AbstractThermoregulationMode`](@ref) of the controller sets which responses
advance together:

| Mode | With each step of core temperature | NicheMapR `TREGMODE` |
|:--|:--|:--|
| [`CoreFirst`](@ref) | nothing else: panting begins when core temperature is at its maximum, sweating when panting is | 1 |
| [`CoreAndPantingFirst`](@ref) | a step of panting | 2 |
| [`CorePantingSweatingFirst`](@ref) | a step of panting and a step of sweating | 3 |

The mode is a hypothesis about the animal. It is a model trait, see
[States, thresholds and traits](states_traits.md).

## An example

The 65 kg animal of [Get started](../get_started.md). `example_thermoregulation_limits` takes the range and
step of each response as keywords:

```@example endotherm_rules
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

function mammal(; limits...)
    physiology_traits = example_heat_exchange_traits()
    behaviour = BehavioralTraits(; thermoregulation = example_thermoregulation_limits(; limits...))
    fibres = insulation_pars(physiology_traits).dorsal
    insulation = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density),
                                     FatLayer(0.0, 901.0u"kg/m^3"))
    Organism(Body(shape_pars(physiology_traits), insulation),
             OrganismTraits(Endotherm(), physiology_traits, behaviour))
end

function respond(animal, air_temperature)
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature))
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    thermoregulate(animal, environment, BiophysicalBehaviour.initial_physiological_state(animal, environment_vars))
end
nothing # hide
```

At four air temperatures:

```@example endotherm_rules
animal = mammal()
temperatures = map(T -> T * u"°C", (10.0, 20.0, 30.0, 45.0))
outputs = [respond(animal, T) for T in temperatures]
row(label, f) = (label, (f(out) for out in outputs)...) # hide
markdown_table(["", ("$(T)" for T in temperatures)...], [ # hide
    row("metabolic rate", out -> out.energy_flows.metabolic_heat_flow), # hide
    row("axis ratio", out -> out.thermoregulation.axis_ratio_b), # hide
    row("flesh conductivity", out -> out.thermoregulation.flesh_conductivity), # hide
    row("core temperature", out -> celsius(out.thermoregulation.core_temperature)), # hide
    row("panting multiplier", out -> out.thermoregulation.pant), # hide
    row("skin wetness", out -> out.thermoregulation.skin_wetness), # hide
    row("skin temperature", out -> celsius(out.thermoregulation.skin_temperature)), # hide
    row("evaporative water loss", out -> out.mass_flows.m_evap), # hide
]) # hide
```

At 10 °C nothing is done and the metabolic rate is above the minimum. At 20 °C the animal has only uncurled.
At 30 °C it is fully vasodilated, its core has warmed to the limit and it has begun to pant. At 45 °C it
is breathing at eight times the resting rate and losing 200 g of water an hour.

The output is the `ThermoregulationOutput` of HeatExchange.jl, with the same four groups: `thermoregulation`,
`morphology`, `energy_flows` and `mass_flows`. The state of every response is in `thermoregulation`.

### The effect of the mode

```@example endotherm_rules
air = 20.0:2.0:48.0
modes = (CoreFirst(), CoreAndPantingFirst(), CorePantingSweatingFirst())
sweeps = map(modes) do mode
    animal = mammal(; thermoregulation_mode = mode, skin_wetness_step = 0.005, skin_wetness_max = 0.5)
    [respond(animal, T * u"°C") for T in air]
end
fig = Figure(size = (760, 520)) # hide
panels = ( # hide
    ("Core temperature (°C)", out -> ustrip(u"°C", out.thermoregulation.core_temperature)), # hide
    ("Panting multiplier", out -> out.thermoregulation.pant), # hide
    ("Skin wetness", out -> out.thermoregulation.skin_wetness), # hide
    ("Water loss (g/h)", out -> ustrip(u"g/hr", out.mass_flows.m_evap)), # hide
) # hide
for (i, (label, f)) in enumerate(panels) # hide
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel = "Air temperature (°C)", ylabel = label) # hide
    for (mode, sweep) in zip(modes, sweeps) # hide
        lines!(ax, air, f.(sweep); linewidth = 2, label = string(nameof(typeof(mode)))) # hide
    end # hide
    i == 1 && axislegend(ax; position = :lt, labelsize = 10) # hide
end # hide
fig # hide
```

With `CoreFirst` the animal spends its tolerance of a warm core before it spends water. With
`CorePantingSweatingFirst` it uses all three together, and at a given air temperature runs cooler and loses
more water.

### The effect of the step

The step of each response is the resolution of the controller. Coarser steps overshoot: the loop stops at the
first state in which the required heat production exceeds the minimum, and the excess is heat that the animal
must then produce.

```@example endotherm_rules
fine = mammal(; pant_step = 0.02)
coarse = mammal(; pant_step = 1.0)
markdown_table(["Panting step", "Panting multiplier at 40 °C", "Metabolic rate", "Water loss"], [ # hide
    (step, out.thermoregulation.pant, out.energy_flows.metabolic_heat_flow, out.mass_flows.m_evap) # hide
    for (step, out) in ((0.02, respond(fine, 40.0u"°C")), (1.0, respond(coarse, 40.0u"°C")))]) # hide
```

## When nothing is left

If every response reaches its limit and the required heat production is still below the minimum, the loop
stops and returns the last state. Its `metabolic_heat_flow` is then below `minimum_heat_flow`: the animal cannot
hold its core temperature in that environment by these means. That is the result, and what follows from it,
heat storage and a rising body temperature, is outside a steady-state model.

## Limits

[`ThermoregulationLimits`](@ref) holds the range of each response as a [`SteppedParameter`](@ref), with
[`InsulationLimits`](@ref) for the dorsal and ventral fur and [`PantingLimits`](@ref) for panting and its cost.
It also holds the weights used by [Thermoregulation by optimisation](optimisation.md), which this controller
ignores. See [Parameters](parameters.md) for every field.

A response is switched off by giving it no range: a `max` equal to its `current` value.

## Bodies of many parts

For an organism on a `CompositeBody` the same loop runs with the responses applied part by part: vasodilation
and sweating to every part, panting to the part that holds the lungs. Fur and posture are not yet changed.

For a body of many parts, posture will be a change of *pose*: raising or lowering wings or ears, bringing the
limbs in to the body or holding them away from it. That changes which surfaces are hidden under joins and what
each part sees of the sky, the ground and its neighbours, through the geometry of BiophysicalGeometry.jl. It is to
replace the change of axis ratio by which a single body curls and uncurls.
See [Bodies of many parts](multipart.md).
