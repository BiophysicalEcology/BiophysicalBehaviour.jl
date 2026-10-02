# Introduction

[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) solves the heat budget of one organism
in one state in one place. Given a body, a set of traits and an environment, it returns the body temperature, or
the metabolic rate and water loss needed to hold a body temperature. It has no opinion about whether that answer
is acceptable to the organism.

Organisms do have an opinion. A lizard whose body temperature would be 48 °C in the open moves into the shade. A
mammal whose heat budget would only balance at a metabolic rate below the minimum it can produce sends blood to
its skin, lets its core temperature rise, pants and sweats. BiophysicalBehaviour.jl computes these responses.

Every one of them has the same form: a change to the organism, or to the environment it is in, followed by
another solution of the heat budget. The package therefore adds no physics. It adds the traits that say what an
organism can change and how far, the targets and thresholds that say when it should, and the controllers that
decide what to change.

```@setup introduction
using Main.FigureHelpers
```

## What changes

| Response | What changes | Function or effector |
|:--|:--|:--|
| seeking or leaving shade | the solar radiation reaching the organism, and the air, ground and sky temperatures around it | [`seek_shade`](@ref), [`avoid_shade`](@ref) |
| retreating underground | the whole environment, replaced by the soil at a chosen depth | [`select_depth`](@ref) |
| climbing | the air temperature and wind speed, taken higher in the profile | [`climb`](@ref) |
| orienting to the sun | the silhouette area that intercepts the direct beam | [`orient_perpendicular`](@ref), [`orient_parallel`](@ref) |
| changing colour | the solar absorptivity of the skin | [`darken`](@ref), [`lighten`](@ref) |
| raising or flattening fur or feathers | the depth of the insulation | [`piloerect`](@ref), [`Piloerect`](@ref) |
| curling up or stretching out | the ratio of the axes of the body, and so its surface area | [`uncurl`](@ref), [`Uncurl`](@ref) |
| vasodilation and vasoconstriction | the thermal conductivity of the flesh | [`vasodilate`](@ref), [`Vasodilate`](@ref) |
| tolerating a warmer body | the core temperature that is defended, or the target body temperature | [`hyperthermia`](@ref), [`increment_target_temperature`](@ref) |
| panting | the volume of air breathed | [`pant`](@ref), [`Pant`](@ref) |
| sweating, licking, wallowing | the fraction of the skin that is wet | [`sweat`](@ref), [`Sweat`](@ref) |

The first three change where the organism is. The rest change one of its parameters. See
[Gradients and control](gradients.md) for why that distinction matters.

## An organism that behaves

An organism in HeatExchange.jl is a body and a set of heat-exchange traits. Here it is a body and an
[`OrganismTraits`](@ref), which holds three things:

| Field | Content |
|:--|:--|
| `thermal_strategy` | [`Ectotherm`](@ref), [`Endotherm`](@ref) or [`Heterotherm`](@ref). Decides which unknown the heat budget is solved for, and which set of responses applies |
| `heat_exchange` | the `HeatExchangeTraits` of HeatExchange.jl, or one of them for each part of a body of many parts |
| `behavior` | a [`BehavioralTraits`](@ref): the thermoregulatory limits and the activity period |

```@example introduction
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

traits = example_organism_traits()
shape = shape_pars(heat_exchange(traits))
mammal = Organism(Body(shape, Naked()), traits)
thermal_strategy(mammal)
```

The thermoregulatory limits are themselves a struct, [`ThermoregulationLimits`](@ref) for an endotherm and
[`EctothermBehavioralLimits`](@ref) for an ectotherm. They hold the range of each thing the organism can change,
as a [`SteppedParameter`](@ref), and the controller that is to use them:

```@example introduction
limits = thermoregulation(mammal)
limits.flesh_conductivity
```

```@example introduction
control_strategy(mammal)
```

## One function

[`thermoregulate`](@ref) is the entry point for all of it. It dispatches first on the thermal strategy of the
organism and then on its control strategy:

| Organism | Controller | What `thermoregulate` does |
|:--|:--|:--|
| [`Ectotherm`](@ref) | [`RuleBasedSequentialControl`](@ref) | finds the position, posture and colour that bring body temperature into the range for activity, for one hour of a microclimate, see [Ectotherm thermoregulation](ectotherm.md) |
| [`Endotherm`](@ref) | [`RuleBasedSequentialControl`](@ref) | applies physiological responses in a fixed order until the heat budget balances at or above the minimum metabolic rate, see [Endotherm thermoregulation by rules](endotherm_rules.md) |
| [`Endotherm`](@ref) | [`IPOPTControl`](@ref) | finds the combination of responses that minimises a cost, subject to the heat budget, see [Thermoregulation by optimisation](optimisation.md) |

An endotherm can also be given a set of available environments, in which case it first chooses where to be, as
an ectotherm does, and then thermoregulates physiologically in that place.

## What is not here yet

This documentation describes steady-state thermoregulation. Changes of pose for bodies of many parts, such as
raising wings or ears or bringing limbs in to the body, are planned, and will replace curling and uncurling as
the model of posture. Heat budgets through time, a library of control
primitives, life stages and the arrest of development are in development on another branch of the package, and
will be documented when they are merged.

## The ecosystem

| Package | Role here |
|:--|:--|
| [BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl) | the body, of one part or many, and the areas, lengths and view factors derived from it |
| [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) | the heat and water budget that every response is evaluated with |
| [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl) | the properties of air and water vapour |
| [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) | the environments an organism can choose between, hour by hour, see [Activity and available environments](environments.md) |
| [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) | metabolic rates and body proportions as functions of size |
| [ThermalPhysiology.jl](https://github.com/BiophysicalEcology/ThermalPhysiology.jl) | the consequences of the body temperatures computed here for performance and survival |
