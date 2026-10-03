# Introduction

[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) solves the heat budget of one organism
in one state in one place. It returns the body temperature, or the metabolic rate and water loss needed to hold
one. It has no opinion about whether that answer is acceptable to the organism.

Organisms do. A lizard whose body temperature would be 48 °C in the open moves into shade. A mammal whose heat
budget would balance only at a metabolic rate below its minimum sends blood to its skin, lets its core warm,
pants and sweats. BiophysicalBehaviour.jl computes these responses.

Each has the same form: a change to the organism, or to the environment it is in, followed by another solution
of the heat budget. The package adds no physics. It adds:

- the traits that say what an organism can change, and how far;
- the targets and thresholds that say when it should, see [States, thresholds and traits](states_traits.md);
- the controllers that decide what to change, see [Behaviour as control](control.md).

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
| tolerating a warmer body | the core temperature defended, or the target body temperature | [`hyperthermia`](@ref), [`increment_target_temperature`](@ref) |
| panting | the volume of air breathed | [`pant`](@ref), [`Pant`](@ref) |
| sweating, licking, wallowing | the fraction of the skin that is wet | [`sweat`](@ref), [`Sweat`](@ref) |

The first three change where the organism is. The rest change one of its parameters. See
[Gradients and control](gradients.md#Three-ways-to-act) for why that matters.

## An organism that behaves

An organism in HeatExchange.jl is a body and a set of heat-exchange traits. Here it is a body and an
[`OrganismTraits`](@ref), which holds:

| Field | Content |
|:--|:--|
| `thermal_strategy` | [`Ectotherm`](@ref), [`Endotherm`](@ref) or [`Heterotherm`](@ref). Decides which unknown the heat budget is solved for, and which responses apply |
| `heat_exchange` | the `HeatExchangeTraits` of HeatExchange.jl, or one for each part of a body of many parts |
| `behavior` | a [`BehavioralTraits`](@ref): the thermoregulatory limits and the activity period |

```@example introduction
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

traits = example_organism_traits()
shape = shape_pars(heat_exchange(traits))
mammal = Organism(Body(shape, Naked()), traits)
thermal_strategy(mammal)
```

The limits are a struct, [`ThermoregulationLimits`](@ref) for an endotherm and
[`EctothermBehavioralLimits`](@ref) for an ectotherm. They hold the range of each thing the organism can
change, as a [`SteppedParameter`](@ref), and the controller that is to use them, see
[Parameters](parameters.md):

```@example introduction
limits = thermoregulation(mammal)
limits.flesh_conductivity
```

```@example introduction
control_strategy(mammal)
```

## One function

[`thermoregulate`](@ref) is the entry point. It dispatches on the thermal strategy of the organism and then on
its control strategy:

| Organism | Controller | What `thermoregulate` does |
|:--|:--|:--|
| [`Ectotherm`](@ref) | [`RuleBasedSequentialControl`](@ref) | finds the position, posture and colour that bring body temperature into the range for activity, for one hour of a microclimate, see [Ectotherm thermoregulation](ectotherm.md) |
| [`Endotherm`](@ref) | [`RuleBasedSequentialControl`](@ref) | applies physiological responses in a fixed order until the heat budget balances at or above the minimum metabolic rate, see [Endotherm thermoregulation by rules](endotherm_rules.md) |
| [`Endotherm`](@ref) | [`IPOPTControl`](@ref) | finds the combination of responses that minimises a cost, subject to the heat budget, see [Thermoregulation by optimisation](optimisation.md) |

An endotherm can also be given a set of available environments. It then first chooses where to be, as an
ectotherm does, and thermoregulates physiologically there, see
[A desert mammal through the year](../tutorials/endotherm_year.md).

## What is not here yet

This documentation describes steady-state thermoregulation. Planned or in development:

- **Pose.** For bodies of many parts, raising wings or ears and bringing limbs in to the body, to replace
  curling and uncurling as the model of posture, see [Bodies of many parts](multipart.md#Thermoregulating).
- **Transients.** Behaviour while body temperature is changing. The transient heat budget itself is to be part
  of HeatExchange.jl.
- **Control primitives, life stages and the arrest of development**, on another branch of this package.
- **Dynamic programming.** Behaviour chosen to maximise expected fitness given reserves of energy and water and
  the risks of each place (Mangel and Clark 1988), see [Gradients and control](gradients.md#Ways-to-decide).

## The ecosystem

| Package | Role here |
|:--|:--|
| [BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl) | the body, of one part or many, and the areas, lengths and view factors derived from it |
| [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) | the heat and water budget that every response is evaluated with |
| [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl) | the properties of air and water vapour |
| [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) | the environments an organism can choose between, hour by hour, see [Activity and available environments](environments.md) |
| [MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) | those environments for any place, from gridded data |
| [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) | metabolic rates and body proportions as functions of size |
| [ThermalPhysiology.jl](https://github.com/BiophysicalEcology/ThermalPhysiology.jl) | the consequences of the body temperatures computed here for performance and survival |

See [Environments and the ecosystem](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/ecosystem)
in the documentation of HeatExchange.jl for how they fit together.
