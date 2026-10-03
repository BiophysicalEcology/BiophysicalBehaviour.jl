# Gradients and control

In [HeatExchange.jl](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/gradients) many terms of
the heat budget are flows down a difference in potential, against a resistance:

```math
\text{flow} = \frac{\text{potential difference}}{\text{resistance}}
```

Heat flows down gradients of temperature, and water vapour down gradients of vapour density. In
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) heat and water move through soil and
air the same way. These are gradients of physical potential: energy or matter flows down them whether or not
anything is watching.

This package is about a gradient of another kind.

## A gradient of information

What an organism does is also decided by a difference: between the state it senses and the state it expects or
prefers. For example, a body temperature above the preferred one, a threshold level of dehydration, or, for an
endotherm, a steady-state heat budget implying a metabolic rate lower than it can produce. That difference drives
behaviour as a temperature difference drives heat. The response is larger the larger the
difference, and continues while the difference is worth reducing.

But nothing flows down it. It is a gradient of *information*: sensed, compared with a reference and acted upon,
each step belonging to the organism, not the physics. In the terms of [Behaviour as control](control.md) it is
the error signal of a control loop.

### Free energy

The word gradient is meant literally. Under the free energy principle (Friston 2010) a living thing must stay
within the small set of states in which it remains what it is. A state outside that set is improbable, and so
surprising, and organisms act to keep surprise low. What they minimise is a bound on it, the *free energy*: at
its simplest, a sum of squared prediction errors, each weighted by its precision. It is not the physical energy
flowing through the organism.

Minimising that prediction error is *homeostasis* (Cannon 1932). Thermoregulation is the homeostasis of body
temperature, and it is what this package computes. A response made ahead of the error is *allostasis*
(Sterling 2012), as in an activity period, or a retreat before the surface becomes lethal. Homeostasis is
feedback, allostasis feed-forward.

| Free energy principle | Here |
|:--|:--|
| expected states | targets and thresholds, see [States, thresholds and traits](states_traits.md) |
| prediction error | the difference between the computed state and the target |
| precision | the weight on an error |
| action | every response in the table below |
| free energy | the objective of [Thermoregulation by optimisation](optimisation.md) |
| homeostasis | what [`thermoregulate`](@ref) does |

Action descends the gradient of free energy as heat descends the gradient of temperature, and in the opposite
sense: physics takes an organism towards equilibrium with its surroundings, and homeostasis holds it away.

Two caveats. The error is not driven to zero: a response costs water, energy, time and safety, and stops when
a further step is not worth it, see [Ways to decide](#Ways-to-decide). And the free energy principle is a
mathematical framing from a physics perspective, not a mechanism: the controllers here are rules applied in
order and a constrained optimisation.

| | Physical gradient | Gradient of information |
|:--|:--|:--|
| between | two places | a state and a target or threshold for it |
| of | temperature, vapour density, water potential | the same quantities, compared with a reference |
| what happens | heat or water flows | the organism acts |
| set by | the environment and the properties of the organism | the traits of the organism |
| computed by | HeatExchange.jl, Microclimate.jl | BiophysicalBehaviour.jl |
| vanishes when | the two places reach the same potential | the state reaches the target, a further step is not worth its cost, or the means are exhausted |

The two are coupled both ways. The physical gradients determine the state. The state, compared with its
reference, determines what the organism does. What it does changes the physical gradients.

## Three ways to act

The heat budget has exchanges, each a difference in potential against a resistance, and sources. An organism
can change any of the three. In the bond-graph terms of
[HeatExchange.jl](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/gradients#The-budget-as-a-network),
it modulates a resistor, moves to other sources of effort, or changes a source of flow.

- A **resistance**, by changing one of its own parameters. It stays put and exchanges heat differently with the
  same surroundings.
- A **gradient**, by changing the potential on either side of the organism–environment divide. Moving replaces
  the environmental potentials. Letting its own body temperature change alters the potential on the organism's
  side.
- A **source, sink or storage**: how much solar radiation it absorbs, how much heat it makes, how much air flows
  through its lungs, and, away from steady state, heat stored as body temperature lags behind its surroundings.

| Response | Acts on | By changing |
|:--|:--|:--|
| raising or flattening fur | resistance | the thickness of the insulating layer |
| vasodilation | resistance | the conductivity of the flesh between core and skin |
| curling up, stretching out | resistance | the areas exposed to the air, the sky and the ground |
| pressing to the ground | resistance | the area in contact with the substrate |
| sweating | resistance | the area of skin from which water evaporates |
| panting | a sink | the air carrying heat and vapour from the lungs: transport by moving air, not a passive resistance |
| changing colour | a source | the fraction of solar radiation absorbed |
| orienting to the sun | a source | the area that intercepts the direct beam |
| seeking shade | a source and a gradient | the solar radiation arriving, and the temperatures of air, ground and sky |
| climbing | gradient | the air temperature and wind speed |
| retreating underground | gradient | every temperature around the body, replaced by that of the soil |
| raising core or target temperature | gradient | the temperature at the organism's end of every flow |

The classification has practical consequences.

- **What it needs.** A change of resistance or source needs nothing from outside the organism, and can be
  computed from one description of the environment. That is how the endotherm controllers work. A change of
  gradient by moving needs a description of the other places the organism could be. That is what
  [`AvailableEnvironments`](@ref) provides, see [Activity and available environments](environments.md).
- **What it costs.** Moving costs time that could be spent feeding, and may expose the animal to predators.
  Evaporation costs water. A warmer core costs performance, and eventually survival. The order in which a
  rule-based controller tries responses, and the weights of an optimising one, are statements about these
  costs.

## Ways to decide

The controllers differ in how they turn the error into action, not in what the error is or in the physics the
action is tested against.

| Controller | The error is | Action |
|:--|:--|:--|
| threshold | compared with a fixed value | switched on or off |
| rules in sequence | tested for its sign | one step of the first response not exhausted, in a fixed order |
| optimisation | a term of an objective | all responses together, minimising the objective subject to the physics |
| dynamic programming (planned) | a cost in expected fitness | the action, through time, that maximises expected fitness given internal state |

In the optimisation the two kinds of gradient are the two parts of the problem. The gradients of information
are the objective: squared differences between state and target, each scaled by its range and weighted, which
is the form of a free energy. The physical gradients are the constraints: the heat budget of each part, as
residuals that must be zero. The solver then follows a third gradient, the derivative of objective and
constraints with respect to each thing the organism can change, by automatic differentiation of the physics.
See [How the optimisation is built](nlp.md).

Neither controller says when a further step is *worth it*: the order of the rules and the weights of the
objective are set by hand. Dynamic state-variable models (Mangel and Clark 1988; Clark and Mangel 2000) make
that the question. The state includes what the organism carries, its energy reserve, gut contents and water,
and each action, such as staying in a burrow, basking or foraging in sun or shade, changes that state and risks
death by predation, heat or starvation. Working backward from the end of a period, dynamic programming finds
the action at each time and state that maximises expected fitness. Thermal error is then one cost among many,
priced in the same currency as food, water and safety. The body temperatures and water losses in each place
come from the heat budget. This is planned for BiophysicalBehaviour.jl, with the energy and water budgets in
[AnimalMapper.jl](https://github.com/BiophysicalEcology/AnimalMapper.jl).
