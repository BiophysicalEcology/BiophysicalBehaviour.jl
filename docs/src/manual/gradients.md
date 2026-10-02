# Gradients and control

In [HeatExchange.jl](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/gradients) almost every
term of the heat budget is a flow down a gradient of a potential, against a resistance:

```math
\text{flow} = \frac{\text{potential difference}}{\text{resistance}}
```

Heat flows down gradients of temperature, water vapour down gradients of vapour density. In
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) heat and water move through soil and
air in the same way, down gradients of temperature and water potential. These are gradients of physical
potential. What flows down them is energy or matter, and it flows whether or not anything is watching.

This package is about a gradient of another kind.

## A gradient of information

An organism responds to the difference between the state it is in and a state it would rather be in: a body
temperature above the preferred one, or a required metabolic rate below the least it can produce. That
difference drives behaviour as a temperature difference drives a flow of heat. The response is larger the
larger the difference, it continues while the difference persists, and it stops when the difference is gone.

But nothing flows down it. It is a gradient of *information*. It must be sensed, compared with a reference and
acted upon, and each of those steps belongs to the organism, not to the physics. In the terms of
[Behaviour as control](control.md) it is the error of a control loop.

| | Physical gradient | Gradient of information |
|:--|:--|:--|
| between | two places | a state and a target or threshold for that state |
| of | temperature, vapour density, water potential | the same quantities, compared with a reference |
| what happens | heat or water flows | the organism acts |
| set by | the environment and the properties of the organism | the traits of the organism |
| computed by | HeatExchange.jl, Microclimate.jl | BiophysicalBehaviour.jl |
| vanishes when | the two places reach the same potential | the state reaches the target, or the means are exhausted |

The two are coupled in one direction through the physics and in the other through behaviour. The physical
gradients determine the state. The state, compared with its reference, determines what the organism does. What
it does changes the physical gradients.

## Two ways to act

The equation for a flow has two terms, and an organism can change either.

It can change a **resistance**, by changing one of its own parameters. It stays where it is and exchanges heat
differently with the same surroundings.

It can change a **gradient**, by changing the potentials at one end or the other. Moving replaces the
potentials around it. Letting its own temperature change moves the potential at its own end.

| Response | Acts on | By changing |
|:--|:--|:--|
| raising or flattening fur | resistance | the thickness of the insulating layer |
| vasodilation | resistance | the conductivity of the flesh between core and skin |
| curling up, stretching out | resistance | the area through which every flow passes |
| pressing to the ground | resistance | the area in contact with the substrate |
| sweating | resistance | the fraction of the skin from which water evaporates |
| panting | resistance | the volume of air carrying heat and vapour from the lungs |
| changing colour | a source | the fraction of solar radiation absorbed |
| orienting to the sun | a source | the area that intercepts the direct beam |
| seeking shade | a source and a gradient | the solar radiation arriving, and the temperatures of the air, ground and sky |
| climbing | gradient | the air temperature and wind speed |
| retreating underground | gradient | every temperature around the body, replaced by that of the soil |
| raising core or target temperature | gradient | the temperature at the organism's end of every flow |

Solar radiation is a source of heat, not a flow down a gradient, so colour and orientation act on neither
term: they change how much of the source is taken up.

The classification has practical consequences. A change of resistance needs nothing from outside the organism,
and can be computed from one description of the environment. That is how the endotherm controllers work: they
take a single environment and change only parameters. A change of gradient by moving needs a description of
the other places the organism could be. That is what [`AvailableEnvironments`](@ref) provides, see
[Activity and available environments](environments.md).

The cost differs too. Changing a gradient by moving costs time that might have been spent feeding, and may
expose the animal to predators. Changing a resistance to evaporation costs water. Changing the gradient by
letting core temperature rise costs performance, and eventually survival. The order in which a rule-based
controller tries responses, and the weights of an optimising one, are statements about those costs.

## Three ways to decide

The controllers of this package differ in how they turn the error into action. They do not differ in what the
error is, or in the physics that the action is tested against.

| Controller | The error is | Action |
|:--|:--|:--|
| threshold | compared with a fixed value | switched on or off |
| rules in sequence | tested for its sign | one step of the first response that is not exhausted, in a fixed order |
| optimisation | a term of an objective | all responses together, to minimise the objective subject to the physics |

In the optimisation the two kinds of gradient appear as the two parts of the problem. The gradients of
information are the objective: squared differences between state and target, each divided by the range over
which it can vary. The physical gradients are the constraints: the heat budget of each part of the body,
written as residuals that must be zero. The solver then follows a third kind of gradient, the derivative of
the objective and constraints with respect to each thing the organism can change, obtained by automatic
differentiation of the physics. See [How the optimisation is built](nlp.md).
