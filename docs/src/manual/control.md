# Behaviour as control

This package could as well have been called BiophysicalControl.jl. What it calls behaviour, the things an
organism does to keep its body temperature, metabolic rate and water loss within bounds, is the regulation of a
physical system towards a reference. That is the subject of control theory, and its terms describe the package
exactly.

Physiology's word for the same thing is *homeostasis* (Cannon 1932): holding the internal state of a body
within the limits in which it can live. A control loop is how it is done, and minimising the error of that loop
is what it consists of, see [Gradients and control](gradients.md#Free-energy).

```@setup control
using Main.FigureHelpers
using CairoMakie
```

```@example control
control_loop_diagram() # hide
```

## The parts of the loop

| Control theory | In an organism | In this package |
|:--|:--|:--|
| **plant**: the system controlled | the heat and water budget of the body | `heat_balance`, `solve_temperature`, `solve_metabolic_rate` and the multi-part solvers of HeatExchange.jl |
| **state**: what the plant is doing | body temperature, or the metabolic rate and water loss needed to hold it | the output of a solve |
| **setpoint** or **reference**: where the state should be | a preferred body temperature, a regulated core temperature, a minimum metabolic rate | `target_temperature`, `active_temperature_min`, `active_temperature_max`, `core_temperature`, `minimum_heat_flow` |
| **error**: state minus setpoint | too hot, too cold, or making less heat than the least it can | the conditions tested in each loop, and the objective of the optimiser |
| **controller**: turns an error into an action | the thermoregulatory system | an [`AbstractControlStrategy`](@ref): [`RuleBasedSequentialControl`](@ref) or [`IPOPTControl`](@ref) |
| **actuator**: acts on the plant | the effectors of physiology, and movement | the [`Effector`](@ref) types and behaviour functions, see [Introduction](introduction.md#What-changes) |
| **saturation**: the limit of an actuator | fur can be raised only so far, a body is only so wet | the `max` and `reference` of each [`SteppedParameter`](@ref) |
| **disturbance**: moves the state from outside | weather | the environment, from Microclimate.jl |

The loop is closed because the result of an action is observed before the next is chosen. Every controller
here solves the heat budget, compares state with reference, acts, and solves again.

## Two plants

The same loop describes an ectotherm and an endotherm, with the plant solved for a different unknown.

For an **ectotherm** the state is body temperature. The reference is a range of temperatures for activity, and
a preferred temperature within it. The actuators are mostly its position: shade, height, depth. The error is a
temperature difference.

For an **endotherm** the core temperature is itself the reference, and is held. The state is the metabolic
rate needed to hold it. There is no error while that rate is at or above the least the animal can produce: the
heat budget gives the amount directly. The error appears when the required rate falls *below* the minimum, and
the animal makes more heat than it can lose. The actuators are those that increase heat loss, and the error is
a difference in power.

This is why the lower critical temperature of an endotherm needs no controller and the upper one does, and why
[`thermoregulate`](@ref) for an endotherm does nothing in the cold but solve the heat budget once.

## Kinds of controller

| Kind | How it acts | Here |
|:--|:--|:--|
| **threshold** (bang-bang) | switches an actuator fully on or off when the state crosses a threshold, as a thermostat does | a lizard is perpendicular to the sun or not, pressed to the ground or not, underground or above it. The thresholds are traits, such as `basking_temperature_min` |
| **stepped** | moves a graded actuator one step at a time until the error is gone or the actuator saturates. Integral control at its simplest | every [`SteppedParameter`](@ref): shade by 3 %, flesh conductivity by 0.1 W/m/K, panting by 0.1 |
| **sequential** | with several actuators, moves them in a fixed order, exhausting each before the next | [`RuleBasedSequentialControl`](@ref), the scheme of NicheMapR, see [Ectotherm thermoregulation](ectotherm.md) and [Endotherm thermoregulation by rules](endotherm_rules.md) |
| **optimal** | states what the organism is trying to achieve and lets a solver find the actions, all actuators at once | [`IPOPTControl`](@ref), see [Thermoregulation by optimisation](optimisation.md) |

Three points on these.

**The step is the resolution of the controller.** A large step overshoots: the animal ends a little cooler, or
wetter, than it needed to be. A small step costs more solutions of the heat budget.

**The order is a hypothesis about the organism**: cheap responses first, those that cost water last. An
[`AbstractThermoregulationMode`](@ref) varies it, by letting panting or sweating begin alongside a rise in core
temperature.

**In optimal control, weights take the place of the order.** The error, and the use of each costly actuator,
are terms of an objective. The heat budget is a set of constraints. The limits of the actuators are bounds.
*Control variables* are those the organism sets: flesh conductivity, skin wetness, panting rate, metabolic heat
production. *State variables* follow from them through the physics: core, skin and surface temperatures.

## Negative and positive feedback

A controller is stable when its action reduces the error that caused it. Shade lowers the body temperature that
prompted the search for it. Sweating increases the heat loss whose shortfall prompted it. That is negative
feedback.

Some responses of an endotherm to heat carry a positive feedback too. A warmer core widens the gradient to the
environment, which increases heat loss, but raises the metabolic rate through the ``Q_{10}`` effect. Panting
loses heat by evaporation and makes heat by muscular work. Both are included: the minimum metabolic rate that
the controller compares against is rescaled each time the core temperature or the panting rate rises,

```math
Q_{\text{min}} = (Q_{\text{basal}} + Q_{\text{pant}})\, Q_{10}^{(T_c - T_{c,\text{ref}})/10}
```

so the reference moves away as the controller approaches it. In an environment hot enough, the gain from each
step is less than its cost, the loop does not close, and every actuator saturates. That is a prediction of the
model: the environment exceeds what the animal can regulate against (Kearney et al. 2021).

## Open and closed loops

Not everything here is feedback. The activity period of an organism, [`Diurnal`](@ref), [`Nocturnal`](@ref) or
[`Crepuscular`](@ref), is open-loop control: it depends on the position of the sun, not on the state of the
animal. So is the reset at the start of each hour, when an ectotherm is returned to its reference position.
[`ResponsiveActivity`](@ref) lets the activity period be any function of the conditions. See
[Activity and available environments](environments.md#When-to-be-active).

## Where the dynamics are

Control theory is mostly about systems with dynamics: a state with inertia, a controller that can oscillate or
lag. The controllers here act on a steady-state plant. Each solution of the heat budget is the state the
organism would settle to, and the loop runs to convergence within one time step, usually an hour. That is a
good approximation for a small animal, whose body temperature follows its surroundings within minutes, and for
the metabolic rate of an endotherm.

For a large ectotherm, or on any time scale shorter than the thermal time constant, the plant has dynamics, and
the controller must act on a body temperature that is still changing. Controllers for transient heat budgets,
including proportional control, are in development.
