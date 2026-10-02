# Behaviour as control

This package could as well have been called BiophysicalControl.jl. What it calls behaviour, the things an
organism does to keep its body temperature, metabolic rate and water loss within bounds, is the regulation of a
physical system towards a reference. That is the subject of control theory, and its terms describe the package
exactly.

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
| **plant**: the system being controlled | the heat and water budget of the body | `heat_balance`, `solve_temperature`, `solve_metabolic_rate` and the multi-part solvers of HeatExchange.jl |
| **state**: what the plant is doing | body temperature, or the metabolic rate and water loss needed to hold it | the output of a solve: `core_temperature`, `metabolic_heat_flow`, the mass flows |
| **setpoint** or **reference**: where the state should be | a preferred body temperature, a regulated core temperature, a minimum metabolic rate | `target_temperature`, `active_temperature_min`, `active_temperature_max`, `core_temperature`, `minimum_heat_flow` |
| **error**: the difference between state and setpoint | too hot, too cold, or producing less heat than the least that can be produced | the conditions tested in each loop, and the objective of the optimiser |
| **controller**: what turns an error into an action | the thermoregulatory system | an [`AbstractControlStrategy`](@ref): [`RuleBasedSequentialControl`](@ref) or [`IPOPTControl`](@ref) |
| **actuator**: what acts on the plant | the effectors of physiology, and movement | the [`Effector`](@ref) types and the behaviour functions, see [Introduction](introduction.md#What-changes) |
| **saturation**: the limit of an actuator | fur can be raised only so far, a body is only so wet | the `max` and `reference` of each [`SteppedParameter`](@ref) |
| **disturbance**: what moves the state from outside | weather | the environment, from Microclimate.jl |

The loop is closed because the result of an action is observed before the next is chosen. Every controller here
solves the heat budget, compares the state with its reference, acts, and solves again.

## Two plants

The same loop describes an ectotherm and an endotherm, with the plant solved for a different unknown.

For an **ectotherm** the state is body temperature. The reference is a range of temperatures at which it is
active, and a preferred temperature within that range. The actuators are mostly its position: shade, height and
depth. The error is a temperature difference.

For an **endotherm** the core temperature is itself the reference, and is held. The state is the metabolic rate
needed to hold it. There is no error while that rate is at or above the least the animal can produce, because
producing more heat is then all that is needed, and the heat budget gives the amount directly. The error
appears when the required rate falls *below* the minimum: the animal makes more heat than it can lose. The
actuators are those that increase heat loss, and the error is a difference in power.

This is why the lower critical temperature of an endotherm needs no controller and the upper one does, and why
[`thermoregulate`](@ref) for an endotherm does nothing in the cold but solve the heat budget once.

## Kinds of controller

### Threshold control

The simplest feedback controller switches an actuator fully on or off when the state crosses a threshold, as a
thermostat does. Control theory calls it bang-bang or on-off control. Several ectotherm responses are of this
kind: a lizard is either perpendicular to the sun or it is not, pressed to the ground or not, underground or
above it. The thresholds are traits, such as `basking_temperature_min` and `emerge_temperature_min`.

### Stepped control

Most actuators here are graded, and are moved one step at a time until the error is gone. A
[`SteppedParameter`](@ref) holds the current value, the limits and the size of a step. Shade is increased by 3 %
of full shade, flesh conductivity by 0.1 W/m/K, the panting multiplier by 0.1. This is integral control in its
most elementary form: the action accumulates for as long as the error persists, and stops when it is gone or
the actuator saturates.

The size of the step is the resolution of the controller. A large step overshoots: the animal ends a little
cooler, or a little wetter, than it needed to be. A small step costs more solutions of the heat budget.

### Sequential control

With several actuators, a rule is needed for which to move. [`RuleBasedSequentialControl`](@ref) moves them in
a fixed order of priority, exhausting each before starting the next: cheap responses first, those that cost
water last. The order is a hypothesis about the organism, and an
[`AbstractThermoregulationMode`](@ref) varies it, by letting panting or sweating begin alongside a rise in core
temperature.

This is the control scheme of NicheMapR, and is described in
[Ectotherm thermoregulation](ectotherm.md) and [Endotherm thermoregulation by rules](endotherm_rules.md).

### Optimal control

The alternative is to state what the organism is trying to achieve and let a solver find the actions.
[`IPOPTControl`](@ref) poses thermoregulation as a constrained optimisation: the error, and the use of each
costly actuator, are terms of an objective to be minimised; the heat budget is a set of constraints that must
hold; the limits of the actuators are bounds. All actuators move at once. Weights on the terms of the objective
take the place of the order of the rules.

The formulation separates the variables as control theory does. *Control variables* are those the organism
sets: flesh conductivity, skin wetness, panting rate, metabolic heat production. *State variables* are those
that follow from them through the physics: core, skin and surface temperatures. See
[Thermoregulation by optimisation](optimisation.md).

## Negative and positive feedback

A controller is stable when its action reduces the error that caused it. Shade lowers the body temperature that
prompted the search for it. Sweating increases the heat loss whose shortfall prompted it. That is negative
feedback.

Some responses of an endotherm to heat carry a positive feedback as well. Raising the core temperature widens
the gradient to the environment, which increases heat loss, but it also raises the metabolic rate through the
``Q_{10}`` effect, which increases the heat to be lost. Panting loses heat by evaporation and makes heat by
muscular work. Both are included: the minimum metabolic rate that the controller compares against is rescaled
each time the core temperature or the panting rate rises,

```math
Q_{\text{min}} = (Q_{\text{basal}} + Q_{\text{pant}})\, Q_{10}^{(T_c - T_{c,\text{ref}})/10}
```

so that the reference itself moves away as the controller approaches it. In an environment hot enough the
gain from each step is less than its cost, the loop does not close, and every actuator saturates. That
condition is a prediction of the model: the environment exceeds what the animal can regulate against
(Kearney et al. 2021).

## Open and closed loops

Not everything here is feedback. The activity period of an organism, [`Diurnal`](@ref), [`Nocturnal`](@ref) or
[`Crepuscular`](@ref), is open-loop control: it depends on the position of the sun, not on the state of the
animal. So is the reset at the start of each hour, when an ectotherm is returned to its reference position
before the loop begins. [`ResponsiveActivity`](@ref) allows the activity period to be any function of the
conditions.

## Where the dynamics are

Control theory is mostly concerned with systems that have dynamics: a state with inertia, a controller that can
oscillate or lag. The controllers documented here act on a steady-state plant. Each solution of the heat budget
is the state the organism would settle to, and the loop runs to convergence within one time step, usually an
hour. That is a good approximation for a small animal, whose body temperature follows its surroundings within
minutes, and for the metabolic rate of an endotherm.

For a large ectotherm, or for anything on a time scale shorter than its thermal time constant, the plant has
dynamics, and the controller must act on a body temperature that is still changing. Transient heat budgets and
controllers for them, including proportional control, are in development and will be documented with them.
