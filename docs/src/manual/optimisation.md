# Thermoregulation by optimisation

[Rules in sequence](endotherm_rules.md) say what an animal does and in what order. The alternative is to say what
the animal is trying to achieve and what each response costs, and let a solver find the responses.
[`IPOPTControl`](@ref) poses thermoregulation as a constrained optimisation, a nonlinear program, and solves it
with the interior-point solver IPOPT, using derivatives of the heat budget from automatic differentiation.

This page describes the problem and how to use it. [How the optimisation is built](nlp.md) describes the
machinery.

!!! warning "In development"
    The optimiser is the newest part of the package. It agrees with the rule-based controller where the two can
    be compared, but can return a poor local optimum under high heat loads, and does not yet report whether
    the solver converged. See [What to check](#What-to-check) and [Present limits](#Present-limits).

```@setup optimisation
using Main.FigureHelpers
using CairoMakie
```

## The problem

In the terms of [Behaviour as control](control.md#Kinds-of-controller), the variables are of two kinds, the
controls ``u`` and states ``x`` of [The formulation](nlp.md#The-formulation).

**Control variables** are what the animal sets:

| Variable | Scope | Bounds |
|:--|:--|:--|
| metabolic heat production | whole organism | from `minimum_heat_flow` to `metabolic_heat_flow_max_multiplier` times it |
| panting rate | whole organism | `panting.pant.reference` to `panting.pant.max` |
| flesh conductivity | each part | `flesh_conductivity.reference` to `flesh_conductivity.max` |
| skin wetness | each part | `skin_wetness.reference` to `skin_wetness.max` |

**State variables** follow from them through the physics:

| Variable | Scope | Bounds |
|:--|:--|:--|
| core temperature | whole organism | `core_temperature.reference` to `core_temperature.max` |
| skin temperature | each part | a little below air temperature to a little above the maximum core temperature |
| surface temperature of the insulation | each part | the same |

The bounds are the limits of [`ThermoregulationLimits`](@ref), the same the rule-based controller uses. An
organism of ``N`` parts has ``3 + 4N`` variables. A single `Body` is the case ``N = 1``.

### Constraints: the physics

The heat budget must hold. For each part, two equations from `part_surface_residuals` of HeatExchange.jl are
zero:

1. the heat arriving at the surface from inside equals the heat leaving it;
2. the skin temperature is consistent with the heat conducted to it from the core.

For the whole organism, heat produced less heat lost in breathing equals the sum over parts of the heat each
passes to its surface:

```math
Q_{\text{gen}} - Q_{\text{resp}} = \sum_{\text{parts}} Q_{\text{gen,net}}
```

and heat production cannot fall below the minimum, scaled for core temperature:

```math
Q_{\text{gen}} \ge Q_{\text{min}}\, Q_{10}^{(T_c - T_{c,\text{ref}})/10}
```

That is ``2N + 2`` constraints, the equations that `solve_metabolic_rate` drives to zero by iteration. Here the
solver satisfies them while choosing the control variables.

### Objective: the costs

What is minimised is a weighted sum of squared departures from where the animal would rather be, each divided
by its range so the terms are comparable:

```math
J = w_c \left(\frac{T_c - T_{c,\text{ref}}}{\Delta T_c}\right)^2
  + w_m \left(\frac{Q_{\text{gen}} - Q_{\text{min}}}{\Delta Q}\right)^2
  + w_p \left(\frac{p - p_{\text{ref}}}{\Delta p}\right)^2
  + w_w \left(\frac{\bar{w} - w_{\text{ref}}}{\Delta w}\right)^2
  + w_k \left(\frac{\bar{k} - k_{\text{ref}}}{\Delta k}\right)^2
  + w_g \left(\frac{(T_c - \bar{T}_s) - \Delta T^*}{\Delta T_c}\right)^2
```

| Weight | Field of `ThermoregulationLimits` | Default | Penalises |
|:--|:--|:--|:--|
| ``w_c`` | `core_temperature_weight` | 1 | core temperature away from the setpoint |
| ``w_m`` | `metabolic_heat_weight` | 0.1 | heat production above the minimum |
| ``w_p`` | `panting_weight` | 1 | panting |
| ``w_w`` | `skin_wetness_weight` | 1 | mean skin wetness above its resting value |
| ``w_k`` | `flesh_conductivity_weight` | 0 | mean flesh conductivity above its resting value |
| ``w_g`` | `gradient_weight` | 0 | a core-to-skin temperature difference away from `target_core_skin_gradient` |

Each term is a [gradient of information](gradients.md#A-gradient-of-information), a difference between a state
and a reference, and the sum has the form of a free energy, see
[Gradients and control](gradients.md#Free-energy). The weights say how much the animal minds each, and take the
place of the order of the rules: a response with a low weight is used freely, one with a high weight last.

## Using it

The controller is a field of the limits. Here is the 65 kg animal of
[Endotherm thermoregulation by rules](endotherm_rules.md) with either controller, and any weights as keywords:

```@example optimisation
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
using Setfield: @set

function mammal(control; weights = (;), limits...)
    physiology_traits = example_heat_exchange_traits()
    regulation = example_thermoregulation_limits(; limits...)
    regulation = @set regulation.control = control
    regulation = ConstructionBase.setproperties(regulation, weights)
    fibres = insulation_pars(physiology_traits).dorsal
    insulation = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density),
                                     FatLayer(0.0, 901.0u"kg/m^3"))
    Organism(Body(shape_pars(physiology_traits), insulation),
             OrganismTraits(Endotherm(), physiology_traits, BehavioralTraits(; thermoregulation = regulation)))
end
import ConstructionBase # hide

environment(air_temperature) = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature)),
)
nothing # hide
```

`thermoregulate` is called as before:

```@example optimisation
optimiser = mammal(IPOPTControl(); skin_wetness_max = 0.5)
warm = environment(30.0u"°C")
init = BiophysicalBehaviour.initial_physiological_state(optimiser, warm.environment_vars)
out = thermoregulate(optimiser, warm, init)
keys(out)
```

The output is that of the multi-part solver, see [Bodies of many parts](multipart.md#Solving): whole-organism
values, and under `parts` a NamedTuple for each part. A single `Body` has the one part `body`.

```@example optimisation
part = out.parts.body
markdown_table(["", "At 30 °C"], [ # hide
    ("metabolic rate", out.metabolic_heat_flow), # hide
    ("core temperature", celsius(out.core_temperature)), # hide
    ("panting multiplier", out.panting_rate), # hide
    ("flesh conductivity", part.flesh_conductivity), # hide
    ("skin wetness", part.skin_wetness), # hide
    ("skin temperature", celsius(part.skin_temperature)), # hide
    ("heat lost in breathing", part.flows.respiration_heat_flow), # hide
    ("heat lost from wet skin", part.flows.skin_evaporation_heat_flow), # hide
]) # hide
```

Where the rules would have raised the core temperature to its limit before panting, the optimiser has left the
core almost where it was and used panting and sweating together, each a little. That is what equal weights on
the three say.

## A sweep, with a warm start

Over conditions that change smoothly, each solution is a good first guess for the next. An
[`IPOPTSolverCache`](@ref) keeps the solution and multipliers of the previous solve and starts the next from
them. It is passed to the method of `thermoregulate` that names the strategy and controller:

```@example optimisation
function sweep(animal, air_temperatures)
    control = control_strategy(animal)
    first_environment = environment(first(air_temperatures))
    init = BiophysicalBehaviour.initial_physiological_state(animal, first_environment.environment_vars)
    cache = IPOPTSolverCache(control, animal, first_environment, init)
    map(air_temperatures) do air_temperature
        thermoregulate(Endotherm(), control, animal, environment(air_temperature), init; cache)
    end
end

air = [T * u"°C" for T in 0.0:2.5:40.0]
optimised = sweep(optimiser, air)
nothing # hide
```

For comparison, the rule-based controller on the same animal, with its change of posture switched off, since
the optimiser does not change geometry:

```@example optimisation
rules = mammal(RuleBasedSequentialControl(); skin_wetness_max = 0.5, skin_wetness_step = 0.005, axis_ratio_max = 1.1)
ruled = map(air) do air_temperature
    surroundings = environment(air_temperature)
    thermoregulate(rules, surroundings, BiophysicalBehaviour.initial_physiological_state(rules, surroundings.environment_vars))
end
fig = Figure(size = (760, 620)) # hide
x = ustrip.(u"°C", air) # hide
panels = ( # hide
    ("Metabolic rate (W)", r -> watts(r.energy_flows.metabolic_heat_flow), r -> watts(r.metabolic_heat_flow)), # hide
    ("Core temperature (°C)", r -> ustrip(u"°C", r.thermoregulation.core_temperature), r -> ustrip(u"°C", r.core_temperature)), # hide
    ("Flesh conductivity (W/m/K)", r -> ustrip(u"W/m/K", r.thermoregulation.flesh_conductivity), r -> ustrip(u"W/m/K", r.parts.body.flesh_conductivity)), # hide
    ("Panting multiplier", r -> r.thermoregulation.pant, r -> r.panting_rate), # hide
    ("Skin wetness", r -> r.thermoregulation.skin_wetness, r -> r.parts.body.skin_wetness), # hide
    ("Skin temperature (°C)", r -> ustrip(u"°C", r.thermoregulation.skin_temperature), r -> ustrip(u"°C", r.skin_temperature)), # hide
) # hide
for (i, (label, f_rules, f_optimised)) in enumerate(panels) # hide
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel = "Air temperature (°C)", ylabel = label) # hide
    lines!(ax, x, f_rules.(ruled); linewidth = 2, color = :grey40, label = "rules") # hide
    lines!(ax, x, f_optimised.(optimised); linewidth = 2, color = :darkorange, label = "optimisation") # hide
    i == 1 && axislegend(ax; position = :rt, labelsize = 11) # hide
end # hide
fig # hide
```

In the cold the two agree: there is nothing to decide, and the heat budget fixes the metabolic rate. Both
vasodilate over the same range of air temperature, the only response without a cost in the default weights.
Above it they part. The rules spend core temperature first, then panting. The optimiser spreads the load.

## What the weights do

At 30 °C, with the weight on panting or on skin wetness reduced tenfold:

```@example optimisation
choices = (
    "equal weights" => (;),
    "panting cheap" => (; panting_weight = 0.1),
    "sweating cheap" => (; skin_wetness_weight = 0.1),
    "core temperature cheap" => (; core_temperature_weight = 0.01),
)
rows = map(choices) do (label, weights)
    animal = mammal(IPOPTControl(); skin_wetness_max = 0.5, weights)
    result = thermoregulate(animal, warm, init)
    (label, celsius(result.core_temperature), result.panting_rate, result.parts.body.skin_wetness,
     result.parts.body.flows.respiration_heat_flow, result.parts.body.flows.skin_evaporation_heat_flow)
end
markdown_table(["Weights", "Core temperature", "Panting", "Skin wetness", "Heat lost in breathing", "Heat lost from skin"], rows) # hide
```

| To represent | Set |
|:--|:--|
| an animal that pants before it sweats: birds, dogs, rabbits | `skin_wetness_weight > panting_weight` |
| an animal that sweats first: humans, horses | `skin_wetness_weight < panting_weight` |
| an animal that lets its core temperature drift, saving water: camels, many desert birds | a low `core_temperature_weight` |
| vasodilation held in reserve | a non-zero `flesh_conductivity_weight` |
| reluctance to raise metabolic rate | a higher `metabolic_heat_weight` |

See [A bird: rules and optimisation](../tutorials/budgerigar.md) for weights chosen for one species.

## What to check

Test the result before using it. `thermoregulate` returns the point at which IPOPT stopped, converged or not.
Each part carries the residuals of its heat budget at that point, and at a solution they are zero:

```@example optimisation
function is_feasible(result; power = 0.01u"W", temperature = 0.01u"K")
    all(result.parts) do part
        flows = part.flows
        abs(flows.residual_energy_balance - flows.residual_internal_conduction) < power &&
            abs(flows.residual_skin_temperature) < temperature
    end
end
all(is_feasible, optimised)
```

A feasible point need not be the best. The problem is not convex, and IPOPT finds a local optimum. Two signs of
a poor one: a metabolic rate far above the minimum in the heat, and a flesh conductivity back at its resting
value when the animal is hot. Starting from a previous solution, as the cache does, and raising
`metabolic_heat_weight`, both help.

## Rules or optimisation

| | Rules in sequence | Optimisation |
|:--|:--|:--|
| The animal is described by | an order of responses and a mode | a cost for each response |
| Responses are used | one at a time, each to its limit | together |
| Fur depth and posture | changed | not changed |
| Result | always one, the same from any start | a local optimum, which can depend on the start |
| Tuning | step sizes, the mode | weights |
| Corresponds to | NicheMapR's `endoR` | nothing in NicheMapR |
| Use when | the order is known, or NicheMapR is to be reproduced | the order is in question, or responses trade off |

## Present limits

- **Fur and posture are fixed.** Piloerection and uncurling rebuild the geometry of the body, and are not yet
  variables. Posture is to become a change of pose of a body of many parts, see
  [Bodies of many parts](multipart.md#Thermoregulating). In the cold the optimiser therefore leaves the coat
  at its resting depth where the rules raise it, see [A bird: rules and optimisation](../tutorials/budgerigar.md).
- **Local optima under high heat loads.** Because the heat lost in breathing rises with metabolic rate, there is
  a second family of solutions in which the animal raises its metabolic rate several-fold and pants the heat
  away. For the animal on this page IPOPT finds them above about 40 °C.
- **No report of convergence.** When no feasible point exists, as in air hotter than the animal can regulate
  against, the returned point does not satisfy the heat budget. Use a check like `is_feasible`.
- **One core.** All parts share one regulated core temperature.
- **One surface per part.** Each part has one skin and one surface temperature, so a difference between the
  sunlit and shaded sides of a part is averaged.
