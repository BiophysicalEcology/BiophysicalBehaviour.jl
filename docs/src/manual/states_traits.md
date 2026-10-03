# States, thresholds and traits

The documentation of HeatExchange.jl sets out why the things it computes
[are not traits](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/units_traits#State-variables-are-not-traits).
Body temperature is a *state variable*. Metabolic rate and water loss are *processes*. Each belongs to an
organism in a place at a time. What can be a trait is a particular value of a state or process at which
something happens to the organism, or at which it acts: a *threshold*.

This package is where thresholds live. A heat budget returns a body temperature of 60 °C without comment. The
model of behaviour says that 60 °C will not do, and what the organism does instead.

```@setup states_traits
using Main.FigureHelpers
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
```

## Four kinds of term

A model of an organism as a dynamical system has four kinds of term (Kearney et al. 2021). Behaviour adds to
each.

| Term | In HeatExchange.jl | Added here |
|:--|:--|:--|
| **state variable** | core, skin and surface temperatures | where the organism is: shade, height, depth; what it is doing: [`Resting`](@ref), [`Basking`](@ref), [`Active`](@ref) |
| **process** | heat flows, metabolic rate, water loss | none: the processes of the heat budget, evaluated in the state behaviour has chosen |
| **environmental variable** | air, ground and sky temperatures, wind, humidity, radiation | the set of environments within reach, see [`AvailableEnvironments`](@ref) |
| **parameter** | the properties of the organism | the range over which each can be changed, and the thresholds at which it is |

## Threshold traits of an ectotherm

For an ectotherm the thresholds are body temperatures, fields of [`EctothermBehavioralLimits`](@ref):

```@example states_traits
limits = example_ectotherm_behavioral_limits()
markdown_table(["Threshold", "Value", "What happens there"], [ # hide
    ("`critical_temperature_min`", celsius(limits.critical_temperature_min), "below it the animal cannot move; it retreats if it can"), # hide
    ("`emerge_temperature_min`", celsius(limits.emerge_temperature_min), "the body temperature at which it can leave its retreat"), # hide
    ("`basking_temperature_min`", celsius(limits.basking_temperature_min), "above it the animal may be out basking, but not yet foraging"), # hide
    ("`active_temperature_min`", celsius(limits.active_temperature_min), "the lower bound of the range for activity"), # hide
    ("`target_temperature`", celsius(limits.target_temperature.reference), "the preferred temperature: above it, the animal begins to respond to heat"), # hide
    ("`active_temperature_max`", celsius(limits.active_temperature_max), "the upper bound of the range for activity"), # hide
    ("`critical_temperature_max`", celsius(limits.critical_temperature_max), "above it the animal cannot move; used when choosing a depth to retreat to"), # hide
]) # hide
```

These are the `CT_min`, `T_RB_min`, `T_B_min`, `T_F_min`, `T_pref`, `T_F_max` and `CT_max` of the NicheMapR
ectotherm model (Kearney and Porter 2020). They are properties of the animal: measured in a thermal gradient or
observed in the field, independent of the weather, with dimensions of temperature.

Together they divide the state variable into bands, and the band decides the consequence:

```@example states_traits
using CairoMakie # hide
fig = Figure(size = (720, 170)) # hide
ax = Axis(fig[1, 1]; xlabel = "Body temperature (°C)", yticksvisible = false, yticklabelsvisible = false, ygridvisible = false) # hide
c(x) = ustrip(u"°C", x) # hide
edges = [0.0, c(limits.critical_temperature_min), c(limits.basking_temperature_min), c(limits.active_temperature_min), c(limits.active_temperature_max), c(limits.critical_temperature_max), 46.0] # hide
names = ["immobile", "resting", "basking", "active", "resting", "immobile"] # hide
colours = [:grey60, RGBf(0.45, 0.55, 0.75), RGBf(0.98, 0.75, 0.35), RGBf(0.80, 0.30, 0.25), RGBf(0.45, 0.55, 0.75), :grey60] # hide
for i in 1:6 # hide
    poly!(ax, Rect2f(edges[i], 0, edges[i + 1] - edges[i], 1); color = (colours[i], 0.8)) # hide
    text!(ax, (edges[i] + edges[i + 1]) / 2, 0.5; text = names[i], align = (:center, :center), fontsize = 12) # hide
end # hide
vlines!(ax, [c(limits.target_temperature.reference)]; color = :black, linestyle = :dash) # hide
text!(ax, c(limits.target_temperature.reference) + 0.3, 0.9; text = "target", align = (:left, :center), fontsize = 11) # hide
xlims!(ax, 0, 46) # hide
ylims!(ax, 0, 1) # hide
fig # hide
```

An hour in which the animal can be active is an hour in which it can feed. The count of such hours is a main
output of a mechanistic niche model, and it comes from two thresholds applied to a computed state.

## Threshold traits of an endotherm

For an endotherm the thresholds are mostly limits on processes and parameters, fields of
[`ThermoregulationLimits`](@ref):

```@example states_traits
limits = example_thermoregulation_limits()
markdown_table(["Threshold", "Value", "Role"], [ # hide
    ("`minimum_heat_flow`", limits.minimum_heat_flow, "the least heat the animal can produce: the reference of the controller"), # hide
    ("`core_temperature.reference`", celsius(limits.core_temperature.reference), "the core temperature that is defended"), # hide
    ("`core_temperature.max`", celsius(limits.core_temperature.max), "the highest core temperature tolerated"), # hide
    ("`flesh_conductivity.max`", limits.flesh_conductivity.max, "flesh conductivity at full vasodilation"), # hide
    ("`axis_ratio_factor.max`", limits.axis_ratio_factor.max, "the most elongated posture"), # hide
    ("`panting.pant.max`", limits.panting.pant.max, "the greatest multiple of resting ventilation"), # hide
    ("`skin_wetness.max`", limits.skin_wetness.max, "the greatest fraction of the skin that can be wet"), # hide
]) # hide
```

The minimum metabolic rate is a threshold on a process: a solution of the heat budget below it is the signal
that the animal must act. The others are the ends of the ranges over which a parameter can be changed.

## Plastic parameters

A parameter in HeatExchange.jl has one value. Many are not constant in a living animal: fur is raised and
flattened, skin is wetted, blood is sent to the skin. A [`SteppedParameter`](@ref) describes a parameter the
organism can change:

```@example states_traits
limits.flesh_conductivity
```

| Field | Meaning |
|:--|:--|
| `current` | the value now. A state, changed by the controller |
| `reference` | the value in the resting or heat-conserving condition |
| `max` | the limit of the response |
| `step` | the increment by which a rule-based controller changes it |

`reference` and `max` are traits. They bound the plasticity of the parameter and can be measured: the flesh
conductivity of a vasoconstricted and a vasodilated limb, the depth of a flattened and a raised coat. `current`
is a state. `step` belongs to the numerical method, not the animal, though it should not be finer than the
animal's own control.

So a parameter of the heat budget becomes, here, a state with trait-valued bounds. In that sense physiological
and behavioural thermoregulation are the same thing: both move a parameter of the physics between limits that
are properties of the organism.

## What is not a threshold trait

- **Critical air temperatures.** The lower and upper critical temperatures of an endotherm, and its
  thermoneutral zone, are *air* temperatures. They belong to the environment, and hold only for the wind,
  humidity, radiation and posture of the measurement (Kearney et al. 2021). Here they are outputs, see
  [A mammal across air temperatures](../tutorials/mammal.md#The-critical-temperatures-are-outputs).
- **Not to be confused with critical body temperatures.** The critical thermal limits of an ectotherm are
  *body* temperatures, and are traits. The difference is whether the value is a property of the state of the
  organism or of the conditions that produce it.
- **Patterns of behaviour.** Hours of activity, the depth of a burrow used, and the shade selected are states,
  however often they are reported as traits. This package computes them.

## The four classes of functional trait

In the classification of Kearney et al. (2021):

| Class | In BiophysicalBehaviour.jl |
|:--|:--|
| **parameter** | the `reference` and `max` of each [`SteppedParameter`](@ref); the cost of panting; the weights of the optimisation, which state the relative cost of each response |
| **threshold** | the body temperatures of [`EctothermBehavioralLimits`](@ref); `minimum_heat_flow` and the maximum core temperature of [`ThermoregulationLimits`](@ref) |
| **model** | the thermal strategy, [`Ectotherm`](@ref) or [`Endotherm`](@ref); the activity period; the capability flags `can_climb`, `can_retreat_underground`, `can_seek_shade`, `can_pant` and the rest; the [`AbstractThermoregulationMode`](@ref); which part of the body holds the lungs |
| **estimation** | preferred temperatures from a thermal gradient, field body temperatures of active animals, evaporative water loss against air temperature in a metabolic chamber: observations from which thresholds and limits are estimated, given the conditions of measurement |

The model traits are types and flags, on which methods are chosen. The parameter and threshold traits are the
data of a limits struct, see [Parameters](parameters.md). The states are what [`thermoregulate`](@ref)
returns.
