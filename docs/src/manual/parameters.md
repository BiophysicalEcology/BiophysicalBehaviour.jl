# Parameters

The behavioural traits of an organism are held beside its heat-exchange traits in an
[`OrganismTraits`](@ref). This page lists them. For the heat-exchange traits see
[Parameters](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/parameters) in the documentation
of HeatExchange.jl.

```@setup parameters
using Main.FigureHelpers
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
```

## The structure

```
Organism
├─ body                         a Body or CompositeBody (BiophysicalGeometry.jl)
└─ traits :: OrganismTraits
   ├─ thermal_strategy          Ectotherm(), Endotherm() or Heterotherm()
   ├─ heat_exchange             HeatExchangeTraits, or a NamedTuple of them by part
   ├─ behavior :: BehavioralTraits
   │  ├─ thermoregulation       EctothermBehavioralLimits or ThermoregulationLimits
   │  │  └─ control             RuleBasedSequentialControl or IPOPTControl
   │  └─ activity_period        Diurnal(), Nocturnal(), Crepuscular(), ...
   ├─ lung_part                 the name of the part that holds the lungs
   └─ couplings                 a HeatCoupling for each join of a CompositeBody
```

| Accessor | Returns |
|:--|:--|
| [`thermal_strategy`](@ref) | the thermal strategy |
| [`heat_exchange`](@ref) | the heat-exchange traits |
| [`behavior`](@ref) | the `BehavioralTraits` |
| [`thermoregulation`](@ref) | the limits struct |
| [`control_strategy`](@ref) | the controller |
| [`activity_period`](@ref) | the activity period |
| [`lung_part`](@ref), [`couplings`](@ref) | the lung part and the couplings |
| [`physiology`](@ref) | the heat-exchange traits by part, with the lung part marked |

Each takes an `Organism` or an `OrganismTraits`.

## A range of a parameter

A [`SteppedParameter`](@ref) is a parameter that the organism can change:

| Field | Meaning |
|:--|:--|
| `current` | the present value |
| `reference` | the resting value. Defaults to `current` |
| `max` | the limit |
| `step` | the increment of one action of the rule-based controller |

## Ectotherms

[`EctothermBehavioralLimits`](@ref), with the defaults of `example_ectotherm_behavioral_limits`, which are those
of the NicheMapR ectotherm model:

```@example parameters
limits = example_ectotherm_behavioral_limits()
stepped(p) = "$(p.reference) to $(p.max), step $(p.step)" # hide
markdown_table(["Field", "Default", "Meaning"], [ # hide
    ("`shade`", stepped(limits.shade), "shade fraction"), # hide
    ("`depth`", "node 1 to the deepest, step 1", "soil node. 1 is the surface"), # hide
    ("`depth_min_underground`", limits.depth_min_underground, "shallowest node of a retreat"), # hide
    ("`height`", "node 1 to the highest, step 1", "height node"), # hide
    ("`absorptivity`", stepped(limits.absorptivity), "solar absorptivity"), # hide
    ("`pant_rate`", stepped(limits.pant_rate), "multiplier on ventilation"), # hide
    ("`target_temperature`", "$(celsius(limits.target_temperature.reference)) to $(celsius(limits.target_temperature.max)), step $(limits.target_temperature.step)", "preferred body temperature, raised towards the maximum for activity"), # hide
    ("`active_temperature_min`", celsius(limits.active_temperature_min), "lower limit for activity"), # hide
    ("`active_temperature_max`", celsius(limits.active_temperature_max), "upper limit for activity"), # hide
    ("`basking_temperature_min`", celsius(limits.basking_temperature_min), "lower limit for basking"), # hide
    ("`emerge_temperature_min`", celsius(limits.emerge_temperature_min), "soil temperature at which the animal can leave its retreat"), # hide
    ("`critical_temperature_min`", celsius(limits.critical_temperature_min), "critical thermal minimum"), # hide
    ("`critical_temperature_max`", celsius(limits.critical_temperature_max), "critical thermal maximum"), # hide
    ("`can_retreat_underground`", limits.can_retreat_underground, ""), # hide
    ("`can_seek_shade`", limits.can_seek_shade, ""), # hide
    ("`can_climb`", limits.can_climb, ""), # hide
    ("`can_solar_orient`", limits.can_solar_orient, "turn broadside to the sun to bask"), # hide
    ("`can_press_to_ground`", limits.can_press_to_ground, ""), # hide
    ("`can_change_absorptivity`", limits.can_change_absorptivity, ""), # hide
    ("`can_pant`", limits.can_pant, ""), # hide
    ("`solve_underground`", limits.solve_underground, "solve the heat budget in the retreat, or take soil temperature"), # hide
    ("`burrow_shade_mode`", limits.burrow_shade_mode, "which microclimate the retreat is in"), # hide
    ("`emerge_signal`", limits.emerge_signal, "rate of change of soil temperature needed to emerge. Zero for none"), # hide
]) # hide
```

Two more fields, `sun_orientation` and `pressed_to_ground`, record the posture chosen in the current hour.

## Endotherms

[`ThermoregulationLimits`](@ref), with the defaults of `example_thermoregulation_limits`:

```@example parameters
limits = example_thermoregulation_limits()
range_of(p) = "$(p.current) to $(p.max), step $(p.step)" # hide
markdown_table(["Field", "Default", "Meaning"], [ # hide
    ("`control`", limits.control, "the controller"), # hide
    ("`minimum_heat_flow`", limits.minimum_heat_flow, "least metabolic heat production"), # hide
    ("`insulation.dorsal`, `.ventral`", "$(limits.insulation.dorsal.reference) to $(limits.insulation.dorsal.max)", "depth of the coat. The step is a fraction of fibre length"), # hide
    ("`axis_ratio_factor`", range_of(limits.axis_ratio_factor), "ratio of the long axis of the body to the short"), # hide
    ("`flesh_conductivity`", range_of(limits.flesh_conductivity), "thermal conductivity of flesh"), # hide
    ("`core_temperature`", "$(celsius(limits.core_temperature.reference)) to $(celsius(limits.core_temperature.max)), step $(limits.core_temperature.step)", "core temperature"), # hide
    ("`panting.pant`", range_of(limits.panting.pant), "multiplier on ventilation"), # hide
    ("`panting.multiplier`", limits.panting.multiplier, "metabolic rate at full panting, as a multiple of the minimum"), # hide
    ("`skin_wetness`", range_of(limits.skin_wetness), "fraction of the skin that is wet"), # hide
]) # hide
```

The remaining fields are used only by [`IPOPTControl`](@ref), see
[Thermoregulation by optimisation](optimisation.md):

```@example parameters
markdown_table(["Field", "Default", "Meaning"], [ # hide
    ("`core_temperature_weight`", limits.core_temperature_weight, "weight on core temperature away from the setpoint"), # hide
    ("`metabolic_heat_weight`", limits.metabolic_heat_weight, "weight on heat production above the minimum"), # hide
    ("`panting_weight`", limits.panting_weight, "weight on panting"), # hide
    ("`skin_wetness_weight`", limits.skin_wetness_weight, "weight on skin wetness"), # hide
    ("`flesh_conductivity_weight`", limits.flesh_conductivity_weight, "weight on flesh conductivity"), # hide
    ("`gradient_weight`", limits.gradient_weight, "weight on the core-to-skin temperature difference"), # hide
    ("`target_core_skin_gradient`", limits.target_core_skin_gradient, "the difference that `gradient_weight` holds to"), # hide
    ("`skin_temperature_undershoot`", limits.skin_temperature_undershoot, "lower bound of skin and surface temperature, below air temperature"), # hide
    ("`skin_temperature_core_overshoot`", limits.skin_temperature_core_overshoot, "upper bound of skin and surface temperature, above the maximum core temperature"), # hide
    ("`metabolic_heat_flow_max_multiplier`", limits.metabolic_heat_flow_max_multiplier, "upper bound of heat production, as a multiple of the minimum"), # hide
    ("`minimum_normalisation_range`", limits.minimum_normalisation_range, "floor on the ranges that normalise the objective"), # hide
]) # hide
```

## Controllers

| Controller | Fields |
|:--|:--|
| [`RuleBasedSequentialControl`](@ref) | `mode`, an [`AbstractThermoregulationMode`](@ref), default `CoreFirst()`; `tolerance`, the fraction below the minimum metabolic rate that is accepted, default 0.005; `max_iterations`, default 1000 |
| [`IPOPTControl`](@ref) | `nlp_strategy`, default `MultipartNLP()`; `smoothing`, default `SmoothBound(1e-5)` |

## Example constructors

| Function | Builds |
|:--|:--|
| `example_ectotherm_behavioral_limits` | [`EctothermBehavioralLimits`](@ref) |
| `example_ectotherm_behavioral_traits` | [`BehavioralTraits`](@ref) of an ectotherm, with an `activity_period` |
| `example_ectotherm_organism_traits` | [`OrganismTraits`](@ref) of an ectotherm. Keywords go to the limits |
| `example_thermoregulation_limits` | [`ThermoregulationLimits`](@ref) |
| `example_behavioral_traits` | [`BehavioralTraits`](@ref) of an endotherm |
| `example_organism_traits` | [`OrganismTraits`](@ref) of an endotherm. Keywords go to `example_heat_exchange_traits` |

A limits struct is immutable. To change one field of an existing one, use `@set` of
[Setfield.jl](https://github.com/jw3126/Setfield.jl), as the package itself does:

```julia
using Setfield: @set
limits = @set limits.flesh_conductivity.max = 2.0u"W/m/K"
```
