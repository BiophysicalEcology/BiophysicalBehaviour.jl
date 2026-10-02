# API

## Thermoregulation

```@docs
thermoregulate
```

## Organisms and traits

```@docs
OrganismTraits
BehavioralTraits
thermal_strategy
heat_exchange
behavior
thermoregulation
activity_period
control_strategy
```

## Thermal strategies

```@docs
AbstractThermalStrategy
Endotherm
Ectotherm
Heterotherm
```

## Controllers

```@docs
AbstractControlStrategy
RuleBasedSequentialControl
IPOPTControl
IPOPTSolverCache
MultipartNLP
AbstractThermoregulationMode
CoreFirst
CoreAndPantingFirst
CorePantingSweatingFirst
```

## Limits

```@docs
SteppedParameter
ThermoregulationLimits
InsulationLimits
PantingLimits
EctothermBehavioralLimits
BurrowShadeMode
MinShadeOnly
AdaptiveBurrowShade
MaxShadeOnly
```

## Activity and state

```@docs
ActivityPeriod
Diurnal
Nocturnal
Crepuscular
CombinedActivity
ResponsiveActivity
is_active
OrganismState
Resting
Basking
Active
```

## Environments

```@docs
AvailableEnvironments
interpolate_environment
solve_body_temperature
```

## Ectotherm behaviours

```@docs
seek_shade
avoid_shade
climb
descend
select_depth
reset_position
increment_target_temperature
darken
lighten
orient_perpendicular
orient_parallel
press_to_ground
```

## Endotherm responses

```@docs
Effector
Piloerect
Uncurl
Vasodilate
Hyperthermia
Pant
Sweat
effect
piloerect
uncurl
vasodilate
hyperthermia
pant
sweat
```

## Bodies of many parts

```@docs
PartSelector
WholeBody
ByName
Compartment
part_names
select_names
map_parts
foldl_parts
set_part
map_part_physiology
couplings
organism_compartment_graph
LungPart
panting_capacity
is_lung_part
unwrap_physiology
broadcast_physiology
physiology
part_physiology
lung_part
lung_physiology
pant_selector
solve_multipart_metabolic_rate
part_surface_setups
precompute_view_partition
ShapeCache
precompute_shape_cache
refresh
```

## Example constructors

```@docs
example_thermoregulation_limits
example_behavioral_traits
example_organism_traits
example_ectotherm_behavioral_limits
example_ectotherm_behavioral_traits
example_ectotherm_organism_traits
```

## Other types

```@docs
AbstractBehavior
AbstractMovementBehavior
AbstractTemperatureRegulation
NullBehavior
BurrowTemperatureRegulation
```

## Internal

```@docs
BiophysicalBehaviour.initial_physiological_state
BiophysicalBehaviour.reset_warm_start!
BiophysicalBehaviour.q10_scale
BiophysicalBehaviour.rebuild_body
```
