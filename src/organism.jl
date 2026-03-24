abstract type AbstractBehaviourParameters end

# =============================================================================
# Control Strategy Types
# =============================================================================

"""
    AbstractControlStrategy

Abstract supertype for thermoregulation control strategies.

Control strategies determine how the thermoregulation loop solves for
heat balance. Different strategies may use different algorithmic approaches.
"""
abstract type AbstractControlStrategy end

"""
    AbstractThermoregulationMode

Abstract supertype for thermoregulation modes used by `RuleBasedSequentialControl`.
"""
abstract type AbstractThermoregulationMode end

"""
    CoreFirst <: AbstractThermoregulationMode

Basic thermoregulation: allow core temperature to rise before panting or sweating.
"""
struct CoreFirst <: AbstractThermoregulationMode end

"""
    CoreAndPantingFirst <: AbstractThermoregulationMode

Allows panting in parallel with core temperature rise (hyperthermia).
"""
struct CoreAndPantingFirst <: AbstractThermoregulationMode end

"""
    CorePantingSweatingFirst <: AbstractThermoregulationMode

Allows both panting and sweating in parallel with core temperature rise.
"""
struct CorePantingSweatingFirst <: AbstractThermoregulationMode end

"""
    RuleBasedSequentialControl{M,T,I} <: AbstractControlStrategy

Rule-based sequential controller (priority-based bang-bang control).

Applies thermoregulation behaviors in a fixed priority order, with each
effector operating in discrete steps until saturation before moving to
the next. The loop iterates until heat balance is achieved within tolerance.

This is the default control strategy, mimicking biological thermoregulation
where organisms engage responses in a prioritized sequence based on
metabolic cost and effectiveness.

# Fields
- `mode::M`: Thermoregulation mode (`CoreFirst`, `CoreAndPantingFirst`, or `CorePantingSweatingFirst`).
- `tolerance::T`: Fraction below Q_minimum allowed
- `max_iterations::I`: Maximum iterations before warning
"""
Base.@kwdef struct RuleBasedSequentialControl{M<:AbstractThermoregulationMode,T,I} <: AbstractControlStrategy
    mode::M = CoreFirst()
    tolerance::T = 0.005
    max_iterations::I = 1000
end

"""
    PDEControl <: AbstractControlStrategy

Partial differential equation-based control strategy.

Uses a PDE formulation to solve the thermoregulation problem, allowing
for spatially-resolved temperature distributions and continuous control
of effectors.

!!! warning
    This control strategy is not yet implemented.
"""
struct PDEControl <: AbstractControlStrategy 
    # Add any reqired settings here
end

# =============================================================================
# Behavior Types
# =============================================================================

"""
    AbstractBehavior

Abstract supertype for organism behaviors.

Behaviors respond to physiological state and external stimuli,
to either control *exposure* to the external stimuli, or modify
physiological parameters to change the *effect* of the external stimuli.

- Response to environmental information that leads to:
    - Changes in metabolic parameters
    - Changes in the environment or position within the environment
"""
abstract type AbstractBehavior end

"""
    AbstractMovementBehavior

Abstract supertype behaviors that modify location based
on inforation from the environment, such as moving up into
cooler air or moving underground into a warmer/cooler retreat.
"""
abstract type AbstractMovementBehavior <: AbstractBehavior end

abstract type AbstractTemperatureRegulation end

struct NullBehavior <: AbstractBehavior end

initialise_state(::NullBehavior) = ()

struct BurrowTemperatureRegulation{T} <: AbstractTemperatureRegulation
    burrowat::T
    emergeat::T
end

abstract type ActivityPeriod end

struct Diurnal <: ActivityPeriod end
struct Nocturnal <: ActivityPeriod end
struct Crepuscular <: ActivityPeriod end

# =============================================================================
# Organism State
# =============================================================================

"""
    OrganismState

Abstract supertype for the instantaneous activity state of an organism.

Concrete subtypes mirror NicheMapR's `ACT` output column:
- [`Resting`](@ref) — underground or thermally unable to be active (ACT = 0)
- [`Basking`](@ref) — above ground, warming up; `T_bask ≤ Tb < T_active_min` (ACT = 1)
- [`Active`](@ref) — above ground, within activity thermal window;
  `T_active_min ≤ Tb ≤ T_active_max` (ACT = 2)
"""
abstract type OrganismState end

"Resting state: underground or outside the thermal window for surface activity."
struct Resting <: OrganismState end

"Basking state: above ground but below `T_active_min`; warming toward activity temperature."
struct Basking <: OrganismState end

"Active state: above ground within the activity thermal window `[T_active_min, T_active_max]`."
struct Active <: OrganismState end

"""
    CombinedActivity{T<:Tuple} <: ActivityPeriod

    CombinedActivity(periods...)

Combines multiple activity periods: active if active in any of them.

# Example
```julia
# Active during both day and dawn/dusk
CombinedActivity(Diurnal(), Crepuscular())
```
"""
struct CombinedActivity{T<:Tuple} <: ActivityPeriod
    periods::T
end
CombinedActivity(periods::ActivityPeriod...) = CombinedActivity(periods)

# =============================================================================
# Thermal Strategy
# =============================================================================

"""
    AbstractThermalStrategy

Abstract supertype for organism thermal regulation strategies.
"""
abstract type AbstractThermalStrategy end

"""
    Endotherm <: AbstractThermalStrategy

Organism that generates internal heat to maintain body temperature.
"""
struct Endotherm <: AbstractThermalStrategy end

"""
    Ectotherm <: AbstractThermalStrategy

Organism that relies on external heat sources for body temperature regulation.
"""
struct Ectotherm <: AbstractThermalStrategy end

"""
    Heterotherm <: AbstractThermalStrategy

Organism that can switch between endothermic and ectothermic strategies.
"""
struct Heterotherm <: AbstractThermalStrategy end

"""
    ResponsiveActivity <: ActivityPeriod

    ResponsiveActivity(isactive)

# Arguments

- `isactive` a `Function` or functor that recieves a `ModelParEnvironment`
    object with the current system state, and decides whether to be "active"
    or "innactive".
"""
struct ResponsiveActivity{F} <: ActivityPeriod
    isactive::F
end

# =============================================================================
# Combined Traits
# =============================================================================

"""
    BehavioralTraits{T,A}

Behavioral traits of an organism.

# Fields
- `thermoregulation::T`: Thermoregulation limits and parameters
- `activity_period::A`: Activity phase (Diurnal, Nocturnal, etc.)
"""
Base.@kwdef struct BehavioralTraits{T,A}
    thermoregulation::T
    activity_period::A = Diurnal()
end

"""
    OrganismTraits{S,P,B} <: AbstractFunctionalTraits

Combined physiological and behavioral traits for an organism.

This allows an `Organism` to carry both its physical/physiological properties
(from HeatExchange) and its behavioral capabilities (from BiophysicalBehaviour)
in a single traits object.

# Fields
- `thermal_strategy::S`: Thermal strategy (Endotherm, Ectotherm, Heterotherm)
- `heat_exchange::P`: Heat exchange traits (HeatExchangeTraits) — morphology and physiology
- `behavior::B`: Behavioral traits (BehavioralTraits)

# Example
```julia
traits = OrganismTraits(
    Endotherm(),
    HeatExchangeTraits(...),
    BehavioralTraits(thermoregulation_limits, Diurnal())
)
organism = Organism(body, traits)
```
"""
struct OrganismTraits{
    S<:AbstractThermalStrategy,
    P<:HeatExchange.AbstractFunctionalTraits,
    B<:BehavioralTraits,
} <: HeatExchange.AbstractFunctionalTraits
    thermal_strategy::S
    heat_exchange::P
    behavior::B
end

# =============================================================================
# Forwarding methods for physiology accessors
# =============================================================================

# Forward all HeatExchange accessor methods to the heat_exchange field
HeatExchange.shapepars(t::OrganismTraits) = HeatExchange.shapepars(t.heat_exchange)
HeatExchange.insulationpars(t::OrganismTraits) = HeatExchange.insulationpars(t.heat_exchange)
function HeatExchange.conductionpars_external(t::OrganismTraits)
    HeatExchange.conductionpars_external(t.heat_exchange)
end
function HeatExchange.conductionpars_internal(t::OrganismTraits)
    HeatExchange.conductionpars_internal(t.heat_exchange)
end
HeatExchange.convectionpars(t::OrganismTraits) = HeatExchange.convectionpars(t.heat_exchange)
HeatExchange.radiationpars(t::OrganismTraits) = HeatExchange.radiationpars(t.heat_exchange)
HeatExchange.evaporationpars(t::OrganismTraits) = HeatExchange.evaporationpars(t.heat_exchange)
HeatExchange.hydraulicpars(t::OrganismTraits) = HeatExchange.hydraulicpars(t.heat_exchange)
HeatExchange.respirationpars(t::OrganismTraits) = HeatExchange.respirationpars(t.heat_exchange)
HeatExchange.metabolismpars(t::OrganismTraits) = HeatExchange.metabolismpars(t.heat_exchange)
HeatExchange.options(t::OrganismTraits) = HeatExchange.options(t.heat_exchange)

# =============================================================================
# OrganismTraits accessors
# =============================================================================

"""
    thermal_strategy(t::OrganismTraits)
    thermal_strategy(o::Organism)

Get the thermal strategy from an OrganismTraits or Organism.
"""
thermal_strategy(t::OrganismTraits) = t.thermal_strategy
thermal_strategy(o::Organism) = thermal_strategy(HeatExchange.traits(o))

"""
    behavior(t::OrganismTraits)
    behavior(o::Organism)

Get the behavioral traits from an OrganismTraits or Organism.
"""
behavior(t::OrganismTraits) = t.behavior
behavior(o::Organism) = behavior(HeatExchange.traits(o))

"""
    heat_exchange(t::OrganismTraits)
    heat_exchange(o::Organism)

Get the `HeatExchangeTraits` (morphology + physiology) from an OrganismTraits or Organism.
"""
heat_exchange(t::OrganismTraits) = t.heat_exchange
heat_exchange(o::Organism) = heat_exchange(HeatExchange.traits(o))

# =============================================================================
# BehavioralTraits accessors
# =============================================================================

"""
    thermoregulation(t::BehavioralTraits)

Get the thermoregulation limits from BehavioralTraits.
"""
thermoregulation(t::BehavioralTraits) = t.thermoregulation

"""
    activity_period(t::BehavioralTraits)

Get the activity period from BehavioralTraits.
"""
activity_period(t::BehavioralTraits) = t.activity_period

# =============================================================================
# OrganismTraits accessors (forward to behavior)
# =============================================================================

"""
    thermoregulation(t::OrganismTraits)
    thermoregulation(o::Organism)

Get the thermoregulation limits from an OrganismTraits or Organism.
"""
thermoregulation(t::OrganismTraits) = thermoregulation(t.behavior)
thermoregulation(o::Organism) = thermoregulation(HeatExchange.traits(o))

"""
    activity_period(t::OrganismTraits)
    activity_period(o::Organism)

Get the activity period from an OrganismTraits or Organism.
"""
activity_period(t::OrganismTraits) = activity_period(t.behavior)
activity_period(o::Organism) = activity_period(HeatExchange.traits(o))

"""
    control_strategy(t::BehavioralTraits)

Get the control strategy from BehavioralTraits.
"""
control_strategy(t::BehavioralTraits) = t.thermoregulation.control

"""
    control_strategy(t::OrganismTraits)
    control_strategy(o::Organism)

Get the control strategy from an OrganismTraits or Organism.
"""
control_strategy(t::OrganismTraits) = control_strategy(t.behavior)
control_strategy(o::Organism) = control_strategy(HeatExchange.traits(o))
