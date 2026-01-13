abstract type AbstractBehaviourParameters end

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
cooler air or moving underground into a warmer/cooler burrow.
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
- `activity::A`: Activity period (Diurnal, Nocturnal, etc.)
"""
Base.@kwdef struct BehavioralTraits{T,A}
    thermoregulation::T
    activity::A = Diurnal()
end

"""
    OrganismTraits{S,P,B} <: AbstractFunctionalTraits

Combined physiological and behavioral traits for an organism.

This allows an `Organism` to carry both its physical/physiological properties
(from HeatExchange) and its behavioral capabilities (from BiophysicalBehaviour)
in a single traits object.

# Fields
- `thermal_strategy::S`: Thermal strategy (Endotherm, Ectotherm, Heterotherm)
- `physiology::P`: Physiological traits (HeatExchangeTraits)
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
    physiology::P
    behavior::B
end

# =============================================================================
# Forwarding methods for physiology accessors
# =============================================================================

# Forward all physiology accessor methods to the physiology field
HeatExchange.shapepars(t::OrganismTraits) = HeatExchange.shapepars(t.physiology)
HeatExchange.insulationpars(t::OrganismTraits) = HeatExchange.insulationpars(t.physiology)
function HeatExchange.conductionpars_external(t::OrganismTraits)
    HeatExchange.conductionpars_external(t.physiology)
end
function HeatExchange.conductionpars_internal(t::OrganismTraits)
    HeatExchange.conductionpars_internal(t.physiology)
end
HeatExchange.convectionpars(t::OrganismTraits) = HeatExchange.convectionpars(t.physiology)
HeatExchange.radiationpars(t::OrganismTraits) = HeatExchange.radiationpars(t.physiology)
HeatExchange.evaporationpars(t::OrganismTraits) = HeatExchange.evaporationpars(t.physiology)
HeatExchange.hydraulicpars(t::OrganismTraits) = HeatExchange.hydraulicpars(t.physiology)
HeatExchange.respirationpars(t::OrganismTraits) = HeatExchange.respirationpars(t.physiology)
HeatExchange.metabolismpars(t::OrganismTraits) = HeatExchange.metabolismpars(t.physiology)

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
    physiology(t::OrganismTraits)
    physiology(o::Organism)

Get the physiological traits from an OrganismTraits or Organism.
"""
physiology(t::OrganismTraits) = t.physiology
physiology(o::Organism) = physiology(HeatExchange.traits(o))

# =============================================================================
# BehavioralTraits accessors
# =============================================================================

"""
    thermoregulation(t::BehavioralTraits)

Get the thermoregulation limits from BehavioralTraits.
"""
thermoregulation(t::BehavioralTraits) = t.thermoregulation

"""
    activity(t::BehavioralTraits)

Get the activity period from BehavioralTraits.
"""
activity(t::BehavioralTraits) = t.activity

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
    activity(t::OrganismTraits)
    activity(o::Organism)

Get the activity period from an OrganismTraits or Organism.
"""
activity(t::OrganismTraits) = activity(t.behavior)
activity(o::Organism) = activity(HeatExchange.traits(o))
