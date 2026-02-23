abstract type AbstractBehaviourParameters end

# =============================================================================
# Thermoregulation Mode Types
# =============================================================================

"""
    AbstractThermoregulationMode

Abstract supertype for thermoregulation modes.

Modes determine which effectors are available during thermoregulation:
- `Core`: Basic thermoregulation only (piloerection, uncurl, vasodilate, hyperthermia)
- `CoreAndPantingFirst`: Adds panting during hyperthermia
- `CorePantingSweatingFirst`: Adds both panting and sweating
"""
abstract type AbstractThermoregulationMode end

"""
    CoreFirst <: AbstractThermoregulationMode

Core thermoregulation comes first in the sequence.
"""
struct CoreFirst <: AbstractThermoregulationMode end

"""
    CoreAndPantingFirst <: AbstractThermoregulationMode

Simultaneous core and panting thermoregulation come first in the sequence.
"""
struct CoreAndPantingFirst <: AbstractThermoregulationMode end

"""
    CorePantingSweatingFirst <: AbstractThermoregulationMode

Simultaneous core, panting and sweating thermoregulation come first in the sequence.
"""
struct CorePantingSweatingFirst <: AbstractThermoregulationMode end

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
    RuleBasedSequentialControl{M,T,I} <: AbstractControlStrategy

Rule-based sequential controller (priority-based bang-bang control).

Applies thermoregulation behaviors in a fixed priority order, with each
effector operating in discrete steps until saturation before moving to
the next. The loop iterates until heat balance is achieved within tolerance.

This is the default control strategy, mimicking biological thermoregulation
where organisms engage responses in a prioritized sequence based on
metabolic cost and effectiveness.

# Fields
- `mode::M`: Thermoregulation mode (`Core`, `CoreAndPantingFirst`, or `CorePantingSweatingFirst`)
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
HeatExchange.options(t::OrganismTraits) = HeatExchange.options(t.physiology)

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
