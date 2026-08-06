abstract type AbstractBehaviourParameters end

# =============================================================================
# Thermoregulation Mode Types
# =============================================================================

"""
    AbstractThermoregulationMode

Abstract supertype for thermoregulation modes.

Modes determine which effectors are available during thermoregulation and in what sequence.
- `Core`: fully sequential, piloerection, uncurl, vasodilate, hyperthermia, pant, sweat
- `CoreAndPantingFirst`: Adds panting during hyperthermia
- `CorePantingSweatingFirst`: Adds both panting and sweating during hyperthermia
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
- `mode::M`: Thermoregulation mode (`CoreFirst`, `CoreAndPantingFirst`, or `CorePantingSweatingFirst`).
- `tolerance::T`: Fraction below minimum_heat_flow allowed
- `max_iterations::I`: Maximum iterations before warning
"""
Base.@kwdef struct RuleBasedSequentialControl{M<:AbstractThermoregulationMode,T,I} <: AbstractControlStrategy
    mode::M = CoreFirst()
    tolerance::T = 0.005
    max_iterations::I = 1000
end

"""
    IPOPTControl <: AbstractControlStrategy

IPOPT-based nonlinear programming control strategy.

Solves the thermoregulation problem as a constrained optimisation:
minimise deviation from the setpoint core temperature subject to heat-balance
equality constraints, with all physiological effectors (flesh_conductivity,
pant, skin_wetness) as continuous decision variables.

# Fields
- `nlp_strategy`: NLP formulation. `MultipartNLP()` (default) solves a genuine
  multi-part organism — one regulated core shared across all parts, per-part surface
  temperatures and effectors, a single whole-organism respiration balance. A single
  `Body` is handled as the one-part case.
- `smoothing`: smoothing policy passed to the heat-balance physics so autodiff (AD)
  sees differentiable kinks. Defaults to `SmoothBound(1.0e-5)`; pass `HardBound()`
  to match the rule-based path's exact `abs`/`max`/`step` behaviour.

Requires `Ipopt.jl`.
"""
Base.@kwdef struct IPOPTControl{S<:HeatExchange.SmoothingStrategy} <: AbstractControlStrategy
    nlp_strategy::HeatExchange.NLPStrategy = MultipartNLP()
    smoothing::S = HeatExchange.SmoothBound(1.0e-5)
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
- [`Basking`](@ref) — above ground, warming up; `basking_temperature_min ≤ core_temperature < active_temperature_min` (ACT = 1)
- [`Active`](@ref) — above ground, within activity thermal window;
  `active_temperature_min ≤ core_temperature ≤ active_temperature_max` (ACT = 2)
"""
abstract type OrganismState end

"Resting state: underground or outside the thermal window for surface activity."
struct Resting <: OrganismState end

"Basking state: above ground but below `active_temperature_min`; warming toward activity temperature."
struct Basking <: OrganismState end

"Active state: above ground within the activity thermal window `[active_temperature_min, active_temperature_max]`."
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
    P,
    B<:BehavioralTraits,
    L,
} <: HeatExchange.AbstractFunctionalTraits
    thermal_strategy::S
    heat_exchange::P
    behavior::B
    lung_part::L
end

# Backwards-compatible constructor: a single lumped `HeatExchangeTraits` is the
# physiology of the sole `:body` part, which hosts the lung (§3.9). Multi-part
# organisms pass a per-part physiology NamedTuple as `heat_exchange` and name the
# lung part explicitly.
# The lung part name is stored as a `Val` (its type parameter `L === Val{name}`), so
# the name lives in the *type*, not as a runtime `Symbol` field. Physiology routing
# (`_wrap_lung`, `pant_selector`, `_lung_mass`) then infers concretely — the sole
# runtime→type step is this one `Val(lung_part)` at construction.
OrganismTraits(thermal_strategy::AbstractThermalStrategy, heat_exchange, behavior::BehavioralTraits;
               lung_part::Symbol=SINGLE_PART_NAME) =
    OrganismTraits(thermal_strategy, heat_exchange, behavior, Val(lung_part))

# =============================================================================
# Forwarding methods for physiology accessors
# =============================================================================

# Whole-organism physiology for the single-valued forwarding accessors. A single
# lumped `HeatExchangeTraits` is returned directly (single-body path — the common,
# hot case, fully type-stable). When per-part physiology is stored, whole-organism
# queries (metabolism, respiration, options) resolve through the lung part.
_whole_physiology(t::OrganismTraits) = _whole_physiology(t.heat_exchange, lung_part(t))
_whole_physiology(heat_exchange, lung_part) = heat_exchange
_whole_physiology(heat_exchange::NamedTuple, lung_part) = unwrap_physiology(heat_exchange[lung_part])

# Forward all HeatExchange accessor methods to the whole-organism physiology
HeatExchange.shape_pars(t::OrganismTraits) = HeatExchange.shape_pars(_whole_physiology(t))
HeatExchange.insulation_pars(t::OrganismTraits) = HeatExchange.insulation_pars(_whole_physiology(t))
function HeatExchange.conduction_pars_external(t::OrganismTraits)
    HeatExchange.conduction_pars_external(_whole_physiology(t))
end
function HeatExchange.conduction_pars_internal(t::OrganismTraits)
    HeatExchange.conduction_pars_internal(_whole_physiology(t))
end
HeatExchange.convection_pars(t::OrganismTraits) = HeatExchange.convection_pars(_whole_physiology(t))
HeatExchange.radiation_pars(t::OrganismTraits) = HeatExchange.radiation_pars(_whole_physiology(t))
HeatExchange.evaporation_pars(t::OrganismTraits) = HeatExchange.evaporation_pars(_whole_physiology(t))
HeatExchange.hydraulic_pars(t::OrganismTraits) = HeatExchange.hydraulic_pars(_whole_physiology(t))
HeatExchange.respiration_pars(t::OrganismTraits) = HeatExchange.respiration_pars(_whole_physiology(t))
HeatExchange.metabolism_pars(t::OrganismTraits) = HeatExchange.metabolism_pars(_whole_physiology(t))
HeatExchange.options(t::OrganismTraits) = HeatExchange.options(_whole_physiology(t))

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
