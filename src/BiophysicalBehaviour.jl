module BiophysicalBehaviour

import ConstructionBase
import FluidProperties

using BiophysicalGeometry
using DataInterpolations
using HeatExchange
using ModelParameters
import OrdinaryDiffEqTsit5
using Unitful
using UnitfulMoles

using BiophysicalGeometry: AbstractBody, shape

using Enzyme
import Ipopt
using LinearAlgebra: dot

using ConstructionBase: getproperties, setproperties
using HeatExchange: heat_balance, zbrent, zbrac, find_zero, Bisection, A42, AlefeldPotraShi, FalsePosition, Order0,
    SmoothingStrategy, SurfaceSolveStrategy, HardBound, SmoothBound, safe_step, safe_clamp, safe_abs
using Setfield: @set
using ComponentArrays

# Organism and traits
export AbstractBehavior,
    AbstractMovementBehavior,
    AbstractTemperatureRegulation,
    NullBehavior,
    BurrowTemperatureRegulation

# Activity periods
export ActivityPeriod,
    Diurnal,
    Nocturnal,
    Crepuscular,
    CombinedActivity,
    ResponsiveActivity

# Organism states
export OrganismState,
    Resting,
    Basking,
    Active

# Thermal strategies
export AbstractThermalStrategy,
    Endotherm,
    Ectotherm,
    Heterotherm

# Control strategies
export AbstractControlStrategy,
    RuleBasedSequentialControl,
    IPOPTControl,
    IPOPTSolverCache

# NLP strategy types (re-exported from HeatExchange for IPOPTControl.nlp_strategy)
export NLPStrategy, WeightedMeanNLP, MultiSidedNLP

# Thermoregulation modes
export AbstractThermoregulationMode,
    CoreFirst,
    CoreAndPantingFirst,
    CorePantingSweatingFirst

# Traits structs
export BehavioralTraits,
    OrganismTraits

# Trait accessors
export thermal_strategy,
    behavior,
    heat_exchange,
    thermoregulation,
    activity_period,
    control_strategy

# Endotherm thermoregulation functions
export piloerect, uncurl, vasodilate, hyperthermia, pant, sweat

export thermoregulate

# Thermoregulation limit structs (shared)
export SteppedParameter,
    InsulationLimits,
    PantingLimits,
    ThermoregulationLimits

# Ectotherm types
export EctothermBehavioralLimits,
    AvailableEnvironments,
    BurrowShadeMode, MinShadeOnly, AdaptiveBurrowShade, MaxShadeOnly

# Ectotherm behaviour functions
export is_active,
    seek_shade,
    avoid_shade,
    climb,
    descend,
    select_depth,
    reset_position,
    interpolate_environment,
    solve_body_temperature,
    darken,
    lighten,
    orient_perpendicular,
    orient_parallel,
    press_to_ground,
    increment_target_temperature

# Example constructors – endotherm (heat exchange examples now in HeatExchange.jl)
export example_thermoregulation_limits,
    example_behavioral_traits,
    example_organism_traits

# Example constructors – ectotherm
export example_ectotherm_behavioral_limits,
    example_ectotherm_behavioral_traits,
    example_ectotherm_organism_traits

# Transient (lumped-capacitance) body-temperature simulation
export EnvironmentForcing, simulate_onelump, simulate_twolump
export TransientBehavioralPhase, SleepPhase, BaskPhase, ForagePhase, CoolPhase, ClimbPhase, BurrowPhase, RefugePhase
export phase_forcing, phase_posture, phase_state, phase_key, phase_condition, next_phase
export simulate_transient_behavior
export simulate_endotherm_activity_cycle

# Arrest: generalized dormancy/diapause/quiescence controllers. See src/arrest/.
# Extension-point types/functions (Abstract*, metric_*/bound_*, signed_gap,
# direction_gate, describe_*/node_label, trigger_conditions/arrest_conditions,
# signal_value/signal_rate) are not exported; use BiophysicalBehaviour.foo.
export NeverController
export BelowBound, AboveBound
export AnyDirection, RisingDirection, FallingDirection
export RawSignal, RawProgress, Accumulate
export FixedBound
export ThresholdController
export ProportionalController
export FunctionController
export AnyController, AllController
export ComposedArrest, AnyArrestModel, AllArrestModel
export initial_controller_state, controller_rate, controller_level, register_callback
export initial_arrest_state, advance_arrest, arrest_level, step_state
export arrest_component
export print_arrest_structure

# Stages: domain-agnostic ordered life-cycle topology. See src/stages/ and
# docs/stages.md. No domain vocabulary (insect, plant, or otherwise) is
# exported or shipped here -- a downstream package supplies concrete
# AbstractStage/AbstractTransition subtypes.
export AbstractStage, AbstractTransition, StageSequence
export stages, transitions
export stage_traits, transition_controller, stage_key, transition_key
export initial_stage_state, advance_stage
export stage_weights, current_stage, stage_value
export stage_conditions
export stage_component
export print_stage_structure

include("organism.jl")
include("endotherm/endotherm_traits.jl")
include("endotherm/thermoregulation/shared.jl")
include("endotherm/thermoregulation/rulebased.jl")
include("endotherm/thermoregulation/ipopt.jl")
include("endotherm/example_variables_and_parameters.jl")
include("ectotherm/ectotherm_traits.jl")
include("ectotherm/thermoregulation.jl")
include("ectotherm/ectothermy.jl")
include("ectotherm/example_variables_and_parameters.jl")
include("endotherm/thermoregulation/behavioural.jl")
include("transient/forcing.jl")
include("transient/simulate.jl")
include("transient/ectotherm/behavioral_driver.jl")
include("transient/endotherm/behavioral_driver.jl")

include("control/conditions.jl")
include("control/metrics.jl")
include("control/bounds.jl")
include("control/threshold.jl")
include("control/proportional.jl")
include("control/function_controller.jl")
include("control/composition.jl")
include("control/step_state.jl")
include("control/structure.jl")

include("arrest/arrest_model.jl")
include("arrest/component.jl")
include("arrest/structure.jl")

include("stages/stage_sequence.jl")
include("stages/lookup.jl")
include("stages/weights.jl")
include("stages/state.jl")
include("stages/conditions.jl")
include("stages/component.jl")
include("stages/structure.jl")

end # module BiophysicalBehaviour
