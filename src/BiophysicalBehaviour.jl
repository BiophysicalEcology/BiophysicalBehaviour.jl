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

# Arrest: generalized dormancy/diapause/quiescence controllers -- any life
# stage, any process. See src/arrest/ for the design notes.
export AbstractArrestController, NeverController
export AbstractComparison, Below, Above
export AbstractDirection, AnyDirection, Rising, Falling
export signal_value, signal_rate
export initial_controller_state, controller_rate, level
export trigger_conditions, register_callback
export signed_gap, direction_gate

export AbstractMetric, RawSignal, RawProgress, Accumulate
export metric_state, metric_rate, metric_value, metric_rate_value
export AbstractBound, FixedBound
export bound_state, bound_value
export ThresholdController
export FunctionController
export AnyController, AllController

export AbstractArrestModel, ComposedArrest, AnyArrestModel, AllArrestModel
export initial_arrest_state, advance_arrest, arrest_level, arrest_conditions, step_state
export arrest_component
export print_arrest_structure
export describe_metric, describe_bound, describe_comparison, describe_direction, node_label

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

include("arrest/controllers.jl")
include("arrest/metrics.jl")
include("arrest/bounds.jl")
include("arrest/threshold.jl")
include("arrest/function_controller.jl")
include("arrest/composition.jl")
include("arrest/arrest_model.jl")
include("arrest/component.jl")
include("arrest/structure.jl")

end # module BiophysicalBehaviour
