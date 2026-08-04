module BiophysicalBehaviour

import ConstructionBase
import FluidProperties

using BiophysicalGeometry
using HeatExchange
using ModelParameters
using Unitful
using UnitfulMoles

using BiophysicalGeometry: AbstractBody, shape

using Enzyme
import Ipopt
using LinearAlgebra: dot

using ConstructionBase: getproperties, setproperties
using HeatExchange: heat_balance, zbrent, zbrac, find_zero, Bisection, A42, AlefeldPotraShi, FalsePosition, Order0, SmoothingStrategy
using Setfield: @set

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

# Part selectors and traversal
export PartSelector,
    WholeBody,
    ByName,
    Compartment,
    part_names,
    select_names,
    map_parts,
    foldl_parts,
    set_part

# Effector tags + interface
export Effector,
    Piloerect,
    Uncurl,
    Vasodilate,
    Hyperthermia,
    Pant,
    Sweat,
    effect

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

include("organism.jl")
include("parts.jl")
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

end # module BiophysicalBehaviour
