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
using Optimization
using OptimizationIpopt
import SciMLBase

using ConstructionBase: getproperties, setproperties
using HeatExchange: heat_balance, zbrent, zbrac, find_zero, Bisection, A42, AlefeldPotraShi, FalsePosition, Order0
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
    IPOPTControl

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

# Example constructors – endotherm
export example_environment_vars,
    example_environment_pars,
    example_ellipsoid_shape_pars,
    example_shape_pars,
    example_insulation_pars,
    example_conduction_pars_external,
    example_conduction_pars_internal,
    example_radiation_pars,
    example_evaporation_pars,
    example_hydraulic_pars,
    example_respiration_pars,
    example_metabolism_pars,
    example_metabolic_rate_options,
    example_thermoregulation_limits,
    example_behavioral_traits,
    example_organism_traits,
    example_heat_exchange_traits

# Example constructors – ectotherm
export example_ectotherm_behavioral_limits,
    example_ectotherm_behavioral_traits,
    example_ectotherm_organism_traits,
    example_ectotherm_conduction_pars_external,
    example_ectotherm_conduction_pars_internal,
    example_ectotherm_radiation_pars,
    example_ectotherm_evaporation_pars,
    example_ectotherm_respiration_pars,
    example_ectotherm_hydraulic_pars,
    example_ectotherm_metabolism_pars,
    example_ectotherm_heat_exchange_traits

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

end # module BiophysicalBehaviour
