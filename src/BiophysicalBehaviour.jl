module BiophysicalBehaviour

import ConstructionBase
import FluidProperties

using BiophysicalGeometry
using HeatExchange
using ModelParameters
using Unitful
using UnitfulMoles

using BiophysicalGeometry: AbstractBody, shape

using ConstructionBase: getproperties, setproperties
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
    ResponsiveActivity

# Thermal strategies
export AbstractThermalStrategy,
    Endotherm,
    Ectotherm,
    Heterotherm

# Control strategies
export AbstractControlStrategy,
    RuleBasedSequentialControl,
    PDEControl

# Thermoregulation modes
export AbstractThermoregulationMode,
    CoreOnly,
    CoreAndPanting,
    CorePantingSweating

# Traits structs
export BehavioralTraits,
    OrganismTraits

# Trait accessors
export thermal_strategy,
    behavior,
    physiology,
    thermoregulation,
    activity,
    control_strategy

# Thermoregulation functions
export piloerect, uncurl, vasodilate, hyperthermia, pant, sweat

export thermoregulate

# Thermoregulation limit structs
export SteppedParameter,
    InsulationLimits,
    PantingLimits,
    ThermoregulationLimits

# Example constructors
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

include("organism.jl")
include("endotherm/endotherm_traits.jl")
include("endotherm/thermoregulation.jl")
include("endotherm/homeothermy.jl")
include("endotherm/example_variables_and_parameters.jl")

end # module BiophysicalBehaviour
