module BiophysicalBehaviour

using Unitful, UnitfulMoles, ModelParameters
using BiophysicalGeometry
using BiophysicalGeometry: AbstractBody, shape
using ConstructionBase
using HeatExchange

export AbstractBehaviouralParameters
export piloerect!, uncurl!, vasodilate!, hyperthermia!, pant!, sweat!
export endotherm_thermoregulation_original
export example_environment_vars, example_environment_pars,
            example_shape_pars, example_insulation_pars, example_conduction_pars_external,
            example_conduction_pars_internal, example_radiation_pars, example_evaporation_pars,
            example_hydraulic_pars, example_respiration_pars, example_insulation_pars, 
            example_metabolism_pars, example_model_pars, example_endotherm_thermoregulation_pars

include("organism.jl")
include("endotherm/thermoregulation.jl")
include("endotherm/homeothermy.jl")
include("endotherm/example_variables_and_parameters.jl")
include("endotherm/endotherm_traits.jl")

end # module BiophysicalBehaviour
