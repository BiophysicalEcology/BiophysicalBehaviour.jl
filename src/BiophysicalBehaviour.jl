module BiophysicalBehaviour

using Unitful, UnitfulMoles, ModelParameters
using BiophysicalGeometry
using BiophysicalGeometry: AbstractBody, shape
using ConstructionBase
using HeatExchange

export piloerect!, uncurl!, vasodilate!, hyperthermia!, pant!, sweat!
export endotherm_thermoregulation_original

include("endotherm/thermoregulation.jl")
include("endotherm/homeothermy.jl")

end # module BiophysicalBehaviour
