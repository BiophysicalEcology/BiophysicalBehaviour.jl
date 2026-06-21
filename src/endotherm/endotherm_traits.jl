"""
    SteppedParameter{T,S}

A parameter that can be adjusted in steps toward a maximum value.

# Fields
- `current::T`: Current value
- `reference::T`: Reference/baseline value (defaults to current)
- `max::T`: Maximum allowed value
- `step::S`: Step size for adjustment
"""
Base.@kwdef struct SteppedParameter{T,S}
    current::T
    reference::T = current
    max::T
    step::S
end

"""
    InsulationLimits{D,S}

Limits for dorsal and ventral insulation depth adjustment (piloerection).

# Fields
- `dorsal::SteppedParameter{D,S}`: Dorsal insulation limits
- `ventral::SteppedParameter{D,S}`: Ventral insulation limits
"""
Base.@kwdef struct InsulationLimits{D,S}
    dorsal::SteppedParameter{D,S}
    ventral::SteppedParameter{D,S}
end

"""
    PantingLimits{P,S,C,M,T}

Limits for panting behavior with associated metabolic costs.

# Fields
- `pant::SteppedParameter{P,S}`: Panting rate limits
- `cost::C`: Current panting cost (W)
- `multiplier::M`: Metabolic cost multiplier at max panting
- `core_temperature_ref::T`: Reference core temperature for Q10 calculation
"""
Base.@kwdef struct PantingLimits{P,S,C,M,T}
    pant::SteppedParameter{P,S}
    cost::C = 0.0u"W"
    multiplier::M = 1.05
    core_temperature_ref::T
end

"""
    ThermoregulationLimits{C,Q,I,Sh,K,Tc,P,Sw} <: AbstractBehaviourParameters

Parameters controlling endotherm thermoregulation behavior.

Contains limits for all adjustable parameters: insulation depth, body shape,
tissue conductivity, core temperature, panting, and skin wetness.

# Fields
- `control::C`: Control strategy (RuleBasedSequentialControl, IPOPTControl, etc.)
- `minimum_heat_flow::Q`: Reference minimum metabolic rate
- `insulation::InsulationLimits`: Piloerection limits (dorsal/ventral)
- `axis_ratio_factor::SteppedParameter`: Body shape adjustment limits (uncurling)
- `flesh_conductivity::SteppedParameter`: Tissue conductivity limits (vasodilation)
- `core_temperature::SteppedParameter`: Core temperature limits (hyperthermia)
- `panting::PantingLimits`: Panting limits and costs
- `skin_wetness::SteppedParameter`: Sweating/cutaneous evaporation limits
- `core_temperature_penalty::Float64`: IPOPT objective penalty for core temperature deviation from setpoint. Default 1.0.
- `metabolic_heat_penalty::Float64`: Regularisation weight on metabolic heat generation. A small value
  (default 0.1) prevents high-panting/high-metabolic_heat_flow degeneracy in cold conditions. In hot
  conditions the Q10 inequality constraint overrides this and forces metabolic_heat_flow up with
  core_temperature, so this value does not impede thermogenesis.
- `panting_penalty::Float64`: IPOPT objective penalty for panting (normalised to [0,1] range).
  Relative to `skin_wetness_penalty` this controls which activates first. Default 1.0.
- `skin_wetness_penalty::Float64`: IPOPT objective penalty for skin wetness (normalised to [0,1] range).
  Set `skin_wetness_penalty > panting_penalty` for panting-first (rabbits, birds);
  `skin_wetness_penalty < panting_penalty` for sweating-first (humans);
  equal for parallel activation. Default 1.0.
- `gradient_penalty::Float64`: IPOPT objective penalty for deviation from `target_core_skin_gradient`.
  Zero (default) disables the term. Non-zero values bias the solution toward maintaining the
  specified core–skin temperature difference, which can activate vasodilation and evaporation
  before absolute core_temperature deviation becomes the primary signal.
- `target_core_skin_gradient::Float64`: Target core_temperature − skin_temperature difference (K).
  Only used when `gradient_penalty > 0`. Typical resting value is ~3 K. Default 2.0.
"""
Base.@kwdef struct ThermoregulationLimits{C<:AbstractControlStrategy,Q,I,Sh,K,Tc,P,Sw} <: AbstractBehaviourParameters
    control::C = RuleBasedSequentialControl()
    minimum_heat_flow::Q
    insulation::I
    axis_ratio_factor::Sh
    flesh_conductivity::K
    core_temperature::Tc
    panting::P
    skin_wetness::Sw
    core_temperature_penalty::Float64  = 1.0
    metabolic_heat_penalty::Float64    = 0.1
    panting_penalty::Float64           = 1.0
    skin_wetness_penalty::Float64      = 1.0
    gradient_penalty::Float64          = 0.0
    target_core_skin_gradient::Float64 = 2.0
end
