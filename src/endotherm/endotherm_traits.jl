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
- `T_core_ref::T`: Reference core temperature for Q10 calculation
"""
Base.@kwdef struct PantingLimits{P,S,C,M,T}
    pant::SteppedParameter{P,S}
    cost::C = 0.0u"W"
    multiplier::M = 1.05
    T_core_ref::T
end

"""
    ThermoregulationLimits{C,Q,I,Sh,K,Tc,P,Sw} <: AbstractBehaviourParameters

Parameters controlling endotherm thermoregulation behavior.

Contains limits for all adjustable parameters: insulation depth, body shape,
tissue conductivity, core temperature, panting, and skin wetness.

# Fields
- `control::C`: Control strategy (RuleBasedSequentialControl, PDEControl, etc.)
- `Q_minimum_ref::Q`: Reference minimum metabolic rate
- `insulation::InsulationLimits`: Piloerection limits (dorsal/ventral)
- `shape_b::SteppedParameter`: Body shape adjustment limits
- `k_flesh::SteppedParameter`: Tissue conductivity limits (vasodilation)
- `T_core::SteppedParameter`: Core temperature limits (hyperthermia)
- `panting::PantingLimits`: Panting limits and costs
- `skin_wetness::SteppedParameter`: Sweating limits
"""
Base.@kwdef struct ThermoregulationLimits{C<:AbstractControlStrategy,Q,I,Sh,K,Tc,P,Sw} <: AbstractBehaviourParameters
    control::C = RuleBasedSequentialControl()
    Q_minimum_ref::Q
    insulation::I
    shape_b::Sh
    k_flesh::K
    T_core::Tc
    panting::P
    skin_wetness::Sw
end
