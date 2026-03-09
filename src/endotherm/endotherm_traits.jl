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
    TorpidLimits{Tc,Ti,Tt}

Parameters controlling heterotherm torpor entry and depth.

Mirrors NicheMapR SOLVENDO.f TORPOR=1 parameters (TC_MIN, TC_INC, TORPTOL).

# Fields
- `T_core_min::Tc`: Minimum allowed core temperature during torpor (TC_MIN, default 19°C).
  The torpor loop will not lower T_core below this floor.
- `T_core_step::Ti`: Core temperature decrement per iteration (TC_INC, default 0.5 K).
- `torptol::Tt`: Tolerance on the Q_gen/Q_basal ratio for the "environment too warm"
  stopping condition (TORPTOL, default 0.05).
"""
Base.@kwdef struct TorpidLimits{Tc,Ti,Tt}
    T_core_min::Tc  = u"K"(19.0u"°C")
    T_core_step::Ti = 0.5u"K"
    torptol::Tt     = 0.05
end

"""
    ThermoregulationLimits{C,Q,I,Sh,K,Tc,P,Sw,To} <: AbstractBehaviourParameters

Parameters controlling endotherm thermoregulation behavior.

Contains limits for all adjustable parameters: insulation depth, body shape,
tissue conductivity, core temperature, panting, skin wetness, and torpor.

# Fields
- `control::C`: Control strategy (RuleBasedSequentialControl, PDEControl, etc.)
- `Q_minimum_ref::Q`: Reference minimum metabolic rate
- `insulation::InsulationLimits`: Piloerection limits (dorsal/ventral)
- `shape_b::SteppedParameter`: Body shape adjustment limits
- `k_flesh::SteppedParameter`: Tissue conductivity limits (vasodilation)
- `T_core::SteppedParameter`: Core temperature limits (hyperthermia)
- `panting::PantingLimits`: Panting limits and costs
- `skin_wetness::SteppedParameter`: Sweating limits
- `torpid::TorpidLimits`: Torpor parameters (heterotherms only; ignored for pure endotherms)
"""
Base.@kwdef struct ThermoregulationLimits{C<:AbstractControlStrategy,Q,I,Sh,K,Tc,P,Sw,To} <: AbstractBehaviourParameters
    control::C = RuleBasedSequentialControl()
    Q_minimum_ref::Q
    insulation::I
    shape_b::Sh
    k_flesh::K
    T_core::Tc
    panting::P
    skin_wetness::Sw
    torpid::To = TorpidLimits()
end
