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
The `*_weight` fields are the weights of the IPOPT objective's `RegulationTarget` terms
(one weighted squared normalised deviation each); their *relative* magnitudes set which
effectors the solver reaches for first.
- `core_temperature_weight::Float64`: weight on the core-temperature-toward-setpoint term. Default 1.0.
- `metabolic_heat_weight::Float64`: weight on the metabolic-heat-toward-minimum term. A small value
  (default 0.1) prevents high-panting/high-metabolic_heat_flow degeneracy in cold conditions. In hot
  conditions the Q10 inequality constraint overrides this and forces metabolic_heat_flow up with
  core_temperature, so this value does not impede thermogenesis.
- `panting_weight::Float64`: weight on the panting term (normalised to [0,1] range).
  Relative to `skin_wetness_weight` this controls which activates first. Default 1.0.
- `skin_wetness_weight::Float64`: weight on the skin-wetness term (normalised to [0,1] range).
  Set `skin_wetness_weight > panting_weight` for panting-first (rabbits, birds);
  `skin_wetness_weight < panting_weight` for sweating-first (humans);
  equal for parallel activation. Default 1.0.
- `flesh_conductivity_weight::Float64`: weight on the flesh-conductivity-toward-reference term
  (vasodilation), normalised to the `[reference, max]` range. Default 0.0.
- `gradient_weight::Float64`: weight on the core–skin gradient term (deviation from `target_core_skin_gradient`).
  Zero (default) disables the term. Non-zero values bias the solution toward maintaining the
  specified core–skin temperature difference, which can activate vasodilation and evaporation
  before absolute core_temperature deviation becomes the primary signal.
- `target_core_skin_gradient`: Target core_temperature − skin_temperature difference,
  a temperature `Quantity`. Only used when `gradient_weight > 0`. Typical resting
  value is ~3 K. Default `2.0u"K"`.
"""
Base.@kwdef struct ThermoregulationLimits{C<:AbstractControlStrategy,Q,I,Sh,K,Tc,P,Sw,Tu} <: AbstractBehaviourParameters
    control::C = RuleBasedSequentialControl()
    minimum_heat_flow::Q
    insulation::I
    axis_ratio_factor::Sh
    flesh_conductivity::K
    core_temperature::Tc
    panting::P
    skin_wetness::Sw
    core_temperature_weight::Float64  = 1.0
    metabolic_heat_weight::Float64    = 0.1
    panting_weight::Float64           = 1.0
    skin_wetness_weight::Float64      = 1.0
    flesh_conductivity_weight::Float64 = 0.0
    gradient_weight::Float64          = 0.0
    # A core→skin temperature difference: carries K, not a bare Float64 (so no
    # consumer has to attach a unit to it).
    target_core_skin_gradient::Tu      = 2.0u"K"
    # Named replacements for the inline NLP magic numbers (§3.7). Consumed by the
    # multi-part NLP template/objective builders (Phase 7.5).
    skin_temperature_undershoot::Tu        = 5.0u"K"
    skin_temperature_core_overshoot::Tu    = 5.0u"K"
    metabolic_heat_flow_max_multiplier::Float64 = 20.0
    minimum_normalisation_range::Float64        = 1e-6
end

