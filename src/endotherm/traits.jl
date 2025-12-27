"""
    EndoThermoregParameters <: AbstractBehavParameters

Behavioural and physiological parameters governing thermoregulatory
responses such as panting, sweating, core-temperature adjustment, and
changes in body shape or tissue conductivity.

# Fields

- `insulation_step` — Incremental reduction in insulation depth from the
  piloerect (maximum fur depth) state.
- `shape_b_step` — Increment by which the body‐shape parameter `shape_b`
  increases per iteration, allowing the animal to uncurl toward
  `shape_b_max`.
- `shape_b_max` — Maximum allowed value of shape_b when adjusting posture.
- `T_core_max` — Maximum core temperature (K).
- `T_core_min` — Minimum core temperature during torpor (K).
- `T_core_step` — Increment by which core temperature is elevated per
  iteration (K).
- `k_flesh_step` — Increment in flesh thermal conductivity (W/m/K).
- `k_flesh_max` — Maximum flesh thermal conductivity (W/m/K).
- `pant_step` — Increment in panting multiplier per iteration.
- `pant_multiplier` — Multiplier applied to basal metabolic rate at
  maximum panting effort.
- `pant_max` — Maximum panting multiplier.
- `skin_wetness_step` — Increment in surface wetness fraction used to
  model sweating behaviour.
- `skin_wetness_max` — Maximum fraction of body surface area that can be
  wetted (0–1).

All parameters use `Param` wrappers where appropriate for unit support
and bounds checking.
"""
Base.@kwdef struct EndoThermoregParameters{
    INS, SBS, SBM, TMX, TMN, TCS, KFS, KFM,
    PTS, PTM, PMX, SWS, SWM
} <: AbstractBehavParameters
    insulation_step::INS   = Param(1.0)
    shape_b_step::SBS      = Param(0.1)
    shape_b_max::SBM       = Param(5.0)
    T_core_max::TMX        = Param(39u"°C" |> u"K")
    T_core_min::TMN        = Param(19u"°C" |> u"K")
    T_core_step::TCS       = Param(0.1u"K")
    k_flesh_step::KFS      = Param(0.1u"W/m/K")
    k_flesh_max::KFM       = Param(2.8u"W/m/K")
    pant_step::PTS         = Param(0.1)
    pant_multiplier::PTM   = Param(1.05)
    pant_max::PMX          = Param(5.0)
    skin_wetness_step::SWS = Param(0.001)
    skin_wetness_max::SWM  = Param(1.0, bounds=(0.0, 1.0))
end