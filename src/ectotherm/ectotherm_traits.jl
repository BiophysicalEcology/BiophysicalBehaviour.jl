"""
    EctothermBehavioralLimits{C,Sh,D,H,T,Tp,Ab,Pa} <: AbstractBehaviourParameters

Parameters controlling ectotherm thermoregulation through position and physiological
state changes.

An ectotherm cannot generate substantial heat internally; instead it thermoregulates
behaviourally by choosing its microhabitat (shade level, height above ground, depth
below ground) and by adjusting its physical state (absorptivity, posture, substrate
contact, panting). This struct encodes the limits on those adjustments and the
thermal thresholds that trigger them.

Ported from NicheMapR ectotherm.f / ectotherm.R (DAYACT/NOCTURN/CREPUS/BURROW/CLIMB
flags, TMINPR/TMAXPR/TBASK/CTMIN/CTMAX/TEMERGE, alpha_max/alpha_min, postur, pct_cond,
panting parameters).

# Fields – position
- `control::C`: Control strategy (default `RuleBasedSequentialControl()`)
- `shade::Sh`: `SteppedParameter` for shade fraction (0–1). `reference`/`current`
  are the starting shade; `max` is maximum available shade; `step` is the increment
  per iteration (SHADEADJUST.f `DSHD`).
- `depth::D`: `SteppedParameter` for soil-node index. `reference` = 1 (surface node);
  `max` is the deepest accessible underground node; `step` = 1.
- `depth_min_underground::Int`: Shallowest accessible underground node (SELDEP.f `MINNODE`,
  NicheMapR `mindepth`). Active animals reset to `depth.reference` (surface = 1);
  retreating underground begins at `depth_min_underground`.
- `height::H`: `SteppedParameter` for atmospheric profile height-node index.
  `reference` = 1 (lowest node); `max` is the highest accessible node; `step` = 1.

# Fields – absorptivity
- `absorptivity::Ab`: `SteppedParameter` for dorsal solar absorptivity (0–1).
  `reference` = alpha_min (lightest, default); `max` = alpha_max (darkest).
  Steps toward `max` when cold (`darken`), back toward `reference` when hot (`lighten`).

# Fields – panting
- `pant_rate::Pa`: `SteppedParameter` for respiratory panting multiplier (≥ 1.0).
  `reference` = 1.0 (baseline breathing); `max` = maximum multiplier. Steps up when
  too hot via `pant`.

# Fields – thermal thresholds
Hierarchy: `T_critical_min` ≤ `T_emerge` ≤ `T_bask` ≤ `T_active_min` ≤ `T_active_max` ≤ `T_critical_max`

- `T_target::Tp`: `SteppedParameter` for the operative target body temperature
  (NicheMapR TPREF). `reference`/`current` = starting target temperature;
  `max` = `T_active_max`; `step` = increment per iteration. When Tb > `T_target.current`
  but < `T_active_max`, the target temperature steps up (hyperthermia tolerance) before
  cooling behaviours are triggered. Reset each hour to `reference`.
- `T_active_min::T`: Lower activity temperature, TMINPR (T_F_min). Below this the animal
  is basking rather than active.
- `T_active_max::T`: Upper activity temperature, TMAXPR (T_F_max). When `T_target`
  has stepped up to this value and Tb is still above it, cooling behaviours begin.
- `T_bask::T`: Minimum basking temperature, TBASK (T_B_min). Thermoregulation loop warms the
  animal to at least `T_bask`; the animal is basking (not active) while `T_bask ≤ Tb < T_active_min`.
- `T_critical_min::T`: Critical thermal minimum, CTMIN (lethal lower limit).
- `T_critical_max::T`: Critical thermal maximum, CTMAX (lethal upper limit).
- `T_emerge::T`: Minimum soil temperature at the previous-hour underground depth to
  allow emergence, TEMERGE (T_RB_min).

# Fields – capability flags
- `can_retreat_underground::Bool`: Whether the animal can retreat underground (BURROW flag).
- `can_climb::Bool`: Whether the animal can climb to elevated height nodes (CLIMB flag).
- `can_seek_shade::Bool`: Whether the animal seeks shade.
- `can_change_absorptivity::Bool`: Whether the animal can darken/lighten (NicheMapR
  `alpha_max ≠ alpha_min`).
- `can_solar_orient::Bool`: Whether the animal can dynamically orient toward/away from the sun
  (NicheMapR `postur`). When `true`, the thermoregulation loop starts each hour at `Intermediate()`
  (neutral/active posture) and can shift to `NormalToSun` (too cold/basking) or `ParallelToSun`
  (too hot). When `false`, `solar_orientation` in `radiation_pars` is used as a fixed orientation
  each hour: `NormalToSun`, `ParallelToSun`, `Intermediate`, or `ZenithAngleVarying` (silhouette
  area computed from sun elevation angle via HeatExchange using `e_vars.zenith_angle`).
- `can_press_to_ground::Bool`: Whether the animal presses against the substrate (tracked
  in output). Conduction fraction is set via organism physiology (`conduction_pars_external`).
- `can_pant::Bool`: Whether the animal can pant for evaporative cooling.

# Fields – posture state
- `sun_orientation::Float64`: Angle (degrees) between body axis and sun direction.
  45.0 = `Intermediate` (neutral/active start); 90.0 = `NormalToSun` (maximum silhouette,
  too-cold/basking response); 0.0 = `ParallelToSun` (minimum silhouette, too-hot response).
  Reset to 45.0 each hour by `reset_position`.

# Fields – conduction state
- `pressed_to_ground::Bool`: `true` if currently pressed to substrate. Reset each hour.
  Conduction fraction is controlled via organism morphology (`conduction_pars_external`).
- `underground_tb_equals_soil::Bool`: When `true` (default), body temperature when underground
  is set equal to the soil temperature at the chosen node (NicheMapR BELOWGROUND.f behaviour).
  When `false`, the full heat balance is solved underground, which gives `Tb > T_soil`
  because metabolic heat has no way to dissipate when all surrounding temperatures equal
  `T_soil`.
- `underground_shaded::Bool`: When `true` (default, matches NicheMapR `shade_burrow=1`), soil
  temperatures in the underground retreat come from the maximum-shade microclimate run (retreat
  is in a shaded location, e.g. under vegetation). When `false`, soil temperatures come from the
  minimum-shade run (exposed, unshaded location). No blending — always one extreme or the other.
"""
Base.@kwdef struct EctothermBehavioralLimits{
    C<:AbstractControlStrategy,Sh,D,H,T,Tp,Ab,Pa
} <: AbstractBehaviourParameters
    control::C          = RuleBasedSequentialControl()
    shade::Sh
    depth::D
    depth_min_underground::Int = 2
    height::H
    absorptivity::Ab    = SteppedParameter(; current=0.9, reference=0.9, max=0.9, step=0.0)
    pant_rate::Pa       = SteppedParameter(; current=1.0, reference=1.0, max=1.0, step=0.1)
    T_target::Tp            # SteppedParameter: rises from reference toward T_active_max each iteration
    T_active_min::T
    T_active_max::T
    T_bask::T
    T_critical_min::T
    T_critical_max::T
    T_emerge::T
    # Capability flags
    can_retreat_underground::Bool   = true
    can_climb::Bool                 = false
    can_seek_shade::Bool            = true
    can_change_absorptivity::Bool   = false
    can_solar_orient::Bool          = true
    can_press_to_ground::Bool       = true
    can_pant::Bool                  = false
    # Posture state (reset each hour)
    sun_orientation::Float64        = 45.0  # degrees; 45 = Intermediate (neutral), 90 = NormalToSun, 0 = ParallelToSun
    # Conduction state (reset each hour)
    pressed_to_ground::Bool         = false
    # Body temperature model when underground
    # true  → Tb = T_soil directly (NicheMapR behaviour; no heat balance solved underground)
    # false → solve full heat balance underground (Tb > T_soil due to metabolic heat)
    underground_tb_equals_soil::Bool     = true
    # Shaded underground flag (NicheMapR shade_burrow)
    # true  → soil temperatures in the underground retreat come from the maximum-shade run
    #         (retreat located under vegetation / shade cover; NicheMapR default)
    # false → soil temperatures come from the minimum-shade run (exposed, unshaded location)
    underground_shaded::Bool             = true
end
