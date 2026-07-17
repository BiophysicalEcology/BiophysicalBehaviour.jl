"""
    BurrowShadeMode

Abstract type for the three NicheMapR `shdburrow` options controlling which soil temperatures
are used for underground retreats.

- `MinShadeOnly()`: always use minimum-shade soil temperatures (NicheMapR `shdburrow=0`)
- `AdaptiveBurrowShade()`: use max-shade soil if min-shade soil at the shallowest accessible
  node is outside `[escape_temperature_min, active_temperature_max]`; otherwise use min-shade (NicheMapR `shdburrow=1`)
- `MaxShadeOnly()`: always use maximum-shade soil temperatures (NicheMapR `shdburrow=2`)
"""
abstract type BurrowShadeMode end
"""Always use minimum-shade soil temperatures for the underground retreat (NicheMapR `shdburrow=0`)."""
struct MinShadeOnly <: BurrowShadeMode end
"""Adaptive: use max-shade soil temperatures when min-shade soil at the shallowest accessible
node is outside `[escape_temperature_min, active_temperature_max]` (NicheMapR `shdburrow=1`)."""
struct AdaptiveBurrowShade <: BurrowShadeMode end
"""Always use maximum-shade soil temperatures for the underground retreat (NicheMapR `shdburrow=2`)."""
struct MaxShadeOnly <: BurrowShadeMode end

"""
    EctothermBehavioralLimits{C,Sh,D,H,T,Tp,Ab,Pa,Ws,Bs} <: AbstractBehaviourParameters

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
- `depth::D`: `SteppedParameter` for soil-node index. `reference` is the animal's normal/foraging
  depth (node 1 = surface, for a typical ground-dwelling ectotherm; a deeper node for a fossorial
  one); `max` is the retreat depth (SELDEP.f `MINNODE`/NicheMapR `mindepth` is now just
  `reference + step`, the search's starting point); `step` = 1.
- `height::H`: `SteppedParameter` for atmospheric profile height-node index. `reference` is the
  animal's normal/foraging height (node 1 = ground, for a typical species; an elevated node for
  an arboreal one); `max` is the retreat height; `step` = 1.

# Fields – absorptivity
- `absorptivity::Ab`: `SteppedParameter` for dorsal solar absorptivity (0–1).
  `reference` = alpha_min (lightest, default); `max` = alpha_max (darkest).
  Steps toward `max` when cold (`darken`), back toward `reference` when hot (`lighten`).

# Fields – panting
- `pant_rate::Pa`: `SteppedParameter` for respiratory panting multiplier (≥ 1.0).
  `reference` = 1.0 (baseline breathing); `max` = maximum multiplier. Steps up when
  too hot via `pant`.

# Fields – thermal thresholds
Hierarchy: `escape_temperature_min` ≤ `emerge_temperature_min` ≤ `basking_temperature_min` ≤ `active_temperature_min` ≤ `active_temperature_max` ≤ `escape_temperature_max`

- `target_temperature::Tp`: `SteppedParameter` for the operative target body temperature
  (NicheMapR TPREF). `reference`/`current` = starting target temperature;
  `max` = `active_temperature_max`; `step` = increment per iteration. When Tb > `target_temperature.current`
  but < `active_temperature_max`, the target temperature steps up (hyperthermia tolerance) before
  cooling behaviours are triggered. Reset each hour to `reference`.
- `active_temperature_min::T`: Lower activity temperature, TMINPR (T_F_min). Below this the animal
  is basking rather than active.
- `active_temperature_max::T`: Upper activity temperature, TMAXPR (T_F_max). When `target_temperature`
  has stepped up to this value and Tb is still above it, cooling behaviours begin.
- `basking_temperature_min::T`: Minimum basking temperature, TBASK (T_B_min). Thermoregulation loop warms
  the animal to at least `basking_temperature_min`; the animal is basking (not active) while
  `basking_temperature_min ≤ Tb < active_temperature_min`.
- `escape_temperature_min::T`: Core temperature below which cold-escape behaviour (climbing
  above a cold air layer, retreating underground) begins.
- `escape_temperature_max::T`: Core temperature above which heat-escape behaviour (underground
  retreat) begins.
- `emerge_temperature_min::T`: Minimum soil temperature at the previous-hour underground depth to
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
- `solve_underground::Bool`: When `true`, the full heat balance is solved underground. When
  `false` (default), body temperature is set equal to soil temperature at the chosen node
  (NicheMapR BELOWGROUND.f behaviour), because metabolic heat has no way to dissipate when
  all surrounding temperatures equal `soil_temperature`.
- `burrow_shade_mode`: A `BurrowShadeMode` value controlling which soil temperatures are used
  for underground retreats. `MinShadeOnly()` = min-shade always; `MaxShadeOnly()` = max-shade
  always; `AdaptiveBurrowShade()` = max-shade when min-shade soil at the shallowest accessible
  node is outside `[escape_temperature_min, active_temperature_max]` (NicheMapR `shdburrow` 0/1/2).
- `emerge_signal`: Minimum soil-temperature change rate (units `K/hr`) required before the
  animal emerges from a retreat deeper than node 2, when `soil_temperature >= T_emerge_min` and
  no activity has occurred yet today (NicheMapR `warmsig`, default `0.0u"K/hr"` = disabled).
  `> 0`: must detect soil warming of at least this rate (diurnal basker waits for morning warm-up).
  `< 0`: must detect soil cooling of at least this rate (nocturnal animal waits for evening cool-down).
  `= 0`: emerge whenever `soil_temperature >= emerge_temperature_min` regardless of trend.
"""
Base.@kwdef struct EctothermBehavioralLimits{
    C<:AbstractControlStrategy,Sh,D,H,T,Tp,Ab,Pa,Ws,Bs<:BurrowShadeMode
} <: AbstractBehaviourParameters
    control::C          = RuleBasedSequentialControl()
    shade::Sh
    depth::D
    height::H
    absorptivity::Ab    = SteppedParameter(; current=0.9, reference=0.9, max=0.9, step=0.0)
    pant_rate::Pa       = SteppedParameter(; current=1.0, reference=1.0, max=1.0, step=0.1)
    target_temperature::Tp          # SteppedParameter: rises from reference toward active_temperature_max each iteration
    active_temperature_min::T
    active_temperature_max::T
    basking_temperature_min::T
    escape_temperature_min::T
    escape_temperature_max::T
    emerge_temperature_min::T
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
    # false → Tb = T_soil directly (NicheMapR behaviour; no heat balance solved underground)
    # true  → solve full heat balance underground (Tb > T_soil due to metabolic heat)
    solve_underground::Bool = false
    # Shaded underground flag (NicheMapR shade_burrow)
    # true  → soil temperatures in the underground retreat come from the maximum-shade run
    #         (retreat located under vegetation / shade cover; NicheMapR default)
    # false → soil temperatures come from the minimum-shade run (exposed, unshaded location)
    burrow_shade_mode::Bs                = MaxShadeOnly()
    emerge_signal::Ws                  = 0.0u"K/hr"
end

# Used by both thermoregulate and simulate_diurnal_behavior.
function _validate_ectotherm_thresholds(limits::EctothermBehavioralLimits)
    limits.escape_temperature_min < limits.emerge_temperature_min ||
        throw(ArgumentError("escape_temperature_min must be below emerge_temperature_min"))
    limits.emerge_temperature_min < limits.basking_temperature_min ||
        throw(ArgumentError("emerge_temperature_min must be below basking_temperature_min"))
    limits.basking_temperature_min < limits.active_temperature_min ||
        throw(ArgumentError("basking_temperature_min must be below active_temperature_min"))
    limits.active_temperature_min < limits.active_temperature_max ||
        throw(ArgumentError("active_temperature_min must be below active_temperature_max"))
    limits.active_temperature_max < limits.escape_temperature_max ||
        throw(ArgumentError("active_temperature_max must be below escape_temperature_max"))
    return nothing
end
