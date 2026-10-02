# Ectotherm Thermoregulation: Logic and Control Flow

**Source files:**
- `src/ectotherm/ectothermy.jl` — outer loop, body-temperature solver, output builder
- `src/ectotherm/thermoregulation.jl` — all individual behaviour functions, depth/environment helpers
- `src/ectotherm/ectotherm_traits.jl` — `EctothermBehavioralLimits`, `BurrowShadeMode`

**NicheMapR correspondence:** `ECTOTHERM.f`, `THERMO.f`, `SELDEP.f`, `ABOVEGROUND.f`, `BELOWGROUND.f`

---

## Overview

`thermoregulate(organism, available_environments, limits, env_pars, step[, previous_depth])`

This is the top-level function called once per hourly time step. It:

1. Dispatches to `thermoregulate(::Ectotherm, ::RuleBasedSequentialControl, ...)` via the organism's
   `thermal_strategy` and `control_strategy`.
2. Runs a sequential priority-ordered loop of behavioural adjustments.
3. Returns a `NamedTuple` describing the final state of the organism at that hour.

The key architectural principle is that **all state is immutable**. Behaviour functions return
new `limits` and `organism_current` values (using `@set` from Setfield.jl). Nothing is mutated
in place; the caller accumulates outputs into arrays.

---

## Inputs

| Argument | Type | Description |
|---|---|---|
| `organism` | `Organism` | Full organism specification (body shape, physiology, activity period, etc.) |
| `available_environments` | `AvailableEnvironments` | Two `MicroResult` objects at min/max shade fractions |
| `limits` | `EctothermBehavioralLimits` | Behavioural bounds and thermal thresholds (see below) |
| `env_pars` | `EnvironmentalPars` | Substrate optical properties etc. (constant across steps) |
| `step` | `Int` | 1-based index into the microclimate time series |
| `previous_depth` | `Int` | Soil-node depth index from the previous time step (1 = surface) |
| `activity_commenced` | `Bool` | Whether the organism has been active (ACT=2) at any earlier step today (WARMSIG logic) |

### `AvailableEnvironments`

Wraps two `Microclimate.MicroResult` objects (one at `min_shade_fraction`, one at `max_shade_fraction`).
All microclimate quantities (temperature, wind, humidity, radiation) are blended linearly between
them using:

```
blend_factor = (shade.current - min_shade_fraction) / (max_shade_fraction - min_shade_fraction)    (clamped 0–1)
```

The `depths` and `heights` vectors record the physical depths and heights corresponding to
soil-node and height-node indices.

### `EctothermBehavioralLimits`

Tracks both **fixed bounds** (what the animal is capable of) and **mutable state** (what it has done so
far this hour). Key fields:

| Field | Meaning | NicheMapR name |
|---|---|---|
| `shade` | `SteppedParameter` — current/min/max/step shade fraction | SHADE / MAXSHD / DSHD |
| `depth` | `SteppedParameter` — current/reference(=1) depth node | BURROW depth node |
| `depth_min_underground` | Shallowest underground node allowed | MINNODE (mindepth) |
| `height` | `SteppedParameter` — current/reference height node | CLIMB height node |
| `absorptivity` | `SteppedParameter` — current/min(ref)/max alpha | alpha_min / alpha_max |
| `pant_rate` | `SteppedParameter` — current/max pant multiplier | pantmax |
| `target_temperature` | `SteppedParameter` — current/max preferred temperature | TPREF / T_F_max |
| `active_temperature_min/max` | Foraging temperature range | TMINPR / TMAXPR (T_F_min / T_F_max) |
| `basking_temperature_min` | Minimum basking temperature | TBASK (T_B_min) |
| `critical_temperature_min/max` | Lethal thermal limits | CTMIN / CTMAX |
| `emerge_temperature_min` | Minimum soil temperature to emerge from underground retreat | TEMERGE (T_RB_min) |
| `can_*` | Boolean capability flags | BURROW / CLIMB / SHADE / postur / panting flags |
| `sun_orientation` | Current posture angle (°): 45=Intermediate, 90=NormalToSun, 0=ParallelToSun | postur |
| `pressed_to_ground` | Current ground-contact state | pct_cond |
| `burrow_shade_mode` | `MinShadeOnly` / `AdaptiveBurrowShade` / `MaxShadeOnly` | shdburrow 0/1/2 |
| `emerge_signal` | Required soil temp change rate before emergence (K/hr) | warmsig |

---

## Phase 1 — Reset to reference state

```julia
limits           = reset_position(limits)
organism_current = organism
```

At the start of every hour, `reset_position` sets all stepped parameters back to their reference values:

- `shade.current` → `shade.reference` (starting shade, typically 0)
- `depth.current` → `depth.reference` (= 1, surface)
- `height.current` → `height.reference` (= 1, lowest node)
- `absorptivity.current` → `absorptivity.max` (darkest = maximum solar gain, ready to lighten if hot)
- `pant_rate.current` → `pant_rate.reference` (= 1.0, no panting)
- `target_temperature.current` → `target_temperature.reference` (starting preferred temperature)
- `sun_orientation` → 45.0 (Intermediate, neutral foraging posture)
- `pressed_to_ground` → false

**Solar orientation reset (when `can_solar_orient = true`):**
The organism's `radiation_pars.solar_orientation` is reset to `Intermediate()` and the silhouette
area to the average of normal/parallel. This allows the loop to shift to `NormalToSun` (too cold)
or `ParallelToSun` (too hot) from a neutral starting point.

When `can_solar_orient = false`, `solar_orientation` is used as set in `radiation_pars` and never
changed: `NormalToSun`, `ParallelToSun`, `Intermediate` (fixed area), or `ZenithAngleVarying`
(silhouette computed dynamically from `e_vars.zenith_angle` by HeatExchange).

**Array-size clamping:**
`height.max` and `depth.max` default to `typemax(Int)` (meaning "use all nodes"). They are
clamped to the actual array sizes from `available_environments` each hour.

---

## Phase 2 — Activity period check

```julia
zenith    = min_shade.solar_radiation.zenith_angle[step]
sunlight = min_shade.global_radiation[step]
active    = is_active(activity_period(organism), zenith, sunlight)
```

`is_active` dispatches on the organism's activity period type:

| Type | Active when |
|---|---|
| `Diurnal` | zenith < 90° **and** sunlight > 0 |
| `Nocturnal` | zenith ≥ 90° **or** sunlight ≤ 0 |
| `Crepuscular` | 85° ≤ zenith ≤ 95° |
| `CombinedActivity` | active in any constituent period |
| `ResponsiveActivity` | user-supplied function |

If `active = false`, the algorithm proceeds to Phase 3a.
If `active = true`, it jumps to Phase 3b (emergence check) or Phase 4 (above-ground thermoregulation).

---

## Phase 3a — Inactive period: burrow as thermal refuge

*Runs only when `active = false`.*

If `can_retreat_underground = true`, the algorithm checks whether any accessible soil node
is warmer than the above-ground air temperature at the reference height (blended between
min/max shade environments). If so, `select_depth` is called to find the best node.

The underground blend factor (which soil temperature profile to use) is determined by
`burrow_shade_mode`:

| `BurrowShadeMode` | Blend factor | NicheMapR shdburrow |
|---|---|---|
| `MinShadeOnly()` | 0.0 (always min-shade soil) | 0 |
| `MaxShadeOnly()` | 1.0 (always max-shade soil) | 2 |
| `AdaptiveBurrowShade()` | 1.0 if min-shade soil at shallowest node is outside [critical_temperature_min, active_temperature_max], else 0.0 | 1 |

After position is determined, `interpolate_environment` builds the `EnvironmentalVars` and
`_build_ectotherm_output` is called and **returned immediately** (no further thermoregulation).

---

## Phase 3b — Emergence check

*Runs only when `active = true` and `previous_depth > depth.reference` (was underground).*

The soil temperature at the previous-step depth is computed. The animal stays underground if either:

**Condition A — Too cold to emerge:**
```
soil_temperature_at_depth < emerge_temperature_min
```

**Condition B — No emergence signal from soil temperature** (mirrors ECTOTHERM.f lines 2218–2244):
```
emerge_signal != 0  AND  previous_depth > 2  AND  !activity_commenced  AND  step > 1
```
`soil_delta = (soil_temperature[step] - soil_temperature[step-1]) / 1hr`  (K/hr; one step = one hour)
- If `emerge_signal > 0`: requires `soil_delta ≥ emerge_signal` (diurnal basker waits for morning warm-up)
- If `emerge_signal < 0`: requires `soil_delta ≤ emerge_signal` (nocturnal animal waits for evening cool-down)

If either condition holds: `select_depth` is called again (so the animal can move to a better node),
`interpolate_environment` is called, and `_build_ectotherm_output` is returned with `active = false`.

If neither condition holds, the animal emerges and proceeds to Phase 4.

---

## Phase 4 — Above-ground behavioural thermoregulation loop

Initial environment and body temperature are computed at the current (reset) position:

```julia
env              = interpolate_environment(available_environments, step, limits, environmental_params)
core_temperature = solve_body_temperature(organism_current, env, environmental_params)
```

Then the iteration loop runs up to `max_iterations` times. Each iteration:

1. Evaluates which thermal branch to enter using direct unit-aware comparisons (all temperatures carry `K` units throughout).
2. Applies exactly **one** behaviour per iteration.
3. Recalculates `env` and `core_temperature` at the end of the iteration (except when retreating underground).

### Body temperature solver

`solve_body_temperature` calls `zbrent` (Brent's method, ported from `ZBRENT.f`) to find the root
of `ectotherm(T, organism, e).heat_balance` on the interval `[273.15 K, 343.15 K]` (0–70°C). If the
solver fails (e.g., no root in bracket), it returns `air_temperature` as a fallback.

### Acceptance (loop exit) condition

```
active_temperature_min ≤ core_temperature ≤ target_temperature.current  →  break
```

**Special case — revert perpendicular→intermediate (NicheMapR THERMO.f phase 2, first sub-case):**
If `can_solar_orient = true` and the organism is fully perpendicular (`sun_orientation = 90.0`)
and `core_temperature` is now in `[active_temperature_min, target_temperature]`, orient back to `Intermediate()` first,
recompute `core_temperature`, then break. This ensures the animal does not remain in a basking
posture once it has warmed to its target.

---

### Too-hot branch (`core_temperature > target_temperature.current`)

Behaviours are tried in strict priority order. Only the **first applicable** one fires per iteration.

**Step 0 — Revert perpendicular→intermediate (THERMO.f phase 2, no-RETURN sub-case):**
If fully perpendicular and now too hot, first revert to `Intermediate()`. Then fall through
to the shade-seeking priority below in the **same** iteration (no `break`, no `continue`).

| Priority | Condition | Behaviour | Effect |
|---|---|---|---|
| 1 | `can_change_absorptivity` and `absorptivity.current > absorptivity.reference` | `lighten` | Decrease dorsal absorptivity by one step toward `absorptivity_min` |
| 2 | `can_seek_shade` and `shade.current < shade.max` | `seek_shade` | Increase shade by `shade.step` |
| 3 | `target_temperature.current < target_temperature.max` | `increment_target_temperature` | Raise tolerance threshold by one step toward `active_temperature_max` |
| 4 | `can_climb` and `height.current < height.max` | `climb` | Move up one height node (cooler, windier air) |
| 5 | `can_pant` and `pant_rate.current < pant_rate.max` | `pant` | Increase pant multiplier by one step |
| 6 | `can_retreat_underground` | `select_depth` then **break** | Burrow underground; exit loop |
| fallback | — | **break** | No options left; accept current Tb |

**Notes:**
- Both `lighten` and `darken` update both `body_absorptivity_ventral` and `body_absorptivity_ventral`.
- `increment_target_temperature` mirrors NicheMapR's TPREF incrementing: the animal tolerates getting
  hotter before triggering shade-seeking. Once `target_temperature.current = active_temperature_max`, shade-seeking
  begins. This means shade-seeking in the model starts only when the animal is above `active_temperature_max`,
  not above the initial `target_temperature`.
- `climb` moves to higher, typically cooler and windier, air layers. Shade is **not** reset when climbing.
- When retreating underground (priority 6), `_underground_blend_factor` is computed, `select_depth`
  finds the best node, `interpolate_environment` is called, then the loop **breaks** immediately
  (no further `Tb` recalculation at the end of the iteration).

---

### Too-cold branch (`core_temperature < basking_temperature_min`)

| Priority | Condition | Behaviour | Effect |
|---|---|---|---|
| 1 | `can_change_absorptivity` and `absorptivity.current < absorptivity.max` | `darken` | Increase dorsal absorptivity by one step toward `absorptivity_max` |
| 2 | `can_solar_orient` and `sun_orientation < 90.0` | `orient_perpendicular` | Orient NormalToSun; maximise silhouette area |
| 3 | `can_press_to_ground` and `!pressed_to_ground` | `press_to_ground` | Record ground contact (conduction fraction from organism physiology) |
| 4 | zenith < 90° (daytime) and `shade.current > shade.reference` | `avoid_shade` | Decrease shade by one step toward minimum |
| 5 | `can_seek_shade` and zenith ≥ 90° (night) and `shade.current < shade.max` | `seek_shade` | At night, seek shade to reduce longwave cooling to cold sky |
| 6 | `can_climb` and `core_temperature < critical_temperature_min` and `height.current < height.max` | `climb` | Emergency climb above critically cold layer |
| 7 | `can_retreat_underground` | conditional `select_depth` then **break** | Retreat if `core_temperature < critical_temperature_min` (emergency) or underground is warmer than air |
| fallback | — | **break** | No options left |

**Notes:**
- `press_to_ground` only sets the state flag; the actual conduction fraction is pre-set in
  `organism.traits.heat_exchange.conduction_pars_external.conduction_fraction`. If set to 0.0,
  pressing to ground has no thermal effect (though it is still recorded in the output).
- Night shade-seeking (priority 5) is the reverse of daytime shade-seeking: the sky is colder than
  vegetation at night, so moving under cover reduces longwave radiative loss.
- Underground retreat (priority 7) is **conditional**: the animal only retreats if `Tb < critical_temperature_min`
  **or** if any accessible soil node is warmer than above-ground air. In contrast to the too-hot case,
  the loop does not `break` immediately after the condition check — it breaks regardless of whether
  it actually moved (both paths converge to `break`).

---

### Basking range (`basking_temperature_min ≤ core_temperature < active_temperature_min`)

```julia
if limits.can_solar_orient && limits.sun_orientation < 90.0 && core_temperature < active_temperature_min
    orient_perpendicular(...)
else
    break
end
```

If the animal can solar orient and is not yet fully perpendicular, it rotates to `NormalToSun()`.
Otherwise it accepts the current state and breaks.

---

### End-of-iteration recalculation

After any behaviour that does **not** break out of the loop:

```julia
env              = interpolate_environment(available_environments, step, limits, environmental_params)
core_temperature = solve_body_temperature(organism_current, env, environmental_params)
```

This recalculates both the microenvironment (blending min/max shade at the new shade/height) and
the steady-state body temperature. The next iteration then re-evaluates which branch to enter.

---

## Phase 5 — Output construction (`_build_ectotherm_output`)

After the loop exits:

1. Determines `is_underground = depth.current > depth.reference`.
2. **Body temperature when underground:** By default (`solve_underground = false`),
   sets `core_temperature = soil_temperature` directly at the chosen node (NicheMapR `BELOWGROUND.f`
   behaviour). The full heat balance would give `core_temperature > soil_temperature` because
   metabolic heat has nowhere to go when all surrounding temperatures equal `soil_temperature`,
   so this override is necessary for physical realism.
3. Runs `ectotherm(core_temperature, organism, e)` once to get the full heat balance output at the final state.
4. Converts `depth_node` and `height_node` to physical depths/heights using `available_environments.depths`
   and `available_environments.heights`.
5. **Activity state classification** (mirrors NicheMapR ACT column):
   - `Resting()` (ACT=0): not in active period, or underground, or `core_temperature < basking_temperature_min` or `core_temperature > active_temperature_max`
   - `Basking()` (ACT=1): `basking_temperature_min ≤ core_temperature < active_temperature_min`
   - `Active()` (ACT=2): `active_temperature_min ≤ core_temperature ≤ active_temperature_max`

---

## Output `NamedTuple`

| Field | Type | Description |
|---|---|---|
| `core_temperature` | `K` | Final body temperature |
| `shade` | `Float64` | Shade fraction chosen (above-ground) or burrow shade fraction (underground) |
| `depth_node` | `Int` | Soil-node index (1 = surface) |
| `height` | length | Physical height (positive above ground, negative depth in m for underground) |
| `absorptivity` | `Float64` | Final dorsal solar absorptivity |
| `sun_orientation` | `Float64` | Final posture: 45=Intermediate, 90=NormalToSun, 0=ParallelToSun |
| `pressed_to_ground` | `Bool` | Whether organism is in ground contact |
| `pant_rate` | `Float64` | Final pant multiplier (1.0 = baseline) |
| `state` | `OrganismState` | `Resting()`, `Basking()`, or `Active()` |
| `ectotherm_out` | `NamedTuple` | Full heat balance output from `HeatExchange.ectotherm()` |

---

## Helper functions

### `interpolate_environment`

Builds an `EnvironmentalVars` for the organism's current position by blending the two microclimate
runs.

**Above-ground:** Blends air temperature, wind speed, sky/ground temperatures, radiation, and
humidity at the current height node. Solar radiation is passed pre-shade (HeatExchange applies
`(1 - shade)` internally). Relative humidity conserves actual vapour pressure at min-shade
reference conditions then recomputes RH at blended temperature (ABOVEGROUND.f WETAIR approach).

**Below-ground:** Sets `air_temperature = sky_temperature = substrate_temperature = ground_temperature = soil_temperature` at the chosen node,
wind = 0.01 m/s, `global_radiation = 0`. Underground blend factor is determined by
`burrow_shade_mode` (binary: 0.0 or 1.0, not a linear blend).

### `select_depth` (SELDEP.f)

Iterates soil nodes from `depth_min_underground` to `depth.max`. Returns the **shallowest** node
where:

```
critical_temperature_min < soil_temperature < critical_temperature_max - (critical_temperature_max - active_temperature_max)/2
```

If no node satisfies the condition, falls back to the deepest available node.

### `solve_body_temperature`

Calls `zbrent` (Brent's method, tolerance 1e-3 K) on `ectotherm(T, organism, e).heat_balance` bracketed
between 273.15 K and 343.15 K. Returns `air_temperature` as fallback if root-finding fails.

### Behaviour functions (all in `thermoregulation.jl`)

Each function accepts `(organism, limits)` (or just `limits` for positional-only changes) and
returns updated `(limits, organism)` or just `limits`. All updates use `@set` (Setfield.jl).

| Function | Changes |
|---|---|
| `reset_position(limits)` | Resets all stepped parameters to reference; called once per hour |
| `lighten(organism, limits)` | Decrease `body_absorptivity_dorsal` and `body_absorptivity_ventral` by `absorptivity.step` |
| `darken(organism, limits)` | Increase `body_absorptivity_dorsal` and `body_absorptivity_ventral` by `absorptivity.step` |
| `seek_shade(limits)` | Increase `shade.current` by `shade.step` |
| `avoid_shade(limits)` | Decrease `shade.current` by `shade.step` |
| `increment_target_temperature(limits)` | Increase `target_temperature.current` by `target_temperature.step` |
| `climb(limits)` | Increase `height.current` by `height.step` |
| `descend(limits)` | Increase `depth.current` by `depth.step` (not called in main loop) |
| `orient_parallel(organism, limits)` | Set `ParallelToSun()`, `sun_orientation=0.0`, update `A_silhouette` |
| `orient_perpendicular(organism, limits)` | Set `NormalToSun()`, `sun_orientation=90.0`, update `A_silhouette` |
| `orient_intermediate(organism, limits)` | Set `Intermediate()`, `sun_orientation=45.0`, `A_silhouette` = mean |
| `press_to_ground(organism, limits)` | Set `pressed_to_ground=true` |
| `pant(organism, limits)` | Increase `pant_rate.current` and `respiration_pars.pant` by one step |
| `select_depth(limits, ...)` | Find and set `depth.current` to best soil node |

---

## Control flow summary (flowchart)

```
Each hour:
  │
  ├─ reset_position() + reset solar_orientation if can_solar_orient
  │
  ├─ is_active(activity_period, zenith, sunlight)?
  │    │
  │    NO ──→ can_retreat_underground?
  │              YES → underground warmer than air? → select_depth()
  │              → interpolate_environment() → _build_output() → RETURN
  │    │
  │    YES
  │    │
  │    ├─ previous_depth > 1? (was underground)
  │    │    YES → soil_temperature < emerge_temperature_min OR emerge_signal not met?
  │    │              YES → select_depth() → interpolate() → _build_output(active=false) → RETURN
  │    │              NO  → continue (emerge)
  │    │
  │    └─ Above-ground thermoregulation loop (max_iterations):
  │         │
  │         ├─ [revert perpendicular → intermediate if core_temperature now in acceptance window] → break
  │         ├─ active_temperature_min ≤ core_temperature ≤ target_temperature? → break
  │         │
  │         ├─ core_temperature > target_temperature? (too hot)
  │         │    [revert perpendicular if applicable, fall through]
  │         │    lighten → seek_shade → increment_target_temperature → climb → pant → select_depth+break → break
  │         │
  │         ├─ core_temperature < basking_temperature_min? (too cold)
  │         │    darken → orient_perpendicular → press_to_ground → avoid_shade →
  │         │    night-seek_shade → climb(CT_min) → select_depth+break → break
  │         │
  │         └─ basking_temperature_min ≤ core_temperature < active_temperature_min? (basking)
  │              orient_perpendicular → break
  │
  └─ _build_ectotherm_output() → RETURN
       (core_temperature = soil_temperature if underground; classify Resting/Basking/Active)
```
