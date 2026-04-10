"""
    AvailableEnvironments{M}

Wraps two `MicroResult` objects (from `Microclimate.solve`) representing the same
location at minimum and maximum shade, allowing the ectotherm code to interpolate
microclimate conditions at any intermediate shade fraction.

The two runs must cover the same time period and use the same spatial resolution.
Atmospheric profiles, soil temperatures, radiation, and humidity are all blended
linearly between the two runs.

# Fields
- `min_shade_result::M`: `MicroResult` run at `min_shade_fraction` shade
- `max_shade_result::M`: `MicroResult` run at `max_shade_fraction` shade
- `min_shade_fraction::Float64`: Shade fraction used for `min_shade_result` (e.g. 0.0)
- `max_shade_fraction::Float64`: Shade fraction used for `max_shade_result` (e.g. 0.9)

# Example
```julia
low_shade_result = Microclimate.solve(problem_at_0pct_shade)
high_shade_result = Microclimate.solve(problem_at_90pct_shade)
available_environments = AvailableEnvironments(low_shade_result, high_shade_result, 0.0, 0.9)
```
"""
struct AvailableEnvironments{M, D, H}
    min_shade_result::M
    max_shade_result::M
    min_shade_fraction::Float64
    max_shade_fraction::Float64
    depths::D
    heights::H
end

# =============================================================================
# Body temperature solver
# =============================================================================

"""
    solve_body_temperature(organism, env_vars, env_pars) → T_core

Find the steady-state body temperature for an ectotherm using Brent's method
(`zbrent` from HeatExchange.jl, ported from NicheMapR `ZBRENT.f`) on the
`ectotherm` heat balance function.

Returns the air temperature as a fallback if root-finding fails.
"""
function solve_body_temperature(organism, env_vars, env_pars, T_bask=nothing, T_active_max=nothing)
    e = (environment_pars=env_pars, environment_vars=env_vars)
    # Fixed bracket covering physiologically plausible ectotherm body temps (0–70°C).
    lo = 273.15
    hi = 343.15
    # Eyes open only when basking or active (NicheMapR SEVAP.f line 119:
    # IF((TC.GE.TBASK).AND.(TC.LE.TMAXPR))). Pre-build closed-eye organism.
    org_closed = if !isnothing(T_bask)
        @set organism.traits.heat_exchange.evaporation_pars.eye_fraction = 0.0
    else
        organism
    end
    try
        T_sol = zbrent(
            T -> begin
                T_K = T * u"K"
                org = (!isnothing(T_bask) && (T_K < T_bask || T_K > T_active_max)) ?
                    org_closed : organism
                ustrip(u"W", ectotherm(T_K, org, e).heat_balance)
            end,
            lo, hi, 1e-3,
        )
        T_sol * u"K"
    catch
        env_vars.air_temperature
    end
end

# =============================================================================
# Top-level dispatch – routes on AvailableEnvironments type
# =============================================================================

"""
    thermoregulate(organism, available_environments, limits, environmental_params, step[, previous_depth])

Run the ectotherm behavioural thermoregulation algorithm for one time step.

Dispatches on the organism's `thermal_strategy` and `control_strategy`.

# Arguments
- `organism`: `Organism` with `Ectotherm` thermal strategy
- `available_environments::AvailableEnvironments`: Wraps min/max shade `MicroResult` objects
- `limits::EctothermBehavioralLimits`: Behavioural bounds (shade, depth, height, thresholds)
- `environmental_params`: `EnvironmentalPars` (substrate absorptivity, emissivity, etc.)
- `step::Int`: Current simulation time step (1-based)
- `previous_depth::Int`: Soil-node index occupied at the previous time step
  (1 = surface; >1 = underground). Controls emergence logic. Default 1.

# Returns
NamedTuple with fields:
- `T_core`: Final body temperature (K)
- `shade`: Shade fraction chosen
- `depth`: Soil-node index (depth.reference = surface, >reference = underground)
- `height`: Height-node index chosen
- `state`: `OrganismState` — `Resting()`, `Basking()`, or `Active()` (mirrors NicheMapR ACT 0/1/2)
- `ectotherm_out`: Full output from the final `ectotherm()` call
"""
function thermoregulate(
    organism::Organism,
    available_environments::AvailableEnvironments,
    limits::EctothermBehavioralLimits,
    environmental_params,
    step::Int,
    previous_depth::Int=1;
    activity_commenced::Bool=false,
)
    thermoregulate(
        thermal_strategy(organism),
        control_strategy(organism),
        organism,
        available_environments,
        limits,
        environmental_params,
        step,
        previous_depth;
        activity_commenced,
    )
end

# =============================================================================
# Core algorithm – Ectotherm + RuleBasedSequentialControl (ECTOTHERM.f)
# =============================================================================

"""
    thermoregulate(::Ectotherm, ::RuleBasedSequentialControl, organism, available_environments, limits, environmental_params, step, previous_depth)

Core ectotherm behavioural thermoregulation loop (ECTOTHERM.f / ectotherm.R logic).

## Behaviour sequence

1. **Reset position** to reference state (surface, no shade, lowest height, default
   absorptivity, perpendicular posture, not pressed to ground, no panting).
2. **Determine activity period** from solar zenith angle and radiation.
3. **Inactive period**: if the animal can retreat underground, call `select_depth` to find
   the optimal soil node; otherwise stay above ground.
4. **Emergence check**: if the animal was underground at the previous step and the soil
   at that depth is still below `T_emerge`, stay underground.
5. **Active period thermoregulation loop** — organism traits modified in-place via `@set`:

   *Too hot* (Tb > T_target.current):
   revert perpendicular→intermediate (same iteration as shade seeking, mirrors
   NicheMapR THERMO.f phase 2 no-RETURN) → `lighten` → `seek_shade` → `climb` →
   `pant` → `increment_T_target` → `retreat_underground`

   *Too cold* (Tb < T_bask):
   `darken` → `orient_perpendicular` → `press_to_ground` → `avoid_shade` →
   `retreat_underground` (if Tb < T_critical_min)

   *Basking* (T_bask ≤ Tb < T_active_min): `orient_perpendicular`

   Each behaviour modifies both `limits` (state flags / stepped parameters) and a
   local `organism_current` (organism traits). The loop exits when Tb is in the
   target range or all options are exhausted.
"""
function thermoregulate(
    ::Ectotherm,
    ::RuleBasedSequentialControl,
    organism::Organism,
    available_environments::AvailableEnvironments,
    limits::EctothermBehavioralLimits,
    environmental_params,
    step::Int,
    previous_depth::Int;
    activity_commenced::Bool=false,
)
    low_shade = available_environments.min_shade_result
    high_shade = available_environments.max_shade_result
    (; max_iterations) = limits.control

    # -------------------------------------------------------------------------
    # 1. Reset position / state for this hour; organism starts at reference state
    # -------------------------------------------------------------------------
    limits           = reset_position(limits)
    organism_current = organism   # fresh reference each hour

    # When solar orientation is enabled, reset to Intermediate (neutral/active) each hour so
    # the loop can shift to NormalToSun (too cold/basking) or ParallelToSun (too hot) as needed.
    # When disabled, solar_orientation in radiation_pars is used as-is each hour:
    # NormalToSun, ParallelToSun, Intermediate → fixed area via HeatExchange dispatch;
    # ZenithAngleVarying → area computed from e_vars.zenith_angle inside HeatExchange.
    if limits.can_solar_orient
        organism_current = @set organism_current.traits.heat_exchange.radiation_pars.solar_orientation =
            Intermediate()
    end

    # Clamp height/depth max to actual environment array sizes so that the
    # default typemax(Int) "use all" sentinel works correctly at runtime.
    limits = @set limits.height.max = min(limits.height.max, length(available_environments.heights))
    limits = @set limits.depth.max  = min(limits.depth.max,  length(available_environments.depths))

    # -------------------------------------------------------------------------
    # 2. Determine activity period
    # -------------------------------------------------------------------------
    zenith    = low_shade.solar_radiation.zenith_angle[step]
    solar_rad = low_shade.global_radiation[step]
    active    = is_active(activity_period(organism), zenith, solar_rad)

    # -------------------------------------------------------------------------
    # 3a. Inactive period → always retreat underground (mirrors NicheMapR burrow=1
    #     logic: a diurnal/nocturnal animal goes underground during its inactive
    #     period unconditionally when can_retreat_underground=true).
    # -------------------------------------------------------------------------
    if !active
        if limits.can_retreat_underground
            underground_bf = _underground_blend_factor(limits.burrow_shade_mode, limits, low_shade, high_shade, step)
            limits = select_depth(limits, low_shade, high_shade, step, limits.shade.current, underground_bf)
        end
        env = interpolate_environment(available_environments, step, limits, environmental_params)
        return _build_ectotherm_output(organism_current, env, environmental_params, limits, active, available_environments, underground_bf)
    end

    # -------------------------------------------------------------------------
    # 3b. Active but was underground last step → check T_emerge condition
    # -------------------------------------------------------------------------
    if previous_depth > limits.depth.reference
        # Underground temperature is binary (shaded or unshaded location),
        # not a linear blend — matches interpolate_environment BELOWGROUND.f logic.
        underground_bf  = _underground_blend_factor(limits.burrow_shade_mode, limits, low_shade, high_shade, step)
        T_soil_at_depth = _blend(
            low_shade.soil_temperature[step, previous_depth],
            high_shade.soil_temperature[step, previous_depth],
            underground_bf,
        )

        # WARMSIG: require a soil temperature change signal before the animal emerges.
        # Mirrors NicheMapR ectotherm.f lines 2218-2244. Applies only when:
        #   Δsoil_signal != 0, depth_node > 2, and no activity has occurred today yet.
        # Δsoil_signal > 0: soil must be warming at >= Δsoil_signal K/hr
        # Δsoil_signal < 0: soil must be cooling at >= |Δsoil_signal| K/hr
        stay_for_signal = false
        if !iszero(limits.Δsoil_signal) && previous_depth > 2 && !activity_commenced && step > 1
            T_soil_prev_step = _blend(
                low_shade.soil_temperature[step-1, previous_depth],
                high_shade.soil_temperature[step-1, previous_depth],
                underground_bf,
            )
            # Each time step is 1 hour so the K change per step equals the K/hr rate.
            soil_delta = (T_soil_at_depth - T_soil_prev_step) / 1u"hr"
            stay_for_signal = (limits.Δsoil_signal > zero(limits.Δsoil_signal) && soil_delta < limits.Δsoil_signal) ||
                              (limits.Δsoil_signal < zero(limits.Δsoil_signal) && soil_delta > limits.Δsoil_signal)
        end

        if T_soil_at_depth < limits.T_emerge || stay_for_signal
            # Stay underground — re-run select_depth to find the best available node for
            # this hour (mirrors NicheMapR calling SELDEP.f each inactive hour). This
            # allows the animal to move deeper when its current node cools below CT_min.
            limits = select_depth(limits, low_shade, high_shade, step, limits.shade.current, underground_bf)
            env    = interpolate_environment(available_environments, step, limits, environmental_params)
            return _build_ectotherm_output(organism_current, env, environmental_params, limits, false, available_environments, underground_bf)
        end
    end

    # -------------------------------------------------------------------------
    # 4. Active above ground: behavioural thermoregulation loop
    # -------------------------------------------------------------------------
    env = interpolate_environment(available_environments, step, limits, environmental_params)
    Tb  = solve_body_temperature(organism_current, env, environmental_params, limits.T_bask, limits.T_active_max)
    underground_bf = 0.0  # updated if the animal retreats underground during the loop

    iteration = 0
    while iteration < max_iterations
        iteration += 1

        # NicheMapR THERMO.f phase 2 (first sub-case): revert perpendicular → intermediate
        # when Tb has risen into the accepted active range [T_active_min, T_target].
        # In NicheMapR THERMO.f the revert fires BEFORE the acceptance check, so the
        # animal always relaxes its posture once it has warmed to T_active_min.
        # After reverting, recalculate Tb (now lower with intermediate area) and break;
        # no further actions are attempted (matches THERMO.f phase 4 no-action → return).
        if limits.can_solar_orient && limits.sun_orientation == 90.0 &&
               Tb >= limits.T_active_min && Tb <= limits.T_target.current
            limits, organism_current = orient_intermediate(organism_current, limits)
            env = interpolate_environment(available_environments, step, limits, environmental_params)
            Tb  = solve_body_temperature(organism_current, env, environmental_params, limits.T_bask, limits.T_active_max)
            break
        end

        # NicheMapR ECTOTHERM.f line 3333: accept if Tb ∈ [T_active_min, T_target].
        if Tb >= limits.T_active_min && Tb <= limits.T_target.current
            break
        end

        if Tb > limits.T_target.current
            # -- Too hot --
            # NicheMapR THERMO.f phase 2 (second sub-case, no RETURN): revert
            # perpendicular → intermediate in the SAME iteration as shade seeking.
            if limits.can_solar_orient && limits.sun_orientation == 90.0
                limits, organism_current = orient_intermediate(organism_current, limits)
                # fall through to shade seeking below
            end

            if limits.can_change_absorptivity &&
                   limits.absorptivity.current > limits.absorptivity.reference
                limits, organism_current = lighten(organism_current, limits)

            elseif limits.can_seek_shade && limits.shade.current < limits.shade.max
                limits = seek_shade(limits)

            elseif limits.T_target.current < limits.T_target.max
                limits = increment_T_target(limits)

            elseif limits.can_climb && limits.height.current < limits.height.max
                limits = climb(limits)

            elseif limits.can_pant && limits.pant_rate.current < limits.pant_rate.max
                limits, organism_current = pant(organism_current, limits)

            elseif limits.can_retreat_underground
                underground_bf = _underground_blend_factor(limits.burrow_shade_mode, limits, low_shade, high_shade, step)
                limits = select_depth(limits, low_shade, high_shade, step, limits.shade.current, underground_bf)
                env    = interpolate_environment(available_environments, step, limits, environmental_params)
                break
            else
                break
            end

        elseif Tb < limits.T_bask
            # -- Too cold: darken → perpendicular → press to ground → avoid shade → retreat_underground --
            if limits.can_change_absorptivity &&
               limits.absorptivity.current < limits.absorptivity.max
                limits, organism_current = darken(organism_current, limits)

            elseif limits.can_solar_orient && limits.sun_orientation < 90.0
                limits, organism_current = orient_perpendicular(organism_current, limits)

            elseif limits.can_press_to_ground && !limits.pressed_to_ground
                limits, organism_current = press_to_ground(organism_current, limits)

            elseif zenith < 90u"°" && limits.shade.current > limits.shade.reference
                limits = avoid_shade(limits)

            elseif limits.can_seek_shade && zenith >= 90u"°" && limits.shade.current < limits.shade.max
                # Night: seek shade to reduce longwave loss to cold sky
                limits = seek_shade(limits)

            elseif limits.can_climb && Tb < limits.T_critical_min && limits.height.current < limits.height.max
                limits = climb(limits)

            elseif limits.can_retreat_underground
                # All warming behaviours are exhausted — always retreat underground.
                # (depth_min_underground ≥ 2 means surface is never a valid burrow node.)
                underground_bf = _underground_blend_factor(limits.burrow_shade_mode, limits, low_shade, high_shade, step)
                limits = select_depth(limits, low_shade, high_shade, step, limits.shade.current, underground_bf)
                env    = interpolate_environment(available_environments, step, limits, environmental_params)
                break
            else
                break
            end

        elseif limits.can_solar_orient && limits.sun_orientation < 90.0 &&
               Tb < limits.T_active_min
            # -- Basking range [T_bask, T_active_min): orient NormalToSun to maximise solar gain --
            limits, organism_current = orient_perpendicular(organism_current, limits)

        else
            break
        end

        env = interpolate_environment(available_environments, step, limits, environmental_params)
        Tb  = solve_body_temperature(organism_current, env, environmental_params, limits.T_bask, limits.T_active_max)
    end

    return _build_ectotherm_output(organism_current, env, environmental_params, limits, active, available_environments, underground_bf)
end

# =============================================================================
# Output builder
# =============================================================================

function _build_ectotherm_output(organism, env_vars, env_pars, limits, in_active_period, available_environments, underground_bf::Float64=0.0)
    e        = (environment_pars=env_pars, environment_vars=env_vars)
    is_underground = limits.depth.current > limits.depth.reference
    # When underground and underground_tb_equals_soil=true (default, matches NicheMapR BELOWGROUND.f):
    # set Tb = T_soil directly. The full heat balance gives Tb > T_soil because Q_metab > 0
    # forces equilibrium above the uniform soil temperature when all surrounding temperatures
    # equal T_soil and heat loss pathways are near zero.
    Tb       = (is_underground && limits.underground_tb_equals_soil) ?
               env_vars.air_temperature :
               solve_body_temperature(organism, env_vars, env_pars, limits.T_bask, limits.T_active_max)
    org_out  = (limits.T_bask <= Tb <= limits.T_active_max) ? organism :
               @set organism.traits.heat_exchange.evaporation_pars.eye_fraction = 0.0
    ecto_out = ectotherm(Tb, org_out, e)
    height = if is_underground
        -uconvert(u"m", available_environments.depths[limits.depth.current])
    else
        available_environments.heights[limits.height.current]
    end
    # Classify organism state (mirrors NicheMapR ACT column: 0=Resting, 1=Basking, 2=Active)
    # NicheMapR ECTOTHERM.f lines 3645-3684 (ACTHR logic):
    #   ACT=0 (Resting):  Tc < T_B_min (TBASK=17.5°C) OR Tc > T_F_max (TMAXPR=34°C)
    #   ACT=1 (Basking):  T_B_min <= Tc < T_F_min (TMINPR=24°C)
    #   ACT=2 (Active):   T_F_min <= Tc <= T_F_max
    state = if !in_active_period || is_underground
        Resting()
    elseif limits.T_active_min <= Tb <= limits.T_active_max
        Active()
    elseif limits.T_bask <= Tb < limits.T_active_min
        Basking()
    else
        Resting()
    end
    return (;
        T_core                = Tb,
        shade                 = is_underground ? (underground_bf > 0.5 ? available_environments.max_shade_fraction : available_environments.min_shade_fraction) : limits.shade.current,
        depth_node            = limits.depth.current,
        height,
        absorptivity          = limits.absorptivity.current,
        sun_orientation       = limits.sun_orientation,
        pressed_to_ground     = limits.pressed_to_ground,
        pant_rate             = limits.pant_rate.current,
        state,
        ectotherm_out         = ecto_out,
    )
end
