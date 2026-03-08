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
min_result = Microclimate.solve(problem_at_0pct_shade)
max_result = Microclimate.solve(problem_at_90pct_shade)
available_environments = AvailableEnvironments(min_result, max_result, 0.0, 0.9)
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

Find the steady-state body temperature for an ectotherm using bisection root-finding
on the `ectotherm` heat balance function from HeatExchange.jl.

Returns the air temperature as a fallback if root-finding fails.
"""
function solve_body_temperature(organism, env_vars, env_pars)
    e = (environment_pars=env_pars, environment_vars=env_vars)
    T_air = env_vars.T_air
    lo = T_air - 40u"K"
    hi = T_air + 100u"K"
    try
        find_zero(
            T_x -> ectotherm(T_x, organism, e).Q_bal,
            (lo, hi),
            Bisection(),
        )
    catch
        T_air
    end
end

# =============================================================================
# Top-level dispatch – routes on AvailableEnvironments type
# =============================================================================

"""
    thermoregulate(organism, available_environments, limits, environmental_params, step[, prev_step_depth])

Run the ectotherm behavioural thermoregulation algorithm for one time step.

Dispatches on the organism's `thermal_strategy` and `control_strategy`.

# Arguments
- `organism`: `Organism` with `Ectotherm` thermal strategy
- `available_environments::AvailableEnvironments`: Wraps min/max shade `MicroResult` objects
- `limits::EctothermBehavioralLimits`: Behavioural bounds (shade, depth, height, thresholds)
- `environmental_params`: `EnvironmentalPars` (substrate absorptivity, emissivity, etc.)
- `step::Int`: Current simulation time step (1-based)
- `prev_step_depth::Int`: Soil-node index occupied at the previous time step
  (1 = surface; >1 = burrowed). Controls emergence logic. Default 1.

# Returns
NamedTuple with fields:
- `T_core`: Final body temperature (K)
- `shade`: Shade fraction chosen
- `depth`: Soil-node index (depth.reference = surface, >reference = burrowed)
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
    prev_step_depth::Int=1,
)
    thermoregulate(
        thermal_strategy(organism),
        control_strategy(organism),
        organism,
        available_environments,
        limits,
        environmental_params,
        step,
        prev_step_depth,
    )
end

# =============================================================================
# Core algorithm – Ectotherm + RuleBasedSequentialControl (ECTOTHERM.f)
# =============================================================================

"""
    thermoregulate(::Ectotherm, ::RuleBasedSequentialControl, organism, available_environments, limits, environmental_params, step, prev_step_depth)

Core ectotherm behavioural thermoregulation loop (ECTOTHERM.f / ectotherm.R logic).

## Behaviour sequence

1. **Reset position** to reference state (surface, no shade, lowest height, default
   absorptivity, perpendicular posture, not pressed to ground, no panting).
2. **Determine activity period** from solar zenith angle and radiation.
3. **Inactive period**: if the animal can burrow, call `select_depth` to find the
   optimal soil node; otherwise stay above ground.
4. **Emergence check**: if the animal was burrowed at the previous step and the soil
   at that depth is still below `T_emerge`, stay burrowed.
5. **Active period thermoregulation loop** — organism traits modified in-place via `@set`:

   *Too hot* (Tb > T_preferred_max):
   `lighten` → `orient_parallel` → `seek_shade` → `climb` → `pant` → `burrow`

   *Too cold* (Tb < T_preferred_min):
   `darken` → `orient_perpendicular` → `press_to_ground` → `avoid_shade` →
   `burrow` (if Tb < T_critical_min)

   Each behaviour modifies both `limits` (state flags / stepped parameters) and a
   local `organism_current` (organism traits). The loop exits when Tb is in the
   preferred range or all options are exhausted.
"""
function thermoregulate(
    ::Ectotherm,
    ::RuleBasedSequentialControl,
    organism::Organism,
    available_environments::AvailableEnvironments,
    limits::EctothermBehavioralLimits,
    environmental_params,
    step::Int,
    prev_step_depth::Int,
)
    min_r = available_environments.min_shade_result
    max_r = available_environments.max_shade_result
    (; tolerance, max_iterations) = limits.control

    # -------------------------------------------------------------------------
    # 1. Reset position / state for this hour; organism starts at reference state
    # -------------------------------------------------------------------------
    limits           = reset_position(limits)
    organism_current = organism   # fresh reference each hour

    # When sun orientation is disabled, fix silhouette at Intermediate() —
    # the geometric mean of NormalToSun and ParallelToSun areas.
    if !limits.can_solar_orient
        organism_current = @set organism_current.traits.heat_exchange.radiation_pars.solar_orientation = Intermediate()
        organism_current = @set organism_current.traits.heat_exchange.radiation_pars.A_silhouette =
            _silhouette_area(organism_current.body, Intermediate())
    end

    # Clamp height/depth max to actual environment array sizes so that the
    # default typemax(Int) "use all" sentinel works correctly at runtime.
    limits = @set limits.height.max = min(limits.height.max, length(available_environments.heights))
    limits = @set limits.depth.max  = min(limits.depth.max,  length(available_environments.depths))

    # -------------------------------------------------------------------------
    # 2. Determine activity period
    # -------------------------------------------------------------------------
    zenith    = min_r.solar_radiation.zenith_angle[step]
    solar_rad = min_r.global_radiation[step]
    active    = is_active(activity_period(organism), zenith, solar_rad)

    # Helper: blend factor ∈ [0,1] from chosen shade
    shade_range = available_environments.max_shade_fraction - available_environments.min_shade_fraction
    @inline blend_factor_for(shade_chosen) = shade_range > 0 ?
        clamp((shade_chosen - available_environments.min_shade_fraction) / shade_range, 0.0, 1.0) :
        0.0

    # -------------------------------------------------------------------------
    # 3a. Inactive period → retreat to optimal burrow depth (or stay above)
    # -------------------------------------------------------------------------
    if !active
        if limits.can_burrow
            burrow_bf = limits.burrow_shaded ? 1.0 : 0.0
            limits = select_depth(limits, min_r, max_r, step, limits.shade.current, burrow_bf)
        end
        env = interpolate_environment(available_environments, step, limits, environmental_params)
        return _build_ectotherm_output(organism_current, env, environmental_params, limits, active, available_environments)
    end

    # -------------------------------------------------------------------------
    # 3b. Active but was burrowed last step → check T_emerge condition
    # -------------------------------------------------------------------------
    if prev_step_depth > limits.depth.reference
        bf          = blend_factor_for(limits.shade.reference)
        T_soil_prev = _blend(
            min_r.soil_temperature[step, prev_step_depth],
            max_r.soil_temperature[step, prev_step_depth],
            bf,
        )
        if T_soil_prev < limits.T_emerge
            limits = @set limits.depth.current = prev_step_depth
            env    = interpolate_environment(available_environments, step, limits, environmental_params)
            return _build_ectotherm_output(organism_current, env, environmental_params, limits, false, available_environments)
        end
    end

    # -------------------------------------------------------------------------
    # 4. Active above ground: behavioural thermoregulation loop
    # -------------------------------------------------------------------------
    env = interpolate_environment(available_environments, step, limits, environmental_params)
    Tb  = solve_body_temperature(organism_current, env, environmental_params)

    iteration = 0
    while iteration < max_iterations
        iteration += 1
        Tb_strip  = ustrip(u"K", Tb)
        Tpref_max = ustrip(u"K", limits.T_preferred_max)
        Tbask     = ustrip(u"K", limits.T_bask)
        Tcrit_min = ustrip(u"K", limits.T_critical_min)

        if Tb_strip > Tpref_max * (1 - tolerance)
            # -- Too hot: lighten → parallel → seek shade → climb → pant → burrow --
            if limits.can_change_absorptivity &&
               limits.absorptivity.current > limits.absorptivity.reference
                limits, organism_current = lighten(organism_current, limits)

            elseif limits.can_solar_orient && limits.sun_orientation > 0.0
                limits, organism_current = orient_parallel(organism_current, limits)

            elseif limits.can_seek_shade && limits.shade.current < limits.shade.max
                limits = seek_shade(limits)

            elseif limits.can_climb && limits.height.current < limits.height.max
                limits = climb(limits)

            elseif limits.can_pant && limits.pant_rate.current < limits.pant_rate.max
                limits, organism_current = pant(organism_current, limits)

            elseif limits.can_burrow
                burrow_bf = limits.burrow_shaded ? 1.0 : 0.0
                limits = select_depth(limits, min_r, max_r, step, limits.shade.current, burrow_bf)
                env    = interpolate_environment(available_environments, step, limits, environmental_params)
                break
            else
                break
            end

        elseif Tb_strip < Tbask * (1 + tolerance)
            # -- Too cold: darken → perpendicular → press to ground → avoid shade → burrow --
            # Target is T_bask (minimum basking temperature); T_preferred_min only used for state.
            if limits.can_change_absorptivity &&
               limits.absorptivity.current < limits.absorptivity.max
                limits, organism_current = darken(organism_current, limits)

            elseif limits.can_solar_orient && limits.sun_orientation < 90.0
                limits, organism_current = orient_perpendicular(organism_current, limits)

            elseif limits.can_press_to_ground && !limits.pressed_to_ground
                limits, organism_current = press_to_ground(organism_current, limits)

            elseif limits.shade.current > limits.shade.reference
                limits = avoid_shade(limits)

            elseif limits.can_burrow && Tb_strip < Tcrit_min
                burrow_bf = limits.burrow_shaded ? 1.0 : 0.0
                limits = select_depth(limits, min_r, max_r, step, limits.shade.current, burrow_bf)
                env    = interpolate_environment(available_environments, step, limits, environmental_params)
                break
            else
                break
            end

        else
            break  # within basking/foraging range [T_bask, T_preferred_max]
        end

        env = interpolate_environment(available_environments, step, limits, environmental_params)
        Tb  = solve_body_temperature(organism_current, env, environmental_params)
    end

    return _build_ectotherm_output(organism_current, env, environmental_params, limits, active, available_environments)
end

# =============================================================================
# Output builder
# =============================================================================

function _build_ectotherm_output(organism, env_vars, env_pars, limits, in_active_period, available_environments)
    e        = (environment_pars=env_pars, environment_vars=env_vars)
    is_burrowed = limits.depth.current > limits.depth.reference
    # When burrowed and burrow_tb_equals_soil=true (default, matches NicheMapR BELOWGROUND.f):
    # set Tb = T_soil directly. The full heat balance gives Tb > T_soil because Q_metab > 0
    # forces equilibrium above the uniform soil temperature when all surrounding temperatures
    # equal T_soil and heat loss pathways are near zero.
    Tb       = (is_burrowed && limits.burrow_tb_equals_soil) ?
               env_vars.T_air :
               solve_body_temperature(organism, env_vars, env_pars)
    ecto_out = ectotherm(Tb, organism, e)
    height = if is_burrowed
        -uconvert(u"m", available_environments.depths[limits.depth.current])
    else
        available_environments.heights[limits.height.current]
    end
    # Classify organism state (mirrors NicheMapR ACT column: 0=Resting, 1=Basking, 2=Active)
    state = if !in_active_period || is_burrowed
        Resting()
    elseif limits.T_preferred_min <= Tb <= limits.T_preferred_max
        Active()
    elseif limits.T_bask <= Tb
        Basking()
    else
        Resting()
    end
    return (;
        T_core                = Tb,
        shade                 = is_burrowed ? (limits.burrow_shaded ? available_environments.max_shade_fraction : available_environments.min_shade_fraction) : limits.shade.current,
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
