# =============================================================================
# Endotherm behavioural habitat-selection layer.
#
# Outer loop selects shade, height, depth, posture, and absorptivity to bring
# operative temperature into the active window before handing off to the
# inner physiological thermoregulate (rule-based or IPOPT) at the chosen
# position.
# =============================================================================

"""
    thermoregulate(::Endotherm, ::RuleBasedSequentialControl, organism,
                   available_environments::AvailableEnvironments,
                   behavioral_limits::EctothermBehavioralLimits,
                   environmental_params, step[, prev_step_depth])

Endotherm thermoregulation with an outer behavioural habitat-selection loop
followed by the inner physiological thermoregulation loop.

## Motivation

Panting and sweating cost water. An endotherm should therefore exhaust cheap
positional options (shade, climbing, retreating underground, posture, absorptivity change)
**before** the physiology resorts to evaporative cooling.

## Behavioural loop trigger

Operative temperature `Te` (= `solve_body_temperature`: the equilibrium body
temperature of a passive body in the current environment) is compared to
`behavioral_limits.active_temperature_min/max`:

- `Te > active_temperature_max` → environment too hot in the open → seek cool microhabitat
- `Te < active_temperature_min` → environment too cold → seek warm microhabitat

## Priority sequences

*Too hot* (water-loss avoidance first, no behavioural panting):
`lighten` → `orient_parallel` → `seek_shade` → `climb` → `retreat_underground`

*Too cold*:
`darken` → `orient_perpendicular` → `press_to_ground` → `avoid_shade` →
`retreat_underground` (only if `Te < critical_temperature_min`)

## After the behavioural loop

The full endotherm physiological thermoregulate (piloerection, vasodilation,
hyperthermia, panting, sweating — `rulebased.jl` or `ipopt.jl`) is run at the
chosen behavioural position.

## Returns

NamedTuple with:
- `Te`: operative temperature (K)
- `shade`, `depth`, `height`: chosen positional state
- `absorptivity`, `posture_perpendicular`, `pressed_to_ground`, `pant_rate`: chosen
  physical state
- `active`: `true` if within the organism's activity period
- `endotherm_out`: full output of the endotherm physiological thermoregulate at the
  chosen position (fields: `thermoregulation`, `morphology`, `energy_fluxes`,
  `mass_fluxes`)
"""
function thermoregulate(
    ::Endotherm,
    ::RuleBasedSequentialControl,
    organism::Organism,
    available_environments::AvailableEnvironments,
    behavioral_limits::EctothermBehavioralLimits,
    environmental_params,
    step::Int,
    prev_step_depth::Int = 1,
)
    min_shade = available_environments.min_shade_result
    max_shade = available_environments.max_shade_result
    (; tolerance, max_iterations) = behavioral_limits.control

    # -------------------------------------------------------------------------
    # 1. Reset behavioural state; organism starts at active reference state
    # -------------------------------------------------------------------------
    behavioral_limits = reset_position(behavioral_limits)
    organism_current  = organism

    # Clamp height/depth max to actual environment sizes.
    behavioral_limits = @set behavioral_limits.height.max = min(behavioral_limits.height.max, length(available_environments.heights))
    behavioral_limits = @set behavioral_limits.depth.max  = min(behavioral_limits.depth.max,  length(available_environments.depths))

    # When sun orientation is disabled, fix silhouette at Intermediate().
    if !behavioral_limits.can_solar_orient
        organism_current = @set organism_current.traits.heat_exchange.radiation_pars.solar_orientation = Intermediate()
        organism_current = @set organism_current.traits.heat_exchange.radiation_pars.silhouette_area =
            _silhouette_area(organism_current.body, Intermediate())
    end

    # -------------------------------------------------------------------------
    # 2. Determine activity period
    # -------------------------------------------------------------------------
    zenith    = min_shade.solar_radiation.zenith_angle[step]
    sunlight = min_shade.global_radiation[step]
    active    = is_active(activity_period(organism), zenith, sunlight)

    # Underground shade factor helper (binary 0/1, matches _underground_blend_factor logic)
    shade_range = available_environments.max_shade_fraction - available_environments.min_shade_fraction
    @inline underground_shade_factor_for(shade_chosen) = shade_range > 0 ?
        clamp((shade_chosen - available_environments.min_shade_fraction) / shade_range, 0.0, 1.0) :
        0.0

    # -------------------------------------------------------------------------
    # 3a. Inactive period → retreat to optimal underground depth (or stay above)
    # -------------------------------------------------------------------------
    if !active
        if behavioral_limits.can_retreat_underground
            underground_shade_factor = underground_shade_factor_for(behavioral_limits.shade.current)
            behavioral_limits = select_depth(behavioral_limits, min_shade, max_shade, step,
                                             underground_shade_factor)
        end
        env = interpolate_environment(available_environments, step, behavioral_limits, environmental_params)
        return _build_endotherm_behavioral_output(
            organism_current, env, environmental_params, behavioral_limits, active, available_environments)
    end

    # -------------------------------------------------------------------------
    # 3b. Active but was underground last step → check emerge_temperature_min condition
    # -------------------------------------------------------------------------
    if prev_step_depth > behavioral_limits.depth.reference
        underground_shade_factor = underground_shade_factor_for(behavioral_limits.shade.reference)
        soil_temperature_prev    = _blend(
            min_shade.soil_temperature[step, prev_step_depth],
            max_shade.soil_temperature[step, prev_step_depth],
            underground_shade_factor,
        )
        if soil_temperature_prev < behavioral_limits.emerge_temperature_min
            behavioral_limits = @set behavioral_limits.depth.current = prev_step_depth
            env = interpolate_environment(available_environments, step, behavioral_limits, environmental_params)
            return _build_endotherm_behavioral_output(
                organism_current, env, environmental_params, behavioral_limits, false, available_environments)
        end
    end

    # -------------------------------------------------------------------------
    # 4. Active: behavioural loop using operative temperature as trigger
    # -------------------------------------------------------------------------
    env = interpolate_environment(available_environments, step, behavioral_limits, environmental_params)
    Te  = solve_body_temperature(organism_current, env, environmental_params)

    iteration = 0
    while iteration < max_iterations
        iteration += 1

        if Te > behavioral_limits.active_temperature_max * (1 - tolerance)
            # Too hot: seek cool microhabitat before panting (water-loss avoidance)
            # lighten → parallel → shade → climb → retreat_underground
            if behavioral_limits.can_change_absorptivity &&
               behavioral_limits.absorptivity.current > behavioral_limits.absorptivity.reference
                behavioral_limits, organism_current = lighten(organism_current, behavioral_limits)

            elseif behavioral_limits.can_solar_orient && behavioral_limits.sun_orientation > 0.0
                behavioral_limits, organism_current = orient_parallel(organism_current, behavioral_limits)

            elseif behavioral_limits.can_seek_shade &&
                   behavioral_limits.shade.current < behavioral_limits.shade.max
                behavioral_limits = seek_shade(behavioral_limits)

            elseif behavioral_limits.can_climb &&
                   behavioral_limits.height.current < behavioral_limits.height.max
                behavioral_limits = @set behavioral_limits.shade.current = behavioral_limits.shade.reference
                behavioral_limits = climb(behavioral_limits)

            elseif behavioral_limits.can_retreat_underground
                underground_shade_factor = underground_shade_factor_for(behavioral_limits.shade.current)
                behavioral_limits = select_depth(behavioral_limits, min_shade, max_shade, step,
                                                 underground_shade_factor)
                break
            else
                break
            end

        elseif Te < behavioral_limits.active_temperature_min * (1 + tolerance)
            # Too cold: seek warm microhabitat
            # darken → perpendicular → press to ground → avoid shade → retreat_underground
            if behavioral_limits.can_change_absorptivity &&
               behavioral_limits.absorptivity.current < behavioral_limits.absorptivity.max
                behavioral_limits, organism_current = darken(organism_current, behavioral_limits)

            elseif behavioral_limits.can_solar_orient && behavioral_limits.sun_orientation < 90.0
                behavioral_limits, organism_current = orient_perpendicular(organism_current, behavioral_limits)

            elseif behavioral_limits.can_press_to_ground && !behavioral_limits.pressed_to_ground
                behavioral_limits, organism_current = press_to_ground(organism_current, behavioral_limits)

            elseif behavioral_limits.shade.current > behavioral_limits.shade.reference
                behavioral_limits = avoid_shade(behavioral_limits)

            elseif behavioral_limits.can_retreat_underground && Te < behavioral_limits.critical_temperature_min
                underground_shade_factor = underground_shade_factor_for(behavioral_limits.shade.current)
                behavioral_limits = select_depth(behavioral_limits, min_shade, max_shade, step,
                                                 underground_shade_factor)
                break
            else
                break
            end

        else
            break  # operative temperature within target range
        end

        env = interpolate_environment(available_environments, step, behavioral_limits, environmental_params)
        Te  = solve_body_temperature(organism_current, env, environmental_params)
    end

    return _build_endotherm_behavioral_output(
        organism_current, env, environmental_params, behavioral_limits, active, available_environments)
end

# =============================================================================
# Output builder — runs physiological thermoregulate at chosen position
# =============================================================================

function _build_endotherm_behavioral_output(organism, env_vars, env_pars, behavioral_limits, active, available_environments)
    Te  = solve_body_temperature(organism, env_vars, env_pars)
    e   = (environment_pars=env_pars, environment_vars=env_vars)

    init = initial_physiological_state(organism, env_vars)
    endotherm_out = thermoregulate(organism, e,
                                   init.metabolic_heat_flow,
                                   init.skin_temperature,
                                   init.insulation_temperature)

    return (;
        Te,
        shade                 = behavioral_limits.shade.current,
        depth                 = available_environments.depths[behavioral_limits.depth.current],
        depth_node            = behavioral_limits.depth.current,
        height                = available_environments.heights[behavioral_limits.height.current],
        height_node           = behavioral_limits.height.current,
        absorptivity          = behavioral_limits.absorptivity.current,
        sun_orientation       = behavioral_limits.sun_orientation,
        pressed_to_ground     = behavioral_limits.pressed_to_ground,
        pant_rate             = behavioral_limits.pant_rate.current,
        active,
        endotherm_out,
    )
end
