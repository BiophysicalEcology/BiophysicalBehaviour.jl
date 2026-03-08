# =============================================================================
# Internal silhouette-area helper
# =============================================================================

# Prefer the direct shape+orientation dispatch (e.g. DesertIguana has
# silhouette_area(shape, ::ParallelToSun) → scalar) over the body-level
# NamedTuple route (silhouette_area(body) → (; normal, parallel)) which
# requires the 3-arg silhouette_area(shape, insulation, body) form that some
# shapes (notably DesertIguana in the installed package) do not implement.
function _silhouette_area(body, o)
    sh = shape(body)
    if hasmethod(silhouette_area, (typeof(sh), typeof(o)))
        return silhouette_area(sh, o)
    end
    areas = silhouette_area(body)   # → (; normal, parallel)
    o isa NormalToSun   ? areas.normal   :
    o isa ParallelToSun ? areas.parallel :
    (areas.normal + areas.parallel) * 0.5
end

# =============================================================================
# Activity period checks
# =============================================================================

"""
    is_active(activity, zenith_angle, solar_rad) → Bool

Return `true` if the animal is in its active period given current solar conditions.

Ported from ECTOTHERM.f (DAYACT/NOCTURN/CREPUS flags). The zenith angle and
solar radiation together identify day, night, and twilight.

- `Diurnal`: active when the sun is above the horizon (solar_rad > 0, zenith < 90°)
- `Nocturnal`: active at night (solar_rad == 0 or zenith ≥ 90°)
- `Crepuscular`: active at twilight (zenith within ±5° of 90°)
- `CombinedActivity`: active if active in any of its constituent periods
- `ResponsiveActivity`: delegates to the user-supplied function
"""
is_active(::Diurnal, zenith, solar_rad) =
    ustrip(u"°", zenith) < 90 && solar_rad > zero(solar_rad)

is_active(::Nocturnal, zenith, solar_rad) =
    ustrip(u"°", zenith) >= 90 || solar_rad <= zero(solar_rad)

is_active(::Crepuscular, zenith, solar_rad) =
    85 <= ustrip(u"°", zenith) <= 95

is_active(a::CombinedActivity, zenith, solar_rad) =
    any(p -> is_active(p, zenith, solar_rad), a.periods)

is_active(a::ResponsiveActivity, zenith, solar_rad) = a.isactive(zenith, solar_rad)

# =============================================================================
# Position helper – blend between two scalar values
# =============================================================================

_blend(min_val, max_val, t) = min_val * (1 - t) + max_val * t

# =============================================================================
# Individual behaviour functions (SHADEADJUST.f, ABOVEGROUND.f)
# =============================================================================

"""
    seek_shade(limits::EctothermBehavioralLimits) → limits

Increase shade by one step toward the maximum (SHADEADJUST.f).

Returns updated `EctothermBehavioralLimits`.
"""
function seek_shade(limits::EctothermBehavioralLimits)
    new_shade = min(limits.shade.current + limits.shade.step, limits.shade.max)
    @set limits.shade.current = new_shade
end

"""
    avoid_shade(limits::EctothermBehavioralLimits) → limits

Decrease shade by one step toward the reference level (sky exposure).

At night, seeking sky (reducing shade) increases longwave cooling; during the day,
seeking shade reduces solar gain. The name reflects the direction of movement rather
than the thermal outcome, which depends on time of day.

Returns updated `EctothermBehavioralLimits`.
"""
function avoid_shade(limits::EctothermBehavioralLimits)
    new_shade = max(limits.shade.reference, limits.shade.current - limits.shade.step)
    @set limits.shade.current = new_shade
end

"""
    climb(limits::EctothermBehavioralLimits) → limits

Increase height by one node step (ABOVEGROUND.f climbing logic).

Returns updated `EctothermBehavioralLimits`.
"""
function climb(limits::EctothermBehavioralLimits)
    new_height = min(limits.height.current + limits.height.step, limits.height.max)
    @set limits.height.current = new_height
end

"""
    descend(limits::EctothermBehavioralLimits) → limits

Increase depth by one node step (move one node deeper into the soil).

Returns updated `EctothermBehavioralLimits`.
"""
function descend(limits::EctothermBehavioralLimits)
    new_depth = min(limits.depth.current + limits.depth.step, limits.depth.max)
    @set limits.depth.current = new_depth
end

"""
    reset_position(limits::EctothermBehavioralLimits) → limits

Reset all behavioural state to reference values for a new hour.

Resets shade, depth, height, absorptivity, pant rate, posture (perpendicular),
and ground-contact state. Called at the start of each hourly time step.
"""
function reset_position(limits::EctothermBehavioralLimits)
    limits = @set limits.shade.current         = limits.shade.reference
    limits = @set limits.depth.current         = limits.depth.reference
    limits = @set limits.height.current        = limits.height.reference
    limits = @set limits.absorptivity.current  = limits.absorptivity.reference
    limits = @set limits.pant_rate.current     = limits.pant_rate.reference
    limits = @set limits.sun_orientation    = 90.0
    limits = @set limits.pressed_to_ground = false
    return limits
end

# =============================================================================
# Absorptivity behaviour (NicheMapR alpha_max / alpha_min)
# =============================================================================

"""
    darken(organism, limits) → (limits, organism)

Increase dorsal solar absorptivity by one step toward `limits.absorptivity.max`.

Called when body temperature is below the preferred minimum to maximise solar
heat gain (NicheMapR `alpha_max`). Returns updated limits and organism.
"""
function darken(organism::Organism, limits::EctothermBehavioralLimits)
    new_alpha = min(limits.absorptivity.current + limits.absorptivity.step,
                   limits.absorptivity.max)
    limits   = @set limits.absorptivity.current = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.α_body_dorsal = new_alpha
    return limits, organism
end

"""
    lighten(organism, limits) → (limits, organism)

Decrease dorsal solar absorptivity by one step toward `limits.absorptivity.reference`.

Called when body temperature is above the preferred maximum to reduce solar
heat gain (NicheMapR `alpha_min`). Returns updated limits and organism.
"""
function lighten(organism::Organism, limits::EctothermBehavioralLimits)
    new_alpha = max(limits.absorptivity.reference,
                   limits.absorptivity.current - limits.absorptivity.step)
    limits   = @set limits.absorptivity.current = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.α_body_dorsal = new_alpha
    return limits, organism
end

# =============================================================================
# Postural behaviour (NicheMapR postur / orient)
# =============================================================================

"""
    orient_parallel(organism, limits) → (limits, organism)

Orient body parallel to the sun (`sun_orientation = 0.0°`) to minimise
direct-beam solar interception.

Sets `solar_orientation = ParallelToSun()` and `A_silhouette` from body geometry
via `silhouette_area(body, ParallelToSun())`. Only call when
`limits.sun_orientation > 0.0` (organism not already fully parallel).
Returns updated limits and organism.
"""
function orient_parallel(organism::Organism, limits::EctothermBehavioralLimits)
    limits   = @set limits.sun_orientation = 0.0
    organism = @set organism.traits.heat_exchange.radiation_pars.solar_orientation = ParallelToSun()
    organism = @set organism.traits.heat_exchange.radiation_pars.A_silhouette =
               _silhouette_area(organism.body, ParallelToSun())
    return limits, organism
end

"""
    orient_perpendicular(organism, limits) → (limits, organism)

Orient body perpendicular to the sun (`sun_orientation = 90.0°`) to maximise
direct-beam solar interception.

Sets `solar_orientation = NormalToSun()` and computes `A_silhouette` from the
body geometry via `silhouette_area(body, NormalToSun())`. Only call when
`limits.sun_orientation < 90.0` (organism not already fully perpendicular).
Returns updated limits and organism.
"""
function orient_perpendicular(organism::Organism, limits::EctothermBehavioralLimits)
    limits   = @set limits.sun_orientation = 90.0
    organism = @set organism.traits.heat_exchange.radiation_pars.solar_orientation = NormalToSun()
    organism = @set organism.traits.heat_exchange.radiation_pars.A_silhouette =
               _silhouette_area(organism.body, NormalToSun())
    return limits, organism
end

# =============================================================================
# Conduction to ground (NicheMapR pct_cond / PCTCOND)
# =============================================================================

"""
    press_to_ground(organism, limits) → (limits, organism)

Record that the organism is pressed against the substrate.

Sets `pressed_to_ground = true` in limits for output tracking. Conduction fraction
is controlled via organism physiology (`conduction_pars_external`), not here.
Returns updated limits and unchanged organism.
"""
function press_to_ground(organism::Organism, limits::EctothermBehavioralLimits)
    limits = @set limits.pressed_to_ground = true
    return limits, organism
end

# =============================================================================
# Panting – ectotherm dispatch (NicheMapR panting flag)
# =============================================================================

"""
    pant(organism, limits::EctothermBehavioralLimits) → (limits, organism)

Increase ectotherm panting rate by one step for evaporative cooling.

Steps the `pant` multiplier in respiration parameters toward `limits.pant_rate.max`.
Unlike the endotherm version, this does not adjust metabolic rate. Returns updated
limits and organism.
"""
function pant(organism::Organism, limits::EctothermBehavioralLimits)
    new_pant = min(limits.pant_rate.current + limits.pant_rate.step,
                  limits.pant_rate.max)
    limits   = @set limits.pant_rate.current = new_pant
    organism = @set organism.traits.heat_exchange.respiration_pars.pant = new_pant
    return limits, organism
end

# =============================================================================
# Depth selection (SELDEP.f)
# =============================================================================

"""
    select_depth(limits, min_result, max_result, step, chosen_shade) → limits

Find the shallowest accessible soil node with a tolerable temperature (SELDEP.f).

Iterates from `limits.depth_min_burrow` to `limits.depth.max`, selecting the first
node where the blended soil temperature lies between `T_critical_min` and a
mid-point threshold between `T_preferred_max` and `T_critical_max`. Falls back to
the deepest node if no node is within tolerance.

# Arguments
- `limits`: Current behavioral limits (depth range and thermal thresholds)
- `min_result`: `MicroResult` from the minimum-shade microclimate run
- `max_result`: `MicroResult` from the maximum-shade microclimate run
- `step`: Current simulation time step (1-based hour index)
- `chosen_shade`: Currently chosen shade fraction (0–1), used to blend soil temps
"""
function select_depth(limits::EctothermBehavioralLimits, min_result, max_result, step, chosen_shade, blend_factor)
    T_max_threshold = limits.T_critical_max -
        (limits.T_critical_max - limits.T_preferred_max) / 2

    n_nodes = size(min_result.soil_temperature, 2)
    depth_min = clamp(limits.depth_min_burrow, 1, n_nodes)
    depth_max = clamp(limits.depth.max,        1, n_nodes)

    for node in depth_min:depth_max
        T_soil = _blend(
            min_result.soil_temperature[step, node],
            max_result.soil_temperature[step, node],
            blend_factor,
        )
        if limits.T_critical_min < T_soil < T_max_threshold
            return @set limits.depth.current = node
        end
    end

    # Fallback: use deepest available node
    return @set limits.depth.current = depth_max
end

# =============================================================================
# Environment interpolation (ABOVEGROUND.f / BELOWGROUND.f)
# =============================================================================

"""
    interpolate_environment(available_environments, step, limits, environmental_params) → EnvironmentalVars

Build an `EnvironmentalVars` for the organism's current position by blending
between the minimum-shade and maximum-shade `MicroResult` objects.

## Above-ground (depth node == depth.reference, i.e. surface node)

Air temperature, humidity, and wind speed are taken from the atmospheric profile
at `limits.height.current` and blended. Solar radiation is reduced by
`limits.shade.current`. Sky and substrate temperatures are blended. The `shade`
field on `EnvironmentalVars` is set to `limits.shade.current` so that HeatExchange
can adjust longwave view factors (ABOVEGROUND.f).

## Below-ground (depth node > depth.reference)

Sets: Q_solar = 0, wind = 0.01 m/s, rh from soil humidity profile, T_air =
T_sky = T_substrate = T_ground = blended soil temperature at the chosen node
(BELOWGROUND.f).

# Arguments
- `available_environments::AvailableEnvironments`: Wraps min/max shade `MicroResult` objects
- `step::Int`: Simulation hour index
- `limits::EctothermBehavioralLimits`: Current position and thermal limits
- `environmental_params`: `EnvironmentalPars` passed through to the output
"""
function interpolate_environment(available_environments, step, limits::EctothermBehavioralLimits, environmental_params)
    min_r = available_environments.min_shade_result
    max_r = available_environments.max_shade_result

    # Blend factor: 0 = fully min-shade environment, 1 = fully max-shade environment
    shade_range = available_environments.max_shade_fraction - available_environments.min_shade_fraction
    blend_factor = if shade_range > 0
        clamp((limits.shade.current - available_environments.min_shade_fraction) / shade_range, 0.0, 1.0)
    else
        0.0
    end

    chosen_shade = limits.shade.current
    burrowed     = limits.depth.current > limits.depth.reference

    # Common fields regardless of position
    zenith_angle = min_r.solar_radiation.zenith_angle[step]
    P_atmos      = _blend(min_r.pressure[step], max_r.pressure[step], blend_factor)

    if burrowed
        # ---- BELOWGROUND environment (BELOWGROUND.f) ----
        # Burrow shade is binary: shaded burrow uses max-shade soil temps,
        # unshaded burrow uses min-shade soil temps (NicheMapR shade_burrow flag).
        burrow_bf = limits.burrow_shaded ? 1.0 : 0.0
        node = limits.depth.current
        T_soil = _blend(
            min_r.soil_temperature[step, node],
            max_r.soil_temperature[step, node],
            burrow_bf,
        )
        k_sub = _blend(
            min_r.soil_thermal_conductivity[step, node],
            max_r.soil_thermal_conductivity[step, node],
            burrow_bf,
        )
        rh = _blend(
            min_r.soil_humidity[step, node],
            max_r.soil_humidity[step, node],
            burrow_bf,
        )
        return EnvironmentalVars(;
            T_air            = T_soil,
            T_air_reference  = T_soil,
            T_sky            = T_soil,
            T_ground         = T_soil,
            T_substrate      = T_soil,
            rh               = clamp(rh, 0.0, 1.0),
            wind_speed       = 0.01u"m/s",
            P_atmos,
            zenith_angle,
            k_substrate      = k_sub,
            global_radiation = zero(min_r.global_radiation[step]),
            diffuse_fraction = 0.0,
            shade            = chosen_shade,
        )
    else
        # ---- ABOVEGROUND environment (ABOVEGROUND.f) ----
        height_node = limits.height.current
        min_profile = min_r.profile[step]
        max_profile = max_r.profile[step]

        T_air = _blend(
            min_profile.air_temperature[height_node],
            max_profile.air_temperature[height_node],
            blend_factor,
        )
        rh = _blend(
            min_profile.relative_humidity[height_node],
            max_profile.relative_humidity[height_node],
            blend_factor,
        )
        wind_speed = _blend(
            min_profile.wind_speed[height_node],
            max_profile.wind_speed[height_node],
            blend_factor,
        )
        T_sky = _blend(
            min_r.sky_temperature[step],
            max_r.sky_temperature[step],
            blend_factor,
        )
        T_ground = _blend(
            min_r.soil_temperature[step, limits.depth.reference],
            max_r.soil_temperature[step, limits.depth.reference],
            blend_factor,
        )
        T_substrate = T_ground
        k_sub = _blend(
            min_r.soil_thermal_conductivity[step, limits.depth.reference],
            max_r.soil_thermal_conductivity[step, limits.depth.reference],
            blend_factor,
        )
        # Solar radiation: stored pre-shade in MicroResult; reduce by chosen shade
        global_radiation = min_r.global_radiation[step] * (1 - chosen_shade)
        diffuse_fraction = _blend(
            min_r.diffuse_fraction[step],
            max_r.diffuse_fraction[step],
            blend_factor,
        )

        return EnvironmentalVars(;
            T_air,
            T_air_reference  = _blend(
                min_r.reference_temperature[step],
                max_r.reference_temperature[step],
                blend_factor,
            ),
            T_sky,
            T_ground,
            T_substrate,
            rh               = clamp(rh, 0.0, 1.0),
            wind_speed,
            P_atmos,
            zenith_angle,
            k_substrate      = k_sub,
            global_radiation,
            diffuse_fraction,
            shade            = chosen_shade,
        )
    end
end
