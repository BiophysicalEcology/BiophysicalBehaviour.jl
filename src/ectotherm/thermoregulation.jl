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
    is_active(activity, zenith_angle, sunlight) → Bool

Return `true` if the animal is in its active period given current solar conditions.

Ported from ECTOTHERM.f (DAYACT/NOCTURN/CREPUS flags). The zenith angle and
solar radiation together identify day, night, and twilight.

- `Diurnal`: active when the sun is above the horizon (sunlight > 0, zenith < 90°)
- `Nocturnal`: active at night (sunlight == 0 or zenith ≥ 90°)
- `Crepuscular`: active at twilight (zenith within ±5° of 90°)
- `CombinedActivity`: active if active in any of its constituent periods
- `ResponsiveActivity`: delegates to the user-supplied function
"""
is_active(::Diurnal, zenith, sunlight) =
    zenith < 90u"°" && sunlight > zero(sunlight)

is_active(::Nocturnal, zenith, sunlight) =
    zenith >= 90u"°" || sunlight <= zero(sunlight)

is_active(::Crepuscular, zenith, sunlight) =
    85u"°" <= zenith <= 95u"°"

is_active(a::CombinedActivity, zenith, sunlight) =
    any(p -> is_active(p, zenith, sunlight), a.periods)

is_active(a::ResponsiveActivity, zenith, sunlight) = a.isactive(zenith, sunlight)

# =============================================================================
# Position helper – blend between two scalar values
# =============================================================================

_blend(min_val, max_val, t) = min_val * (1 - t) + max_val * t

# =============================================================================
# Underground shade blend factor (NicheMapR shdburrow 0/1/2)
# =============================================================================

_underground_blend_factor(::MinShadeOnly, limits, min_shade, max_shade, step) = 0.0
_underground_blend_factor(::MaxShadeOnly, limits, min_shade, max_shade, step) = 1.0

function _underground_blend_factor(::AdaptiveBurrowShade, limits, min_shade, max_shade, step)
    n_nodes          = size(min_shade.soil_temperature, 2)
    min_node         = clamp(limits.depth_min_underground, 1, n_nodes)
    soil_temperature = min_shade.soil_temperature[step, min_node]
    too_hot          = soil_temperature > limits.active_temperature_max
    too_cold         = soil_temperature < limits.escape_temperature_min
    return (too_hot || too_cold) ? 1.0 : 0.0
end

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
    limits = @set limits.absorptivity.current  = limits.absorptivity.max
    limits = @set limits.pant_rate.current     = limits.pant_rate.reference
    limits = @set limits.target_temperature.current  = limits.target_temperature.reference
    limits = @set limits.sun_orientation   = 45.0
    limits = @set limits.pressed_to_ground = false
    return limits
end

"""
    increment_target_temperature(limits::EctothermBehavioralLimits) → limits

Increment the target body temperature by one step toward `active_temperature_max`.

Matches NicheMapR's TPREF incrementing (ECTOTHERM.f ENB>0 branch): the organism
tolerates a progressively warmer body temperature before triggering cooling
behaviours. Once `target_temperature.current` reaches `active_temperature_max`, shade-seeking and
other cooling responses begin.

Returns updated `EctothermBehavioralLimits`.
"""
function increment_target_temperature(limits::EctothermBehavioralLimits)
    new_T = min(limits.target_temperature.current + limits.target_temperature.step,
                limits.target_temperature.max)
    @set limits.target_temperature.current = new_T
end

# =============================================================================
# Absorptivity behaviour (NicheMapR alpha_max / alpha_min)
# =============================================================================

"""
    darken(organism, limits) → (limits, organism)

Increase dorsal and ventral solar absorptivity by one step toward `limits.absorptivity.max`.

Called when body temperature is below the target minimum to maximise solar
heat gain (NicheMapR `alpha_max`). Returns updated limits and organism.
"""
function darken(organism::Organism, limits::EctothermBehavioralLimits)
    new_alpha = min(limits.absorptivity.current + limits.absorptivity.step,
                   limits.absorptivity.max)
    limits   = @set limits.absorptivity.current = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.body_absorptivity_dorsal = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.body_absorptivity_ventral = new_alpha
    return limits, organism
end

"""
    lighten(organism, limits) → (limits, organism)

Decrease dorsal solar absorptivity by one step toward `limits.absorptivity.reference`.

Called when body temperature is above the target maximum to reduce solar
heat gain (NicheMapR `alpha_min`). Returns updated limits and organism.
"""
function lighten(organism::Organism, limits::EctothermBehavioralLimits)
    new_alpha = max(limits.absorptivity.reference,
                   limits.absorptivity.current - limits.absorptivity.step)
    limits   = @set limits.absorptivity.current = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.body_absorptivity_dorsal = new_alpha
    organism = @set organism.traits.heat_exchange.radiation_pars.body_absorptivity_ventral = new_alpha
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
    organism = @set organism.traits.heat_exchange.radiation_pars.silhouette_area =
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
    organism = @set organism.traits.heat_exchange.radiation_pars.silhouette_area =
               _silhouette_area(organism.body, NormalToSun())
    return limits, organism
end

"""
    orient_intermediate(organism, limits) → (limits, organism)

Revert solar orientation to the neutral/foraging posture (`sun_orientation = 45.0°`).

Sets `solar_orientation = Intermediate()` and `A_silhouette` to the average of
normal and parallel areas. Called when body temperature re-enters the active range
after basking (perpendicular → intermediate) or when Tb drops back to target_temperature
after parallel cooling (parallel → intermediate).

Mirrors NicheMapR THERMO.f lines 119-128: revert when TC ≥ TMINPR (from perpendicular)
or TC ≤ TPREF (from parallel).
Returns updated limits and organism.
"""
function orient_intermediate(organism::Organism, limits::EctothermBehavioralLimits)
    A_intermediate = (_silhouette_area(organism.body, NormalToSun()) +
                      _silhouette_area(organism.body, ParallelToSun())) * 0.5
    limits   = @set limits.sun_orientation = 45.0
    organism = @set organism.traits.heat_exchange.radiation_pars.solar_orientation = Intermediate()
    organism = @set organism.traits.heat_exchange.radiation_pars.silhouette_area = A_intermediate
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
    select_depth(limits, min_shade, max_shade, step, underground_shade_factor) → limits

Find the shallowest accessible soil node with a tolerable temperature (SELDEP.f).

Iterates from `limits.depth_min_underground` to `limits.depth.max`, selecting the first
node where the blended soil temperature lies between `escape_temperature_min` and
`escape_temperature_max`. Falls back to the deepest node if no node is within tolerance.

# Arguments
- `limits`: Current behavioral limits (depth range and thermal thresholds)
- `min_shade`: `MicroResult` from the minimum-shade microclimate run
- `max_shade`: `MicroResult` from the maximum-shade microclimate run
- `step`: Current simulation time step (1-based hour index)
- `underground_shade_factor`: Blend factor (0=min-shade, 1=max-shade) for soil temperatures
"""
function select_depth(limits::EctothermBehavioralLimits, min_shade, max_shade, step, underground_shade_factor)
    depth_selection_threshold = limits.escape_temperature_max

    n_nodes   = size(min_shade.soil_temperature, 2)
    depth_min = clamp(limits.depth_min_underground, 1, n_nodes)
    depth_max = clamp(limits.depth.max,             1, n_nodes)

    for node in depth_min:depth_max
        soil_temperature = _blend(
            min_shade.soil_temperature[step, node],
            max_shade.soil_temperature[step, node],
            underground_shade_factor,
        )
        if limits.escape_temperature_min < soil_temperature < depth_selection_threshold
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

Air temperature, wind speed, sky/substrate temperatures are blended between
min- and max-shade results. Solar radiation is reduced by `limits.shade.current`.
The `shade` field on `EnvironmentalVars` is set to `limits.shade.current` so that
HeatExchange can adjust longwave view factors (ABOVEGROUND.f).

Relative humidity is computed by conserving the actual vapour pressure at the
min-shade reference temperature (matching ABOVEGROUND.f WETAIR approach): the
actual vapour pressure at `(T_air_min_shade, rh_min_shade)` is preserved, then
divided by the saturation vapour pressure at the blended `T_air` to give the
corrected RH.

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
    min_shade = available_environments.min_shade_result
    max_shade = available_environments.max_shade_result

    # Shade blend factor: 0 = fully min-shade environment, 1 = fully max-shade environment.
    # Interpolates shade.current between min_shade_fraction and max_shade_fraction.
    shade_range        = available_environments.max_shade_fraction - available_environments.min_shade_fraction
    shade_blend_factor = shade_range > 0 ?
        clamp((limits.shade.current - available_environments.min_shade_fraction) / shade_range, 0.0, 1.0) :
        0.0

    current_shade  = limits.shade.current
    is_underground = limits.depth.current > limits.depth.reference

    # Common fields regardless of position
    zenith_angle        = min_shade.solar_radiation.zenith_angle[step]
    atmospheric_pressure = min_shade.pressure[step]

    if is_underground
        # ---- BELOWGROUND environment (BELOWGROUND.f) ----
        # Underground shade is binary: shaded retreat uses max-shade soil temps,
        # unshaded retreat uses min-shade soil temps (NicheMapR shade_burrow flag).
        underground_shade_factor = _underground_blend_factor(limits.burrow_shade_mode, limits, min_shade, max_shade, step)
        node             = limits.depth.current
        soil_temperature = _blend(
            min_shade.soil_temperature[step, node],
            max_shade.soil_temperature[step, node],
            underground_shade_factor,
        )
        substrate_conductivity = _blend(
            min_shade.soil_thermal_conductivity[step, node],
            max_shade.soil_thermal_conductivity[step, node],
            underground_shade_factor,
        )
        # Conserve actual vapour pressure at the chosen underground conditions,
        # then compute RH at soil_temperature (consistent with ABOVEGROUND.f WETAIR approach).
        rh_ref  = _blend(min_shade.soil_humidity[step, node], max_shade.soil_humidity[step, node], underground_shade_factor)
        P_v_ref = FluidProperties.vapour_pressure(soil_temperature) * rh_ref
        rh      = clamp(ustrip(u"Pa", P_v_ref) /
                        ustrip(u"Pa", FluidProperties.vapour_pressure(soil_temperature)), 0.0, 1.0)
        return EnvironmentalVars(;
            air_temperature           = soil_temperature,
            reference_air_temperature = soil_temperature,
            sky_temperature           = soil_temperature,
            ground_temperature        = soil_temperature,
            substrate_temperature     = soil_temperature,
            relative_humidity         = clamp(rh, 0.0, 1.0),
            wind_speed                = 0.01u"m/s",
            atmospheric_pressure,
            zenith_angle,
            substrate_conductivity,
            global_radiation          = zero(min_shade.global_radiation[step]),
            diffuse_fraction          = 0.0,
            shade                     = current_shade,
        )
    else
        # ---- ABOVEGROUND environment (ABOVEGROUND.f) ----
        height_node = limits.height.current

        air_temperature = _blend(
            min_shade.profile.air_temperature[step, height_node],
            max_shade.profile.air_temperature[step, height_node],
            shade_blend_factor,
        )
        # Conserve actual vapour pressure at min-shade reference conditions,
        # then compute RH at the blended temperature (ABOVEGROUND.f WETAIR logic).
        rh_ref  = min_shade.profile.relative_humidity[step, height_node]
        P_v_ref = FluidProperties.vapour_pressure(min_shade.profile.air_temperature[step, height_node]) * rh_ref
        rh      = clamp(ustrip(u"Pa", P_v_ref) /
                        ustrip(u"Pa", FluidProperties.vapour_pressure(air_temperature)), 0.0, 1.0)
        wind_speed = _blend(
            min_shade.profile.wind_speed[step, height_node],
            max_shade.profile.wind_speed[step, height_node],
            shade_blend_factor,
        )
        sky_temperature = _blend(
            min_shade.sky_temperature[step],
            max_shade.sky_temperature[step],
            shade_blend_factor,
        )
        ground_temperature = _blend(
            min_shade.soil_temperature[step, limits.depth.reference],
            max_shade.soil_temperature[step, limits.depth.reference],
            shade_blend_factor,
        )
        substrate_conductivity = _blend(
            min_shade.soil_thermal_conductivity[step, limits.depth.reference],
            max_shade.soil_thermal_conductivity[step, limits.depth.reference],
            shade_blend_factor,
        )
        # Solar radiation: stored pre-shade in MicroResult. HeatExchange applies
        # (1 - shade) internally, so pass the raw pre-shade value here.
        global_radiation = min_shade.global_radiation[step]
        diffuse_fraction = _blend(
            min_shade.diffuse_fraction[step],
            max_shade.diffuse_fraction[step],
            shade_blend_factor,
        )

        return EnvironmentalVars(;
            air_temperature,
            reference_air_temperature = _blend(
                min_shade.reference_temperature[step],
                max_shade.reference_temperature[step],
                shade_blend_factor,
            ),
            sky_temperature,
            ground_temperature,
            substrate_temperature = ground_temperature,
            relative_humidity     = clamp(rh, 0.0, 1.0),
            wind_speed,
            atmospheric_pressure,
            zenith_angle,
            substrate_conductivity,
            global_radiation,
            diffuse_fraction,
            shade = current_shade,
        )
    end
end
