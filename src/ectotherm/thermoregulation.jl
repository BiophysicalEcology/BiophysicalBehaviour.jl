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
    zenith < 90u"°" && solar_rad > zero(solar_rad)

is_active(::Nocturnal, zenith, solar_rad) =
    zenith >= 90u"°" || solar_rad <= zero(solar_rad)

is_active(::Crepuscular, zenith, solar_rad) =
    85u"°" <= zenith <= 95u"°"

is_active(a::CombinedActivity, zenith, solar_rad) =
    any(p -> is_active(p, zenith, solar_rad), a.periods)

is_active(a::ResponsiveActivity, zenith, solar_rad) = a.isactive(zenith, solar_rad)

# =============================================================================
# Position helper – blend between two scalar values
# =============================================================================

_blend(min_val, max_val, t) = min_val * (1 - t) + max_val * t

# =============================================================================
# Underground shade blend factor (NicheMapR shdburrow 0/1/2)
# =============================================================================

_underground_blend_factor(::MinShadeOnly, limits, low_shade, high_shade, step) = 0.0
_underground_blend_factor(::MaxShadeOnly, limits, low_shade, high_shade, step) = 1.0

function _underground_blend_factor(::AdaptiveBurrowShade, limits, low_shade, high_shade, step)
    n_nodes  = size(low_shade.soil_temperature, 2)
    min_node = clamp(limits.depth_min_underground, 1, n_nodes)
    T_soil   = low_shade.soil_temperature[step, min_node]
    too_hot  = T_soil > limits.T_active_max
    too_cold = T_soil < limits.T_critical_min
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
    limits = @set limits.T_target.current   = limits.T_target.reference
    limits = @set limits.sun_orientation    = 45.0
    limits = @set limits.pressed_to_ground = false
    return limits
end

"""
    increment_T_target(limits::EctothermBehavioralLimits) → limits

Increment the target body temperature by one step toward `T_active_max`.

Matches NicheMapR's TPREF incrementing (ECTOTHERM.f ENB>0 branch): the organism
tolerates a progressively warmer body temperature before triggering cooling
behaviours. Once `T_target.current` reaches `T_active_max`, shade-seeking and
other cooling responses begin.

Returns updated `EctothermBehavioralLimits`.
"""
function increment_T_target(limits::EctothermBehavioralLimits)
    new_T = min(limits.T_target.current + limits.T_target.step,
                limits.T_target.max)
    @set limits.T_target.current = new_T
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
after basking (perpendicular → intermediate) or when Tb drops back to T_target
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
    select_depth(limits, low_shadeesult, high_shadeesult, step, chosen_shade) → limits

Find the shallowest accessible soil node with a tolerable temperature (SELDEP.f).

Iterates from `limits.depth_min_underground` to `limits.depth.max`, selecting the first
node where the blended soil temperature lies between `T_critical_min` and a
mid-point threshold between `T_active_max` and `T_critical_max`. Falls back to
the deepest node if no node is within tolerance.

# Arguments
- `limits`: Current behavioral limits (depth range and thermal thresholds)
- `low_shadeesult`: `MicroResult` from the minimum-shade microclimate run
- `high_shadeesult`: `MicroResult` from the maximum-shade microclimate run
- `step`: Current simulation time step (1-based hour index)
- `chosen_shade`: Currently chosen shade fraction (0–1), used to blend soil temps
"""
function select_depth(limits::EctothermBehavioralLimits, low_shadeesult, high_shadeesult, step, chosen_shade, blend_factor)
    T_max_threshold = limits.T_critical_max -
        (limits.T_critical_max - limits.T_active_max) / 2

    n_nodes = size(low_shadeesult.soil_temperature, 2)
    depth_min = clamp(limits.depth_min_underground,  1, n_nodes)
    depth_max = clamp(limits.depth.max,        1, n_nodes)

    for node in depth_min:depth_max
        T_soil = _blend(
            low_shadeesult.soil_temperature[step, node],
            high_shadeesult.soil_temperature[step, node],
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
    low_shade = available_environments.min_shade_result
    high_shade = available_environments.max_shade_result

    # Blend factor: 0 = fully min-shade environment, 1 = fully max-shade environment.
    # Matches NicheMapR ABOVEGROUND.f: blend = SHADE / MAXSHD.
    # The min-shade run is treated as the zero-shade reference regardless of its
    # actual shade fraction, so blend(20%) = 20/90 = 0.222 (not 0.0).
    blend_factor = available_environments.max_shade_fraction > 0 ?
        clamp(limits.shade.current / available_environments.max_shade_fraction, 0.0, 1.0) :
        0.0

    chosen_shade = limits.shade.current
    is_underground = limits.depth.current > limits.depth.reference

    # Common fields regardless of position
    zenith_angle = low_shade.solar_radiation.zenith_angle[step]
    P_atmos      = low_shade.pressure[step]

    if is_underground
        # ---- BELOWGROUND environment (BELOWGROUND.f) ----
        # Underground shade is binary: shaded retreat uses max-shade soil temps,
        # unshaded retreat uses min-shade soil temps (NicheMapR shade_burrow flag).
        underground_bf = _underground_blend_factor(limits.burrow_shade_mode, limits, low_shade, high_shade, step)
        node = limits.depth.current
        T_soil = _blend(
            low_shade.soil_temperature[step, node],
            high_shade.soil_temperature[step, node],
            underground_bf,
        )
        k_sub = _blend(
            low_shade.soil_thermal_conductivity[step, node],
            high_shade.soil_thermal_conductivity[step, node],
            underground_bf,
        )
        # Conserve actual vapour pressure at the chosen underground conditions,
        # then compute RH at T_soil (consistent with ABOVEGROUND.f WETAIR approach).
        rh_ref = _blend(low_shade.soil_humidity[step, node], high_shade.soil_humidity[step, node], underground_bf)
        P_v_ref = FluidProperties.vapour_pressure(T_soil) * rh_ref
        rh = clamp(ustrip(u"Pa", P_v_ref) /
                   ustrip(u"Pa", FluidProperties.vapour_pressure(T_soil)), 0.0, 1.0)
        return EnvironmentalVars(;
            air_temperature           = T_soil,
            reference_air_temperature = T_soil,
            sky_temperature           = T_soil,
            ground_temperature        = T_soil,
            substrate_temperature     = T_soil,
            relative_humidity         = clamp(rh, 0.0, 1.0),
            wind_speed                = 0.01u"m/s",
            atmospheric_pressure      = P_atmos,
            zenith_angle,
            substrate_conductivity    = k_sub,
            global_radiation          = zero(low_shade.global_radiation[step]),
            diffuse_fraction          = 0.0,
            shade                     = chosen_shade,
        )
    else
        # ---- ABOVEGROUND environment (ABOVEGROUND.f) ----
        height_node = limits.height.current

        T_air = _blend(
            low_shade.profile.air_temperature[step, height_node],
            high_shade.profile.air_temperature[step, height_node],
            blend_factor,
        )
        # Conserve actual vapour pressure at min-shade reference conditions,
        # then compute RH at the blended temperature (ABOVEGROUND.f WETAIR logic).
        rh_ref  = low_shade.profile.relative_humidity[step, height_node]
        P_v_ref = FluidProperties.vapour_pressure(low_shade.profile.air_temperature[step, height_node]) * rh_ref
        rh      = clamp(ustrip(u"Pa", P_v_ref) /
                        ustrip(u"Pa", FluidProperties.vapour_pressure(T_air)), 0.0, 1.0)
        wind_speed = _blend(
            low_shade.profile.wind_speed[step, height_node],
            high_shade.profile.wind_speed[step, height_node],
            blend_factor,
        )
        T_sky = _blend(
            low_shade.sky_temperature[step],
            high_shade.sky_temperature[step],
            blend_factor,
        )
        T_ground = _blend(
            low_shade.soil_temperature[step, limits.depth.reference],
            high_shade.soil_temperature[step, limits.depth.reference],
            blend_factor,
        )
        T_substrate = T_ground
        k_sub = _blend(
            low_shade.soil_thermal_conductivity[step, limits.depth.reference],
            high_shade.soil_thermal_conductivity[step, limits.depth.reference],
            blend_factor,
        )
        # Solar radiation: stored pre-shade in MicroResult. HeatExchange applies
        # (1 - shade) internally, so pass the raw pre-shade value here.
        global_radiation = low_shade.global_radiation[step]
        diffuse_fraction = _blend(
            low_shade.diffuse_fraction[step],
            high_shade.diffuse_fraction[step],
            blend_factor,
        )

        return EnvironmentalVars(;
            air_temperature = T_air,
            reference_air_temperature = _blend(
                low_shade.reference_temperature[step],
                high_shade.reference_temperature[step],
                blend_factor,
            ),
            sky_temperature       = T_sky,
            ground_temperature    = T_ground,
            substrate_temperature = T_substrate,
            relative_humidity     = clamp(rh, 0.0, 1.0),
            wind_speed,
            atmospheric_pressure  = P_atmos,
            zenith_angle,
            substrate_conductivity = k_sub,
            global_radiation,
            diffuse_fraction,
            shade                  = chosen_shade,
        )
    end
end
