const _ENVIRONMENTAL_VARS_FIELDS = (
    :air_temperature, :reference_air_temperature, :sky_temperature, :ground_temperature,
    :substrate_temperature, :bush_temperature, :vegetation_temperature, :relative_humidity,
    :wind_speed, :atmospheric_pressure, :zenith_angle, :substrate_conductivity,
    :global_radiation, :diffuse_fraction, :shade,
)

"""
    EnvironmentForcing(times, vars::HeatExchange.EnvironmentalVarsVec)
    EnvironmentForcing(times, vars)(t; smoothing=HardBound())

Linear-interpolation wrapper around a time series of environmental variables,
callable at arbitrary `t` to produce a scalar environment NamedTuple (same field
names as `HeatExchange.EnvironmentalVars`) for `HeatExchange.onelump`/`twolump`.
Mirrors NicheMapR's `approxfun`-based time-varying forcing.
"""
struct EnvironmentForcing{I}
    interpolators::I
end

function EnvironmentForcing(times, vars::HeatExchange.EnvironmentalVarsVec)
    # Stripped to plain seconds: DataInterpolations' fast path for Float64-valued
    # fields (shade, diffuse_fraction, ...) mishandles a Unitful time grid.
    times_s = ustrip.(u"s", times)
    interpolators = NamedTuple(
        field => DataInterpolations.LinearInterpolation(
            getproperty(vars, field), times_s; extrapolation=DataInterpolations.ExtrapolationType.Constant,
        )
        for field in _ENVIRONMENTAL_VARS_FIELDS
    )
    return EnvironmentForcing(interpolators)
end

function (forcing::EnvironmentForcing)(t; smoothing::HeatExchange.SmoothingStrategy=HeatExchange.HardBound())
    (; interpolators) = forcing
    t_s = ustrip(u"s", t)
    wind_speed = HeatExchange.safe_relu(smoothing, interpolators.wind_speed(t_s))
    global_radiation = HeatExchange.safe_relu(smoothing, interpolators.global_radiation(t_s))
    zenith_angle = HeatExchange.safe_clamp(smoothing, interpolators.zenith_angle(t_s), 0.0u"°", 90.0u"°")
    return HeatExchange.EnvironmentalVars(;
        air_temperature=interpolators.air_temperature(t_s),
        reference_air_temperature=interpolators.reference_air_temperature(t_s),
        sky_temperature=interpolators.sky_temperature(t_s),
        ground_temperature=interpolators.ground_temperature(t_s),
        substrate_temperature=interpolators.substrate_temperature(t_s),
        bush_temperature=interpolators.bush_temperature(t_s),
        vegetation_temperature=interpolators.vegetation_temperature(t_s),
        relative_humidity=interpolators.relative_humidity(t_s),
        wind_speed,
        atmospheric_pressure=interpolators.atmospheric_pressure(t_s),
        zenith_angle,
        substrate_conductivity=interpolators.substrate_conductivity(t_s),
        global_radiation,
        diffuse_fraction=interpolators.diffuse_fraction(t_s),
        shade=interpolators.shade(t_s),
    )
end
