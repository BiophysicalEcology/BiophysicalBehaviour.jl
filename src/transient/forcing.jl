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
    interpolators = NamedTuple(
        field => DataInterpolations.LinearInterpolation(
            getproperty(vars, field), times; extrapolation=DataInterpolations.ExtrapolationType.Constant,
        )
        for field in _ENVIRONMENTAL_VARS_FIELDS
    )
    return EnvironmentForcing(interpolators)
end

function (forcing::EnvironmentForcing)(t; smoothing::HeatExchange.SmoothingStrategy=HeatExchange.HardBound())
    (; interpolators) = forcing
    wind_speed = HeatExchange.safe_relu(smoothing, interpolators.wind_speed(t))
    global_radiation = HeatExchange.safe_relu(smoothing, interpolators.global_radiation(t))
    zenith_angle = HeatExchange.safe_clamp(smoothing, interpolators.zenith_angle(t), 0.0u"°", 90.0u"°")
    return HeatExchange.EnvironmentalVars(;
        air_temperature=interpolators.air_temperature(t),
        reference_air_temperature=interpolators.reference_air_temperature(t),
        sky_temperature=interpolators.sky_temperature(t),
        ground_temperature=interpolators.ground_temperature(t),
        substrate_temperature=interpolators.substrate_temperature(t),
        bush_temperature=interpolators.bush_temperature(t),
        vegetation_temperature=interpolators.vegetation_temperature(t),
        relative_humidity=interpolators.relative_humidity(t),
        wind_speed,
        atmospheric_pressure=interpolators.atmospheric_pressure(t),
        zenith_angle,
        substrate_conductivity=interpolators.substrate_conductivity(t),
        global_radiation,
        diffuse_fraction=interpolators.diffuse_fraction(t),
        shade=interpolators.shade(t),
    )
end
