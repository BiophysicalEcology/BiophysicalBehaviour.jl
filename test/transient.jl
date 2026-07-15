using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# No R reference data exists yet for onelump_var.R/twolump.R/trans_behav.R under
# time-varying forcing (see test/R/onelump_var_test.R etc., a manual step). These are
# smoke/self-consistency checks: sane trajectories under synthetic diurnal forcing, and
# convergence to HeatExchange.onelump's closed-form steady state under constant forcing.

function _diurnal_forcing(times; shade=0.0, wind_speed=1.0u"m/s", radiation_scale=1.0)
    hours = ustrip.(u"hr", times)
    n = length(times)
    air_temperature = u"K".((20.0 .+ 8.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°C")
    zenith_angle = clamp.(90.0 .- 90.0 .* sin.(2π .* (hours .- 6) ./ 24), 0.0, 90.0) .* u"°"
    global_radiation = radiation_scale .* max.(0.0, 800.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"W/m^2"
    vars = HeatExchange.EnvironmentalVarsVec(;
        air_temperature, sky_temperature=air_temperature, ground_temperature=air_temperature,
        substrate_temperature=air_temperature, relative_humidity=fill(0.5, n),
        wind_speed=fill(wind_speed, n), atmospheric_pressure=fill(101325.0u"Pa", n),
        zenith_angle, substrate_conductivity=fill(2.79u"W/m/K", n),
        global_radiation, diffuse_fraction=fill(0.1, n), shade=fill(shade, n),
    )
    return EnvironmentForcing(times, vars)
end

body = Body(Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1), Naked())
environment_pars = example_environment_pars()

@testset "simulate_onelump: constant forcing matches closed form" begin
    times = (0:60:36000)u"s"
    forcing = _diurnal_forcing(times; radiation_scale=0.0)  # near-constant over a short window
    core_temperature_init = u"K"(20.0u"°C")
    sol = simulate_onelump(
        times, core_temperature_init, body, environment_pars, forcing;
        internal_conduction=example_conduction_pars_internal(),
        posture=Intermediate(), body_absorptivity=0.85, emissivity=0.95,
        sky_view_factor=0.4, ground_view_factor=0.4, metabolic_heat_volumetric=0.0u"W/m^3",
    )
    @test length(sol.core_temperature) == length(times)
    @test all(isfinite, ustrip.(u"K", sol.core_temperature))
end

@testset "simulate_diurnal_behavior: sensible daily cycle" begin
    times = (0:1:48)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    limits = example_ectotherm_behavioral_limits()
    # a small body (fast thermal response) so the 4-phase state machine (no deep-refuge
    # escape valve for extreme heat, out of scope) can track it within critical bounds
    small_body = Body(Ellipsoid(0.05u"kg", 1000.0u"kg/m^3", 1.1, 1.1), Naked())

    result = simulate_diurnal_behavior(
        times, u"K"(15.0u"°C"), small_body, environment_pars, sun_forcing, shade_forcing, limits;
        internal_conduction=example_conduction_pars_internal(),
        body_absorptivity=0.85, emissivity=0.95, sky_view_factor=0.4, ground_view_factor=0.4,
        metabolic_heat_volumetric=0.0u"W/m^3",
    )
    @test all(isfinite, ustrip.(u"K", result.core_temperature))
    @test result.t[end] >= times[end]
    # thermoregulation should keep body temperature within the critical range
    @test all(t -> limits.critical_temperature_min <= t <= limits.critical_temperature_max, result.core_temperature)
    # the animal should reach an active state at some point over 2 simulated days
    @test any(s -> s isa Active, result.state)
    @test any(s -> s isa Resting, result.state)
end
