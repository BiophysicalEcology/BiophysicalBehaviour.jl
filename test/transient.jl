using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# No R reference data exists yet for onelump_var.R/twolump.R/trans_behav.R under
# time-varying forcing (see test/R/onelump_var_test.R etc., a manual step). These are
# self-consistency checks: sane trajectories under synthetic and diurnal forcing.

function _diurnal_forcing(times; shade=0.0, wind_speed=1.0u"m/s", radiation_scale=1.0, mean_temp_c=20.0, amplitude_c=8.0)
    hours = ustrip.(u"hr", times)
    n = length(times)
    air_temperature = u"K".((mean_temp_c .+ amplitude_c .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°C")
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

shape_pars = Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
body = Body(shape_pars, Naked())
traits = example_ectotherm_heat_exchange_traits(; shape_pars)
organism = Organism(body, traits)
environment_pars = example_environment_pars()

@testset "simulate_onelump: finite trajectory under near-constant forcing" begin
    times = (0:60:36000)u"s"
    forcing = _diurnal_forcing(times; radiation_scale=0.0)  # near-constant over a short window
    core_temperature_init = u"K"(20.0u"°C")
    sol = simulate_onelump(times, core_temperature_init, organism, environment_pars, forcing)
    @test length(sol.core_temperature) == length(times)
    @test all(isfinite, ustrip.(u"K", sol.core_temperature))
end

@testset "simulate_transient_behavior: sensible daily cycle" begin
    times = (0:1:48)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    limits = example_ectotherm_behavioral_limits()
    # a small body (fast thermal response) so the 4-phase state machine (no deep-refuge
    # escape valve for extreme heat, out of scope) can track it within critical bounds
    small_shape_pars = Ellipsoid(0.05u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
    small_organism = Organism(Body(small_shape_pars, Naked()), example_ectotherm_heat_exchange_traits(; shape_pars=small_shape_pars))

    result = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits,
    )
    @test all(isfinite, ustrip.(u"K", result.core_temperature))
    @test result.t[end] >= times[end]
    # thermoregulation should keep body temperature within the critical range
    @test all(t -> limits.escape_temperature_min <= t <= limits.escape_temperature_max, result.core_temperature)
    # the animal should reach an active state at some point over 2 simulated days
    @test any(s -> s isa Active, result.state)
    @test any(s -> s isa Resting, result.state)
end

# Shared small-body setup for the phase-broadening tests below (fast thermal response, same
# reasoning as the base test above).
small_shape_pars = Ellipsoid(0.05u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
small_organism = Organism(Body(small_shape_pars, Naked()), example_ectotherm_heat_exchange_traits(; shape_pars=small_shape_pars))

@testset "simulate_transient_behavior: burrow substitutes for sleep" begin
    times = (0:1:48)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    underground_forcing = _diurnal_forcing(times; shade=1.0, radiation_scale=0.0, wind_speed=0.01u"m/s")
    limits = example_ectotherm_behavioral_limits()  # can_retreat_underground=true by default

    result = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits;
        underground_forcing,
    )
    @test result.t[end] >= times[end]  # didn't stall
    @test any(p -> p isa BurrowPhase, result.phase)
    @test !any(p -> p isa SleepPhase, result.phase)  # full substitution, not a competing choice
end

@testset "simulate_transient_behavior: climb substitutes for cool" begin
    times = (0:1:48)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0, mean_temp_c=40.0)  # hot enough to force an escape phase
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    climb_forcing = _diurnal_forcing(times; shade=0.0, wind_speed=3.0u"m/s")
    limits = example_ectotherm_behavioral_limits(; can_climb=true)

    result = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits;
        climb_forcing,
    )
    @test result.t[end] >= times[end]  # didn't stall
    @test any(p -> p isa ClimbPhase, result.phase)
    @test !any(p -> p isa CoolPhase, result.phase)  # full substitution, not a competing choice
end

@testset "simulate_transient_behavior: initial_phase" begin
    # Start at hour 10 (clearly daytime, away from the dawn/dusk boundary) so the requested
    # initial phase isn't immediately re-decided by the pre-check loop.
    times = (10:1:34)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    limits = example_ectotherm_behavioral_limits()

    result = simulate_transient_behavior(
        times, u"K"(20.0u"°C"), organism, environment_pars, sun_forcing, shade_forcing, limits;
        initial_phase=BaskPhase(),
    )
    @test result.phase[1] isa BaskPhase

    @test_throws ArgumentError simulate_transient_behavior(
        times, u"K"(30.0u"°C"), organism, environment_pars, sun_forcing, shade_forcing, limits;
        initial_phase=ClimbPhase(),  # can_climb=false, no climb_forcing
    )
end

@testset "simulate_transient_behavior: sleep_forcing" begin
    times = (0:1:24)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")
    limits = example_ectotherm_behavioral_limits()

    default_sleep = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits,
    )
    unshaded_sleep = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits;
        sleep_forcing=sun_forcing,
    )
    # Adaptive step counts differ between runs - compare a summary statistic, not the raw vectors.
    @test extrema(default_sleep.core_temperature) != extrema(unshaded_sleep.core_temperature)
end

@testset "simulate_transient_behavior: threshold validation" begin
    times = (0:1:24)u"hr" .|> u"s"
    sun_forcing = _diurnal_forcing(times; shade=0.0)
    shade_forcing = _diurnal_forcing(times; shade=0.9)
    # basking_temperature_min within active_min_hysteresis (default 0.15K) of active_temperature_min
    # (default 24°C) - passes the base ordering check (still strictly less) but reintroduces the
    # shared-threshold bug class the active_min_hysteresis margin is meant to prevent.
    bad_limits = example_ectotherm_behavioral_limits(; basking_temperature_min=u"K"(24.0u"°C") - 0.05u"K")
    @test_throws ArgumentError simulate_transient_behavior(
        times, u"K"(20.0u"°C"), organism, environment_pars, sun_forcing, shade_forcing, bad_limits,
    )
end

@testset "simulate_transient_behavior: nocturnal activity_period" begin
    times = (0:1:48)u"hr" .|> u"s"
    # _diurnal_forcing clamps zenith to [0,90], pinning the Nocturnal signal at exactly 0 all
    # night - use an unclamped zenith (continuing past 90 toward 180 at solar midnight) instead.
    hours = ustrip.(u"hr", times)
    n = length(times)
    air_temperature = u"K".((20.0 .+ 8.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°C")
    zenith_angle = (90.0 .- 90.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°"
    unclamped(; shade, wind_speed, radiation_scale=1.0) = EnvironmentForcing(times, HeatExchange.EnvironmentalVarsVec(;
        air_temperature, sky_temperature=air_temperature, ground_temperature=air_temperature,
        substrate_temperature=air_temperature, relative_humidity=fill(0.5, n),
        wind_speed=fill(wind_speed, n), atmospheric_pressure=fill(101325.0u"Pa", n),
        zenith_angle, substrate_conductivity=fill(2.79u"W/m/K", n),
        global_radiation=radiation_scale .* max.(0.0, 800.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"W/m^2",
        diffuse_fraction=fill(0.1, n), shade=fill(shade, n),
    ))
    sun_forcing = unclamped(; shade=0.0, wind_speed=1.0u"m/s")
    shade_forcing = unclamped(; shade=0.9, wind_speed=0.5u"m/s", radiation_scale=0.05)
    # Thresholds lowered to fall within night air temperature (~12-20°C) - a nocturnal
    # ectotherm has no metabolic heat, so it can only warm toward ambient, not above it.
    limits = example_ectotherm_behavioral_limits(;
        emerge_temperature_min=u"K"(10.0u"°C"), basking_temperature_min=u"K"(14.0u"°C"),
        active_temperature_min=u"K"(18.0u"°C"),
    )

    result = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits;
        activity_period=Nocturnal(),
    )
    @test all(isfinite, ustrip.(u"K", result.core_temperature))
    @test result.t[end] >= times[end]
    @test any(s -> s isa Active, result.state)
    # Active only near/within the night window for a nocturnal animal - allow a few degrees
    # past the boundary for the default activity_hysteresis margin plus solver step spacing.
    active_zenith = [ustrip(u"°", sun_forcing(result.t[i]).zenith_angle)
                     for i in eachindex(result.state) if result.state[i] isa Active]
    @test all(>=(85.0), active_zenith)
end

@testset "simulate_transient_behavior: crepuscular activity_period" begin
    times = (0:1:48)u"hr" .|> u"s"
    # Same unclamped-zenith rationale as the nocturnal test above.
    hours = ustrip.(u"hr", times)
    n = length(times)
    air_temperature = u"K".((20.0 .+ 8.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°C")
    zenith_angle = (90.0 .- 90.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"°"
    unclamped(; shade, wind_speed, radiation_scale=1.0) = EnvironmentForcing(times, HeatExchange.EnvironmentalVarsVec(;
        air_temperature, sky_temperature=air_temperature, ground_temperature=air_temperature,
        substrate_temperature=air_temperature, relative_humidity=fill(0.5, n),
        wind_speed=fill(wind_speed, n), atmospheric_pressure=fill(101325.0u"Pa", n),
        zenith_angle, substrate_conductivity=fill(2.79u"W/m/K", n),
        global_radiation=radiation_scale .* max.(0.0, 800.0 .* sin.(2π .* (hours .- 6) ./ 24)) .* u"W/m^2",
        diffuse_fraction=fill(0.1, n), shade=fill(shade, n),
    ))
    sun_forcing = unclamped(; shade=0.0, wind_speed=1.0u"m/s")
    shade_forcing = unclamped(; shade=0.9, wind_speed=0.5u"m/s", radiation_scale=0.05)
    # Crepuscular's active window (zenith 85-95°) is brief - thresholds lowered so it's
    # actually reachable within it, same reasoning as the nocturnal test above.
    limits = example_ectotherm_behavioral_limits(;
        emerge_temperature_min=u"K"(10.0u"°C"), basking_temperature_min=u"K"(15.0u"°C"),
        active_temperature_min=u"K"(18.0u"°C"),
    )

    result = simulate_transient_behavior(
        times, u"K"(15.0u"°C"), small_organism, environment_pars, sun_forcing, shade_forcing, limits;
        activity_period=Crepuscular(),
    )
    @test all(isfinite, ustrip.(u"K", result.core_temperature))
    @test result.t[end] >= times[end]
    @test any(s -> s isa Active, result.state)
    # Active only within the twilight band (zenith 85-95°) for a crepuscular animal.
    active_zenith = [ustrip(u"°", sun_forcing(result.t[i]).zenith_angle)
                     for i in eachindex(result.state) if result.state[i] isa Active]
    @test all(z -> 80.0 <= z <= 100.0, active_zenith)
end
