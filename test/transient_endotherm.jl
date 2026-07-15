using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# No R reference exists for this (transient endotherm body temperature is new, not a
# NicheMapR port) - core physics self-consistency (steady-state convergence, sign/thermal-mass
# sanity) is already covered in HeatExchange.jl/test/endotherm_onelump.jl. These are
# wrapper-level smoke tests: sane trajectories under time-varying forcing.

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

environment_pars = example_environment_pars()

@testset "simulate_endotherm_onelump: Naked, diurnal forcing" begin
    shape_pars = example_shape_pars(; mass=0.5u"kg")
    body = Body(shape_pars, Naked())
    traits = example_heat_exchange_traits(;
        shape_pars,
        metabolism_pars=example_metabolism_pars(; core_temperature=u"K"(38.0u"°C"), model=Kleiber()),
    )
    organism = Organism(body, traits)

    times = (0:1:48)u"hr" .|> u"s"
    forcing = _diurnal_forcing(times)
    result = simulate_endotherm_onelump(times, u"K"(38.0u"°C"), organism, environment_pars, forcing)

    @test length(result.t) == length(result.core_temperature) == length(result.core_temperature_rate) ==
        length(result.skin_temperature) == length(result.diagnostics)
    @test all(isfinite, ustrip.(u"K", result.core_temperature))
    @test all(isfinite, ustrip.(u"K/s", result.core_temperature_rate))
end

@testset "simulate_endotherm_onelump: Insulated, diurnal forcing" begin
    shape_pars = example_shape_pars(; mass=0.0337u"kg")
    insulation_pars = example_insulation_pars()
    conduction_pars_internal = example_conduction_pars_internal(; fat_fraction=0.05)
    radiation_pars = example_radiation_pars()
    fat = FatLayer(conduction_pars_internal.fat_fraction, conduction_pars_internal.fat_density)
    pven = radiation_pars.ventral_fraction
    mean_depth = insulation_pars.dorsal.depth * (1 - pven) + insulation_pars.ventral.depth * pven
    mean_diam = insulation_pars.dorsal.diameter * (1 - pven) + insulation_pars.ventral.diameter * pven
    mean_density = insulation_pars.dorsal.density * (1 - pven) + insulation_pars.ventral.density * pven
    fur = FibrousLayer(mean_depth, mean_diam, mean_density)
    body = Body(shape_pars, CompositeInsulation(fur, fat))

    mhf_bmr = metabolic_rate(McKechnieWolf(), shape_pars.mass)
    metabolism_pars = example_metabolism_pars(;
        core_temperature=u"K"(38.0u"°C"), model=McKechnieWolf(), metabolic_heat_flow=mhf_bmr,
    )
    traits = example_heat_exchange_traits(;
        shape_pars, insulation_pars, conduction_pars_internal, radiation_pars, metabolism_pars,
    )
    organism = Organism(body, traits)

    times = (0:1:48)u"hr" .|> u"s"
    forcing = _diurnal_forcing(times; shade=0.9, radiation_scale=0.05, wind_speed=0.5u"m/s")

    @testset "resting (BMR)" begin
        result = simulate_endotherm_onelump(
            times, u"K"(38.0u"°C"), organism, environment_pars, forcing; metabolic_heat_flow=mhf_bmr,
        )
        @test all(isfinite, ustrip.(u"K", result.core_temperature))
        @test all(isfinite, ustrip.(u"K", result.skin_temperature))
    end

    @testset "sustained activity heats the core" begin
        # elevated fixed metabolic rate (e.g. sustained flight) should drive net warming
        active_mhf = mhf_bmr * 4.0
        result = simulate_endotherm_onelump(
            times, u"K"(38.0u"°C"), organism, environment_pars, forcing; metabolic_heat_flow=active_mhf,
        )
        @test all(isfinite, ustrip.(u"K", result.core_temperature))
        @test maximum(result.core_temperature) > u"K"(38.0u"°C")
        # core_temperature crosses a chosen reference threshold at a sensible, inspectable point
        critical = u"K"(42.0u"°C")
        crossing = findfirst(>=(critical), result.core_temperature)
        @test !isnothing(crossing)
        @test result.core_temperature_rate[crossing] > 0u"K/s"
    end
end
