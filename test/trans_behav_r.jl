using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful
using DataFrames, CSV
using Test

# R-vs-Julia comparison against a real trans_behav.R run (test/R/trans_behav.R), using
# actual NicheMapR micro_global forcing rather than a synthetic approximation.

datadir = joinpath(@__DIR__, "data", "trans_behav")
params = DataFrame(CSV.File(joinpath(datadir, "trans_behav_params.csv")))[1, :]
metout = DataFrame(CSV.File(joinpath(datadir, "trans_behav_metout.csv")))
shadmet = DataFrame(CSV.File(joinpath(datadir, "trans_behav_shadmet.csv")))
soil = DataFrame(CSV.File(joinpath(datadir, "trans_behav_soil.csv")))
shadsoil = DataFrame(CSV.File(joinpath(datadir, "trans_behav_shadsoil.csv")))
day_results = DataFrame(CSV.File(joinpath(datadir, "trans_behav_day_results.csv")))
sum_stats = DataFrame(CSV.File(joinpath(datadir, "trans_behav_sum_stats.csv")))[1, :]

# SOLR is unattenuated in both met/shadmet; shade is applied via the `shade` fraction.
function forcing_from(met, soil, pressure, shade_fraction)
    n = nrow(met)
    times = met.TIME .* u"minute" .|> u"s"
    air_temperature = u"K".(met.TALOC .* u"°C")
    EnvironmentForcing(
        times,
        HeatExchange.EnvironmentalVarsVec(;
            air_temperature, sky_temperature=u"K".(met.TSKYC .* u"°C"),
            ground_temperature=u"K".(soil.D0cm .* u"°C"), substrate_temperature=u"K".(soil.D0cm .* u"°C"),
            relative_humidity=met.RHLOC ./ 100, wind_speed=met.VLOC .* u"m/s",
            atmospheric_pressure=fill(pressure * u"Pa", n), zenith_angle=met.ZEN .* u"°",
            substrate_conductivity=fill(2.79u"W/m/K", n), global_radiation=met.SOLR .* u"W/m^2",
            diffuse_fraction=fill(0.1, n), shade=fill(shade_fraction, n),
        ),
    )
end

sun_forcing = forcing_from(metout, soil, params.press, 0.0)
shade_forcing = forcing_from(shadmet, shadsoil, params.press, params.shade / 100)

# R's shape_b/shape_c are minor:major axis ratios (b axis:a axis); Julia's
# axis_ratio_b/axis_ratio_c are major:minor (a:b) — invert to match.
body = Body(Ellipsoid(params.Ww_g * u"g", params.rho_body * u"kg/m^3", 1 / params.shape_b, 1 / params.shape_c), Naked())
environment_pars = example_environment_pars(; ground_albedo=params.alpha_sub)
internal_conduction = example_conduction_pars_internal(;
    flesh_conductivity=params.k_flesh * u"W/m/K", flesh_specific_heat=params.c_body * u"J/kg/K",
)
limits = example_ectotherm_behavioral_limits(;
    active_temperature_min=u"K"(Float64(params.T_F_min) * u"°C"), active_temperature_max=u"K"(Float64(params.T_F_max) * u"°C"),
    basking_temperature_min=u"K"(Float64(params.T_B_min) * u"°C"), critical_temperature_max=u"K"(Float64(params.CT_max) * u"°C"),
)

times = metout.TIME .* u"minute" .|> u"s"
core_temperature_init = u"K"(day_results.Tb[1] * u"°C")

result = simulate_diurnal_behavior(
    times, core_temperature_init, body, environment_pars, sun_forcing, shade_forcing, limits;
    internal_conduction, body_absorptivity=params.alpha, emissivity=0.95,
    sky_view_factor=0.4, ground_view_factor=0.4, metabolic_heat_volumetric=params.q * u"W/m^3",
)

julia_max_Tb = maximum(ustrip.(u"°C", result.core_temperature))
julia_min_Tb = minimum(ustrip.(u"°C", result.core_temperature))

@testset "trans_behav: R comparison" begin
    # both models cap the thermoregulating animal close to T_F_max, and both show the
    # non-trivial gap between thermoregulating and non-thermoregulating extremes that is
    # the whole point of the behavioral loop — loose tolerances since the two use
    # different simplified state machines and different ODE solvers.
    @test julia_max_Tb ≈ sum_stats.max_Tb atol = 2.0
    @test isfinite(julia_min_Tb)
    @test julia_max_Tb <= params.CT_max
    @test any(s -> s isa Active, result.state)
    @test any(s -> s isa Resting, result.state)
end

# Uncomment to plot the Julia trajectory against R's own day_results output:
# using Plots
# plot(uconvert.(u"hr", day_results.time .* u"s"), day_results.Tb .* u"°C"; label="Tb (R)", color=:black)
# plot!(uconvert.(u"hr", result.t), uconvert.(u"°C", result.core_temperature); label="Tb (Julia)", color=:orange)
# plot!(uconvert.(u"hr", day_results.time .* u"s"), day_results.Tb_open .* u"°C"; label="Tb_open (R)", color=:grey, linestyle=:dash)
# plot!(uconvert.(u"hr", day_results.time .* u"s"), day_results.T_air_shd .* u"°C"; label="local air temp (shade, R)", color=:blue)
# hline!([params.T_F_max * u"°C"]; label="T_F_max", color=:red, linestyle=:dash)
# hline!([params.T_F_min * u"°C"]; label="T_F_min", color=:lightblue, linestyle=:dash)
