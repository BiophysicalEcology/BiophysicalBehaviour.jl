# Diurnal thermoregulatory behavior under a transient heat budget, mimicking trans_behav.R's
# own example (500 g desert ectotherm shuttling between sun and shade), using real
# NicheMapR forcing (test/R/trans_behav.R). `simulate_diurnal_behavior` dispatches on the
# initial state: a plain temperature runs one-lump, a `(; core_temperature,
# shell_temperature)` NamedTuple runs two-lump.
using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful
using DataFrames, CSV

datadir = joinpath(dirname(pathof(BiophysicalBehaviour)), "..", "test", "data", "trans_behav")
params = DataFrame(CSV.File(joinpath(datadir, "trans_behav_params.csv")))[1, :];
metout = DataFrame(CSV.File(joinpath(datadir, "trans_behav_metout.csv")));
shadmet = DataFrame(CSV.File(joinpath(datadir, "trans_behav_shadmet.csv")));
soil = DataFrame(CSV.File(joinpath(datadir, "trans_behav_soil.csv")));
shadsoil = DataFrame(CSV.File(joinpath(datadir, "trans_behav_shadsoil.csv")));

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

sun_forcing = forcing_from(metout, soil, params.press, 0.0);
shade_forcing = forcing_from(shadmet, shadsoil, params.press, params.shade / 100);

# Animal parameters
body_mass = 500.0 * u"g"
body_density = 1000.0 * u"kg/m^3"
axis_ratio_b = 5.0
axis_ratio_c = 5.0
flesh_conductivity = 0.5 * u"W/m/K"
flesh_specific_heat = 3073.0 * u"J/kg/K"

body_absorptivity=0.85
emissivity=0.95
sky_view_factor=0.4
ground_view_factor=0.4
metabolic_heat_volumetric=0.0 * u"W/m^3"

active_temperature_min = u"K"(Float64(33.0) * u"°C")
active_temperature_max = u"K"(Float64(43.0) * u"°C")
basking_temperature_min = u"K"(Float64(18.0) * u"°C")
critical_temperature_max = u"K"(Float64(48.0) * u"°C")

body = Body(Ellipsoid(body_mass, body_density, axis_ratio_b, axis_ratio_c), Naked());
environment_pars = example_environment_pars(; ground_albedo=params.alpha_sub)
internal_conduction = example_conduction_pars_internal(;
    flesh_conductivity, flesh_specific_heat,
    fat_conductivity=flesh_conductivity, fat_specific_heat=flesh_specific_heat,
)
limits = example_ectotherm_behavioral_limits(;
    active_temperature_min, active_temperature_max,
    basking_temperature_min, critical_temperature_max,
)

times = metout.TIME .* u"minute" .|> u"s"
air_temperature = u"K".(metout.TALOC .* u"°C")
core_temperature_init = air_temperature[1]  # ~air temperature at midnight
kw = (;
    internal_conduction, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
)

result = simulate_diurnal_behavior(
    times, core_temperature_init, body, environment_pars, sun_forcing, shade_forcing, limits; kw...,
)

# non-thermoregulating baseline: stays in the open (full sun) all day
open_solution = simulate_ectotherm_onelump(
    times, core_temperature_init, body, environment_pars, sun_forcing; posture=Intermediate(), kw...,
)

# bout statistics: contiguous stretches of Active state
function activity_bouts(t, state)
    bouts = typeof(1.0u"s")[]
    bout_start = 0.0u"s"
    for i in eachindex(state)
        if state[i] isa Active && (i == 1 || !(state[i - 1] isa Active))
            bout_start = t[i]
        elseif state[i] isa Active && (i == lastindex(state) || !(state[i + 1] isa Active))
            push!(bouts, t[i] - bout_start)
        end
    end
    return bouts
end

bouts = activity_bouts(result.t, result.state)
println("number of foraging bouts: ", length(bouts))
println("total activity time: ", round(u"minute", sum(bouts; init=0.0u"s")))
println("longest bout: ", round(u"minute", maximum(bouts; init=0.0u"s")))
println("body temperature range (thermoregulating): ", extrema(u"°C".(result.core_temperature)))
println("body temperature range (non-thermoregulating, open): ", extrema(u"°C".(open_solution.core_temperature)))

# two-lump: same scenario, core and shell tracked separately
shell_thickness = 0.002u"m"
state_init = (core_temperature=core_temperature_init, shell_temperature=core_temperature_init)

result_twolump = simulate_diurnal_behavior(
    times, state_init, body, environment_pars, sun_forcing, shade_forcing, limits; kw..., shell_thickness,
)

println()
println("two-lump core temperature range: ", extrema(u"°C".(result_twolump.core_temperature)))
println("two-lump shell temperature range: ", extrema(u"°C".(result_twolump.shell_temperature)))

# Uncomment to reproduce trans_behav.R's plot (Tb_open in grey, air temp in blue,
# thermoregulating Tb in orange, T_F_min/T_F_max/CT_max threshold lines):
using Plots
plot(uconvert.(u"hr", open_solution.t), uconvert.(u"°C", open_solution.core_temperature); label="Tb, non-thermoregulating", color=:grey, legend=:topleft, ylim=(15, 50), xlabel="time (hr)", ylabel="temperature (°C)")
plot!(uconvert.(u"hr", result.t), uconvert.(u"°C", shade_forcing.(result.t) .|> e -> e.air_temperature); label="local air temp (shade)", color=:blue)
plot!(uconvert.(u"hr", result.t), uconvert.(u"°C", result.core_temperature); label="Tb, thermoregulating", color=:orange, linewidth=2)
hline!([uconvert(u"°C", limits.active_temperature_max)]; label="T_F_max", color=:red, linestyle=:dash)
hline!([uconvert(u"°C", limits.active_temperature_min)]; label="T_F_min", color=:lightblue, linestyle=:dash)
hline!([uconvert(u"°C", limits.critical_temperature_max)]; label="CT_max", color=:red)

# Uncomment for the two-lump core-vs-shell trajectory:
# plot(uconvert.(u"hr", result_twolump.t), uconvert.(u"°C", result_twolump.core_temperature); label="core", color=:orange, legend=:topleft)
# plot!(uconvert.(u"hr", result_twolump.t), uconvert.(u"°C", result_twolump.shell_temperature); label="shell", color=:orange, linestyle=:dash)
