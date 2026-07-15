# Transient (thermal-mass) endotherm core temperature: a budgerigar-like bird takes off from
# a warm perch into flight. Metabolic rate jumps to flight level and wind speed jumps to flight
# airspeed simultaneously at takeoff; effectors (insulation, panting, skin wetness, flesh
# conductivity) are held fixed for the flight (see docs/transient_body_temperature.md).
# `simulate_endotherm_onelump` dispatches on `insulation(body(organism))`, same as
# `HeatExchange.endotherm_onelump`.
using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Unitful

# --- budgerigar-like organism (parameters from examples/budgerigar.jl) ---
shape_pars = Ellipsoid(33.7u"g", 1000.0u"kg/m^3", 3.0, 3.0)
insulation_pars = InsulationParameters(;
    dorsal=FibreProperties(; diameter=30.0u"μm", length=23.1u"mm", density=5000e4u"1/m^2", depth=5.9u"mm", reflectance=0.248, conductivity=0.209u"W/m/K"),
    ventral=FibreProperties(; diameter=30.0u"μm", length=22.7u"mm", density=5000e4u"1/m^2", depth=5.7u"mm", reflectance=0.351, conductivity=0.209u"W/m/K"),
    depth_compressed=5.7u"mm", longwave_depth_fraction=1.0,
)
ventral_fraction = 0.5
internal_conduction = InternalConductionParameters(;
    fat_fraction=0.05, flesh_conductivity=0.9u"W/m/K", fat_conductivity=0.23u"W/m/K", fat_density=901.0u"kg/m^3",
)
fat = FatLayer(internal_conduction.fat_fraction, internal_conduction.fat_density)
mean_depth = insulation_pars.dorsal.depth * (1 - ventral_fraction) + insulation_pars.ventral.depth * ventral_fraction
mean_diameter = insulation_pars.dorsal.diameter * (1 - ventral_fraction) + insulation_pars.ventral.diameter * ventral_fraction
mean_density = insulation_pars.dorsal.density * (1 - ventral_fraction) + insulation_pars.ventral.density * ventral_fraction
fur = FibrousLayer(mean_depth, mean_diameter, mean_density)
body = Body(shape_pars, CompositeInsulation(fur, fat))

resting_metabolic_heat_flow = metabolic_rate(McKechnieWolf(), shape_pars.mass)
metabolism_pars = example_metabolism_pars(;
    core_temperature=u"K"(38.0u"°C"), model=McKechnieWolf(), metabolic_heat_flow=resting_metabolic_heat_flow,
)
radiation_pars = example_radiation_pars(; ventral_fraction)
traits = example_heat_exchange_traits(;
    shape_pars, insulation_pars, conduction_pars_internal=internal_conduction, radiation_pars, metabolism_pars,
)
organism = Organism(body, traits)

# --- warm environment, perched (resting) vs airborne (flight) wind speed ---
environment_pars = example_environment_pars()
perch_wind_speed = 0.5u"m/s"
flight_wind_speed = 9.0u"m/s"  # ~ cruising airspeed of a small parrot
air_temperature = u"K"(32.0u"°C")  # warm day

perch_environment_vars = example_environment_vars(; air_temperature, wind_speed=perch_wind_speed)

# resting core temperature before takeoff: steady state at BMR, perched in the shade
resting = solve_metabolic_rate(
    organism, (; environment_pars, environment_vars=perch_environment_vars),
    u"K"(35.0u"°C"), air_temperature,
)
core_temperature_init = resting.thermoregulation.core_temperature
println("resting core temperature before takeoff: ", u"°C"(core_temperature_init))

# --- flight: metabolic rate and wind speed both jump at takeoff (t = 0) ---
flight_metabolic_heat_flow = resting_metabolic_heat_flow * 8.0  # elevated for sustained flight
times = (0:1.0:180)u"minute" .|> u"s"
flight_forcing = EnvironmentForcing(
    [first(times), last(times)],
    HeatExchange.EnvironmentalVarsVec(;
        air_temperature=fill(air_temperature, 2), sky_temperature=fill(air_temperature, 2),
        ground_temperature=fill(air_temperature, 2), substrate_temperature=fill(air_temperature, 2),
        relative_humidity=fill(0.4, 2), wind_speed=fill(flight_wind_speed, 2),
        atmospheric_pressure=fill(101325.0u"Pa", 2), zenith_angle=fill(30.0u"°", 2),
        substrate_conductivity=fill(2.79u"W/m/K", 2), global_radiation=fill(0.0u"W/m^2", 2),
        diffuse_fraction=fill(0.1, 2), shade=fill(0.0, 2),
    ),
);

result = simulate_endotherm_onelump(
    times, core_temperature_init, organism, environment_pars, flight_forcing;
    metabolic_heat_flow=flight_metabolic_heat_flow,
);

println("flight metabolic rate: ", round(u"W", flight_metabolic_heat_flow; digits=2),
    " (", round(flight_metabolic_heat_flow / resting_metabolic_heat_flow; digits=1), "x resting)")
println("core temperature range during flight: ", extrema(u"°C".(result.core_temperature)))

# --- how long can flight be sustained before reaching a critical body temperature? ---
critical_temperature = u"K"(43.0u"°C")  # matches examples/budgerigar.jl's core_temperature_max
crossing = findfirst(>=(critical_temperature), result.core_temperature)
if isnothing(crossing)
    println("critical temperature ($(u"°C"(critical_temperature))) not reached within ",
        round(u"minute", times[end]), " of flight")
else
    println("critical temperature ($(u"°C"(critical_temperature))) reached after ",
        round(u"minute", result.t[crossing]), " of sustained flight")
end

# Uncomment to plot the core temperature trajectory during flight:
# using Plots
# plot(uconvert.(u"minute", result.t), uconvert.(u"°C", result.core_temperature); label="core temperature", color=:red, linewidth=2)
# hline!([uconvert(u"°C", critical_temperature)]; label="critical temperature", color=:red, linestyle=:dash)
# hline!([uconvert(u"°C", core_temperature_init)]; label="resting (pre-flight)", color=:grey, linestyle=:dash)

# ============================================================================
# Diurnal activity/rest cycling with real microclimate forcing and height choice
# ============================================================================
# `simulate_endotherm_activity_cycle` repeats the takeoff scenario above over a full day of
# real NicheMapR forcing (test/R/trans_behav.R's output, already saved for the ectotherm
# examples). Flying exposes the bird to reference-height conditions (TAREF/VREF); resting
# uses local-height conditions (TALOC/VLOC) — the "height choice" `active_forcing`/
# `resting_forcing` are built from independently, matching how the ectotherm behavioral
# driver's sun/shade forcing pair works. Effective wind speed while flying is
# `max(movement_speed, reference-height ambient wind)`, not movement speed alone.
using DataFrames, CSV

datadir = joinpath(dirname(pathof(BiophysicalBehaviour)), "..", "test", "data", "trans_behav")
metout = DataFrame(CSV.File(joinpath(datadir, "trans_behav_metout.csv")))
soil = DataFrame(CSV.File(joinpath(datadir, "trans_behav_soil.csv")))
params = DataFrame(CSV.File(joinpath(datadir, "trans_behav_params.csv")))[1, :]

diurnal_times = metout.TIME .* u"minute" .|> u"s"
movement_speed = 9.0u"m/s"  # ~ cruising airspeed of a small parrot

resting_forcing = EnvironmentForcing(
    diurnal_times,
    HeatExchange.EnvironmentalVarsVec(;
        air_temperature=u"K".(metout.TALOC .* u"°C"), sky_temperature=u"K".(metout.TSKYC .* u"°C"),
        ground_temperature=u"K".(soil.D0cm .* u"°C"), substrate_temperature=u"K".(soil.D0cm .* u"°C"),
        relative_humidity=metout.RHLOC ./ 100, wind_speed=metout.VLOC .* u"m/s",
        atmospheric_pressure=fill(params.press * u"Pa", nrow(metout)), zenith_angle=metout.ZEN .* u"°",
        substrate_conductivity=fill(2.79u"W/m/K", nrow(metout)), global_radiation=metout.SOLR .* u"W/m^2",
        diffuse_fraction=fill(0.1, nrow(metout)), shade=fill(0.0, nrow(metout)),
    ),
)
active_forcing = EnvironmentForcing(
    diurnal_times,
    HeatExchange.EnvironmentalVarsVec(;
        air_temperature=u"K".(metout.TAREF .* u"°C"), sky_temperature=u"K".(metout.TSKYC .* u"°C"),
        ground_temperature=u"K".(soil.D0cm .* u"°C"), substrate_temperature=u"K".(soil.D0cm .* u"°C"),
        relative_humidity=metout.RH ./ 100, wind_speed=max.(ustrip.(u"m/s", movement_speed), metout.VREF) .* u"m/s",
        atmospheric_pressure=fill(params.press * u"Pa", nrow(metout)), zenith_angle=metout.ZEN .* u"°",
        substrate_conductivity=fill(2.79u"W/m/K", nrow(metout)), global_radiation=metout.SOLR .* u"W/m^2",
        diffuse_fraction=fill(0.1, nrow(metout)), shade=fill(0.0, nrow(metout)),
    ),
)

diurnal_core_temperature_init = u"K"(metout.TALOC[1] * u"°C")  # ~ air temperature at midnight
diurnal_result = simulate_endotherm_activity_cycle(
    diurnal_times, diurnal_core_temperature_init, organism, environment_pars, active_forcing, resting_forcing;
    active_metabolic_heat_flow=flight_metabolic_heat_flow,
    active_temperature_max=critical_temperature, resume_temperature=u"K"(39.0u"°C"),
)

function _activity_bouts(t, state)
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

flight_bouts = _activity_bouts(diurnal_result.t, diurnal_result.state)
println("\nflight bouts over the day: ", length(flight_bouts))
println("total flight time: ", round(u"minute", sum(flight_bouts; init=0.0u"s")))
println("longest bout: ", round(u"minute", maximum(flight_bouts; init=0.0u"s")))
println("core temperature range over the day: ", extrema(u"°C".(diurnal_result.core_temperature)))

# Uncomment to plot the diurnal core temperature trajectory:
# using Plots
# plot(uconvert.(u"hr", diurnal_result.t), uconvert.(u"°C", diurnal_result.core_temperature); label="core temperature", color=:red, legend=:topleft)
# plot!(uconvert.(u"hr", diurnal_times), u"°C".(active_forcing.(diurnal_times) .|> e -> e.air_temperature); label="reference-height air temp", color=:blue, linestyle=:dash)
# hline!([uconvert(u"°C", critical_temperature)]; label="active_temperature_max", color=:red, linestyle=:dash)
