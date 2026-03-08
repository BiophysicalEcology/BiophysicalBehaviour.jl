# Desert iguana thermoregulation example
#
# Simulates 2 summer days at Alice Springs, Australia (lat -23.7°, lon 133.9°).
# Microclimate is solved at 0% and 90% shade; the ectotherm thermoregulation
# loop then chooses shade, height and burrow depth each hour.

using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using Plots
using DataFrames
using CSV
using Statistics
using Flatten: flatten

# ── Location ──────────────────────────────────────────────────────────────
latitude  = -23.7u"°"
longitude = 133.9u"°"
elevation = 579.0u"m"

# ── Time ──────────────────────────────────────────────────────────────────
days  = [15, 46]              # mid-Jan, mid-Feb (Southern Hemisphere summer)
hours = collect(0.0:1:23.0)
nsteps = length(days) * length(hours)

# ── Shade ─────────────────────────────────────────────────────────────────
minimum_shade = 0.0
maximum_shade = 0.9

# ── Soil depths and atmospheric heights ───────────────────────────────────
depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [1.0, 200.0]u"cm"    # 1 cm node (organism height) + 2 m reference

# ── Terrain ───────────────────────────────────────────────────────────────
micro_terrain = MicroTerrain(;
    elevation,
    roughness_height = 0.004u"m",
    karman_constant  = 0.4,
    dyer_constant    = 16.0,
    viewfactor       = 1.0,
)

solar_terrain = SolarTerrain(;
    latitude,
    longitude,
    elevation,
    slope                = 0.0u"°",
    aspect               = 0.0u"°",
    albedo               = 0.2,
    horizon_angles       = fill(0.0u"°", 24),
    atmospheric_pressure = atmospheric_pressure(elevation),
)

# ── Soil thermal properties ───────────────────────────────────────────────
soil_thermal_model = example_soil_thermal_parameters(;
    soil_bulk_density    = 1.5u"Mg/m^3",
    soil_mineral_density = 2.56u"Mg/m^3",
)

# Note: example_soil_moisture_model expects dimensionless bulk_density and
# mineral_density (in Mg/m³ without units — it appends u"Mg/m^3" internally)
soil_moisture_model = example_soil_moisture_model(depths;
    bulk_density    = 1.5,
    mineral_density = 2.56,
)

# ── Hourly skeleton (pressure supplied; weather interpolated from minmax) ─
environment_hourly = HourlyTimeseries(;
    pressure              = fill(atmospheric_pressure(elevation), nsteps),
    reference_temperature = nothing,
    reference_humidity    = nothing,
    reference_wind_speed  = nothing,
    global_radiation      = nothing,
    cloud_cover           = nothing,
    rainfall              = nothing,
    zenith_angle          = nothing,
    longwave_radiation    = nothing,
)

# ── Monthly min/max weather (warm, dry summer) ────────────────────────────
environment_minmax = MonthlyMinMaxEnvironment(;
    reference_temperature_min = fill(25.0u"°C", length(days)),
    reference_temperature_max = fill(42.0u"°C", length(days)),
    reference_wind_min        = fill(0.5u"m/s", length(days)),
    reference_wind_max        = fill(3.0u"m/s", length(days)),
    reference_humidity_min    = fill(0.10,      length(days)),
    reference_humidity_max    = fill(0.30,      length(days)),
    cloud_min                 = fill(0.0,       length(days)),
    cloud_max                 = fill(0.1,       length(days)),
    minima_times              = [0, 0, 1, 1],
    maxima_times              = [1, 1, 0, 0],
)
deep_soil_temperature = mean(mean([u"K".(environment_minmax.reference_temperature_min), 
                        u"K".(environment_minmax.reference_temperature_max)]))
initial_soil_temperature = fill(deep_soil_temperature, length(depths))

# ── Daily environment helper ──────────────────────────────────────────────
function make_daily_env(shade_frac)
    DailyTimeseries(;
        shade                 = fill(shade_frac,    length(days)),
        soil_wetness          = fill(0.0,           length(days)),
        surface_emissivity    = fill(0.96,          length(days)),
        cloud_emissivity      = fill(0.96,          length(days)),
        rainfall              = fill(0.0u"kg/m^2",  length(days)),
        deep_soil_temperature = fill(deep_soil_temperature,     length(days)),
        leaf_area_index       = fill(0.1,           length(days)),
    )
end

function make_problem(daily_env)
    MicroProblem(;
        latitude,
        days,
        hours,
        depths,
        heights,
        solar_terrain,
        micro_terrain,
        soil_thermal_model,
        soil_moisture_model,
        environment_minmax,
        environment_daily  = daily_env,
        environment_hourly,
        iterate_day             = 3,
        runmoist                = false,
        initial_soil_temperature = nothing,
        initial_soil_moisture   = fill(0.2, length(depths)),
    )
end

# ── Solve microclimate at 0% shade and 90% shade ──────────────────────────
println("Solving microclimate (0% shade)...")
min_result = Microclimate.solve(make_problem(make_daily_env(minimum_shade)))

println("Solving microclimate (90% shade)...")
max_result = Microclimate.solve(make_problem(make_daily_env(maximum_shade)))


# ── Comparison mode ───────────────────────────────────────────────────────
# true  = feed NicheMapR microclimate to Julia thermoregulation
#         (isolates thermoregulation algorithm differences from microclimate differences)
# false = use Julia microclimate (full end-to-end Julia vs NicheMapR comparison)
use_nmr_microclimate = true

# ── Organism ──────────────────────────────────────────────────────────────
# 40 g desert iguana (DesertIguana shape: allometric surface areas, NicheMapR defaults)
organism_traits = example_ectotherm_organism_traits(
    activity_period     = Diurnal(),
    can_climb           = false,
    can_burrow          = true,
    can_seek_shade      = true,
    can_solar_orient    = false,
    can_press_to_ground = true,
    burrow_shaded       = false,
    alpha_min           = 0.9,   # matches NicheMapR alpha_min / alpha_max = 0.85
    alpha_max           = 0.9,
    heat_exchange = example_heat_exchange_traits(;
        #metabolism_pars = example_ectotherm_metabolism_pars(model = nothing),
        #evaporation_pars = example_ectotherm_evaporation_pars(skin_wetness = 0.0, eye_fraction = 0.0),
        radiation_pars = example_ectotherm_radiation_pars(α_body_dorsal = 0.9, α_body_ventral = 0.9),
    )
)
body     = Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked())
organism = Organism(body, organism_traits)

limits   = thermoregulation(organism)
env_pars = example_environment_pars(; elevation)

# ── Available environments ─────────────────────────────────────────────────
if use_nmr_microclimate
    # Read NicheMapR microclimate CSVs (288 rows = 12 months × 24 h).
    # Fields not in the R CSVs (soil_thermal_conductivity, soil_humidity,
    # diffuse_fraction) fall back to the Julia MicroResult.
    r_metout_df   = CSV.read(joinpath(@__DIR__, "..", "test", "data", "ectotherm", "metout.csv"),   DataFrame)
    r_soil_df     = CSV.read(joinpath(@__DIR__, "..", "test", "data", "ectotherm", "soil.csv"),     DataFrame)
    r_shadmet_df  = CSV.read(joinpath(@__DIR__, "..", "test", "data", "ectotherm", "shadmet.csv"),  DataFrame)
    r_shadsoil_df = CSV.read(joinpath(@__DIR__, "..", "test", "data", "ectotherm", "shadsoil.csv"), DataFrame)

    r_met0   = r_metout_df[1:nsteps, :]
    r_soil0  = r_soil_df[1:nsteps, :]
    r_met90  = r_shadmet_df[1:nsteps, :]
    r_soil90 = r_shadsoil_df[1:nsteps, :]

    function build_nmr_micro(met_df, soil_df; fallback)
        n = nrow(met_df)
        soil_cols = ["D0cm", "D2.5cm", "D5cm", "D10cm", "D15cm",
                     "D20cm", "D30cm", "D50cm", "D100cm", "D200cm"]
        soil_T = hcat([(soil_df[!, col] .+ 273.15) .* u"K" for col in soil_cols]...)
        profiles = [(
            air_temperature   = [(met_df.TALOC[i] + 273.15)u"K",
                                 (met_df.TAREF[i] + 273.15)u"K"],
            relative_humidity = [met_df.RHLOC[i] / 100.0,
                                 met_df.RH[i]    / 100.0],
            wind_speed        = [met_df.VLOC[i]  * u"m/s",
                                 met_df.VREF[i]  * u"m/s"],
        ) for i in 1:n]
        (
            solar_radiation           = (zenith_angle = met_df.ZEN .* u"°",),
            pressure                  = fill(atmospheric_pressure(elevation), n),
            profile                   = profiles,
            sky_temperature           = (met_df.TSKYC .+ 273.15) .* u"K",
            soil_temperature          = soil_T,
            soil_thermal_conductivity = fallback.soil_thermal_conductivity,
            soil_humidity             = fallback.soil_humidity,
            global_radiation          = met_df.SOLR .* u"W/m^2",
            diffuse_fraction          = fill(0.1, n),  # NicheMapR PDIF default
            reference_temperature     = (met_df.TAREF .+ 273.15) .* u"K",
        )
    end

    nmr_min = build_nmr_micro(r_met0,  r_soil0;  fallback = min_result)
    nmr_max = build_nmr_micro(r_met90, r_soil90; fallback = max_result)
    available_environments = AvailableEnvironments(nmr_min, nmr_max, minimum_shade, maximum_shade, depths, heights)
    run_label = "Julia / NMR μC"
else
    available_environments = AvailableEnvironments(min_result, max_result, minimum_shade, maximum_shade, depths, heights)
    run_label = "Julia / Julia μC"
end

# ── Thermoregulation loop (one call per hour) ─────────────────────────────
results         = NamedTuple[]
prev_depth_node = limits.depth.reference
for step in 1:nsteps
    out = thermoregulate(organism, available_environments, limits, env_pars, step, prev_depth_node)
    prev_depth_node = out.depth_node
    push!(results, out)
end

function flatten_namedtuple(nt, prefix="")
    out = Pair{Symbol,Any}[]
    for (k, v) in pairs(nt)
        name = prefix == "" ? string(k) : string(prefix, "_", k)
        if v isa NamedTuple
            append!(out, flatten_namedtuple(v, name))
        else
            push!(out, Symbol(name) => v)
        end
    end
    return out
end
flat = [NamedTuple(flatten_namedtuple(r)) for r in results]
df = DataFrame(flat)

# ── Extract outputs ───────────────────────────────────────────────────────
T_body = [r.T_core  for r in results]
shade  = [r.shade   for r in results]
height = [r.height  for r in results]
state  = [r.state   for r in results]
orientation = [r.sun_orientation  for r in results]

T_air = [min_result.profile[i].air_temperature[1] for i in 1:nsteps]

# ── Load NicheMapR (R) results for comparison ─────────────────────────────
# R script: test/R/ectotherm_thermoregulate.R writes environ.csv with 12 months
# of hourly output. Rows 1:nsteps correspond to doy 15 (Jan 15) and doy 46
# (Feb 15), matching this Julia simulation exactly.
r_path    = joinpath(@__DIR__, "..", "test", "data", "ectotherm", "environ.csv")
r_environ = CSV.read(r_path, DataFrame)
r         = r_environ[1:nsteps, :]        # first 48 hours only

# Unit conversions to match Julia outputs:
#   TC   (°C)           → already in °C; Julia T_body is in K → convert below
#   SHADE (%)           → divide by 100 to get fraction (like Julia `shade`)
#   DEP  (cm, neg=down) → divide by 100 → m (Julia `height` is neg when burrowed)
#   ACT  (0=Resting, 1=Basking, 2=Active) → Julia `state` is OrganismState
r_T_body_C = r.TC
r_shade    = r.SHADE ./ 100
r_height_m = r.DEP   ./ 100
r_act      = r.ACT   # 0/1/2
julia_act  = [s isa Active ? 2 : s isa Basking ? 1 : 0 for s in state]

T_body_C = ustrip.(u"°C", T_body)

# ── Summary statistics ─────────────────────────────────────────────────────
rmse_T = sqrt(mean((T_body_C .- r_T_body_C).^2))
bias_T = mean(T_body_C .- r_T_body_C)
println("\n── Body temperature ($run_label vs NicheMapR) ──")
println("  $run_label mean  = $(round(mean(T_body_C),    digits=2)) °C")
println("  NicheMapR mean   = $(round(mean(r_T_body_C), digits=2)) °C")
println("  RMSE             = $(round(rmse_T, digits=3)) °C")
println("  Bias (Jl−NMR)    = $(round(bias_T, digits=3)) °C")

state_match = sum(julia_act .== r_act)
println("\n── Activity state agreement (Resting/Basking/Active) ──")
println("  Matching hours = $state_match / $nsteps")
println("  Julia:     Resting=$(sum(julia_act.==0)), Basking=$(sum(julia_act.==1)), Active=$(sum(julia_act.==2))")
println("  NicheMapR: Resting=$(sum(r_act.==0)), Basking=$(sum(r_act.==1)), Active=$(sum(r_act.==2))")

# ── Plot ──────────────────────────────────────────────────────────────────
t         = 1:nsteps
Tpref_min = ustrip(u"°C", limits.T_preferred_min)
Tpref_max = ustrip(u"°C", limits.T_preferred_max)
Tcrit_min = ustrip(u"°C", limits.T_critical_min)
Tcrit_max = ustrip(u"°C", limits.T_critical_max)

p1 = plot(t, T_body_C;
    ylabel = "°C", label = "Tb ($run_label)", lw = 2, color = :red,
    title  = "Body temperature")
plot!(p1, t, r_T_body_C;
    label = "Tb (NicheMapR)", lw = 2, color = :darkred, linestyle = :dash)
plot!(p1, t, u"°C".(T_air);
    label = "T_air (1 cm)", lw = 1, color = :steelblue, linestyle = :dash)
hline!(p1, [Tpref_min, Tpref_max];
    label = "Tpref range", color = :orange, linestyle = :dash)
hline!(p1, [Tcrit_min, Tcrit_max];
    label = "Tcrit range", color = :grey, linestyle = :dot)

p2 = plot(t, shade .* 100;
    ylabel = "%", ylim = (0, 100),
    label = "shade ($run_label)", color = :green, title = "Shade selection")
plot!(p2, t, r_shade .* 100;
    label = "shade (NicheMapR)", color = :darkgreen, linestyle = :dash)
scatter!(p2, t[julia_act.==1], fill(65.0, sum(julia_act.==1));
    ms = 2, color = :orange, label = "basking ($run_label)")
scatter!(p2, t[julia_act.==0], fill(60.0, sum(julia_act.==0));
    ms = 2, color = :blue, label = "resting ($run_label)")
scatter!(p2, t[r_act.==1], fill(45.0, sum(r_act.==1));
    ms = 2, color = :goldenrod, label = "basking (NicheMapR)")
scatter!(p2, t[r_act.==0], fill(40.0, sum(r_act.==0));
    ms = 2, color = :purple, label = "resting (NicheMapR)")

# Height comparison: Julia uses organism height in m (positive=above, neg=burrowed);
# NicheMapR DEP is depth below ground (neg=burrowed, 0=surface, no above-ground height).
p3 = plot(t, ustrip.(u"m", height);
    ylabel = "m", label = "height ($run_label)",
    color = :brown, title = "Height / depth (+ above ground, − burrowed)")
plot!(p3, t, r_height_m;
    label = "depth (NicheMapR)", color = :chocolate, linestyle = :dash)
hline!(p3, [0.0]; color = :black, lw = 1, label = "ground")

plot(p1, p2, p3;
    layout      = (3, 1),
    xlabel      = "time step (h)",
    size        = (900, 800),
    left_margin = 5Plots.mm,
)
