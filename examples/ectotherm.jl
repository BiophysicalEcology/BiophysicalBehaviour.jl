# Desert iguana thermoregulation example
#
# Simulates 2 summer days at Alice Springs, Australia (lat -23.7°, lon 133.9°).
# Microclimate is solved at 0% and 90% shade; the ectotherm thermoregulation
# loop then chooses shade, height and underground depth each hour.

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
# NicheMapR convention: min-shade microclimate run is treated as the 0%-shade
# reference in ABOVEGROUND.f (blend = SHADE/MAXSHD).  Setting min_shade_fraction=0.0
# replicates this — blend_factor = (shade - 0) / (0.9 - 0) = shade / maxshade.
# The organism's behavioural floor (shade_min below) is set to 0.2 = NicheMapR minshade=20.
minimum_shade = 0.0   # AvailableEnvironments fraction (0 = treat min-shade run as 0% reference)
minimum_shade_behaviour = 0.0  # organism's behavioural floor (NicheMapR minshade = 20 %)
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
low_shade_result = Microclimate.solve(make_problem(make_daily_env(minimum_shade)))

println("Solving microclimate (90% shade)...")
high_shade_result = Microclimate.solve(make_problem(make_daily_env(maximum_shade)))


# ── Comparison mode ───────────────────────────────────────────────────────
# true  = feed NicheMapR microclimate to Julia thermoregulation
#         (isolates thermoregulation algorithm differences from microclimate differences)
# false = use Julia microclimate (full end-to-end Julia vs NicheMapR comparison)
use_nmr_microclimate = true

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

# ── Organism ──────────────────────────────────────────────────────────────
# 40 g desert iguana (DesertIguana shape: allometric surface areas, NicheMapR defaults)
organism_traits = example_ectotherm_organism_traits(
    activity_period     = CombinedActivity(Diurnal(), Nocturnal(), Crepuscular()),
    can_climb           = false,
    can_retreat_underground = false,
    can_seek_shade      = false,
    can_solar_orient    = true,
    can_press_to_ground = false,
    underground_shaded  = false,
    shade_min           = minimum_shade_behaviour,  # organism can't be below 20% (NicheMapR minshade=20)
    shade_max           = maximum_shade,
    alpha_min           = 0.9,
    alpha_max           = 0.9,
    heat_exchange = example_ectotherm_heat_exchange_traits(;
        conduction_pars_external = example_ectotherm_conduction_pars_external(conduction_fraction = 0.0), # pct_cond=0 in R test
        evaporation_pars = example_ectotherm_evaporation_pars(eye_fraction = 0.0), # NicheMapR live=2: WEYES=0 (eyes only open in live=1)
        radiation_pars = example_ectotherm_radiation_pars(α_body_dorsal = 0.9, α_body_ventral = 0.9,
        solar_orientation = NormalToSun(), ϵ_body_dorsal = 0.95, ϵ_body_ventral = 0.95),
    )
)
body     = Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked())
organism = Organism(body, organism_traits)

limits   = thermoregulation(organism)
env_pars = example_environment_pars(; elevation)

# ── Available environments ─────────────────────────────────────────────────
if use_nmr_microclimate
    # Read NicheMapR microclimate CSVs (288 rows = 12 months × 24 h).
    # metout/soil = NicheMapR run at minshade=20%; shadmet/shadsoil = at maxshade=90%.
    # min_shade_fraction=0.0 in AvailableEnvironments treats the 20%-shade run as the
    # 0%-shade reference, matching NicheMapR ABOVEGROUND.f blend = SHADE/MAXSHD.
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

    nmr_min = build_nmr_micro(r_met0,  r_soil0;  fallback = low_shade_result)
    nmr_max = build_nmr_micro(r_met90, r_soil90; fallback = high_shade_result)
    available_environments = AvailableEnvironments(nmr_min, nmr_max, minimum_shade, maximum_shade, depths, heights)
    run_label = "Julia / NMR μC"
else
    available_environments = AvailableEnvironments(low_shade_result, high_shade_result, minimum_shade, maximum_shade, depths, heights)
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

T_air = [low_shade_result.profile[i].air_temperature[1] for i in 1:nsteps]

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
#   DEP  (cm, neg=down) → divide by 100 → m; DEP==0 (surface) → 0.01 m to match
#                          Julia heights[1] = 1 cm above ground
#   ACT  (0=Resting, 1=Basking, 2=Active) → Julia `state` is OrganismState
r_T_body_C = r.TC
r_shade    = r.SHADE ./ 100
r_height_m = ifelse.(r.DEP .== 0, ustrip(u"m", heights[1]), r.DEP ./ 100)
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

# ── Heat flux comparison ──────────────────────────────────────────────────────
r_enbal_path = joinpath(@__DIR__, "..", "test", "data", "ectotherm", "enbal.csv")
r_enbal = CSV.read(r_enbal_path, DataFrame)
r_eb    = r_enbal[1:nsteps, :]

# Julia heat fluxes (W). Sign convention — same in both Julia and NicheMapR (FUN.f):
#   gains: Q_solar, Q_ir_in, Q_metab  (positive)
#   losses: Q_ir_out, Q_conv, Q_cond, Q_resp, Q_evap  (positive magnitudes)
#   balance: Q_bal = ENB = QIN − QOUT  (≈ 0 at equilibrium; negative when cooling)
# Q_conv can be negative when air is warmer than skin (convective gain).
# NicheMapR keeps QRESP separate from QEVAP (cutaneous only) — same as Julia.
jl_Q_solar  = [ustrip(u"W", r.ectotherm_out.enbal.Q_solar)  for r in results]
jl_Q_ir_in  = [ustrip(u"W", r.ectotherm_out.enbal.Q_ir_in)  for r in results]
jl_Q_ir_out = [ustrip(u"W", r.ectotherm_out.enbal.Q_ir_out) for r in results]
jl_Q_conv   = [ustrip(u"W", r.ectotherm_out.enbal.Q_conv)   for r in results]
jl_Q_cond   = [ustrip(u"W", r.ectotherm_out.enbal.Q_cond)   for r in results]
jl_Q_metab  = [ustrip(u"W", r.ectotherm_out.enbal.Q_metab)  for r in results]
jl_Q_evap   = [ustrip(u"W", r.ectotherm_out.enbal.Q_evap)   for r in results]
jl_Q_resp   = [ustrip(u"W", r.ectotherm_out.enbal.Q_resp)   for r in results]
jl_Q_bal    = [ustrip(u"W", r.ectotherm_out.enbal.Q_bal)    for r in results]

# (label, Julia W, NicheMapR W) — no sign flips needed; both use same convention
flux_pairs = [
    ("Solar",        jl_Q_solar,   r_eb.QSOL),
    ("IR in",        jl_Q_ir_in,   r_eb.QIRIN),
    ("IR out",       jl_Q_ir_out,  r_eb.QIROUT),
    ("Convection",   jl_Q_conv,    r_eb.QCONV),
    ("Conduction",   jl_Q_cond,    r_eb.QCOND),
    ("Metabolic",    jl_Q_metab,   r_eb.QMET),
    ("Respiratory",  jl_Q_resp,    r_eb.QRESP),
    ("Evap",         jl_Q_evap,    r_eb.QEVAP),
    ("Net balance",  jl_Q_bal,     r_eb.ENB),
]

println("\n── Heat fluxes ($run_label vs NicheMapR) ──────────────────────────")
println("  $(rpad("Component", 14))  $(rpad("Jl mean (W)", 13))  $(rpad("NMR mean (W)", 13))  $(rpad("RMSE (W)", 10))  Bias (W)")
for (name, jl, nmr) in flux_pairs
    rmse_f = sqrt(mean((jl .- nmr).^2))
    bias_f = mean(jl .- nmr)
    println("  $(rpad(name, 14))  $(rpad(string(round(mean(jl),  digits=4)), 13))  " *
            "$(rpad(string(round(mean(nmr), digits=4)), 13))  " *
            "$(rpad(string(round(rmse_f,   digits=4)), 10))  $(round(bias_f, digits=4))")
end

# ── Plot ──────────────────────────────────────────────────────────────────
show_legend = false   # set false to hide all plot legends

t         = 1:nsteps
T_target_min = ustrip(u"°C", limits.T_active_min)
T_target_max = ustrip(u"°C", limits.T_active_max)
Tcrit_min = ustrip(u"°C", limits.T_critical_min)
Tcrit_max = ustrip(u"°C", limits.T_critical_max)
T_emerge = ustrip(u"°C", limits.T_emerge)
T_bask = ustrip(u"°C", limits.T_bask)

lbl(s) = show_legend ? s : ""   # suppress label strings when legend hidden

p1 = plot(t, T_body_C;
    label = lbl("Tb ($run_label)"), lw = 2, color = :red,
    title  = "Body temperature", legend = show_legend)
scatter!(p1, t, r_T_body_C;
    label = lbl("Tb (NicheMapR)"), ms = 2, color = :darkred)
plot!(p1, t, u"°C".(T_air);
    label = lbl("T_air (1 cm)"), lw = 1, color = :steelblue, linestyle = :dash)
hline!(p1, [T_target_min, T_target_max];
    label = lbl("T_target range"), color = :orange, linestyle = :dash)
hline!(p1, [Tcrit_min, Tcrit_max];
    label = lbl("Tcrit range"), color = :grey, linestyle = :dot)
hline!(p1, [T_emerge, T_bask];
    label = lbl("Tcrit range"), color = :gold, linestyle = :dot)

p2 = plot(t, shade .* 100;
    ylabel = "%", ylim = (0, 100),
    label = lbl("shade ($run_label)"), color = :green, title = "Shade selection",
    legend = show_legend)
scatter!(p2, t, r_shade .* 100;
    label = lbl("shade (NicheMapR)"), ms = 2, color = :darkgreen)
scatter!(p2, t[julia_act.==1], fill(65.0, sum(julia_act.==1));
    ms = 3, marker = :xcross, color = :orange, label = lbl("basking ($run_label)"))
scatter!(p2, t[julia_act.==0], fill(60.0, sum(julia_act.==0));
    ms = 3, marker = :xcross, color = :blue, label = lbl("resting ($run_label)"))
scatter!(p2, t[r_act.==1], fill(45.0, sum(r_act.==1));
    ms = 3, marker = :xcross, color = :goldenrod, label = lbl("basking (NicheMapR)"))
scatter!(p2, t[r_act.==0], fill(40.0, sum(r_act.==0));
    ms = 3, marker = :xcross, color = :purple, label = lbl("resting (NicheMapR)"))

# Height comparison: Julia uses organism height in m (positive=above, neg=underground);
# NicheMapR DEP is depth below ground (neg=underground, 0=surface, no above-ground height).
p3 = plot(t, ustrip.(u"m", height);
    ylabel = "m", label = lbl("height ($run_label)"),
    color = :brown, title = "Height / depth (+ above ground, − underground)",
    legend = show_legend)
scatter!(p3, t, r_height_m;
    label = lbl("depth (NicheMapR)"), ms = 2, color = :chocolate)
hline!(p3, [0.0]; color = :black, lw = 1, label = lbl("ground"))

for p in (p1, p2, p3)
    vline!(p, collect(t); color = :grey, lw = 0.5, linestyle = :dash, label = "")
end

plot(p1, p2, p3;
    layout      = (3, 1),
    xlabel      = "time step (h)",
    size        = (900, 800),
    left_margin = 5Plots.mm,
)

# # ── Heat flux plot ────────────────────────────────────────────────────────────
# flux_plots = map(flux_pairs) do (name, jl, nmr)
#     p = plot(t, jl;
#         lw=2, color=:steelblue, label=lbl(run_label),
#         ylabel="W", title=name, legend=show_legend)
#     scatter!(p, t, nmr; ms=2, color=:darkred, label=lbl("NicheMapR"))
#     hline!(p, [0.0]; color=:grey, lw=0.5, label="")
#     vline!(p, collect(t); color=:grey, lw=0.5, linestyle=:dash, label="")
#     p
# end

# plot(flux_plots...;
#     layout      = (5, 2),
#     xlabel      = "time step (h)",
#     size        = (1100, 1200),
#     left_margin = 5Plots.mm,
# )
