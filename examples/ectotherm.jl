# Desert iguana (Dipsosaurus dorsalis) thermoregulation — Palm Springs, California
#
# Mechanistic simulation based on:
#   Porter, W.P., Mitchell, J.W., Beckman, W.A. & DeWitt, C.B. (1973).
#   Behavioral implications of mechanistic ecology. Oecologia, 13, 1–54.
#
# Simulates the middle day of each month at Palm Springs, CA (33.8°N, 116.5°W, 130 m).
# Microclimate is solved at 0% and 90% shade; the thermoregulation loop then chooses
# shade, height above ground (up to 200 cm, simulating climbing into creosote bushes)
# and underground retreat depth each hour.

using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using Plots
using DataFrames
using Statistics

# ── Location ──────────────────────────────────────────────────────────────
latitude  = 33.8u"°"
longitude = -116.5u"°"
elevation = 130.0u"m"

# ── Time ──────────────────────────────────────────────────────────────────
days   = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]  # 15th of each month
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
hours  = collect(0.0:1:23.0)
ndays  = length(days)
nsteps = ndays * 24

# ── Shade ─────────────────────────────────────────────────────────────────
minimum_shade = 0.0
maximum_shade = 0.9

# ── Soil depths and atmospheric heights ───────────────────────────────────
depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [1.0, 50.0, 100.0, 150.0, 200.0]u"cm"   # ground + climbing nodes up to 200 cm (creosote bushes)

# ── Terrain ───────────────────────────────────────────────────────────────
micro_terrain = MicroTerrain(;
    elevation,
    roughness_height = 0.0005u"m",   # z0 = 0.05 cm (Porter et al. 1973, p. 5)
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
    albedo               = 0.33,     # quartz sand solar reflectance (Porter et al. 1973, p. 7–8)
    horizon_angles       = fill(0.0u"°", 24),
    atmospheric_pressure = atmospheric_pressure(elevation),
)

# ── Soil thermal properties ───────────────────────────────────────────────
soil_thermal_model = example_soil_thermal_parameters(;
    soil_bulk_density    = 1.5u"Mg/m^3",
    soil_mineral_density = 2.56u"Mg/m^3",
)

soil_moisture_model = example_soil_moisture_model(depths;
    bulk_density    = 1.5,
    mineral_density = 2.56,
)

# ── Hourly skeleton ────────────────────────────────────────────────────────
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

# ── Monthly climatology — Palm Springs, CA ─────────────────────────────────
# Temperature and wind from CRU historical climate averages (approx. Palm Springs AP);
# humidity and cloud are typical desert values.  All values are for the 15th
# of each calendar month (Jan–Dec). Assuming clear skies.
# humidity: relative humidity (%) ÷ 100 to give a fraction as required by Microclimate
environment_minmax = MonthlyMinMaxEnvironment(;
    reference_temperature_min = [3.7, 6.0, 8.2, 11.2, 15.0, 19.1, 22.6, 22.2, 18.8, 13.5,  7.5,  3.5]u"°C",
    reference_temperature_max = [19.3,21.8,23.6, 27.0, 31.0, 35.9, 39.2, 38.4, 35.4, 30.3, 23.9, 19.7]u"°C",
    reference_wind_min        = [0.24,0.26,0.29, 0.31, 0.31, 0.31, 0.31, 0.29, 0.28, 0.25, 0.24, 0.23]u"m/s",
    reference_wind_max        = [2.40,2.62,2.91, 3.06, 3.06, 3.06, 3.06, 2.91, 2.76, 2.47, 2.40, 2.33]u"m/s",
    reference_humidity_min    = [34.0,33.7,33.4, 31.1, 31.4, 31.1, 33.4, 35.2, 34.2, 33.0, 33.1, 34.0] ./ 100,
    reference_humidity_max    = [95.6,94.3,89.4, 83.3, 82.7, 83.2, 86.2, 89.1, 90.8, 92.2, 94.8, 99.3] ./ 100,
    cloud_min                 = fill(0.0, ndays),
    cloud_max                 = fill(0.0, ndays),
    minima_times              = [0, 0, 1, 1],   # temp/wind min at sunrise; humidity/cloud min at solar noon
    maxima_times              = [1, 1, 0, 0],   # temp/wind max after solar noon; humidity/cloud max at sunrise
)

deep_soil_temperature = mean([mean(u"K".(environment_minmax.reference_temperature_min)),
                               mean(u"K".(environment_minmax.reference_temperature_max))])

function make_daily_env(shade_frac)
    DailyTimeseries(;
        shade                 = fill(shade_frac,    ndays),
        soil_wetness          = fill(0.0,           ndays),
        surface_emissivity    = fill(0.96,          ndays),
        cloud_emissivity      = fill(0.96,          ndays),
        rainfall              = fill(0.0u"kg/m^2",  ndays),
        deep_soil_temperature = fill(deep_soil_temperature, ndays),
        leaf_area_index       = fill(0.1,           ndays),
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
low_shade_result  = Microclimate.solve(make_problem(make_daily_env(minimum_shade)))

println("Solving microclimate (90% shade)...")
high_shade_result = Microclimate.solve(make_problem(make_daily_env(maximum_shade)))

available_environments = AvailableEnvironments(
    low_shade_result, high_shade_result, minimum_shade, maximum_shade, depths, heights
)

# ── Organism ──────────────────────────────────────────────────────────────
# 40 g desert iguana (Dipsosaurus dorsalis) — parameters from Porter et al. (1973)
#   alpha_max = 0.80 (maximum absorption at Tb < T_active_min, p. 10)
#   alpha_min = 0.60 (minimum absorption at Tb > T_active_max, p. 10)
#   T_active_min = 38°C (minimum activity temperature, p. 20)
#   T_active_max = 43°C (maximum voluntary Tb, p. 20)
#   T_target     = 38.5°C (preferred temperature from DeWitt 1967, cited p. 19)
#   can_climb = true (paper explicitly models climbing into creosote bushes to 200 cm)
body            = Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked())
organism_traits = example_ectotherm_organism_traits(
    activity_period         = CombinedActivity(Diurnal(), Crepuscular()), #Diurnal(),
    T_target                = u"K"(38.5u"°C"),
    T_active_min            = u"K"(38.0u"°C"),
    T_active_max            = u"K"(43.0u"°C"),
    T_bask                  = u"K"(34.0u"°C"),
    T_emerge                = u"K"(15.0u"°C"),
    T_critical_min          = u"K"(3.0u"°C"),
    T_critical_max          = u"K"(44.0u"°C"),    
    can_climb               = true,
    can_retreat_underground = true,
    depth_min_underground   = 3,
    burrow_shade_mode       = MinShadeOnly(),
    warm_signal             = 0.0u"K/hr", #TODO better name
    can_seek_shade          = false,
    shade_min               = minimum_shade,
    shade_max               = maximum_shade,
    can_solar_orient        = true,
    can_press_to_ground     = false,    
    can_change_absorptivity = true,
    alpha_min               = 0.6,
    alpha_max               = 0.8,
    alpha_step              = 0.003,    
    can_pant                = false,
    pant_max                = 3.0,
    heat_exchange = example_ectotherm_heat_exchange_traits(;
        conduction_pars_external = example_ectotherm_conduction_pars_external(
            conduction_fraction = 0.0),
        evaporation_pars = example_ectotherm_evaporation_pars(
            eye_fraction = 0.0003, skin_wetness = 0.001),
        radiation_pars = example_ectotherm_radiation_pars(
            α_body_dorsal = 0.8, α_body_ventral = 0.8,
            solar_orientation = Intermediate(),
            ϵ_body_dorsal = 0.95, ϵ_body_ventral = 0.95),
        respiration_pars = example_ectotherm_respiration_pars(mouth_fraction = 0.0),
    )
)
organism = Organism(body, organism_traits)

limits   = thermoregulation(organism)
env_pars = example_environment_pars(; 
    elevation, 
    α_ground = solar_terrain.albedo,
)

# ── Thermoregulation loop (one call per hour) ─────────────────────────────
results         = NamedTuple[]
prev_depth_node = limits.depth.reference # TODO better name, "previous_depth"
activity_today  = false # TODO better name, "activity_commenced"
for step in 1:nsteps
    if (step - 1) % 24 == 0
        activity_today = false
    end
    out = thermoregulate(
        organism, available_environments, limits, env_pars, step, prev_depth_node;
        activity_today
    )
    prev_depth_node = out.depth_node
    activity_today  = activity_today || out.state isa Active || out.state isa Basking
    push!(results, out)
end

# ── Extract outputs ───────────────────────────────────────────────────────
T_body     = [r.T_core         for r in results]
shade      = [r.shade          for r in results]
height     = [r.height         for r in results]
depth_node = [r.depth_node     for r in results]
absorptivity = [r.absorptivity for r in results]
sun_orientation = [r.sun_orientation for r in results]
state      = [r.state          for r in results]
T_air      = [low_shade_result.profile[i].air_temperature[1] for i in 1:nsteps]

# Combined position: climbing = +cm, active at surface = 0, underground = −cm.
# heights[1] is the ground node (1 cm); treat it as 0 so only genuine climbing
# (heights[2]+) appears as positive.
_ground_ht = ustrip(u"cm", heights[1])
pos_cm = [depth_node[i] > 1 ?
    -ustrip(u"cm", depths[depth_node[i]]) :
    (h = ustrip(u"cm", height[i]); h > _ground_ht ? h : 0.0)
    for i in 1:nsteps]

T_body_C = ustrip.(u"°C", T_body)
act      = [s isa Active ? 2 : s isa Basking ? 1 : 0 for s in state]

println("\n── Annual activity summary ──")
println("  Resting=$(sum(act.==0)), Basking=$(sum(act.==1)), Active=$(sum(act.==2))")

# ── Heat fluxes ───────────────────────────────────────────────────────────
jl_Q_solar  = [ustrip(u"W", r.ectotherm_out.enbal.Q_solar)  for r in results]
jl_Q_ir_in  = [ustrip(u"W", r.ectotherm_out.enbal.Q_ir_in)  for r in results]
jl_Q_ir_out = [ustrip(u"W", r.ectotherm_out.enbal.Q_ir_out) for r in results]
jl_Q_conv   = [ustrip(u"W", r.ectotherm_out.enbal.Q_conv)   for r in results]
jl_Q_cond   = [ustrip(u"W", r.ectotherm_out.enbal.Q_cond)   for r in results]
jl_Q_metab  = [ustrip(u"W", r.ectotherm_out.enbal.Q_metab)  for r in results]
jl_Q_evap   = [ustrip(u"W", r.ectotherm_out.enbal.Q_evap)   for r in results]
jl_Q_resp   = [ustrip(u"W", r.ectotherm_out.enbal.Q_resp)   for r in results]
jl_Q_bal    = [ustrip(u"W", r.ectotherm_out.enbal.Q_bal)    for r in results]

# TODO - these data are incorrect
# ── Reference data from Porter et al. (1973) — March and July only ────────
# Approximate body temperatures (°C) digitised from Fig. 8 of Porter et al. (1973).
porter1973_jul_Tb = [   # hours 0–23, July 15
    32.5, 32.0, 31.5, 31.0, 31.0, 31.5,   # 00–05 underground (soil cooling)
    33.0, 37.0,                             # 06–07 emerging / basking
    40.5, 42.0, 42.0,                       # 08–10 active
    41.0, 40.5, 40.0, 39.5, 39.0, 39.5,   # 11–16 underground (midday retreat)
    41.0, 41.0,                             # 17–18 late re-emergence
    38.0, 36.0, 35.0, 34.0, 33.0,         # 19–23 retreat / overnight
]
# TODO where are these data from?
porter1973_mar_Tb = [   # hours 0–23, March 15
    18.0, 17.5, 17.0, 17.0, 17.0, 17.5,   # 00–05 underground
    19.0, 21.0, 26.0,                       # 06–08 pre-emergence soil warming
    34.0, 40.0, 41.5, 42.0, 41.5,          # 09–13 active
    41.0, 40.0, 37.0,                       # 14–16 late active / retreating
    28.0, 22.0, 20.0, 19.0, 18.5, 18.0, 18.0,  # 17–23 underground overnight
]
# TODO check these
# Activity states from paper (1 = basking, 2 = active, 0 = resting/underground)
porter1973_jul_act = [0,0,0,0,0,0, 0,1, 2,2,2, 0,0,0,0,0,0, 2,2, 0,0,0,0,0]
porter1973_mar_act = [0,0,0,0,0,0, 0,0,0, 2,2,2,2,2, 2,2,0, 0,0,0,0,0,0,0]

# ── Split outputs by month ────────────────────────────────────────────────
month_ranges = [(m-1)*24+1 : m*24 for m in 1:ndays]
month_Tb     = [T_body_C[r]              for r in month_ranges]
month_act    = [act[r]                    for r in month_ranges]
month_sh     = [shade[r] .* 100          for r in month_ranges]
month_ht     = [ustrip.(u"m", height[r]) for r in month_ranges]
month_pos    = [pos_cm[r]               for r in month_ranges]
month_Ta     = [ustrip.(u"°C", T_air[r]) for r in month_ranges]

T_target_min = ustrip(u"°C", limits.T_active_min)
T_target_max = ustrip(u"°C", limits.T_active_max)
T_emerge_C   = ustrip(u"°C", limits.T_emerge)
T_bask_C     = ustrip(u"°C", limits.T_bask)

# ── Fig. 1 – Body temperature for all 12 months (4×3 grid) ───────────────
# Porter (1973) Tb overlay shown only for March (month 3) and July (month 7)
panels_Tb = map(1:ndays) do m
    p = plot(hours, month_Tb[m];
        lw = 2, color = :red, label = "",
        title = months[m] * " (day $(days[m]))",
        ylabel = "°C", ylim = (0, 55), titlefontsize = 9)
    plot!(p, hours, month_Ta[m]; lw = 1, color = :steelblue, linestyle = :dash, label = "")
    hline!(p, [T_target_min, T_target_max]; color = :orange, linestyle = :dash, lw = 1, label = "")
    if m == 3   # March — Porter 1973 comparison
        plot!(p, hours, porter1973_mar_Tb; lw = 2, color = :darkred, linestyle = :dash, label = "")
    elseif m == 7   # July — Porter 1973 comparison
        plot!(p, hours, porter1973_jul_Tb; lw = 2, color = :darkred, linestyle = :dash, label = "")
    end
    p
end

display(plot(panels_Tb...;
    layout      = (4, 3),
    xlabel      = "hour",
    size        = (1200, 1000),
    left_margin = 4Plots.mm,
    plot_title  = "Body temperature — Dipsosaurus dorsalis, Palm Springs\n(red solid=Julia, dark-red dashed=Porter 1973 [Mar/Jul only], blue dashed=T_air)",
))

# ── Fig. 2 – Annual heatmaps 2×2 (body temp, activity, shade, position) ──
act_matrix = zeros(Int,     24, ndays)
sh_matrix  = zeros(Float64, 24, ndays)
pos_matrix = zeros(Float64, 24, ndays)
tb_matrix  = zeros(Float64, 24, ndays)
for m in 1:ndays
    act_matrix[:, m] = month_act[m]
    sh_matrix[:, m]  = month_sh[m]
    pos_matrix[:, m] = month_pos[m]
    tb_matrix[:, m]  = month_Tb[m]
end

p_tb_hm = heatmap(months, hours, tb_matrix;
    color = cgrad(:RdYlBu, rev=true), clims = (0, 50),
    colorbar_title = "°C",
    title = "Body temperature (°C)", ylabel = "hour", xlabel = "month")

p_act_hm = heatmap(months, hours, act_matrix;
    color = cgrad([:steelblue, :orange, :firebrick], [0, 0.5, 1]),
    clims = (0, 2), colorbar_title = "0=rest 1=bask 2=active",
    title = "Activity state", ylabel = "hour", xlabel = "month",
    yflip = false)

p_sh_hm = heatmap(months, hours, sh_matrix;
    color = :Greens, clims = (0, 100),
    colorbar_title = "shade %",
    title = "Shade used (%)", ylabel = "hour", xlabel = "month")

# zero (surface) = normalised 0.200; next underground node = −2.5 cm = 0.190;
# next climbing node = +50 cm = 0.400.  Narrow green band [0.198, 0.202] catches
# only exactly 0 (after pos_cm maps ground node → 0).
pos_cmap = cgrad(
    [:saddlebrown, :tan, :limegreen, :limegreen, :skyblue, :steelblue],
    [0.0, 0.185, 0.198, 0.202, 0.22, 1.0],
)
p_ht_hm = heatmap(months, hours, pos_matrix;
    color = pos_cmap, clims = (-50, 200),
    colorbar_title = "cm (+ = above, − = below)",
    title = "Position (above ground +cm / underground −cm)", ylabel = "hour", xlabel = "month")

display(plot(p_tb_hm, p_act_hm, p_sh_hm, p_ht_hm;
    layout      = (2, 2),
    size        = (1200, 800),
    left_margin = 6Plots.mm,
))

# ── Fig. 3 – Monthly activity budget (Porter et al. 1973 Fig. 11 style) ───
using StatsPlots

n_rest_v   = [sum(month_act[m] .== 0) for m in 1:ndays]
n_bask_v   = [sum(month_act[m] .== 1) for m in 1:ndays]
n_active_v = [sum(month_act[m] .== 2) for m in 1:ndays]

display(groupedbar(hcat(n_rest_v, n_bask_v, n_active_v);
    bar_position = :stack,
    xticks       = (1:ndays, months),
    label        = ["rest" "bask" "active"],
    color        = [:steelblue :orange :firebrick],
    title        = "Monthly activity budget — Dipsosaurus dorsalis, Palm Springs",
    ylabel       = "hours per day",
    ylim         = (0, 24),
    size         = (800, 400),
    left_margin  = 5Plots.mm,
))

# ── Fig. 4 – Position per month (4×3 grid) ────────────────────────────────
# Positive = height above ground (cm); negative = depth below ground (cm).
# Ground level (0) shown as a dashed reference line.
depth_cm_max = 50.0
height_cm_max = maximum(ustrip.(u"cm", heights))

panels_pos = map(1:ndays) do m
    p = plot(hours, month_pos[m];
        lw = 2, color = :sienna, label = "",
        title = months[m] * " (day $(days[m]))",
        ylabel = "cm", ylim = (-depth_cm_max, height_cm_max),
        titlefontsize = 9)
    hline!(p, [0.0]; color = :black, lw = 1, linestyle = :dash, label = "")
    hspan!(p, [-depth_cm_max, 0]; fillalpha = 0.06, color = :brown, label = "")
    hspan!(p, [0, height_cm_max]; fillalpha = 0.06, color = :skyblue, label = "")
    p
end

display(plot(panels_pos...;
    layout      = (4, 3),
    xlabel      = "hour",
    size        = (1200, 1000),
    left_margin = 4Plots.mm,
    plot_title  = "Position — above ground (+cm, blue) / underground (−cm, brown)",
))
