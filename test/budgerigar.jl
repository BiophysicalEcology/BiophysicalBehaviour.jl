using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using FluidProperties
using ConstructionBase
using ModelParameters
using Unitful, UnitfulMoles
using Test
using DataFrames, CSV
using Plots

testdir = realpath(joinpath(dirname(pathof(BiophysicalBehaviour)), "../test"))

# budgerigar observations
Weathers1976Fig1 = DataFrame(CSV.File("$testdir/data/budgerigar/Weathers1976Fig1.csv"))
Weathers1976Fig2 = DataFrame(CSV.File("$testdir/data/budgerigar/Weathers1976Fig2.csv"))
Weathers1976Fig3 = DataFrame(CSV.File("$testdir/data/budgerigar/Weathers1976Fig3.csv"))
Weathers1976Fig5 = DataFrame(CSV.File("$testdir/data/budgerigar/Weathers1976Fig5.csv"))

# budgerigar NicheMapR endoR predictions
treg = DataFrame(CSV.File("$testdir/data/budgerigar/treg.csv"))[:, 2:end]
morph = DataFrame(CSV.File("$testdir/data/budgerigar/morph.csv"))[:, 2:end]
enbal = DataFrame(CSV.File("$testdir/data/budgerigar/enbal.csv"))[:, 2:end]
masbal = DataFrame(CSV.File("$testdir/data/budgerigar/masbal.csv"))[:, 2:end]

# budgerigar parameters
shape_pars = example_ellipsoid_shape_pars(; mass = 33.7u"g", shape_b = 1.1, shape_c = 1.1)
insulation_pars = example_insulation_pars(; 
                    fibre_diameter_dorsal = 30.0u"μm",
                    fibre_diameter_ventral = 30.0u"μm",
                    fibre_length_dorsal = 23.1u"mm",
                    fibre_length_ventral = 22.7u"mm",
                    insulation_depth_dorsal = 5.9u"mm",
                    insulation_depth_ventral = 5.7u"mm",
                    insulation_depth_compressed = 5.7u"mm",
                    fibre_density_dorsal = 5000e+04u"1/m^2",
                    fibre_density_ventral = 5000e+04u"1/m^2",
                    insulation_reflectance_dorsal=0.248,
                    insulation_reflectance_ventral=0.351,
                    )
radiation_pars = example_radiation_pars()
Q_metabolism = metabolic_rate(McKechnieWolf(), shape_pars.mass)
metabolism_pars = example_metabolism_pars(; T_core = (38.0 + 273.15)u"K", q10 = 2, Q_metabolism)
respiration_pars = example_respiration_pars(; fO2_extract=0.25, Δ_breath=5.0u"K")
evaporation_pars = example_evaporation_pars(; skin_wetness = 0.005)

# set up geometry
conduction_pars_internal = example_conduction_pars_internal()
fat = Fat(conduction_pars_internal.fat_fraction, conduction_pars_internal.ρ_fat)
mean_insulation_depth = insulation_pars.insulation_depth_dorsal * (1 - radiation_pars.ventral_fraction) + 
    insulation_pars.insulation_depth_ventral * radiation_pars.ventral_fraction
mean_fibre_diameter = insulation_pars.fibre_diameter_dorsal * (1 - radiation_pars.ventral_fraction) + 
    insulation_pars.fibre_diameter_ventral * radiation_pars.ventral_fraction
mean_fibre_density = insulation_pars.fibre_density_dorsal * (1 - radiation_pars.ventral_fraction) + 
    insulation_pars.fibre_density_ventral * radiation_pars.ventral_fraction
fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

traits = Traits(
    shape_pars,
    insulation_pars,
    example_conduction_pars_external(),
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    example_hydraulic_pars(),
    respiration_pars,
    metabolism_pars,
)

organism = Organism(geometry, traits)

# environmental conditions

air_temperatures = (collect(0.0:1.0:50).+273.15)u"K"
P_atmos = 101325.0u"Pa"
ρ_vapour = wet_air_properties(40.0u"°C", 0.3, P_atmos).ρ_vap
saturated_ρ_vapours = DataFrame(wet_air_properties.(air_temperatures, 1, P_atmos)).ρ_vap
experimental_relative_humdities = ρ_vapour ./ saturated_ρ_vapours
experimental_relative_humdities[experimental_relative_humdities .> 1.0] .= 1.0
experimental_relative_humdities[air_temperatures .< 30.0u"°C"] .= 0.15

environment_vars = example_environment_vars(;
                    T_air = air_temperatures[1],
                    rh = experimental_relative_humdities[1])
environment_pars = example_environment_pars()
environment = (; environment_pars, environment_vars)

# update q10s
q10s = fill(1.0, length(air_temperatures))
q10s[air_temperatures .> thermoregulation_pars.T_core_max] .= metabolism_pars.q10
metabolism_pars = example_metabolism_pars(; T_core = (38.0 + 273.15)u"K", q10 = q10s[1], Q_metabolism)

# initial conditions
T_skin = metabolism_pars.T_core - 3.0u"K"
T_insulation = environment_vars.T_air
Q_minimum = metabolism_pars.Q_metabolism
Q_gen = 0.0u"W"

thermoregulation_pars = example_endotherm_thermoregulation_pars(;
    thermoregulation_mode = 3,
    tolerance = 0.005,
    max_iterations = 1000,

    Q_minimum,
    Q_minimum_ref = Q_minimum,

    insulation_depth_dorsal = insulation_pars.insulation_depth_dorsal,
    insulation_depth_ventral = insulation_pars.insulation_depth_ventral,
    insulation_depth_dorsal_max = insulation_pars.insulation_depth_dorsal,
    insulation_depth_ventral_max = insulation_pars.insulation_depth_ventral,
    insulation_depth_dorsal_ref = insulation_pars.insulation_depth_dorsal,
    insulation_depth_ventral_ref = insulation_pars.insulation_depth_ventral,
    insulation_step = 0.0,

    shape_b = 1.1,
    shape_b_step = 0.1,
    shape_b_max = 5.0,

    k_flesh = 0.9u"W/m/K",
    k_flesh_step = 0.1u"W/m/K",
    k_flesh_max = 2.8u"W/m/K",

    T_core = metabolism_pars.T_core,
    T_core_step = 0.1u"K",
    T_core_max = (43.0 + 273.15)u"K",
    T_core_ref = metabolism_pars.T_core,

    pant = 1.0,
    pant_step = 0.01,
    pant_max = 15.0,
    pant_cost = 0.0u"W",
    pant_multiplier = 1.0,

    skin_wetness = evaporation_pars.skin_wetness,
    skin_wetness_step = 0.0025,
    skin_wetness_max = 0.05,
)

model_pars = example_model_pars()

endotherm_out = endotherm_thermoregulation_original(
    Q_gen,
    T_skin,
    T_insulation,
    organism,
    thermoregulation_pars,
    environment,
    model_pars
)
thermoregulation = endotherm_out.thermoregulation
morphology = endotherm_out.morphology
energy_fluxes = endotherm_out.energy_fluxes
mass_fluxes = endotherm_out.mass_fluxes

# now run across all temperatures

results = NamedTuple[]

@assert length(air_temperatures) ==
        length(experimental_relative_humdities) ==
        length(q10s)

for (T_air, rh, q10) in zip(
        air_temperatures,
        experimental_relative_humdities,
        q10s,
    )

    # --- Environment ---
    environment_vars = example_environment_vars(;
        T_air,
        rh,
        wind_speed=0.1u"m/s",
        P_atmos=101325.0u"Pa",
        zenith_angle=20.0u"°",
        k_substrate=2.79u"W/m/K",
        global_radiation=0.0u"W/m^2",
        diffuse_fraction=0,
        shade=0,
    )

    environment = (; environment_pars, environment_vars)

    # --- Metabolism (Q10 changes here) ---
    metabolism_pars = example_metabolism_pars(
        T_core = (38.0 + 273.15)u"K",
        Q_metabolism = Q_minimum,
        q10 = q10,
    )
    traits = Traits(
        shape_pars,
        insulation_pars,
        example_conduction_pars_external(),
        conduction_pars_internal,
        radiation_pars,
        ConvectionParameters(),
        evaporation_pars,
        example_hydraulic_pars(),
        respiration_pars,
        metabolism_pars,
    )

    organism = Organism(geometry, traits)

    #--- Initial conditions (reset every run!) ---
    T_skin = metabolism_pars.T_core - 3.0u"K"
    T_insulation = environment_vars.T_air
    Q_gen = 0.0u"W"

    # --- Thermoregulation ---
    endotherm_out = endotherm_thermoregulation_original(
        Q_gen,
        T_skin,
        T_insulation,
        organism,
        thermoregulation_pars,
        environment,
        model_pars,
    )

    tr = endotherm_out.thermoregulation
    ef = endotherm_out.energy_fluxes
    mf = endotherm_out.mass_fluxes

    push!(results, (
        T_air = T_air,
        rh = rh,
        q10 = q10,

        Q_gen = ef.Q_gen,
        T_core = tr.T_core,
        T_skin_dorsal = tr.T_skin_dorsal,
        T_skin_ventral = tr.T_skin_ventral,
        T_insulation_dorsal = tr.T_insulation_dorsal,
        T_insulation_ventral = tr.T_insulation_ventral,

        pant = tr.pant,
        skin_wetness = tr.skin_wetness,

        V_air = mf.V_air,
        H2O_total = mf.m_evap,
        H2O_resp = mf.m_resp,
        H2O_cut = mf.m_sweat,
    ))
end

predicted = DataFrame(results)


gr()

p1 = plot(
    air_temperatures, predicted.Q_gen,
    lw = 2,
    xlabel = "air temperature",
    ylabel = "metabolic rate",
    ylim = (0.2, 1.2),
    label = "predicted",
)

plot!(
    p1,
    (enbal.TA .+ 273.15)u"K",
    (enbal.QGEN)u"W",
    lw = 2,    
    color = :red,
    label = "NicheMapR",
)

scatter!(
    p1,
    (Weathers1976Fig1.Tair.+273.15)u"K",
    u"W".(HeatExchange.O2_to_Joules(Typical(), 
        (Weathers1976Fig1.mlO2gh * ustrip(u"g", shape_pars.mass))u"ml/hr", respiration_pars.rq)),
    color = :red,
    ms = 4,
    label = "observed",
)

plot!(
    p1,
    legend = :topright,
)

p2 = plot(
    air_temperatures, predicted.H2O_total,
    lw = 2,
    xlabel = "air temperature",
    ylabel = "water loss",
    ylim = (0, 1.5),
    label = "total (pred)",
)

plot!(
    p2,
    air_temperatures, u"g/hr".(predicted.H2O_resp),
    linestyle = :dash,
    label = "respiratory",
)

plot!(
    p2,
    air_temperatures, predicted.H2O_cut,
    linestyle = :dash,
    color = :blue,
    label = "cutaneous",
)

plot!(
    p2,
    (masbal.TA .+ 273.15)u"K",
    (masbal.H2OResp_g .+ masbal.H2OCut_g)u"g/hr",
    lw = 1,    
    color = :red,
    label = "NicheMapR",
)
plot!(
    p2,
    (masbal.TA .+ 273.15)u"K",
    (masbal.H2OResp_g)u"g/hr",
    lw = 1,    
    color = :red,
    label = "NicheMapR",
)
plot!(
    p2,
    (masbal.TA .+ 273.15)u"K",
    (masbal.H2OCut_g)u"g/hr",
    lw = 1,    
    color = :red,
    label = "NicheMapR",
)
scatter!(
    p2,
    (Weathers1976Fig3.Tair.+273.15)u"K",
    (Weathers1976Fig3.mgH2Ogh)u"mg/g/hr" .* shape_pars.mass,
    color = :red,
    ms = 4,
    label = "total (obs)",
)

plot!(
    p2,
    legend = :topleft,
)

p3 = plot(
    air_temperatures, predicted.T_insulation_dorsal,
    color = :grey,
    lw = 2,
    xlabel = "air temperature",
    ylabel = "feather, skin and core temperature",
    ylim = (10+273.15, 50+273.15),
    label = "feathers dorsal",
)

plot!(p3, air_temperatures, predicted.T_insulation_ventral, color = :grey, linestyle = :dash, label = "feathers ventral")
plot!(p3, air_temperatures, predicted.T_skin_dorsal, color = :orange, label = "skin dorsal")
plot!(p3, air_temperatures, predicted.T_skin_ventral, color = :orange, linestyle = :dash, label = "skin ventral")
plot!(p3, air_temperatures, predicted.T_core, color = :red, lw = 2, label = "core (pred)")

plot!(
    p3,
    (treg.TA .+ 273.15)u"K",
    (treg.TFA_D .+ 273.15)u"K",
    lw = 1,    
    color = :red,
)
plot!(
    p3,
    (treg.TA .+ 273.15)u"K",
    (treg.TFA_V .+ 273.15)u"K",
    lw = 1,    
    color = :red,
)
plot!(
    p3,
    (treg.TA .+ 273.15)u"K",
    (treg.TSKIN_D .+ 273.15)u"K",
    lw = 1,    
    color = :red,
)
plot!(
    p3,
    (treg.TA .+ 273.15)u"K",
    (treg.TSKIN_V .+ 273.15)u"K",
    lw = 1,    
    color = :red,
)
plot!(
    p3,
    (treg.TA .+ 273.15)u"K",
    (treg.TC .+ 273.15)u"K",
    lw = 1,    
    color = :red,
)
scatter!(
    p3,
    (Weathers1976Fig3.Tair.+273.15)u"K",
    (Weathers1976Fig2.Tb.+273.15)u"K",
    color = :red,
    ms = 4,
    label = "core (obs)",
)

plot!(
    p3,
    legend = :bottomright,
)
plot!(
    p3,
    legend = :none,
)

p4 = plot(
    air_temperatures,
    u"ml/minute".(predicted.V_air),
    lw = 2,
    xlim = (-5, 50) .+ 273.15,
    ylim = (0, 250),
    xlabel = "air temperature",
    ylabel = "air",
    label = "predicted",
)

plot!(
    p4,
    (masbal.TA .+ 273.15)u"K",
    (masbal.AIR_L)u"L/hr",
    lw = 1,    
    color = :red,
    label = "NicheMapR",
)
scatter!(
    p4,
    (Weathers1976Fig5.Tair.+273.15)u"K",
    (Weathers1976Fig5.breaths_min .* (13.2 .* ustrip(u"kg", shape_pars.mass) .^ 1.08) .*
        ((Weathers1976Fig5.Tair .+ 273.15) ./ 273.15))u"ml/minute",
    color = :red,
    ms = 4,
    label = "observed",
)

plot!(
    p4,
    legend = :topleft,
)

plot(p1, p2, p3, p4, layout = (2, 2))
# plot(p1)
# plot(p2)
# plot(p3)
# plot(p4)