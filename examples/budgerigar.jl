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
shape_pars = example_ellipsoid_shape_pars(; mass = 33.7u"g", shape_coefficient_b = 1.1, shape_coefficient_c = 1.1)
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
metabolic_heat_flow = metabolic_rate(McKechnieWolf(), shape_pars.mass)
metabolism_pars = example_metabolism_pars(; core_temperature = (38.0 + 273.15)u"K", q10 = 2, metabolic_heat_flow)
respiration_pars = example_respiration_pars(; oxygen_extraction_efficiency=0.25, exhaled_temperature_offset=5.0u"K")
evaporation_pars = example_evaporation_pars(; skin_wetness = 0.005)

# set up geometry
conduction_pars_internal = example_conduction_pars_internal()
fat = Fat(conduction_pars_internal.fat_fraction, conduction_pars_internal.fat_density)
mean_insulation_depth = insulation_pars.dorsal.depth * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.ventral.depth * radiation_pars.ventral_fraction
mean_fibre_diameter = insulation_pars.dorsal.diameter * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.ventral.diameter * radiation_pars.ventral_fraction
mean_fibre_density = insulation_pars.dorsal.density * (1 - radiation_pars.ventral_fraction) +
    insulation_pars.ventral.density * radiation_pars.ventral_fraction
fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

# environmental conditions
air_temperatures = (collect(0.0:1.0:50).+273.15)u"K"
atmospheric_pressure = 101325.0u"Pa"
ρ_vapour = wet_air_properties(40.0u"°C", 0.3, atmospheric_pressure).vapour_density
saturated_ρ_vapours = DataFrame(wet_air_properties.(air_temperatures, 1.0, atmospheric_pressure)).vapour_density
experimental_relative_humdities = ρ_vapour ./ saturated_ρ_vapours
experimental_relative_humdities[experimental_relative_humdities .> 1.0] .= 1.0
experimental_relative_humdities[air_temperatures .< 30.0u"°C"] .= 0.15

environment_vars = example_environment_vars(;
                    air_temperature = air_temperatures[1],
                    relative_humidity = experimental_relative_humdities[1])
environment_pars = example_environment_pars()

# initial conditions
skin_temperature = metabolism_pars.core_temperature - 3.0u"K"
insulation_temperature = environment_vars.air_temperature
minimum_metabolic_heat_flow = metabolism_pars.metabolic_heat_flow
generated_heat_flow = 0.0u"W"

# Thermoregulation limits
core_temperature_ref = metabolism_pars.core_temperature
core_temperature_max = (43.0 + 273.15)u"K"

# update q10s
q10s = fill(1.0, length(air_temperatures))
q10s[air_temperatures .>= core_temperature_max] .= metabolism_pars.q10

# Helper function to create organism with current parameters
function create_organism(shape_pars, insulation_pars, conduction_pars_internal, radiation_pars,
                         evaporation_pars, respiration_pars, metabolism_pars, geometry)
    physiology_traits = HeatExchangeTraits(
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
        example_metabolic_rate_options(),
    )

    thermoregulation_limits = ThermoregulationLimits(;
        control=RuleBasedSequentialControl(;
            mode=CorePantingSweatingFirst(),
            tolerance=0.005,
            max_iterations=1000,
        ),
        minimum_metabolic_heat_flow_ref=minimum_metabolic_heat_flow,
        insulation=InsulationLimits(;
            dorsal=SteppedParameter(;
                current=insulation_pars.dorsal.depth,
                reference=insulation_pars.dorsal.depth,
                max=insulation_pars.dorsal.depth,
                step=0.0,
            ),
            ventral=SteppedParameter(;
                current=insulation_pars.ventral.depth,
                reference=insulation_pars.ventral.depth,
                max=insulation_pars.ventral.depth,
                step=0.0,
            ),
        ),
        shape_coefficient_b=SteppedParameter(;
            current=1.1,
            max=5.0,
            step=0.1,
        ),
        flesh_conductivity=SteppedParameter(;
            current=0.9u"W/m/K",
            max=2.8u"W/m/K",
            step=0.1u"W/m/K",
        ),
        core_temperature=SteppedParameter(;
            current=core_temperature_ref,
            reference=core_temperature_ref,
            max=core_temperature_max,
            step=0.1u"K",
        ),
        panting=PantingLimits(;
            panting_rate=SteppedParameter(;
                current=1.0,
                max=15.0,
                step=0.01,
            ),
            cost=0.0u"W",
            multiplier=1.0,
            core_temperature_ref=core_temperature_ref,
        ),
        skin_wetness=SteppedParameter(;
            current=evaporation_pars.skin_wetness,
            max=0.05,
            step=0.0025,
        ),
    )

    behavioral_traits = BehavioralTraits(;
        thermoregulation=thermoregulation_limits,
        activity=Diurnal(),
    )
    organism_traits = OrganismTraits(Endotherm(), physiology_traits, behavioral_traits)

    return Organism(geometry, organism_traits)
end

# Initial run
metabolism_pars_init = example_metabolism_pars(; core_temperature = (38.0 + 273.15)u"K", q10 = q10s[1], metabolic_flux)
organism = create_organism(shape_pars, insulation_pars, conduction_pars_internal, radiation_pars,
                           evaporation_pars, respiration_pars, metabolism_pars_init, geometry)
environment = (; environment_pars, environment_vars)

endotherm_out = thermoregulate(
    organism,
    environment,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
thermoregulation = endotherm_out.thermoregulation
morphology = endotherm_out.morphology
energy_flows = endotherm_out.energy_flows
mass_flows = endotherm_out.mass_flows

# now run across all temperatures

results = NamedTuple[]

@assert length(air_temperatures) ==
        length(experimental_relative_humdities) ==
        length(q10s)

for (air_temp, rel_humidity, q10) in zip(
        air_temperatures,
        experimental_relative_humdities,
        q10s,
    )

    # --- Environment ---
    environment_vars = example_environment_vars(;
        air_temperature=air_temp,
        relative_humidity=rel_humidity,
        wind_speed=0.1u"m/s",
        atmospheric_pressure=101325.0u"Pa",
        zenith_angle=20.0u"°",
        substrate_conductivity=2.79u"W/m/K",
        global_radiation=0.0u"W/m^2",
        diffuse_fraction=0,
        shade=0,
    )

    environment = (; environment_pars, environment_vars)

    # --- Metabolism (Q10 changes here) ---
    metabolism_pars = example_metabolism_pars(
        core_temperature = (38.0 + 273.15)u"K",
        metabolic_heat_flow = minimum_metabolic_heat_flow,
        q10 = q10,
    )

    organism = create_organism(shape_pars, insulation_pars, conduction_pars_internal, radiation_pars,
                               evaporation_pars, respiration_pars, metabolism_pars, geometry)

    #--- Initial conditions (reset every run!) ---
    skin_temperature = metabolism_pars.core_temperature - 3.0u"K"
    insulation_temperature = environment_vars.air_temperature
    generated_heat_flow = 0.0u"W"

    # --- Thermoregulation ---
    endotherm_out = thermoregulate(
        organism,
        environment,
        generated_heat_flow,
        skin_temperature,
        insulation_temperature,
    )

    tr = endotherm_out.thermoregulation
    ef = endotherm_out.energy_flows
    mf = endotherm_out.mass_flows

    push!(results, (
        air_temperature = air_temp,
        relative_humidity = rel_humidity,
        q10 = q10,

        generated_heat_flow = ef.generated_heat_flow,
        core_temperature = tr.core_temperature,
        skin_temperature_dorsal = tr.skin_temperature_dorsal,
        skin_temperature_ventral = tr.skin_temperature_ventral,
        insulation_temperature_dorsal = tr.insulation_temperature_dorsal,
        insulation_temperature_ventral = tr.insulation_temperature_ventral,

        pant = tr.pant,
        skin_wetness = tr.skin_wetness,

        V_air = mf.air_flow,
        H2O_total = mf.m_evap,
        H2O_resp = mf.respiration_mass,
        H2O_cut = mf.m_sweat,
    ))
end

predicted = DataFrame(results)


gr()
default(guidefontsize=8, titlefontsize=10)
plot_NicheMapR_output = false

p1 = plot(
    u"°C".(air_temperatures), predicted.generated_heat_flow,
    lw = 2,
    xlabel = "air temperature",
    title = "metabolic rate",
    ylim = (0.2, 1.2),
    label = "predicted",
)

if plot_NicheMapR_output
plot!(
    p1,
    (enbal.TA)u"°C",
    (enbal.QGEN)u"W",
    lw = 2,
    color = :red,
    label = "NicheMapR",
)
end

scatter!(
    p1,
    (Weathers1976Fig1.Tair.+273.15)u"K",
    u"W".(HeatExchange.O2_to_Joules(Typical(),
        (Weathers1976Fig1.mlO2gh * ustrip(u"g", shape_pars.mass))u"ml/hr", respiration_pars.respiratory_quotient)),
    color = :red,
    ms = 4,
    label = "observed",
)

plot!(
    p1,
    legend = :topright,
    legendfontsize=6,
)

p2 = plot(
    u"°C".(air_temperatures), predicted.H2O_total,
    lw = 2,
    xlabel = "air temperature",
    title = "water loss",
    ylim = (0, 1.5),
    label = "total (pred)",
)

plot!(
    p2,
    u"°C".(air_temperatures), u"g/hr".(predicted.H2O_resp),
    linestyle = :dash,
    label = "respiratory",
)

plot!(
    p2,
    u"°C".(air_temperatures), predicted.H2O_cut,
    linestyle = :dash,
    color = :blue,
    label = "cutaneous",
)

if plot_NicheMapR_output
plot!(
    p2,
    (masbal.TA)u"°C",
    (masbal.H2OResp_g .+ masbal.H2OCut_g)u"g/hr",
    lw = 1,
    color = :red,
    label = "NicheMapR",
)

plot!(
    p2,
    (masbal.TA)u"°C",
    (masbal.H2OResp_g)u"g/hr",
    lw = 1,
    color = :red,
    label = "NicheMapR",
)
plot!(
    p2,
    (masbal.TA)u"°C",
    (masbal.H2OCut_g)u"g/hr",
    lw = 1,
    color = :red,
    label = "NicheMapR",
)
end

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
    legendfontsize=6,
)

p3 = plot(
    u"°C".(air_temperatures), u"°C".(predicted.insulation_temperature_dorsal),
    color = :grey,
    lw = 2,
    xlabel = "air temperature",
    title = "feather, skin and core temperature",
    ylim = (10, 50),
    label = "feathers dorsal",
)

plot!(p3, u"°C".(air_temperatures), u"°C".(predicted.insulation_temperature_ventral), color = :grey, linestyle = :dash, label = "feathers ventral")
plot!(p3, u"°C".(air_temperatures), u"°C".(predicted.skin_temperature_dorsal), color = :orange, label = "skin dorsal")
plot!(p3, u"°C".(air_temperatures), u"°C".(predicted.skin_temperature_ventral), color = :orange, linestyle = :dash, label = "skin ventral")
plot!(p3, u"°C".(air_temperatures), u"°C".(predicted.core_temperature), color = :red, lw = 2, label = "core (pred)")

if plot_NicheMapR_output
plot!(
    p3,
    (treg.TA)u"°C",
    (treg.TFA_D)u"°C",
    lw = 1,
    color = :red,
)
plot!(
    p3,
    (treg.TA)u"°C",
    (treg.TFA_V)u"°C",
    lw = 1,
    color = :red,
)
plot!(
    p3,
    (treg.TA)u"°C",
    (treg.TSKIN_D)u"°C",
    lw = 1,
    color = :red,
)
plot!(
    p3,
    (treg.TA)u"°C",
    (treg.TSKIN_V)u"°C",
    lw = 1,
    color = :red,
)
plot!(
    p3,
    (treg.TA)u"°C",
    (treg.TC)u"°C",
    lw = 1,
    color = :red,
)
end

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
    legendfontsize=6,
)
if plot_NicheMapR_output
plot!(
    p3,
    legend = :none,
)
end

p4 = plot(
    u"°C".(air_temperatures),
    u"ml/minute".(predicted.V_air),
    lw = 2,
    xlim = (-5, 50),
    ylim = (0, 250),
    xlabel = "air temperature",
    title = "ventilation rate",
    label = "predicted",
)

if plot_NicheMapR_output
plot!(
    p4,
    (masbal.TA)u"°C",
    (masbal.AIR_L)u"L/hr",
    lw = 1,
    color = :red,
    label = "NicheMapR",
)
end
scatter!(
    p4,
    (Weathers1976Fig5.Tair)u"°C",
    (Weathers1976Fig5.breaths_min .* (13.2 .* ustrip(u"kg", shape_pars.mass) .^ 1.08) .*
        ((Weathers1976Fig5.Tair .+ 273.15) ./ 273.15))u"ml/minute",
    color = :red,
    ms = 4,
    label = "observed",
)

plot!(
    p4,
    legend = :topleft,
    legendfontsize=6,
)

plot(p1, p2, p3, p4, layout = (2, 2))
# plot(p1, legendfontsize=6)
# plot(p2)
# plot(p3)
# plot(p4)
