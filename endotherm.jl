using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using ConstructionBase
using ModelParameters
using Unitful, UnitfulMoles
using Test
using DataFrames, CSV

testdir = realpath(joinpath(dirname(pathof(BiophysicalBehaviour)), "../test"))

endo_input_names = Symbol.(DataFrame(CSV.File("$testdir/data/endoR_input_names.csv"))[:, 2])
furmult = 1
shape = 4
endo_input_vec = DataFrame(CSV.File("$testdir/data/endoR_input_" * string(shape) * "_" * string(furmult) * ".csv"))[:, 2]
treg_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_treg_" * string(shape) * "_" * string(furmult) * ".csv"))[:, 2:end]))
morph_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_morph_" * string(shape) * "_" * string(furmult) * ".csv"))[:, 2:end]))
enbal_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_enbal_" * string(shape) * "_" * string(furmult) * ".csv"))[:, 2:end]))
masbal_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_masbal_" * string(shape) * "_" * string(furmult) * ".csv"))[:, 2:end]))

endo_input = (; zip(endo_input_names, endo_input_vec)...)

# define shape
if endo_input.SHAPE == 1
    shape_pars = Cylinder((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3",
        (endo_input.SHAPE_B)) # define shape
end
if endo_input.SHAPE == 2
    shape_pars = Sphere((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3")
end
if endo_input.SHAPE == 3
    shape_pars = Plate((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3",
        (endo_input.SHAPE_B), (endo_input.SHAPE_C))
end
if endo_input.SHAPE == 4
    shape_pars = Ellipsoid((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3",
        (endo_input.SHAPE_B), (endo_input.SHAPE_C))
end

fat = Fat(endo_input.FATPCT / 100.0, (endo_input.FATDEN)u"kg/m^3")
mean_insulation_depth = (endo_input.ZFURD * (1 - endo_input.PVEN) + endo_input.ZFURV * endo_input.PVEN)u"m"
mean_fibre_diameter = (endo_input.DHAIRD * (1 - endo_input.PVEN) + endo_input.DHAIRV * endo_input.PVEN)u"m"
mean_fibre_density = (endo_input.RHOD * (1 - endo_input.PVEN) + endo_input.RHOV * endo_input.PVEN)u"1/m^2"
fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

environment_vars = EnvironmentalVars(;
    T_air=u"K"((endo_input.TA)u"°C"),
    T_air_reference=u"K"((endo_input.TAREF)u"°C"),
    T_sky=u"K"((endo_input.TSKY)u"°C"),
    T_ground=u"K"((endo_input.TGRD)u"°C"),
    T_substrate=u"K"((endo_input.TCONDSB)u"°C"),
    T_bush=u"K"((endo_input.TBUSH)u"°C"),
    T_vegetation=u"K"((endo_input.TAREF)u"°C"),
    rh=endo_input.RH / 100.0,
    wind_speed=(endo_input.VEL)u"m/s",
    P_atmos=(endo_input.BP)u"Pa",
    zenith_angle=(endo_input.Z)u"°",
    k_substrate=(endo_input.KSUB)u"W/m/K",
    global_radiation=(endo_input.QSOLR)u"W/m^2",
    diffuse_fraction=endo_input.PDIF,
    shade=endo_input.SHADE / 100,
)

environment_pars = EnvironmentalPars(;
    α_ground=endo_input.ABSSB,
    ϵ_ground=1.0,
    ϵ_sky=1.0,
    elevation=(endo_input.ELEV)u"m",
    fluid=endo_input.FLTYPE,
    fN2=endo_input.N2GAS / 100.0,
    fO2=endo_input.O2GAS / 100.0,
    fCO2=endo_input.CO2GAS / 100.0,
    convection_enhancement=endo_input.CONV_ENHANCE,
)

conduction_pars_external = ExternalConductionParameters(;
    conduction_fraction=endo_input.PCOND,
)

conduction_pars_internal = InternalConductionParameters(;
    fat_fraction=endo_input.FATPCT / 100.0,
    k_flesh=(endo_input.AK1)u"W/m/K",
    k_fat=(endo_input.AK2)u"W/m/K",
    ρ_fat=(endo_input.FATDEN)u"kg/m^3",
)

if endo_input.ORIENT == 0.0
    solar_orientation = Intermediate()
elseif endo_input.ORIENT == 1.0
    solar_orientation = NormalToSun()
else
    solar_orientation = ParallelToSun()
end
radiation_pars = RadiationParameters(;
    α_body_dorsal=(1 - endo_input.REFLD),
    α_body_ventral=(1 - endo_input.REFLV),
    ϵ_body_dorsal=endo_input.EMISAN,
    ϵ_body_ventral=endo_input.EMISAN,
    F_sky=endo_input.FSKREF,
    F_ground=endo_input.FGDREF,
    F_bush=endo_input.FABUSH,
    ventral_fraction=endo_input.PVEN,
    solar_orientation,
)

evaporation_pars = EvaporationParameters(;
    skin_wetness=endo_input.PCTWET / 100.0,
    insulation_wetness=endo_input.FURWET / 100.0,
    eye_fraction=endo_input.PCTEYES / 100.0,
    bare_skin_fraction=endo_input.PCTBAREVAP / 100.0,
    insulation_fraction=1.0,
)

hydraulic_pars = HydraulicParameters(;
)

respiration_pars = RespirationParameters(;
    fO2_extract=endo_input.EXTREF / 100.0,
    pant=endo_input.PANT,
    rq=endo_input.RQ,
    Δ_breath=(endo_input.DELTAR)u"K",
    rh_exit=endo_input.RELXIT / 100.0,
)

metabolism_pars = MetabolismParameters(;
    T_core=u"K"((endo_input.TC)u"°C"),
    Q_metabolism=(endo_input.QBASAL)u"W",
    q10=endo_input.Q10,
    model=Kleiber(),
)

if endo_input.FURTHRMK == 0.0
    insulation_conductivity = nothing
else
    insulation_conductivity = (endo_input.FURTHRMK)u"W/m/K"
end
insulation_pars = InsulationParameters(;
    insulation_conductivity_dorsal=insulation_conductivity,
    insulation_conductivity_ventral=insulation_conductivity,
    fibre_diameter_dorsal=(endo_input.DHAIRD)u"m",
    fibre_diameter_ventral=(endo_input.DHAIRV)u"m",
    fibre_length_dorsal=(endo_input.LHAIRD)u"m",
    fibre_length_ventral=(endo_input.LHAIRV)u"m",
    insulation_depth_dorsal=(endo_input.ZFURD)u"m",
    insulation_depth_ventral=(endo_input.ZFURV)u"m",
    fibre_density_dorsal=(endo_input.RHOD)u"1/m^2",
    fibre_density_ventral=(endo_input.RHOV)u"1/m^2",
    insulation_reflectance_dorsal=endo_input.REFLD,
    insulation_reflectance_ventral=endo_input.REFLV,
    insulation_depth_compressed=(endo_input.ZFURCOMP)u"m",
    fibre_conductivity=(endo_input.KHAIR)u"W/m/K",
    longwave_depth_fraction=endo_input.XR,
)
traits = Traits(
    insulation_pars,
    conduction_pars_external,
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    hydraulic_pars,
    respiration_pars,
    metabolism_pars
)

mammal = Organism(geometry, traits)
environment = (; environment_pars, environment_vars)

model_pars = EndoModelPars(
    respire=Bool(endo_input.RESPIRE),
    simulsol_tolerance=(endo_input.DIFTOL)u"K",
    resp_tolerance=endo_input.BRENTOL,
)

# initial conditions
T_skin = u"K"((endo_input.TS)u"°C")
T_insulation = u"K"((endo_input.TFA)u"°C")
Q_minimum = (endo_input.QBASAL)u"W"
Q_gen = 0.0u"W"
Q_minimum_ref = Q_minimum

T_core_step = (endo_input.TC_INC)u"K"
T_core_max = u"K"((endo_input.TC_MAX)u"°C")
T_core_ref = metabolism_pars.T_core
T_core = T_core_ref
q10 = metabolism_pars.q10

insulation_depth_dorsal_max = (endo_input.ZFURD_MAX)u"m"
insulation_depth_ventral_max = (endo_input.ZFURV_MAX)u"m"
insulation_depth_dorsal_ref = insulation_pars.insulation_depth_dorsal
insulation_depth_ventral_ref = insulation_pars.insulation_depth_ventral
fibre_length_dorsal = insulation_pars.fibre_length_dorsal
fibre_length_ventral = insulation_pars.fibre_length_ventral
fibre_diameter_dorsal = insulation_pars.fibre_diameter_dorsal
fibre_diameter_ventral = insulation_pars.fibre_diameter_ventral
fibre_density_dorsal = insulation_pars.fibre_density_dorsal
fibre_density_ventral = insulation_pars.fibre_density_ventral
insulation_step = endo_input.PZFUR

k_flesh = (endo_input.AK1)u"W/m/K"
k_flesh_step = (endo_input.AK1_INC)u"W/m/K"
k_flesh_max = (endo_input.AK1_MAX)u"W/m/K"
shape_b = endo_input.SHAPE_B
shape_b_step = endo_input.UNCURL
shape_b_max = endo_input.SHAPE_B_MAX

pant = respiration_pars.pant
pant_step = endo_input.PANT_INC
pant_max = endo_input.PANT_MAX
pant_multiplier = endo_input.PANT_MULT
pant_cost = 0.0u"W" # initialise

skin_wetness = evaporation_pars.skin_wetness
skin_wetness_step = endo_input.PCTWET_INC
skin_wetness_max = endo_input.PCTWET_MAX

thermoregulation_mode = endo_input.TREGMODE
thermoregulate = true

# start thermoregulation conduction_pars_external

if insulation_step > 0.0 # start with erect insulation
    insulation_depth_dorsal = insulation_depth_dorsal_max
    insulation_depth_ventral = insulation_depth_ventral_max
    insulation_pars = ConstructionBase.setproperties(
        insulation_pars,
        insulation_depth_dorsal=insulation_depth_dorsal,
        insulation_depth_ventral=insulation_depth_ventral,
    )
    mean_insulation_depth = insulation_depth_dorsal * (1 - radiation_pars.ventral_fraction) +
                            insulation_depth_ventral * radiation_pars.ventral_fraction
    mean_fibre_diameter = fibre_diameter_dorsal * (1 - radiation_pars.ventral_fraction) +
                          fibre_diameter_ventral * radiation_pars.ventral_fraction
    mean_fibre_density = fibre_density_dorsal * (1 - radiation_pars.ventral_fraction) +
                         fibre_density_ventral * radiation_pars.ventral_fraction
    fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
    geometry = Body(shape_pars, CompositeInsulation(fur, fat))
    traits = Traits(
        insulation_pars,
        conduction_pars_external,
        conduction_pars_internal,
        radiation_pars,
        ConvectionParameters(),
        evaporation_pars,
        hydraulic_pars,
        respiration_pars,
        metabolism_pars
    )
    mammal = Organism(geometry, traits)
end

endotherm_out = solve_metabolic_rate(T_skin, T_insulation, mammal, environment, model_pars)

T_skin = endotherm_out.thermoregulation.T_skin
T_insulation = endotherm_out.thermoregulation.T_insulation
Q_gen = endotherm_out.energy_fluxes.Q_gen

@time while Q_gen < Q_minimum * 0.995

    if thermoregulate
        if (insulation_depth_dorsal > insulation_depth_dorsal_ref) &&
           (insulation_depth_ventral > insulation_depth_ventral_ref)
            insulation_depth_dorsal = max(insulation_depth_dorsal_ref,
                insulation_depth_dorsal - insulation_step * fibre_length_dorsal)
            insulation_depth_ventral = max(insulation_depth_ventral_ref,
                insulation_depth_ventral - insulation_step * fibre_length_ventral)
            insulation_pars = ConstructionBase.setproperties(
                insulation_pars,
                insulation_depth_dorsal = insulation_depth_dorsal,
                insulation_depth_ventral = insulation_depth_ventral,
            )                
            mean_insulation_depth = insulation_depth_dorsal * (1 - radiation_pars.ventral_fraction) +
                                    insulation_depth_ventral * radiation_pars.ventral_fraction
            mean_fibre_diameter = fibre_diameter_dorsal * (1 - radiation_pars.ventral_fraction) +
                                  fibre_diameter_ventral * radiation_pars.ventral_fraction
            mean_fibre_density = fibre_density_dorsal * (1 - radiation_pars.ventral_fraction) +
                                 fibre_density_ventral * radiation_pars.ventral_fraction
            fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
            geometry = Body(shape_pars, CompositeInsulation(fur, fat))
        else
            insulation_depth_dorsal = insulation_depth_dorsal_ref
            insulation_depth_ventral = insulation_depth_ventral_ref
            insulation_pars = ConstructionBase.setproperties(
                insulation_pars,
                insulation_depth_dorsal = insulation_depth_dorsal,
                insulation_depth_ventral = insulation_depth_ventral,
            )
            mean_insulation_depth = insulation_depth_dorsal * (1 - radiation_pars.ventral_fraction) +
                                    insulation_depth_ventral * radiation_pars.ventral_fraction
            mean_fibre_diameter = fibre_diameter_dorsal * (1 - radiation_pars.ventral_fraction) +
                                  fibre_diameter_ventral * radiation_pars.ventral_fraction
            mean_fibre_density = fibre_density_dorsal * (1 - radiation_pars.ventral_fraction) +
                                 fibre_density_ventral * radiation_pars.ventral_fraction
            fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
            geometry = Body(shape_pars, CompositeInsulation(fur, fat))
            if shape_b < shape_b_max
                shape_b += shape_b_step
                shape_pars = ConstructionBase.setproperties(
                    shape_pars, b = shape_b, c = shape_b
                )
                geometry = Body(shape_pars, CompositeInsulation(fur, fat))
            else
                shape_b = shape_b_max
                shape_pars = ConstructionBase.setproperties(
                    shape_pars, b = shape_b, c = shape_b
                )
                geometry = Body(shape_pars, CompositeInsulation(fur, fat))
                if k_flesh < k_flesh_max
                    k_flesh += k_flesh_step
                    conduction_pars_internal = ConstructionBase.setproperties(
                        conduction_pars_internal, k_flesh=k_flesh
                    )
                else
                    k_flesh = k_flesh_max
                    conduction_pars_internal = ConstructionBase.setproperties(
                        conduction_pars_internal, k_flesh=k_flesh
                    )

                    if T_core < T_core_max
                        T_core += T_core_step
                        metabolism_pars = ConstructionBase.setproperties(
                            metabolism_pars, T_core=T_core
                        )
                        q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                        if (thermoregulation_mode >= 2) && (pant < pant_max)
                            pant += pant_step
                            respiration_pars = ConstructionBase.setproperties(
                                respiration_pars, pant=pant
                            )
                            pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                        (pant_multiplier - 1) * Q_minimum_ref

                            if thermoregulation_mode == 3
                                skin_wetness += skin_wetness_step
                                evaporation_pars = ConstructionBase.setproperties(
                                    evaporation_pars, skin_wetness=skin_wetness
                                )
                                if skin_wetness > skin_wetness_max
                                    skin_wetness = skin_wetness_max
                                    evaporation_pars = ConstructionBase.setproperties(
                                        evaporation_pars, skin_wetness=skin_wetness
                                    )
                                end
                            end
                        end

                        Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                    else
                        T_core = T_core_max
                        metabolism_pars = ConstructionBase.setproperties(
                            metabolism_pars, T_core=T_core
                        )
                        q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                        if pant < pant_max
                            pant += pant_step
                            respiration_pars = ConstructionBase.setproperties(
                                respiration_pars, pant=pant
                            )
                            pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                        (pant_multiplier - 1) * Q_minimum_ref
                            Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                            if thermoregulation_mode == 3
                                skin_wetness += skin_wetness_step
                                evaporation_pars = ConstructionBase.setproperties(
                                    evaporation_pars, skin_wetness=skin_wetness
                                )
                                if skin_wetness > skin_wetness_max
                                    skin_wetness = skin_wetness_max
                                    evaporation_pars = ConstructionBase.setproperties(
                                        evaporation_pars, skin_wetness=skin_wetness
                                    )
                                end
                            end
                        else
                            pant = pant_max
                            respiration_pars = ConstructionBase.setproperties(
                                respiration_pars, pant=pant
                            )
                            pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                        (pant_multiplier - 1) * Q_minimum_ref

                            Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                            skin_wetness += skin_wetness_step
                            evaporation_pars = ConstructionBase.setproperties(
                                evaporation_pars, skin_wetness=skin_wetness
                            )
                            if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                                skin_wetness = skin_wetness_max
                                evaporation_pars = ConstructionBase.setproperties(
                                    evaporation_pars, skin_wetness=skin_wetness
                                )
                                return
                            end
                        end
                    end
                end
            end
        end
        traits = Traits(
            insulation_pars,
            conduction_pars_external,
            conduction_pars_internal,
            radiation_pars,
            ConvectionParameters(),
            evaporation_pars,
            hydraulic_pars,
            respiration_pars,
            metabolism_pars
        )
        mammal = Organism(geometry, traits)
        endotherm_out = solve_metabolic_rate(T_skin, T_insulation, mammal, environment, model_pars)

        T_skin = endotherm_out.thermoregulation.T_skin
        T_insulation = endotherm_out.thermoregulation.T_insulation
        Q_gen = endotherm_out.energy_fluxes.Q_gen

    else
        break
    end
end
endotherm_out.morphology
endotherm_out.energy_fluxes
endotherm_out.thermoregulation