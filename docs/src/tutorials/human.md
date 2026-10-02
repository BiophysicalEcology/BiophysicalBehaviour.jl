# A human that thermoregulates

HomoTherm (Kearney et al. 2026) is the human model of NicheMapR: a head, a trunk, two arms and two legs, each
with its own heat budget, inside a loop that vasodilates, lets the core warm and sweats. The tutorial
[A human of many parts](https://biophysicalecology.github.io/HeatExchange.jl/dev/tutorials/human) in the
documentation of HeatExchange.jl builds the same person without the loop, and compares the two in the cold,
stopping at 18 °C where HomoTherm begins to thermoregulate. This tutorial adds the loop and carries the
comparison into the heat.

The two loops are not the same, and the tutorial ends by saying where they part and why.

```@setup human
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
import BiophysicalGeometry: Sphere, Top, Bottom
```

## The person

The proportions are those of HomoTherm, from
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl): a 70 kg person with a
density of 1050 kg/m³. The head is an ellipsoid and the other parts are cylinders.

```@example human
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
import BiologicalScaling

proportions = BiologicalScaling.body_part_proportions(BiologicalScaling.Human())
kinds = (:head, :trunk, :arm, :leg)
part_mass = NamedTuple{kinds}(Tuple(proportions.mass_fraction .* 70.0u"kg"))
ratio = NamedTuple{kinds}(Tuple(proportions.aspect_ratio))
density = 1050.0u"kg/m^3"
shape(kind) = kind == :head ? Ellipsoid(part_mass[kind], density, ratio[kind], ratio[kind]) :
                              Cylinder(part_mass[kind], density, ratio[kind])
nothing # hide
```

Each kind of part has its own physiology, the defaults of HomoTherm for a person at rest. The trunk, arms and
legs are under 6 mm of clothing, and the head under a mean of 5 mm of hair. The face is bare.

```@example human
flesh_conductivity = map(k -> k * u"W/m/K", (head = 1.1, trunk = 0.9, arm = 0.5, leg = 0.5))
fat_fraction = (head = 0.035, trunk = 0.252, arm = 0.07, leg = 0.161)
sky_view = (head = 0.50, trunk = 0.42, arm = 0.35, leg = 0.35)
ground_view = (head = 0.38, trunk = 0.42, arm = 0.35, leg = 0.35)
bare_skin = (head = 0.6, trunk = 0.0, arm = 0.0, leg = 0.0)

covering(depth, fibre_diameter) = example_insulation_pars(;
    fibre_diameter_dorsal = fibre_diameter, fibre_diameter_ventral = fibre_diameter,
    fibre_length_dorsal = 50.0u"mm", fibre_length_ventral = 50.0u"mm",
    insulation_depth_dorsal = depth, insulation_depth_ventral = depth, insulation_depth_compressed = depth,
    fibre_density_dorsal = 3.0e8u"m^-2", fibre_density_ventral = 3.0e8u"m^-2",
    insulation_reflectance_dorsal = 0.3, insulation_reflectance_ventral = 0.3,
)
clothing, hair = covering(6.0u"mm", 1.0u"μm"), covering(5.0u"mm", 75.0u"μm")
coat = (head = hair, trunk = clothing, arm = clothing, leg = clothing)

core = u"K"(36.8u"°C")
resting = 105.0u"W"
part_traits(kind; fat = fat_fraction[kind]) = example_heat_exchange_traits(;
    shape_pars = shape(kind),
    insulation_pars = coat[kind],
    conduction_pars_external = example_conduction_pars_external(; conduction_fraction = 0.0),
    conduction_pars_internal = example_conduction_pars_internal(;
        fat_fraction = fat, flesh_conductivity = flesh_conductivity[kind], fat_density = density),
    radiation_pars = example_radiation_pars(;
        body_emissivity_dorsal = 0.98, body_emissivity_ventral = 0.98,
        sky_view_factor = sky_view[kind], ground_view_factor = ground_view[kind]),
    evaporation_pars = example_evaporation_pars(; skin_wetness = 0.01, bare_skin_fraction = bare_skin[kind]),
    respiration_pars = example_respiration_pars(; oxygen_extraction_efficiency = 0.25, exhaled_relative_humidity = 0.92),
    metabolism_pars = example_metabolism_pars(; core_temperature = core, metabolic_heat_flow = resting, q10 = 2.0),
)
part_body(kind; fat = fat_fraction[kind]) = Body(shape(kind), CompositeInsulation(
    FibrousLayer(coat[kind].dorsal.depth, coat[kind].dorsal.diameter, coat[kind].dorsal.density),
    FatLayer(fat, density)))
nothing # hide
```

The parts are joined where a person's are: the head on top of the trunk, the arms at the shoulders, the legs
beneath.

```@example human
trunk, head, arm, leg = part_body(:trunk), part_body(:head), part_body(:arm), part_body(:leg)
trunk_length = trunk.geometry.length.length_skin
r_arm, r_leg = skin_radius(arm), skin_radius(leg)
shoulder(side) = Attachment(Lateral(trunk_length - r_arm, side), Disc(r_arm))
hip(side) = Attachment(EndA(1.1r_leg, side), Disc(r_leg))
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))

human = CompositeBody(;
    parts = (; trunk, head, arm_left = arm, arm_right = arm, leg_left = leg, leg_right = leg),
    joins = (
        Join(trunk = Attachment(EndB(0.0u"m", 0.0), Disc(r_arm)), head = Attachment(PoleB(), Disc(r_arm))),
        Join(trunk = shoulder(0.0), arm_left = Attachment(Lateral(r_arm, π), Disc(r_arm)); twist = π),
        Join(trunk = shoulder(π), arm_right = Attachment(Lateral(r_arm, 0.0), Disc(r_arm)); twist = π),
        Join(trunk = hip(0.0), leg_left = leg_top),
        Join(trunk = hip(π), leg_right = leg_top),
    ),
)
kind_of = (trunk = :trunk, head = :head, arm_left = :arm, arm_right = :arm, leg_left = :leg, leg_right = :leg)
names = keys(human.parts)
composite_views(human; views = (:oblique, :front, :side), titles = ["", "side", "front"], size = (700, 340)) # hide
```

## What the person can do

HomoTherm's person vasodilates, to a flesh conductivity of 5 W/m/K, lets the core warm to 38 °C, and sweats
until all of the skin is wet. It does not pant. The same limits here, with the same step sizes:

```@example human
limits = ThermoregulationLimits(;
    control = RuleBasedSequentialControl(; mode = CorePantingSweatingFirst()),
    minimum_heat_flow = resting,
    insulation = InsulationLimits(;
        dorsal = SteppedParameter(; current = 6.0u"mm", max = 6.0u"mm", step = 0.0),
        ventral = SteppedParameter(; current = 6.0u"mm", max = 6.0u"mm", step = 0.0)),
    axis_ratio_factor = SteppedParameter(; current = 1.9, max = 1.9, step = 0.1),
    flesh_conductivity = SteppedParameter(; current = 0.9u"W/m/K", max = 5.0u"W/m/K", step = 0.05u"W/m/K"),
    core_temperature = SteppedParameter(; current = core, reference = core, max = u"K"(38.0u"°C"), step = 0.05u"K"),
    panting = PantingLimits(; pant = SteppedParameter(; current = 1.0, max = 1.0, step = 0.1),
                            multiplier = 1.0, core_temperature_ref = core),
    skin_wetness = SteppedParameter(; current = 0.01, max = 1.0, step = 0.005),
)
person = Organism(human, OrganismTraits(Endotherm(), map(part_traits, kind_of),
                                        BehavioralTraits(; thermoregulation = limits); lung_part = :trunk))
nothing # hide
```

With panting given no range, the mode `CorePantingSweatingFirst` lets sweating begin with the rise in core
temperature, as it does in a person.

## The room

Still air at 0.1 m/s and 50 % relative humidity, with the walls at air temperature.

```@example human
function respond(person, air_temperature)
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature),
                                                  relative_humidity = 0.5, wind_speed = 0.1u"m/s")
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    init = (; metabolic_heat_flow = resting, skin_temperature = u"K"(30.0u"°C"),
            insulation_temperature = environment_vars.air_temperature + 10.0u"K")
    thermoregulate(person, environment, init)
end

air = [T * u"°C" for T in -10.0:2.0:34.0]
responses = [respond(person, T) for T in air]
all(part.success for response in responses for part in response.parts)
```

The output of the multi-part controller is that of the multi-part solver, see
[Bodies of many parts](../manual/multipart.md). The whole-body quantities that HomoTherm reports follow from
it. Mean skin temperature is weighted by the area of each part. The core temperature reached is recovered from
the lung temperature, which is the mean of the core and the mean skin temperature. The water evaporated from
the skin is the latent heat lost there divided by the latent heat of vaporisation.

```@example human
areas = [ustrip(u"m^2", total_area(part)) for part in human.parts]
mean_skin(response) = sum(ustrip(u"°C", part.skin_temperature) * area for (part, area) in zip(response.parts, areas)) / sum(areas)
core_reached(response) = ustrip(u"°C", 2 * response.lung_temperature - response.skin_temperature)
function skin_water(response, air_temperature)
    latent_heat = HeatExchange.enthalpy_of_vaporisation(u"K"(air_temperature))
    ustrip(u"g/hr", sum(part.flows.skin_evaporation for part in response.parts) / latent_heat)
end
nothing # hide
```

## Compared with HomoTherm

The reference is `HomoTherm(TA, VEL = 0.1, RH = 50)` of NicheMapR 3.3.3 with its defaults, written by
`docs/src/data/nichemapr_reference.R`.

```@example human
using CSV, DataFrames
reference = CSV.read(joinpath(@__DIR__, "..", "data", "homotherm_whole.csv"), DataFrame)
x = ustrip.(u"°C", air) # hide
fig = Figure(size = (760, 600)) # hide
panels = ( # hide
    ("Metabolic rate (W)", [watts(r.metabolic_heat_flow) for r in responses], reference.metabolic_W), # hide
    ("Mean skin temperature (°C)", mean_skin.(responses), reference.skin_C), # hide
    ("Core temperature (°C)", core_reached.(responses), reference.core_C), # hide
    ("Water evaporated from skin (g/h)", skin_water.(responses, air), 1000 .* reference.cutaneous_water_L_h), # hide
) # hide
for (i, (label, here, there)) in enumerate(panels) # hide
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel = i > 2 ? "Air temperature (°C)" : "", ylabel = label) # hide
    lines!(ax, reference.air_temperature_C, there; linewidth = 2, color = :grey55, label = "HomoTherm") # hide
    scatterlines!(ax, x, here; color = :black, markersize = 7, label = "BiophysicalBehaviour.jl") # hide
    i == 1 && hlines!(ax, [105.0]; color = :black, linestyle = :dash) # hide
    i == 1 && axislegend(ax; position = :rt, labelsize = 10) # hide
end # hide
fig # hide
```

The dashed line is the resting metabolic rate. The grey lines continue to 46 °C, beyond where this model is
run.

```@example human
compared = [T in reference.air_temperature_C for T in x]
at(T) = only(reference[reference.air_temperature_C .== T, :])
rows = map(x[1:3:end], responses[1:3:end], air[1:3:end]) do T, response, air_temperature
    r = at(T)
    (T, round(watts(response.metabolic_heat_flow); digits = 1), round(r.metabolic_W; digits = 1),
     round(mean_skin(response); digits = 1), round(r.skin_C; digits = 1),
     round(skin_water(response, air_temperature); digits = 1), round(1000 * r.cutaneous_water_L_h; digits = 1))
end
markdown_table(["Air (°C)", "Metabolic rate (W)", "HomoTherm", "Skin (°C)", "HomoTherm", "Skin water (g/h)", "HomoTherm"], rows) # hide
```

**In the cold**, up to 18 °C, neither model does anything but make heat, and the comparison is that of the
heat budgets: the metabolic rates agree within about 2 %, as in the HeatExchange.jl tutorial. Here the area
hidden where the parts join is computed from the joins, where HomoTherm uses fixed fractions.

**From 20 °C**, both hold the metabolic rate near the resting 105 W. Vasodilation warms the skin, steeply,
over a few degrees of air temperature. Then sweating begins, and the water evaporated from the skin rises
almost linearly with air temperature, in both, to well over 100 g/h at 34 °C. Below 20 °C the small loss of
water through unsweating skin is about half that of HomoTherm.

The clearest difference is the core. Here it rises to its limit of 38 °C by an air temperature of 28 °C, where
HomoTherm holds 36.8 °C until 32 °C and lets it drift up slowly after that. The metabolic rate here is a few
watts higher through the warm range for that reason, the ``Q_{10}`` effect of a warmer core, and the skin is
about 1 °C warmer. All three follow from the differences below.

## The person at 30 °C

```@example human
warm = respond(person, 30.0u"°C")
skin = NamedTuple{names}(map(part -> ustrip(u"°C", part.skin_temperature), warm.parts))
temperature_views(human, skin; views = (:oblique, :front), titles = ["", "skin"]) # hide
```

```@example human
markdown_table(["Part", "Skin", "Clothing or hair surface", "Heat passed to the surface", "Heat lost from wet skin"], [ # hide
    (string(name), celsius(part.skin_temperature), celsius(part.insulation_temperature), part.net_metabolic, part.flows.skin_evaporation) # hide
    for (name, part) in zip(names, warm.parts)]) # hide
```

Most of the heat that the trunk and limbs pass to their skin leaves as evaporation. The value for the head is
negative, a gain of latent heat at its surface, which has not been examined further.

## Where the two loops differ

| | HomoTherm | Here |
|:--|:--|:--|
| core temperature | one for each part: 36.8, 36.8, 36.5 and 36.7 °C at rest | one, that of the trunk |
| flesh conductivity | one for each part, each raised from its own resting value | each part has its own at rest. The first step of vasodilation sets every part to one value |
| vasodilation and fat | as flesh conductivity rises above 0.5 W/m/K, the fat layer is thinned by a tenth at each step, standing for blood that bypasses it | the fat layer is unchanged |
| what starts sweating | skin temperature: above 35 °C, wetness and core temperature rise together | the heat budget: sweating rises with core temperature once vasodilation is complete |
| rate of vasodilation | faster when mean skin temperature is between 32 and 35 °C, slower above | one step size |
| ceiling on sweating | a maximum sweat rate, 0.75 L/h for each m² | all of the skin wet |
| ceiling on core temperature | can pass 38 °C if nothing else is left | stops at 38 °C |
| each part | a dorsal and a ventral side, averaged | one surface |
| lungs | at a temperature half way through the flesh of the trunk; air exhaled cooler than the lungs in cool air | at the mean of core and mean skin temperature; air exhaled at lung temperature |

HomoTherm's loop was written for a person, with rules that follow what is known of human skin temperature and
sweating. The loop here is the general one of [Endotherm thermoregulation by rules](../manual/endotherm_rules.md),
applied to a body of parts. That the two agree as far as they do says that the order of the responses, and the
physics under them, carry most of the result.

## Where this model stops

Above 34 °C the two diverge. HomoTherm continues to 46 °C in this room, with almost all of the skin wet and
the core a little over 38 °C. Here, at 36 °C the surface solve of one part fails, and from 38 °C the loop ends
with every response at its limit and the heat budget still not balanced at the resting rate: by these
responses, within these limits, the person cannot lose enough heat.

### Without the fat

One difference in the table can be tested directly. In HomoTherm a vasodilated part loses its insulating fat,
since the blood carries heat through it. The nearest thing here is to begin again, at the point where the two
diverge, with a person who has no fat layer at all:

```@example human
lean_body(kind) = part_body(kind; fat = 0.0)
lean = let trunk = lean_body(:trunk), head = lean_body(:head), arm = lean_body(:arm), leg = lean_body(:leg)
    CompositeBody(; parts = (; trunk, head, arm_left = arm, arm_right = arm, leg_left = leg, leg_right = leg),
                  joins = human.joins)
end
lean_person = Organism(lean, OrganismTraits(Endotherm(), map(kind -> part_traits(kind; fat = 0.0), kind_of),
                                            BehavioralTraits(; thermoregulation = limits); lung_part = :trunk))
hot = [T * u"°C" for T in 30.0:2.0:36.0]
lean_responses = [respond(lean_person, T) for T in hot]
rows = map(hot, lean_responses) do T, response
    r = at(ustrip(u"°C", T))
    (ustrip(u"°C", T), all(part.success for part in response.parts),
     round(watts(response.metabolic_heat_flow); digits = 1), round(r.metabolic_W; digits = 1),
     round(mean_skin(response); digits = 1), round(r.skin_C; digits = 1),
     round(core_reached(response); digits = 1), round(r.core_C; digits = 1),
     round(skin_water(response, T); digits = 1), round(1000 * r.cutaneous_water_L_h; digits = 1))
end
markdown_table(["Air (°C)", "Converged", "Metabolic rate (W)", "HomoTherm", "Skin (°C)", "HomoTherm", "Core (°C)", "HomoTherm", "Skin water (g/h)", "HomoTherm"], rows) # hide
```

Without fat the person balances the heat budget at 36 °C, where the person with fat does not. Run further, the
same test balances to within a few watts at 40 °C and fails beyond it, with one failure of the surface solve at
38 °C on the way. Those runs are slow, because the loop goes to its limit of iterations, and are not repeated
on this page.

So fat is part of the answer and not all of it. The lean person reaches 36 °C with a core at the 38 °C limit and
a skin about 2 °C warmer than HomoTherm's, and still evaporates less water. At the limit, with all the skin
wet, about 0.18 kg/h leaves the skin here, where HomoTherm evaporates 0.22 to 0.28 kg/h from a cooler skin
between 40 and 46 °C. The remaining difference is therefore in how much water a wet skin under clothing can
evaporate, which is a question for the heat budget of HeatExchange.jl and not for the controller.

What a human-specific controller would add to this one is then clear: vasodilation that thins the fat layer of
the part, the skin-temperature trigger for sweating, and a core temperature for each part.
