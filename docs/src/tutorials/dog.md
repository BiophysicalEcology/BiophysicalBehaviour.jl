# A dog of many parts

A dog is not an ellipsoid. Its legs are thin, with much surface for their mass, and its lungs are in its trunk.
This tutorial builds a dog of six parts, trunk, head and four legs, each with its own heat budget, and lets it
thermoregulate. It is the behavioural counterpart of
[A human of many parts](https://biophysicalecology.github.io/HeatExchange.jl/dev/tutorials/human) in the
documentation of HeatExchange.jl.

```@setup dog
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
import BiophysicalGeometry: Sphere, Top, Bottom
```

## The body

Each part is a cylinder with the same 2 mm pelt and no fat: trunk 18 kg, head 2 kg, each leg 1 kg, 24 kg in all.

```@example dog
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

density = 1000.0u"kg/m^3"
pelt_pars = example_insulation_pars()
fibres = pelt_pars.dorsal
pelt = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3"))
furred(shape) = Body(shape, pelt)

trunk_shape = Cylinder(18.0u"kg", density, 3.0)
head_shape = Cylinder(2.0u"kg", density, 1.5)
leg_shape = Cylinder(1.0u"kg", density, 5.0)

trunk = furred(trunk_shape)
trunk_length = trunk.geometry.length.length_skin
neck, hip = Disc(0.02u"m"), Disc(0.015u"m")
leg_join(name, along, around) = Join(;
    trunk = Attachment(Lateral(along * trunk_length, around), hip),
    (name => Attachment(EndB(0.0u"m", 0.0), hip),)...,
)

dog_body = CompositeBody(;
    parts = (; trunk, head = furred(head_shape),
             leg_fl = furred(leg_shape), leg_fr = furred(leg_shape),
             leg_bl = furred(leg_shape), leg_br = furred(leg_shape)),
    joins = (
        Join(trunk = Attachment(EndA(0.0u"m", 0.0), neck), head = Attachment(EndB(0.0u"m", 0.0), neck)),
        leg_join(:leg_fl, 0.25, π / 2 + 0.3), leg_join(:leg_fr, 0.25, π / 2 - 0.3),
        leg_join(:leg_bl, 0.75, π / 2 + 0.3), leg_join(:leg_br, 0.75, π / 2 - 0.3),
    ),
)
composite_views(dog_body; views = (:oblique, :side, :front), titles = ["", "side", "front"], size = (720, 300)) # hide
```

The joins say where each part meets the trunk and over what area. That area is hidden from the air.

```@example dog
body_graph(dog_body) # hide
```

```@example dog
names = part_names(dog_body)
markdown_table(["Part", "Mass", "Surface area", "Area per unit mass"], [ # hide
    (string(name), mass(part.shape), total_area(part), uconvert(u"cm^2/kg", total_area(part) / mass(part.shape))) # hide
    for (name, part) in pairs(dog_body.parts)]) # hide
```

A leg has three times the surface of the trunk per kilogram.

## The physiology

Each part has its own `HeatExchangeTraits`, here differing only in shape. The trunk holds the lungs. All joins
are the default `SharedCore`, so the dog has one core temperature, 37 °C, and one minimum metabolic rate,
77.6 W. See [Bodies of many parts](../manual/multipart.md) and, for the couplings,
[Bodies of many parts](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/multipart) in the
documentation of HeatExchange.jl.

```@example dog
part_traits(shape) = example_heat_exchange_traits(;
    shape_pars = shape,
    insulation_pars = pelt_pars,
    conduction_pars_external = example_conduction_pars_external(; conduction_fraction = 0.0),
)
physiology_traits = (;
    trunk = part_traits(trunk_shape), head = part_traits(head_shape),
    leg_fl = part_traits(leg_shape), leg_fr = part_traits(leg_shape),
    leg_bl = part_traits(leg_shape), leg_br = part_traits(leg_shape),
)
dog(limits) = Organism(dog_body,
    OrganismTraits(Endotherm(), physiology_traits, BehavioralTraits(; thermoregulation = limits); lung_part = :trunk))
ruled_dog = dog(example_thermoregulation_limits(; skin_wetness_max = 0.3))
lung_part(ruled_dog), pant_selector(ruled_dog), organism_compartment_graph(ruled_dog)
```

## One air temperature

The first guess of the fur surface temperature is put half way between skin and air. With the default guess,
air temperature, the trunk's surface solve does not converge in the cold, see
[Bodies of many parts](../manual/multipart.md#Solving).

```@example dog
function respond(animal, air_temperature)
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature))
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    skin = metabolism_pars(animal).core_temperature - 4.0u"K"
    init = (; metabolic_heat_flow = thermoregulation(animal).minimum_heat_flow, skin_temperature = skin,
            insulation_temperature = (skin + environment_vars.air_temperature) / 2)
    thermoregulate(animal, environment, init)
end

cold = respond(ruled_dog, 5.0u"°C")
markdown_table(["Part", "Skin", "Fur surface", "Heat passed to the surface", "Per unit mass", "Converged"], [ # hide
    (string(name), celsius(part.skin_temperature), celsius(part.insulation_temperature), part.net_metabolic, # hide
     uconvert(u"W/kg", part.net_metabolic / mass(body_part.shape)), part.success) # hide
    for (name, part, body_part) in zip(names, cold.parts, dog_body.parts)]) # hide
```

```@example dog
(metabolic_rate = cold.metabolic_heat_flow, lung_temperature = celsius(cold.lung_temperature))
```

At 5 °C the dog needs about twice its minimum metabolic rate. Its legs are a sixth of its mass and lose two
fifths of the heat that leaves through its surface.

```@example dog
skin = NamedTuple{names}(map(part -> ustrip(u"°C", part.skin_temperature), cold.parts))
surface = NamedTuple{names}(map(part -> ustrip(u"°C", part.insulation_temperature), cold.parts))
temperature_views(dog_body, skin; views = (:oblique, :side), titles = ["skin", ""]) # hide
```

```@example dog
temperature_views(dog_body, surface; views = (:oblique, :side), titles = ["fur surface", ""], label = "Fur surface temperature (°C)") # hide
```

With a shared core, every part has blood at 37 °C at its centre, and its skin temperature is set by how far the
skin is from the centre. The trunk, the thickest part, has the coolest skin. A real dog saves heat by letting
its legs cool, with countercurrent exchange in the limbs. That needs each leg to have its own core temperature,
joined by a `ConductiveCoupling`, which the controllers do not yet support.

## Across air temperatures

```@example dog
air = 0.0:2.5:45.0
responses = [respond(ruled_dog, T * u"°C") for T in air]
all(part.success for response in responses for part in response.parts)
```

```@example dog
fig = Figure(size = (760, 560)) # hide
ax = Axis(fig[1, 1]; ylabel = "Metabolic rate (W)") # hide
lines!(ax, air, [watts(r.metabolic_heat_flow) for r in responses]; color = :black, linewidth = 2) # hide
hlines!(ax, [77.6]; color = :grey, linestyle = :dash) # hide
ax = Axis(fig[1, 2]; ylabel = "Skin temperature (°C)") # hide
for (i, name) in enumerate(names[1:3]) # hide
    lines!(ax, air, [ustrip(u"°C", r.parts[i].skin_temperature) for r in responses]; linewidth = 2, label = i == 3 ? "legs" : string(name)) # hide
end # hide
axislegend(ax; position = :rb, labelsize = 10) # hide
ax = Axis(fig[2, 1]; xlabel = "Air temperature (°C)", ylabel = "Heat passed to the surface (W)") # hide
lines!(ax, air, [watts(r.parts[1].net_metabolic) for r in responses]; linewidth = 2, label = "trunk") # hide
lines!(ax, air, [watts(r.parts[2].net_metabolic) for r in responses]; linewidth = 2, label = "head") # hide
lines!(ax, air, [watts(sum(p.net_metabolic for p in r.parts[3:6])) for r in responses]; linewidth = 2, label = "four legs") # hide
hlines!(ax, [0.0]; color = :grey, linestyle = :dash) # hide
axislegend(ax; position = :rt, labelsize = 10) # hide
ax = Axis(fig[2, 2]; xlabel = "Air temperature (°C)", ylabel = "Heat lost in breathing (W)") # hide
lines!(ax, air, [watts(r.metabolic_heat_flow - r.net_metabolic_total) for r in responses]; color = :black, linewidth = 2) # hide
fig # hide
```

The dashed line is the minimum metabolic rate. Below about 20 °C the dog makes extra heat. Above it the
controller acts, as in [Endotherm thermoregulation by rules](../manual/endotherm_rules.md): vasodilation in
every part (the step up in skin temperature), then a rise in core temperature, then panting.

The panels on the right and below show where the heat goes. As the air warms, the heat each part can pass to
its surface falls, and once the air is warmer than the skin it is negative: the surface gains heat from the air.
All the metabolic heat, and that gain, then leaves through the lungs. Panting applies to the lung part alone,
so the trunk's traits set its cost.

The dog cannot change posture. A real one curls up in the cold and sprawls in the heat. For a body of parts
that is a change of pose, which hides or exposes surface at the joins and changes what each part sees. It is
planned, in place of the change of axis ratio a single body uses, see
[Bodies of many parts](../manual/multipart.md#Thermoregulating).

## The same dog, by optimisation

With [`IPOPTControl`](@ref), flesh conductivity and skin wetness are variables of each part, 27 variables in
all, see [The nonlinear program](../manual/nlp.md).

```@example dog
using Setfield: @set
limits = example_thermoregulation_limits(; skin_wetness_max = 0.3)
optimising_dog = dog(@set limits.control = IPOPTControl())
warm = respond(optimising_dog, 30.0u"°C")
markdown_table(["Part", "Flesh conductivity", "Skin wetness", "Skin temperature", "Heat lost from wet skin"], [ # hide
    (string(name), part.flesh_conductivity, part.skin_wetness, celsius(part.skin_temperature), part.flows.skin_evaporation_heat_flow) # hide
    for (name, part) in pairs(warm.parts)]) # hide
```

```@example dog
(metabolic_rate = warm.metabolic_heat_flow, core_temperature = celsius(warm.core_temperature), panting = warm.panting_rate)
```

The optimiser treats the parts differently. The cost of skin wetness is on its mean over the parts, so the
water goes where it does most good. The rules move every part together. Real dogs sweat only through their
paws, so a model of a dog should lower `skin_wetness_max` and let panting do the work. The point is the
mechanism: with parts, *where* on the body a response is made can be asked.

## Parts that shade each other

In the sun, each part intercepts the direct beam by its silhouette, and the parts hide sky and ground from each
other. [`precompute_view_partition`](@ref) computes both for a position of the sun:

```@example dog
sunny = example_environment_vars(; air_temperature = u"K"(20.0u"°C"), global_radiation = 800.0u"W/m^2", zenith_angle = 30.0u"°")
view = precompute_view_partition(ruled_dog, sunny)
markdown_table(["Part", "Sky", "Ground", "Other parts", "Sunlit silhouette"], [ # hide
    (string(name), v.sky, v.ground, sum(values(v.neighbours)), uconvert(u"cm^2", v.lit_silhouette)) for (name, v) in pairs(view)]) # hide
```

No part sees a full hemisphere of sky and ground: other parts fill a fifth of the view of a leg or the head and
a tenth of the trunk's. Passed to [`solve_multipart_metabolic_rate`](@ref), the partition replaces each part's
view factors:

```@example dog
environment = (; environment_pars = example_environment_pars(), environment_vars = sunny)
guess = (u"K"(33.0u"°C"), u"K"(30.0u"°C"))
alone = solve_multipart_metabolic_rate(ruled_dog, environment, guess...)
together = solve_multipart_metabolic_rate(ruled_dog, environment, guess...; view)
markdown_table(["Part", "Solar heat, parts alone", "Solar heat, with shading", "Fur surface, alone", "Fur surface, with shading"], [ # hide
    (string(name), a.flows.solar, b.flows.solar, celsius(a.insulation_temperature), celsius(b.insulation_temperature)) # hide
    for (name, a, b) in zip(names, alone.parts, together.parts)]) # hide
```

With the sun 30° from overhead, the head is partly shaded by the trunk and absorbs little more than half what it
would alone, while the sunward legs gain.
