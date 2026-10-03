# Bodies of many parts

HeatExchange.jl can solve the heat budget of a body of several parts: each with its own surface, joined through
a shared or a conducting core. Its documentation
([Bodies of many parts](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/multipart)) builds the
inputs of that solve by hand. This package builds them from an organism, and thermoregulates the result.

```@setup multipart
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
import BiophysicalGeometry: Sphere, Top, Bottom
```

## An organism of many parts

The body is a `CompositeBody` of
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl): named parts, each a
`Body`, and the joins between them. Its documentation has an
[interactive builder](https://biophysicalecology.github.io/BiophysicalGeometry.jl/dev/builder) for bodies of
this kind. Here is a trunk with a head, both furred cylinders:

```@example multipart
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

fibres = example_insulation_pars().dorsal
coat = CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3"))
trunk_shape = Cylinder(20.0u"kg", 1000.0u"kg/m^3", 3.0)
head_shape = Cylinder(2.0u"kg", 1000.0u"kg/m^3", 1.5)
neck = Disc(0.03u"m")
animal_body = CompositeBody(;
    parts = (; trunk = Body(trunk_shape, coat), head = Body(head_shape, coat)),
    joins = (Join(trunk = Attachment(EndA(0.0u"m", 0.0), neck), head = Attachment(EndB(0.0u"m", 0.0), neck)),),
)
composite_views(animal_body; views = (:oblique, :side), titles = ["", "side"], size = (560, 260)) # hide
```

The physiology is a NamedTuple with the same names, one `HeatExchangeTraits` for each part, and
[`OrganismTraits`](@ref) is told which part holds the lungs:

```@example multipart
part_traits(shape; kw...) = example_heat_exchange_traits(;
    shape_pars = shape,
    conduction_pars_external = example_conduction_pars_external(; conduction_fraction = 0.0),
    kw...,
)
part_physiology_traits = (; trunk = part_traits(trunk_shape), head = part_traits(head_shape))
behaviour = BehavioralTraits(; thermoregulation = example_thermoregulation_limits(; minimum_heat_flow = 40.0u"W"))
traits = OrganismTraits(Endotherm(), part_physiology_traits, behaviour; lung_part = :trunk)
animal = Organism(animal_body, traits)
part_names(animal_body), lung_part(animal)
```

A single `HeatExchangeTraits` given to a `CompositeBody` is applied to every part. A plain `Body` is treated
throughout as a composite of one part, named `:body`, so there is one code path.

## Whole-organism and per-part quantities

| Quantity | Belongs to | Taken from |
|:--|:--|:--|
| shape, insulation, radiation, evaporation, conduction to the ground | each part | that part's traits |
| flesh conductivity, skin wetness | each part | that part's traits |
| core temperature, minimum metabolic rate, ``Q_{10}`` | the organism | the traits of the lung part |
| ventilation, panting, oxygen extraction | the organism | the traits of the lung part |

Respiration happens in one place. [`physiology`](@ref) returns the per-part traits with the lung part wrapped in
a [`LungPart`](@ref), on which the respiratory calculations dispatch:

```@example multipart
map(is_lung_part, physiology(animal))
```

## Selecting parts

A [`PartSelector`](@ref) names the parts something applies to:

| Selector | Selects |
|:--|:--|
| [`WholeBody`](@ref)`()` | every part |
| [`ByName`](@ref)`(:head, :trunk)` | the parts named |
| [`Compartment`](@ref)`(:trunk)` | the parts of a thermal compartment |

```@example multipart
select_names(WholeBody(), animal_body), select_names(pant_selector(animal), animal_body)
```

[`map_parts`](@ref), [`foldl_parts`](@ref) and [`set_part`](@ref) apply a function to the selected parts of a
body, and [`map_part_physiology`](@ref) to the selected parts' traits. This is how a response is applied to part
of an animal. Here the skin of the head alone is wetted:

```@example multipart
using Setfield: @set
wet_headed = map_part_physiology(ByName(:head), animal) do part
    @set part.evaporation_pars.skin_wetness = 0.5
end
map(part -> evaporation_pars(part).skin_wetness, physiology(wet_headed))
```

## Compartments

How the cores of two joined parts are related is a `HeatCoupling` of HeatExchange.jl, one for each join:

| Coupling | Meaning |
|:--|:--|
| `SharedCore()` | blood mixes the two cores to one temperature. The parts form one compartment |
| `ConductiveCoupling()` | each has its own core temperature, and heat is conducted between them through the join |

They are given as the `couplings` keyword of [`OrganismTraits`](@ref), a tuple in the order of the joins. With
none given, every join is a `SharedCore`, and the animal has one regulated core:

```@example multipart
organism_compartment_graph(animal)
```

!!! note "One core for now"
    Thermoregulation is at present supported for a single shared core. Compartments with their own core
    temperatures can be declared and are solved by `solve_multipart_metabolic_rate`, but that path is still in
    development and not yet used by the controllers.

## Solving

[`solve_multipart_metabolic_rate`](@ref) is the multi-part counterpart of `solve_metabolic_rate`. From the
organism and environment it builds the setup of each part ([`part_surface_setups`](@ref)), with the area hidden
under each join removed, and passes them to `solve_coupled_metabolic_rate` of HeatExchange.jl, which solves
each surface and closes the respiration balance once, in the lungs:

```@example multipart
environment = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(; air_temperature = u"K"(15.0u"°C")),
)
skin_guess, surface_guess = u"K"(32.0u"°C"), u"K"(25.0u"°C")
out = solve_multipart_metabolic_rate(animal, environment, skin_guess, surface_guess)
names = part_names(animal_body) # hide
markdown_table(["Part", "Skin", "Fur surface", "Heat passed to the surface", "Converged"], [ # hide
    (string(name), celsius(part.skin_temperature), celsius(part.insulation_temperature), part.net_metabolic, part.success) # hide
    for (name, part) in zip(names, out.parts)]) # hide
```

```@example multipart
(out.metabolic_heat_flow, out.net_metabolic_total, celsius(out.lung_temperature))
```

| Field of the output | Content |
|:--|:--|
| `metabolic_heat_flow` | the heat production that balances the whole organism |
| `parts` | for each part, in order: `skin_temperature`, `insulation_temperature`, `net_metabolic`, `flows`, `success` |
| `net_metabolic_total` | the sum of the heat passed to the surfaces of the parts |
| `skin_temperature`, `insulation_temperature` | means over the parts |
| `lung_temperature`, `respiration_out` | the respiration balance |

Check the `success` of every part. The surface solve of a large, well-insulated part can fail from a poor first
guess, in particular a fur surface at air temperature in the cold. A guess between skin and air temperature, as
above, is more robust.

### Geometry that does not change

The areas, silhouette and characteristic dimension of each part depend only on its shape and coat. A
[`ShapeCache`](@ref) from [`precompute_shape_cache`](@ref) holds them, and [`refresh`](@ref) rebuilds it only
after a response that changes geometry. The rule-based controller makes one before its loop.

### Parts that see each other

A part does not see a full hemisphere of sky and ground: its neighbours are in the way.
[`precompute_view_partition`](@ref) computes, for the present pose and sun, how much of each part's view is sky,
ground and each other part, and the silhouette of each part the sun reaches:

```@example multipart
view = precompute_view_partition(animal, environment.environment_vars)
markdown_table(["Part", "Sky", "Ground", "Other parts", "Sunlit silhouette"], [ # hide
    (string(name), v.sky, v.ground, sum(values(v.neighbours)), v.lit_silhouette) for (name, v) in pairs(view)]) # hide
```

Passed as the `view` keyword of `solve_multipart_metabolic_rate`, it replaces the view factors of each part's
traits and adds longwave exchange between the parts, see
[Parts that see each other](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/multipart#Parts-that-see-each-other)
in HeatExchange.jl.

## Thermoregulating

[`thermoregulate`](@ref) takes the organism as it would any other. With the rule-based controller, responses are
applied by scope:

| Response | Applied to |
|:--|:--|
| vasodilation | every part |
| rise in core temperature | the organism, through the lung part |
| panting | the lung part |
| sweating | every part |

Raising and flattening the coat and changing posture are not yet applied to bodies of many parts.

!!! note "Posture as pose"
    For a body of many parts, posture will be a change of *pose*: raising or lowering wings or ears, bringing
    the limbs in to the body or holding them away. That changes which surfaces are hidden under joins and what
    each part sees of the sky, the ground and its neighbours, through the geometry of BiophysicalGeometry.jl.
    It is to replace the change of axis ratio by which a single body curls and uncurls.

```@example multipart
function respond(animal, air_temperature)
    environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature))
    surroundings = (; environment_pars = example_environment_pars(), environment_vars)
    core = metabolism_pars(animal).core_temperature
    skin = core - 4.0u"K"
    init = (; metabolic_heat_flow = 0.0u"W", skin_temperature = skin,
            insulation_temperature = (skin + environment_vars.air_temperature) / 2)
    thermoregulate(animal, surroundings, init)
end

air = 10.0:2.5:40.0
responses = [respond(animal, T * u"°C") for T in air]
fig, ax = figure_axis("Air temperature (°C)", "Skin temperature (°C)"; size = (700, 380)) # hide
for (i, name) in enumerate(names) # hide
    lines!(ax, air, [ustrip(u"°C", r.parts[i].skin_temperature) for r in responses]; linewidth = 2, label = string(name)) # hide
end # hide
ax2 = Axis(fig[1, 2]; xlabel = "Air temperature (°C)", ylabel = "Metabolic rate (W)") # hide
lines!(ax2, air, [watts(r.metabolic_heat_flow) for r in responses]; linewidth = 2, color = :black) # hide
axislegend(ax; position = :rb) # hide
fig # hide
```

The trunk has the cooler skin in the cold. It is the thicker part, and with a shared core its skin is further
from the warm centre. See [A dog of many parts](../tutorials/dog.md) for a body with legs, and
[A human that thermoregulates](../tutorials/human.md) for a person.

With [`IPOPTControl`](@ref), flesh conductivity and skin wetness are variables of each part, so the optimiser
can vasodilate or wet one part and not another, see [Thermoregulation by optimisation](optimisation.md).
