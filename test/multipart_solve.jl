# BB-side multi-part solve integration (Phase 4 / 6-7).
#
# Two gates:
#  1. A single `Body` organism, driven through `solve_multipart_metabolic_rate`,
#     reproduces the dorsal/ventral `solve_metabolic_rate` baseline in the
#     symmetric case — the Phase 7 equivalence gate, now from the organism API.
#  2. A two-part organism solves, routes respiration/lung mass through the
#     `lung_part`, and returns one surface result per part.

using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Unitful, UnitfulMoles
using Test

const ρ = 1000.0u"kg/m^3"
const b = 3.0

# Shared, symmetric physiology (dorsal == ventral, no ground contact, no solar).
ins_pars = example_insulation_pars()
fibre    = ins_pars.dorsal
fur      = FibrousLayer(fibre.depth, fibre.diameter, fibre.density)
fat      = FatLayer(0.0, 901.0u"kg/m^3")
internal = example_conduction_pars_internal()
external = example_conduction_pars_external(; conduction_fraction = 0.0)
evap     = example_evaporation_pars()
rad_pars = example_radiation_pars()          # symmetric: sky = ground = 0.5
resp     = example_respiration_pars()
metab    = example_metabolism_pars()
core     = metab.core_temperature

env_vars = example_environment_vars()         # uniform temps, no solar, shade = 0
env_pars = example_environment_pars()
environment = (; environment_pars = env_pars, environment_vars = env_vars)

behav = BehavioralTraits(; thermoregulation = example_thermoregulation_limits())

heat_exchange_of(shape) = example_heat_exchange_traits(;
    shape_pars               = shape,
    insulation_pars          = ins_pars,
    conduction_pars_external = external,
    conduction_pars_internal = internal,
    radiation_pars           = rad_pars,
    evaporation_pars         = evap,
    respiration_pars         = resp,
    metabolism_pars          = metab,
)

skin0  = core - 5u"K"
insul0 = env_vars.air_temperature + 2u"K"

@testset "single Body ≡ dorsal/ventral baseline (from organism API)" begin
    shape   = Cylinder(1.0u"kg", ρ, b)
    body    = Body(shape, CompositeInsulation(fur, fat))
    traits  = OrganismTraits(Endotherm(), heat_exchange_of(shape), behav)
    organism = Organism(body, traits)

    baseline  = solve_metabolic_rate(organism, environment, skin0, insul0)
    multipart = solve_multipart_metabolic_rate(organism, environment, skin0, insul0)

    # One part is built for the single body, and it hosts the lung.
    @test length(multipart.parts) == 1

    @test multipart.metabolic_heat_flow ≈ baseline.energy_flows.metabolic_heat_flow rtol = 1e-3
    @test multipart.skin_temperature ≈ baseline.thermoregulation.skin_temperature rtol = 1e-4
    @test multipart.insulation_temperature ≈ baseline.thermoregulation.insulation_temperature rtol = 1e-4
end

@testset "two-part organism solves, lung routed to torso" begin
    torso_shape = Cylinder(0.6u"kg", ρ, b)
    head_shape  = Cylinder(0.4u"kg", ρ, b)
    torso = Body(torso_shape, CompositeInsulation(fur, fat))
    head  = Body(head_shape,  CompositeInsulation(fur, fat))
    r_join = 0.02u"m"
    dog = CompositeBody(;
        parts = (; torso, head),
        joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_join)),
                      head  = Attachment(EndB(0.0u"m", 0.0), Disc(r_join))),),
    )

    phys = (; torso = heat_exchange_of(torso_shape), head = heat_exchange_of(head_shape))
    traits = OrganismTraits(Endotherm(), phys, behav; lung_part = :torso)
    organism = Organism(dog, traits)

    result = solve_multipart_metabolic_rate(organism, environment, skin0, insul0)

    # One surface solve per part, in body order.
    @test length(result.parts) == 2
    @test all(part -> part.success, result.parts)

    # Respiration is closed at the lung (torso), not the head.
    @test BiophysicalBehaviour._lung_mass(dog, :torso) == 0.6u"kg"

    # A regulated endotherm produces positive metabolic heat and holds skin below core.
    @test result.metabolic_heat_flow > 0.0u"W"
    @test result.skin_temperature < core
    @test result.lung_temperature < core
end
