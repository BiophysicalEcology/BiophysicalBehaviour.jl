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
using Setfield: @set
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

@testset "genuine per-part ≡ dorsal/ventral at the idealised point (Phase 8 gate)" begin
    # The one point where the dorsal/ventral model and the genuine per-part model must
    # coincide: an animal standing upright at midday in symmetric surroundings — sun
    # overhead (zenith 0), uniform radiant environment, symmetric insulation and view
    # factors (ventral_fraction = 0.5, sky = ground = 0.5, dorsal == ventral), and
    # *no ground contact* (conduction_fraction = 0, i.e. upright, not lying down).
    # There the two paths agree on the regulated quantities to numerical tolerance and
    # both conserve energy. Ground contact is deliberately excluded: it is exactly
    # where the old model's conflation bites (ventral solar carries a (1 − 2·cond)
    # factor, ground conduction is weighted by the ground *view factor*), so the two
    # models diverge there by design — that divergence is the reason dorsal/ventral is
    # being retired, not a regression. See memory phase7-done-phase8-scope.
    shape = Cylinder(1.0u"kg", ρ, b)
    body  = Body(shape, CompositeInsulation(fur, fat))

    ideal_env_vars = @set env_vars.zenith_angle = 0.0u"°"          # sun overhead
    ideal_env_vars = @set ideal_env_vars.global_radiation = 500.0u"W/m^2"
    heat_exchange = example_heat_exchange_traits(;
        shape_pars               = shape,
        insulation_pars          = ins_pars,
        conduction_pars_external = external,          # conduction_fraction = 0 (upright)
        conduction_pars_internal = internal,
        radiation_pars           = rad_pars,
        evaporation_pars         = evap,
        respiration_pars         = resp,
        metabolism_pars          = metab,
    )
    org = Organism(body, OrganismTraits(Endotherm(), heat_exchange, behav))
    ideal_env = (; environment_pars = env_pars, environment_vars = ideal_env_vars)
    init = (; metabolic_heat_flow = 0.0u"W", skin_temperature = skin0, insulation_temperature = insul0)
    ctrl = control_strategy(org)

    old = thermoregulate(Endotherm(), ctrl, org, ideal_env, init;
                         inner_solve = BiophysicalBehaviour._multisided_inner)
    new = thermoregulate(Endotherm(), ctrl, org, ideal_env, init;
                         inner_solve = BiophysicalBehaviour._multipart_inner)

    # Regulated quantities agree at the assumption point (essentially exact).
    @test new.energy_flows.metabolic_heat_flow ≈ old.energy_flows.metabolic_heat_flow rtol = 1e-3
    @test new.thermoregulation.core_temperature ≈ old.thermoregulation.core_temperature rtol = 1e-4
    @test new.thermoregulation.skin_temperature ≈ old.thermoregulation.skin_temperature rtol = 1e-4
    @test new.thermoregulation.insulation_temperature ≈ old.thermoregulation.insulation_temperature rtol = 1e-4

    # Both paths conserve energy at the converged state (net balance ≈ 0).
    @test abs(ustrip(u"W", new.energy_flows.heat_balance)) < 1e-3
    @test abs(ustrip(u"W", old.energy_flows.heat_balance)) < 1e-3
end

@testset "occlusion view partition wired through the solve (Phase 8)" begin
    # A dorsal-up two-half body under an overhead sun: the ventral (belly) half is
    # shadowed by the dorsal half. With the precomputed occlusion view the ventral half
    # receives no direct beam and its sky/ground split emerges from geometry; the naive
    # per-part path double-counts the sun on both halves. `precompute_view_partition`
    # reads the body's (biological, dorsal-up) pose — the modeller orients it.
    M = 1.0u"kg"
    sun_vars = @set env_vars.global_radiation = 500.0u"W/m^2"
    sun_vars = @set sun_vars.zenith_angle = 0.0u"°"
    sun_env  = (; environment_pars = env_pars, environment_vars = sun_vars)
    root = BiophysicalGeometry.Pose((0.0u"m", 0.0u"m", 0.0u"m"),
                                    BiophysicalGeometry.rotation_axis_angle((1.0, 0.0, 0.0), π / 2))
    dv = CompositeBody(;
        parts = (; dorsal  = Body(HalfCylinder(M / 2, ρ, b), CompositeInsulation(fur, fat)),
                   ventral = Body(HalfCylinder(M / 2, ρ, b), CompositeInsulation(fur, fat))),
        joins = (Join(dorsal  = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),),
        root_pose = root)
    half_metab = @set metab.metabolic_heat_flow = metab.metabolic_heat_flow / 2
    mkp(sky, cf) = example_heat_exchange_traits(; shape_pars = HalfCylinder(M / 2, ρ, b),
        insulation_pars = ins_pars,
        conduction_pars_external = example_conduction_pars_external(; conduction_fraction = cf),
        conduction_pars_internal = internal,
        radiation_pars = example_radiation_pars(; sky_view_factor = sky),
        evaporation_pars = evap, respiration_pars = resp, metabolism_pars = half_metab)
    org = Organism(dv, OrganismTraits(Endotherm(),
        (; dorsal = mkp(1.0, 0.0), ventral = mkp(0.0, 0.4)), behav; lung_part = :dorsal))

    view = precompute_view_partition(org, sun_vars; ndirections = 200, resolution = 64)

    # The view sees the dorsal/ventral split and the ventral shadow (no direct beam).
    @test view.dorsal.sky > view.dorsal.ground
    @test view.ventral.ground > view.ventral.sky
    @test view.ventral.lit_silhouette < 0.05 * view.dorsal.lit_silhouette

    naive = solve_multipart_metabolic_rate(org, sun_env, skin0, insul0)
    occl  = solve_multipart_metabolic_rate(org, sun_env, skin0, insul0; view)
    # occlusion removes the ventral half's phantom direct beam
    @test occl.parts[2].flows.solar < naive.parts[2].flows.solar
    @test all(p -> p.success, occl.parts)

    # Inter-part surface exchange (Phase 8, step 3). Each half's occluded solid angle
    # exchanges longwave with the *other* half's surface (the reciprocal neighbour
    # fractions) rather than the sky/ground behind it, so the view partitions each sum
    # to one — no solid angle, hence no energy, is lost — and the two halves' surface
    # temperatures are coupled.
    @test view.dorsal.sky  + view.dorsal.ground  + view.dorsal.neighbours.ventral ≈ 1 atol = 1e-6
    @test view.ventral.sky + view.ventral.ground + view.ventral.neighbours.dorsal ≈ 1 atol = 1e-6
    @test view.dorsal.neighbours.ventral ≈ view.ventral.neighbours.dorsal rtol = 0.05

    # An otherwise-identical view with the neighbour term stripped (same sky/ground/lit)
    # isolates the coupling: `occl` has it, `decoupled` does not.
    view0 = NamedTuple{keys(view)}(map(keys(view)) do n
        (; view[n].sky, view[n].ground, neighbours = NamedTuple(), view[n].lit_silhouette)
    end)
    decoupled = solve_multipart_metabolic_rate(org, sun_env, skin0, insul0; view = view0)

    gap(out) = out.parts[1].insulation_temperature - out.parts[2].insulation_temperature
    # Coupling shrinks the warm-dorsal ↔ cool-ventral surface-temperature gap: the warm
    # half sheds heat to the cool half (dorsal cools, ventral warms).
    @test gap(occl) < gap(decoupled)
    @test occl.parts[1].insulation_temperature < decoupled.parts[1].insulation_temperature
    @test occl.parts[2].insulation_temperature > decoupled.parts[2].insulation_temperature
    # The warm half now radiates more (also to the cool sibling); the cool half less.
    @test occl.parts[1].flows.longwave > decoupled.parts[1].flows.longwave
    @test occl.parts[2].flows.longwave < decoupled.parts[2].flows.longwave
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

# Build a reusable two-part (torso lung + head) organism.
function two_part_organism()
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
    return Organism(dog, traits)
end

@testset "panting affects only the lung part (§3.9 exit criterion)" begin
    organism = two_part_organism()
    limits = thermoregulation(organism)

    before = physiology(organism)
    _, _, panted = BiophysicalBehaviour._pant_multipart(organism, limits.panting)
    after = physiology(panted)

    # Pant rate rises on the torso (lung) and is untouched on the head.
    @test HeatExchange.respiration_pars(after.torso).pant > HeatExchange.respiration_pars(before.torso).pant
    @test HeatExchange.respiration_pars(after.head).pant == HeatExchange.respiration_pars(before.head).pant

    # The head's cutaneous physiology (skin wetness, flesh conductivity) is unchanged;
    # panting is respiratory, at the lung, not a whole-body cutaneous response.
    @test HeatExchange.evaporation_pars(after.head).skin_wetness ==
          HeatExchange.evaporation_pars(before.head).skin_wetness
    @test HeatExchange.conduction_pars_internal(after.head).flesh_conductivity ==
          HeatExchange.conduction_pars_internal(before.head).flesh_conductivity
end

@testset "per-part vasodilation and sweating span all parts" begin
    organism = two_part_organism()
    limits = thermoregulation(organism)

    _, dilated = BiophysicalBehaviour._vasodilate_multipart(organism, limits.flesh_conductivity)
    dphys = physiology(dilated)
    step = limits.flesh_conductivity.step
    @test HeatExchange.conduction_pars_internal(dphys.torso).flesh_conductivity ==
          limits.flesh_conductivity.current + step
    @test HeatExchange.conduction_pars_internal(dphys.head).flesh_conductivity ==
          limits.flesh_conductivity.current + step

    _, swept = BiophysicalBehaviour._sweat_multipart(organism, limits.skin_wetness)
    sphys = physiology(swept)
    wet = limits.skin_wetness.current + limits.skin_wetness.step
    @test HeatExchange.evaporation_pars(sphys.torso).skin_wetness == wet
    @test HeatExchange.evaporation_pars(sphys.head).skin_wetness == wet
end

@testset "multi-part rule-based thermoregulate runs the effector ladder" begin
    organism = two_part_organism()
    init = BiophysicalBehaviour.initial_physiological_state(organism, env_vars)

    result = thermoregulate(organism, environment, init)

    # The loop returns the multipart solve result and reaches (or exhausts toward)
    # the metabolic minimum without erroring.
    @test length(result.parts) == 2
    @test result.metabolic_heat_flow > 0.0u"W"
    @test all(part -> part.success, result.parts)
end

@testset "ThermoregulationLimits: magic-number replacement fields default" begin
    # The four magic-number replacement fields default on ThermoregulationLimits.
    limits = example_thermoregulation_limits()
    @test limits.skin_temperature_undershoot == 5.0u"K"
    @test limits.skin_temperature_core_overshoot == 5.0u"K"
    @test limits.metabolic_heat_flow_max_multiplier == 20.0
    @test limits.minimum_normalisation_range == 1e-6
end

@testset "shape cache: cached solve ≡ uncached, refresh dispatch (Phase 5)" begin
    organism = two_part_organism()

    cache = precompute_shape_cache(organism)
    @test cache isa ShapeCache
    @test keys(cache.parts) == (:torso, :head)
    @test haskey(cache.parts.torso, :total_area)
    @test haskey(cache.parts.torso, :silhouette_area)
    @test haskey(cache.parts.torso, :characteristic_dim)

    uncached = solve_multipart_metabolic_rate(organism, environment, skin0, insul0)
    cached   = solve_multipart_metabolic_rate(organism, environment, skin0, insul0; cache)

    # The cache holds exactly the values the on-the-fly path computes, so results
    # are bitwise identical.
    @test cached.metabolic_heat_flow == uncached.metabolic_heat_flow
    @test cached.skin_temperature == uncached.skin_temperature
    @test all(map((c, u) -> c.net_metabolic == u.net_metabolic, cached.parts, uncached.parts))

    # Physiological effectors leave geometry untouched: refresh returns the same
    # cache object (a no-op the compiler elides).
    @test refresh(cache, organism, Vasodilate()) === cache
    @test refresh(cache, organism, Hyperthermia()) === cache
    @test refresh(cache, organism, Pant()) === cache
    @test refresh(cache, organism, Sweat()) === cache

    # Geometry effectors rebuild the cache (values equal here — geometry unchanged).
    rebuilt = refresh(cache, organism, Piloerect())
    @test rebuilt isa ShapeCache
    @test rebuilt.parts.torso.total_area == cache.parts.torso.total_area
end

@testset "map_part_physiology selector: ByName touches only its target" begin
    organism = two_part_organism()
    before = physiology(organism)

    # Bump only the head's flesh conductivity via a ByName selector.
    bumped = map_part_physiology(ByName(:head), organism) do phys
        @set phys.conduction_pars_internal.flesh_conductivity = 9.9u"W/m/K"
    end
    after = physiology(bumped)

    @test HeatExchange.conduction_pars_internal(after.head).flesh_conductivity == 9.9u"W/m/K"
    @test HeatExchange.conduction_pars_internal(after.torso).flesh_conductivity ==
          HeatExchange.conduction_pars_internal(before.torso).flesh_conductivity
    # The lung index and part set are preserved through the update.
    @test lung_part(bumped) == :torso
    @test keys(physiology(bumped)) == (:torso, :head)
end

@testset "compartment graph: SharedCore contracts, ConductiveCoupling splits" begin
    # Single Body → one compartment (the regulated core).
    shape = Cylinder(1.0u"kg", ρ, b)
    body  = Body(shape, CompositeInsulation(fur, fat))
    single = Organism(body, OrganismTraits(Endotherm(), heat_exchange_of(shape), behav))
    @test couplings(single) == ()
    @test HeatExchange.num_compartments(organism_compartment_graph(single)) == 1

    # Two parts, default couplings → all SharedCore → one shared core.
    shared = two_part_organism()
    @test couplings(shared) == ()
    gshared = organism_compartment_graph(shared)
    @test HeatExchange.num_compartments(gshared) == 1
    @test HeatExchange.compartment_of(gshared, :torso) == HeatExchange.compartment_of(gshared, :head)

    # Two parts joined by a ConductiveCoupling → head floats in its own compartment.
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
    conductive = Organism(dog, OrganismTraits(Endotherm(), phys, behav;
        lung_part = :torso, couplings = (ConductiveCoupling(),)))
    @test couplings(conductive) == (ConductiveCoupling(),)
    gcond = organism_compartment_graph(conductive)
    @test HeatExchange.num_compartments(gcond) == 2
    @test HeatExchange.compartment_of(gcond, :torso) != HeatExchange.compartment_of(gcond, :head)
end

@testset "two HalfCylinders + SharedCore ≡ dorsal/ventral (Phase 7 gate)" begin
    M = 1.0u"kg"
    full_shape = Cylinder(M, ρ, b)
    full = Body(full_shape, CompositeInsulation(fur, fat))
    baseline = solve_metabolic_rate(
        Organism(full, OrganismTraits(Endotherm(), heat_exchange_of(full_shape), behav)),
        environment, skin0, insul0)

    # Two equal-mass half-cylinders joined flat-to-flat reconstruct the full cylinder;
    # SharedCore (the default empty couplings) contracts them into one regulated core.
    hshape = HalfCylinder(M / 2, ρ, b)
    dv = CompositeBody(;
        parts = (; dorsal = Body(hshape, CompositeInsulation(fur, fat)),
                   ventral = Body(hshape, CompositeInsulation(fur, fat))),
        joins = (Join(dorsal  = Attachment(Flat(), FullCover()),
                      ventral = Attachment(Flat(), FullCover())),))
    phys = (; dorsal = heat_exchange_of(hshape), ventral = heat_exchange_of(hshape))
    org = Organism(dv, OrganismTraits(Endotherm(), phys, behav; lung_part = :dorsal))
    @test HeatExchange.num_compartments(organism_compartment_graph(org)) == 1

    result = solve_multipart_metabolic_rate(org, environment, skin0, insul0)

    # The two halves reconstruct the full cylinder's exposed surface exactly; the only
    # residual is the half-cylinder characteristic dimension (2^(1/3) smaller than the
    # full cylinder's) feeding convection, so agreement is ~0.5%, not bitwise.
    @test result.metabolic_heat_flow ≈ baseline.energy_flows.metabolic_heat_flow rtol = 1e-2
    @test result.skin_temperature ≈ baseline.thermoregulation.skin_temperature rtol = 1e-3
    @test result.insulation_temperature ≈ baseline.thermoregulation.insulation_temperature rtol = 1e-3
end

@testset "ConductiveCoupling: floating compartment + coupling limits" begin
    torso_shape = Cylinder(0.6u"kg", ρ, b)
    head_shape  = Cylinder(0.4u"kg", ρ, b)
    r_join = 0.02u"m"
    mkdog() = CompositeBody(;
        parts = (; torso = Body(torso_shape, CompositeInsulation(fur, fat)),
                   head  = Body(head_shape,  CompositeInsulation(fur, fat))),
        joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_join)),
                      head  = Attachment(EndB(0.0u"m", 0.0), Disc(r_join))),),
    )
    phys = (; torso = heat_exchange_of(torso_shape), head = heat_exchange_of(head_shape))
    solve_with(coupling) = begin
        o = Organism(mkdog(), OrganismTraits(Endotherm(), phys, behav;
            lung_part = :torso, couplings = (coupling,)))
        (o, solve_multipart_metabolic_rate(o, environment, skin0, insul0))
    end

    # Derived series-resistance coupling: the head floats in its own compartment.
    o_d, res_d = solve_with(ConductiveCoupling())
    reg = HeatExchange.compartment_of(organism_compartment_graph(o_d), :torso)
    flo = HeatExchange.compartment_of(organism_compartment_graph(o_d), :head)

    @test length(res_d.compartment_core_temperatures) == 2
    @test res_d.parts[2].success
    # The regulated (torso) compartment is pinned at the setpoint core.
    @test res_d.compartment_core_temperatures[reg] ≈ core
    # The floating head core sits above its own skin (it is still a heat source into
    # its surface) — a genuine intermediate node, not pinned.
    head_core_d = res_d.compartment_core_temperatures[flo]
    @test res_d.parts[2].skin_temperature < head_core_d

    # Limit 1 — a very tight interface drags the head core to the regulated core
    # (approaching a SharedCore contraction).
    _, res_tight = solve_with(ConductiveCoupling(1.0e6u"W/m^2/K"))
    @test res_tight.compartment_core_temperatures[flo] ≈ core rtol = 1e-3

    # Limit 2 — a very loose interface leaves the head at its isolated equilibrium,
    # further from the setpoint than the derived coupling.
    _, res_loose = solve_with(ConductiveCoupling(1.0e-4u"W/m^2/K"))
    @test abs(res_loose.compartment_core_temperatures[flo] - core) ≥ abs(head_core_d - core)
end

@testset "six-part dog: many-part solve + thermoregulation ladder" begin
    torso_shape = Cylinder(18.0u"kg", ρ, b)
    head_shape  = Cylinder(2.0u"kg", ρ, 1.5)
    leg_shape   = Cylinder(1.0u"kg", ρ, 5.0)
    torso = Body(torso_shape, CompositeInsulation(fur, fat))
    head  = Body(head_shape,  CompositeInsulation(fur, fat))
    mk_leg() = Body(leg_shape, CompositeInsulation(fur, fat))
    r_head = 0.02u"m"
    r_leg  = 0.015u"m"
    L = torso.geometry.length.length_skin       # axial coordinate for Lateral leg mounts

    dog = CompositeBody(;
        parts = (; torso, head,
                   leg_fl = mk_leg(), leg_fr = mk_leg(),
                   leg_bl = mk_leg(), leg_br = mk_leg()),
        joins = (
            Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head)),
                 head  = Attachment(EndB(0.0u"m", 0.0), Disc(r_head))),
            Join(torso = Attachment(Lateral(0.25L, π/2 + 0.3), Disc(r_leg)),
                 leg_fl = Attachment(EndB(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.25L, π/2 - 0.3), Disc(r_leg)),
                 leg_fr = Attachment(EndB(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.75L, π/2 + 0.3), Disc(r_leg)),
                 leg_bl = Attachment(EndB(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.75L, π/2 - 0.3), Disc(r_leg)),
                 leg_br = Attachment(EndB(0.0u"m", 0.0), Disc(r_leg))),
        ),
    )

    phys = (;
        torso  = heat_exchange_of(torso_shape),
        head   = heat_exchange_of(head_shape),
        leg_fl = heat_exchange_of(leg_shape), leg_fr = heat_exchange_of(leg_shape),
        leg_bl = heat_exchange_of(leg_shape), leg_br = heat_exchange_of(leg_shape),
    )
    traits = OrganismTraits(Endotherm(), phys, behav; lung_part = :torso)
    organism = Organism(dog, traits)

    # Direct solve: one surface result per part, lung mass from the torso.
    solved = solve_multipart_metabolic_rate(organism, environment, skin0, insul0)
    @test length(solved.parts) == 6
    @test all(part -> part.success, solved.parts)
    @test BiophysicalBehaviour._lung_mass(dog, :torso) == 18.0u"kg"

    # Full rule-based ladder runs on the six-part organism.
    init = BiophysicalBehaviour.initial_physiological_state(organism, env_vars)
    regulated = thermoregulate(organism, environment, init)
    @test length(regulated.parts) == 6
    @test regulated.metabolic_heat_flow > 0.0u"W"
end
