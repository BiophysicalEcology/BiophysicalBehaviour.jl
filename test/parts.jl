using BiophysicalBehaviour
using BiophysicalGeometry
using Unitful
using Test

const ρ = 1000.0u"kg/m^3"

# A small composite: torso + head + two legs.
function build_dog()
    torso = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    leg   = Body(Cylinder(1u"kg", ρ, 5.0), Naked())
    head  = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L = torso.geometry.length.length_skin
    r_leg = skin_radius(leg)
    r_head = 0.04u"m"
    CompositeBody(;
        parts = (; torso, head, leg_fl = leg, leg_fr = leg),
        joins = (
            Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head)),
                 head = Attachment(PoleA(), Disc(r_head))),
            Join(torso = Attachment(Lateral(0.2*L, π/2), Disc(r_leg)),
                 leg_fl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
            Join(torso = Attachment(Lateral(0.8*L, π/2), Disc(r_leg)),
                 leg_fr = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
        ),
    )
end

@testset "part_names / select_names" begin
    dog = build_dog()
    @test part_names(dog) == (:torso, :head, :leg_fl, :leg_fr)

    body = Body(Sphere(2u"kg", ρ), Naked())
    @test part_names(body) == (:body,)

    @test select_names(WholeBody(), dog) == (:torso, :head, :leg_fl, :leg_fr)
    @test select_names(ByName(:head), dog) == (:head,)
    @test select_names(ByName(:leg_fl, :leg_fr), dog) == (:leg_fl, :leg_fr)
    @test select_names(Compartment(:torso), dog) == (:torso,)
    @test select_names(WholeBody(), body) == (:body,)
end

@testset "map_parts round-trips" begin
    dog = build_dog()
    # Identity over the whole body reproduces geometry.
    same = map_parts((n, p) -> p, WholeBody(), dog)
    @test total_area(same) ≈ total_area(dog)
    @test part_names(same) == part_names(dog)

    # Selecting a single part touches only that part. Swap the head for a
    # bigger head; torso/legs identical, head changed.
    big_head = Body(Ellipsoid(3u"kg", ρ, 1.5, 1.0), Naked())
    changed = map_parts((n, p) -> big_head, ByName(:head), dog)
    @test changed.parts.head === big_head
    @test changed.parts.torso === dog.parts.torso
    @test changed.parts.leg_fl === dog.parts.leg_fl

    # Single Body: WholeBody applies f; empty selection returns it unchanged.
    body = Body(Sphere(2u"kg", ρ), Naked())
    bigger = Body(Sphere(3u"kg", ρ), Naked())
    @test map_parts((n, p) -> bigger, WholeBody(), body) === bigger
    @test map_parts((n, p) -> bigger, ByName(:nope), body) === body
end

@testset "foldl_parts" begin
    dog = build_dog()
    # Sum flesh volume over selected parts.
    total = foldl_parts((acc, n, p) -> acc + flesh_volume(p), 0.0u"m^3", WholeBody(), dog)
    @test total ≈ sum(flesh_volume(p) for p in values(dog.parts))

    legs = foldl_parts((acc, n, p) -> acc + 1, 0, ByName(:leg_fl, :leg_fr), dog)
    @test legs == 2

    # Collect selected names in order.
    names = foldl_parts((acc, n, p) -> (acc..., n), (), WholeBody(), dog)
    @test names == (:torso, :head, :leg_fl, :leg_fr)
end

@testset "set_part" begin
    dog = build_dog()
    big_head = Body(Ellipsoid(3u"kg", ρ, 1.5, 1.0), Naked())
    updated = set_part(dog, :head, big_head)
    @test updated.parts.head === big_head
    @test updated.parts.torso === dog.parts.torso
    @test part_names(updated) == part_names(dog)      # ordering preserved

    # Val form is equivalent.
    updated2 = set_part(dog, Val(:head), big_head)
    @test updated2.parts.head === big_head

    # Single Body: replacing :body yields the new body.
    body = Body(Sphere(2u"kg", ρ), Naked())
    bigger = Body(Sphere(3u"kg", ρ), Naked())
    @test set_part(body, :body, bigger) === bigger
end

@testset "traversal type stability" begin
    dog = build_dog()
    # The lightweight primitives that don't rebuild the body are fully
    # inferred (this is what the Phase 5 effector hot loop relies on).
    @inferred select_names(WholeBody(), dog)
    @inferred foldl_parts((acc, n, p) -> acc + flesh_volume(p), 0.0u"m^3", WholeBody(), dog)

    # On a plain Body, map_parts/set_part don't reconstruct, so they infer.
    body = Body(Sphere(2u"kg", ρ), Naked())
    bigger = Body(Sphere(3u"kg", ρ), Naked())
    @inferred map_parts((n, p) -> bigger, WholeBody(), body)
    @inferred set_part(body, Val(:body), bigger)

    # NOTE: map_parts / set_part over a CompositeBody go through the
    # CompositeBody constructor, which re-solves poses via a merge-fold and
    # is not currently inference-stable. Phase 5's zero-alloc effector loop
    # must update parts/caches without a full body reconstruction rather than
    # rely on rebuilding the composite each iteration.
    @test_broken (@inferred(map_parts((n, p) -> p, WholeBody(), dog)); true)
end
