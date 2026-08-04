using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Unitful
using Test

const ρ = 1000.0u"kg/m^3"

# Distinct per-part physiology: tell parts apart by their shape mass.
physiology_of(mass) = example_heat_exchange_traits(; shape_pars = Cylinder(mass, ρ, 3.0))

@testset "LungPart wrapper" begin
    phys = physiology_of(5.0u"kg")
    lung = LungPart(phys; panting_capacity = 3.0)

    @test lung isa LungPart
    @test is_lung_part(lung)
    @test !is_lung_part(phys)
    @test panting_capacity(lung) == 3.0
    @test panting_capacity(phys) === nothing          # non-lung parts have no capacity
    @test unwrap_physiology(lung) === phys
    @test unwrap_physiology(phys) === phys             # no-op on a plain part

    # HeatExchange physiology accessors forward transparently through the wrapper
    @test HeatExchange.shape_pars(lung) === HeatExchange.shape_pars(phys)
    @test HeatExchange.metabolism_pars(lung) === HeatExchange.metabolism_pars(phys)
    @test HeatExchange.respiration_pars(lung) === HeatExchange.respiration_pars(phys)
end

@testset "broadcast_physiology — single Body" begin
    body = Body(Cylinder(5.0u"kg", ρ, 3.0), Naked())
    phys = physiology_of(5.0u"kg")
    per_part = broadcast_physiology(body, phys)

    @test keys(per_part) == (:body,)
    @test is_lung_part(per_part.body)                  # sole part hosts the lung
    @test unwrap_physiology(per_part.body) === phys
end

@testset "broadcast_physiology — per-part NamedTuple, distinct physiology" begin
    torso = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    head  = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    L = torso.geometry.length.length_skin
    r_head = 0.04u"m"
    dog = CompositeBody(;
        parts = (; torso, head),
        joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head)),
                      head = Attachment(PoleA(), Disc(r_head))),),
    )

    phys = (; torso = physiology_of(20.0u"kg"), head = physiology_of(2.0u"kg"))
    per_part = broadcast_physiology(dog, phys; lung_part = :torso, panting_capacity = 5.0)

    @test keys(per_part) == (:torso, :head)
    # distinct physiology preserved per part
    @test HeatExchange.shape_pars(unwrap_physiology(per_part.torso)).mass == 20.0u"kg"
    @test HeatExchange.shape_pars(per_part.head).mass == 2.0u"kg"
    # only the lung part is wrapped, and only it carries panting capacity
    @test is_lung_part(per_part.torso)
    @test !is_lung_part(per_part.head)
    @test panting_capacity(per_part.torso) == 5.0
    @test panting_capacity(per_part.head) === nothing
end

@testset "Organism accessors — single Body (backwards compatible)" begin
    body = Body(Cylinder(5.0u"kg", ρ, 3.0), Naked())
    het = physiology_of(5.0u"kg")
    behav = BehavioralTraits(; thermoregulation = example_thermoregulation_limits())
    # The 3-argument constructor still works and defaults lung_part to :body
    traits = OrganismTraits(Endotherm(), het, behav)
    organism = Organism(body, traits)

    @test lung_part(organism) == :body
    @test keys(physiology(organism)) == (:body,)
    @test is_lung_part(part_physiology(organism, :body))   # sole part hosts the lung
    @test lung_physiology(organism) === het                # unwrapped for whole-organism queries
    @test pant_selector(organism) === ByName{(:body,)}()

    # Whole-organism forwarding accessors are unchanged for the single-body path
    @test HeatExchange.metabolism_pars(traits) === HeatExchange.metabolism_pars(het)
    @test HeatExchange.insulation_pars(traits) === HeatExchange.insulation_pars(het)
end

@testset "Organism accessors — multi-part, lung routed to torso" begin
    torso = Body(Cylinder(20u"kg", ρ, 3.0), Naked())
    head  = Body(Ellipsoid(2u"kg", ρ, 1.5, 1.0), Naked())
    r_head = 0.04u"m"
    dog = CompositeBody(;
        parts = (; torso, head),
        joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(r_head)),
                      head = Attachment(PoleA(), Disc(r_head))),),
    )
    phys = (; torso = physiology_of(20.0u"kg"), head = physiology_of(2.0u"kg"))
    behav = BehavioralTraits(; thermoregulation = example_thermoregulation_limits())
    traits = OrganismTraits(Endotherm(), phys, behav; lung_part = :torso)
    organism = Organism(dog, traits)

    @test lung_part(organism) == :torso
    @test keys(physiology(organism)) == (:torso, :head)
    @test is_lung_part(part_physiology(organism, :torso))
    @test !is_lung_part(part_physiology(organism, :head))
    # Pant is routed to the lung part only
    @test pant_selector(organism) === ByName{(:torso,)}()
    # Whole-organism respiration/metabolism resolve through the lung part
    @test HeatExchange.shape_pars(lung_physiology(organism)).mass == 20.0u"kg"
    @test HeatExchange.metabolism_pars(traits) === HeatExchange.metabolism_pars(phys.torso)
end
