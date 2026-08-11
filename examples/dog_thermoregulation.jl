# dog_thermoregulation.jl
#
# Multi-part endotherm thermoregulation end-to-end (plan Phase 9).
#
# Builds a six-part dog — torso (hosting the lung), head, and four legs, each a
# fur-insulated cylinder — as a CompositeBody, gives each part its own
# physiology, and runs the rule-based sequential controller. All parts share one
# regulated core (the all-SharedCore compartment); respiration is closed at the
# torso lung and panting affects only that part (§3.9).
#
# Run from the repository root:
#   julia --project=. examples/dog_thermoregulation.jl

using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Unitful, UnitfulMoles

const ρ = 1000.0u"kg/m^3"

# ── Insulation shared across parts (a uniform pelt) ─────────────────────────
ins_pars = example_insulation_pars()
fibre    = ins_pars.dorsal
fur      = FibrousLayer(fibre.depth, fibre.diameter, fibre.density)
fat      = FatLayer(0.0, 901.0u"kg/m^3")

# Per-part physiology keyed by anatomical name — distinguished here by shape mass.
heat_exchange_of(shape) = example_heat_exchange_traits(;
    shape_pars               = shape,
    insulation_pars          = ins_pars,
    conduction_pars_external = example_conduction_pars_external(; conduction_fraction = 0.0),
    conduction_pars_internal = example_conduction_pars_internal(),
    radiation_pars           = example_radiation_pars(),
    evaporation_pars         = example_evaporation_pars(),
    respiration_pars         = example_respiration_pars(),
    metabolism_pars          = example_metabolism_pars(),
)

# ── Geometry: torso + head + four legs ──────────────────────────────────────
torso_shape = Cylinder(18.0u"kg", ρ, 3.0)
head_shape  = Cylinder(2.0u"kg", ρ, 1.5)
leg_shape   = Cylinder(1.0u"kg", ρ, 5.0)

insulated(shape) = Body(shape, CompositeInsulation(fur, fat))
torso = insulated(torso_shape)
L     = torso.geometry.length.length_skin
r_head, r_leg = 0.02u"m", 0.015u"m"

dog = CompositeBody(;
    parts = (; torso,
               head   = insulated(head_shape),
               leg_fl = insulated(leg_shape), leg_fr = insulated(leg_shape),
               leg_bl = insulated(leg_shape), leg_br = insulated(leg_shape)),
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

physiology_nt = (;
    torso  = heat_exchange_of(torso_shape),
    head   = heat_exchange_of(head_shape),
    leg_fl = heat_exchange_of(leg_shape), leg_fr = heat_exchange_of(leg_shape),
    leg_bl = heat_exchange_of(leg_shape), leg_br = heat_exchange_of(leg_shape),
)
traits = OrganismTraits(
    Endotherm(), physiology_nt,
    BehavioralTraits(; thermoregulation = example_thermoregulation_limits());
    lung_part = :torso,
)
organism = Organism(dog, traits)

# ── Environment + solve ─────────────────────────────────────────────────────
environment = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(),
)
init = BiophysicalBehaviour.initial_physiological_state(organism, environment.environment_vars)

result = thermoregulate(organism, environment, init)

println("Six-part dog — lung part: ", lung_part(organism))
println("  metabolic heat flow : ", result.metabolic_heat_flow)
println("  mean skin temperature: ", result.skin_temperature)
println("  lung temperature     : ", result.lung_temperature)
println("  per-part net metabolic heat:")
for (name, part) in pairs(NamedTuple{part_names(dog)}(result.parts))
    println("    ", rpad(name, 7), " skin=", round(u"K", part.skin_temperature; digits = 2),
            "  net_metabolic=", round(u"W", part.net_metabolic; digits = 3))
end
