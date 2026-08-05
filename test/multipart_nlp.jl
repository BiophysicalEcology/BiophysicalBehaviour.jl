# Multi-part NLP (`IPOPTControl` + `MultipartNLP`) — Phase 7.5.
#
# Gates:
#  1. The multi-part NLP residual function vanishes at the iterative coupled
#     solve's solution (physics-consistency: the NLP constraints and the
#     iterative solver agree at a converged point).
#  2. A 1-part organism solves via IPOPT to a physically sane, heat-balanced
#     state (residuals ≈ 0 at the returned solution, Q10 satisfied).
#  3. A two-part organism solves; every part is heat-balanced and physically
#     sane, respiration is closed once whole-organism.

using BiophysicalBehaviour
using BiophysicalGeometry
using HeatExchange
using Unitful, UnitfulMoles
using Test
import Flatten
const BB = BiophysicalBehaviour

const ρ = 1000.0u"kg/m^3"
const b = 3.0

ins_pars = example_insulation_pars()
fibre    = ins_pars.dorsal
fur      = FibrousLayer(fibre.depth, fibre.diameter, fibre.density)
fat      = FatLayer(0.0, 901.0u"kg/m^3")
internal = example_conduction_pars_internal()
external = example_conduction_pars_external(; conduction_fraction = 0.0)
evap     = example_evaporation_pars()
rad_pars = example_radiation_pars()
resp     = example_respiration_pars()
metab    = example_metabolism_pars()
core     = metab.core_temperature

env_vars = example_environment_vars()
env_pars = example_environment_pars()
environment = (; environment_pars = env_pars, environment_vars = env_vars)
behav = BehavioralTraits(; thermoregulation = example_thermoregulation_limits())

heat_exchange_of(shape) = example_heat_exchange_traits(;
    shape_pars = shape, insulation_pars = ins_pars,
    conduction_pars_external = external, conduction_pars_internal = internal,
    radiation_pars = rad_pars, evaporation_pars = evap,
    respiration_pars = resp, metabolism_pars = metab)

skin0  = core - 5u"K"
insul0 = env_vars.air_temperature + 2u"K"

control = IPOPTControl(; nlp_strategy = MultipartNLP(), smoothing = HeatExchange.SmoothBound(1e-5))
limits  = example_thermoregulation_limits()

# Rebuild the flat decision-variable vector from a solution by flattening a
# variable structure of the same shape the packer uses — no index arithmetic.
function pack_variables(packed, core_temperature, metabolic_heat_flow, panting_rate, part_leaves)
    # Build a Unitful structure of the template's shape and let `flatten(_, Real)`
    # rip out the Float64 magnitudes (in the template's canonical units) — the
    # inverse of the packer's `reconstruct`. `log(ustrip(u"W", _))` is the targeted
    # strip the log-transform variable genuinely needs; no other unit stripping.
    structure = BB._variable_structure(packed.part_names, core_temperature,
        log(ustrip(u"W", metabolic_heat_flow)), Float64(panting_rate), part_leaves)
    return collect(Float64, Flatten.flatten(structure, Real))
end

function pack_solution(packed, out)
    part_leaves = ntuple(length(out.parts)) do i
        p = values(out.parts)[i]
        BB._part_leaf(p.skin_temperature, p.insulation_temperature, p.flesh_conductivity, p.skin_wetness)
    end
    return pack_variables(packed, out.core_temperature, out.metabolic_heat_flow, out.panting_rate, part_leaves)
end

function residuals_at(out, organism)
    packed = HeatExchange.nlp_pack(MultipartNLP(), organism, environment, skin0, insul0;
                                   smoothing = control.smoothing)
    x = pack_solution(packed, out)
    _, ncon = BB._problem_size(packed)
    params = (; nlp_packed = packed,
                minimum_heat_flow = ustrip(u"W", limits.minimum_heat_flow),
                q10 = Float64(metab.q10),
                setpoint_temperature = ustrip(u"K", core))
    r = zeros(ncon)
    BB._heat_balance_residuals!(packed, r, x, params)
    return r, params
end

@testset "MultipartNLP residuals vanish at the iterative coupled solution" begin
    shape  = Cylinder(1.0u"kg", ρ, b)
    body   = Body(shape, CompositeInsulation(fur, fat))
    org    = Organism(body, OrganismTraits(Endotherm(), heat_exchange_of(shape), behav))

    ref    = solve_multipart_metabolic_rate(org, environment, skin0, insul0)
    packed = HeatExchange.nlp_pack(MultipartNLP(), org, environment, skin0, insul0)

    N = length(packed.setups)
    _, ncon = BB._problem_size(packed)
    part_leaves = ntuple(N) do i
        p = values(ref.parts)[i]
        BB._part_leaf(p.skin_temperature, p.insulation_temperature,
                      internal.flesh_conductivity, Float64(evap.skin_wetness))
    end
    x = pack_variables(packed, core, ref.metabolic_heat_flow, resp.pant, part_leaves)
    params = (; nlp_packed = packed, minimum_heat_flow = ustrip(u"W", metab.metabolic_heat_flow),
                q10 = Float64(metab.q10), setpoint_temperature = ustrip(u"K", core))
    r = zeros(ncon)
    BB._heat_balance_residuals!(packed, r, x, params)

    # Equality residuals (surface balances, skin temps, whole-organism balance) vanish.
    @test maximum(abs, r[1:2N+1]) < 1e-1
end

@testset "MultipartNLP — 1-part organism solves to a heat-balanced state" begin
    shape  = Cylinder(1.0u"kg", ρ, b)
    body   = Body(shape, CompositeInsulation(fur, fat))
    org    = Organism(body, OrganismTraits(Endotherm(), heat_exchange_of(shape), behav))
    init   = BB.initial_physiological_state(org, env_vars)

    out = thermoregulate(Endotherm(), control, org, environment, init)

    @test length(out.parts) == 1
    # Physical sanity: core between setpoint and its ceiling, skin between air and core.
    @test core <= out.core_temperature <= limits.core_temperature.max + 1e-6u"K"
    @test env_vars.air_temperature < out.parts[1].skin_temperature < out.core_temperature
    @test out.metabolic_heat_flow > 0.0u"W"

    # Heat-balance residuals vanish at the returned NLP solution.
    r, params = residuals_at(out, org)
    N = length(out.parts)
    @test maximum(abs, r[1:2N+1]) < 1e-1
    # Q10 inequality satisfied (metabolic ≥ minimum · q10^((core−setpoint)/10)).
    @test r[2N+2] > -1e-3
end

@testset "MultipartNLP — two-part organism solves, respiration closed once" begin
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
    org  = Organism(dog, OrganismTraits(Endotherm(), phys, behav; lung_part = :torso))
    init = BB.initial_physiological_state(org, env_vars)

    out = thermoregulate(Endotherm(), control, org, environment, init)

    @test length(out.parts) == 2
    @test keys(out.parts) == (:torso, :head)
    for p in values(out.parts)
        @test env_vars.air_temperature < p.skin_temperature < out.core_temperature
        @test p.insulation_temperature < p.skin_temperature
    end
    @test out.metabolic_heat_flow > 0.0u"W"

    # Every part is heat-balanced and the whole-organism respiration balance closes.
    r, params = residuals_at(out, org)
    N = length(out.parts)
    @test maximum(abs, r[1:2N+1]) < 1e-1
    @test r[2N+2] > -1e-3
end

@testset "MultipartNLP — IPOPTSolverCache warm-start reproduces the fresh solve" begin
    torso_shape = Cylinder(0.6u"kg", ρ, b)
    head_shape  = Cylinder(0.4u"kg", ρ, b)
    torso = Body(torso_shape, CompositeInsulation(fur, fat))
    head  = Body(head_shape,  CompositeInsulation(fur, fat))
    dog = CompositeBody(;
        parts = (; torso, head),
        joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(0.02u"m")),
                      head  = Attachment(EndB(0.0u"m", 0.0), Disc(0.02u"m"))),),
    )
    phys = (; torso = heat_exchange_of(torso_shape), head = heat_exchange_of(head_shape))
    org  = Organism(dog, OrganismTraits(Endotherm(), phys, behav; lung_part = :torso))
    init = BB.initial_physiological_state(org, env_vars)

    fresh = thermoregulate(Endotherm(), control, org, environment, init)
    cache = IPOPTSolverCache(control, org, environment, init)
    warm1 = thermoregulate(Endotherm(), control, org, environment, init; cache)  # cold in cache
    warm2 = thermoregulate(Endotherm(), control, org, environment, init; cache)  # primal+dual warm-started

    @test warm1.metabolic_heat_flow ≈ fresh.metabolic_heat_flow rtol = 1e-3
    @test warm2.metabolic_heat_flow ≈ fresh.metabolic_heat_flow rtol = 1e-3
    @test warm2.core_temperature ≈ fresh.core_temperature rtol = 1e-4
end
