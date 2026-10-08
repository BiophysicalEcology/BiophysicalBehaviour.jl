# A grid of thermoregulation runs, to pin the current results while the solve is rewritten.
# `generate_reference.jl` writes them to data/reference_thermoregulation.csv.
#
# The includer defines `cylinder(mass, axis_ratio_b)` and `half_cylinder(mass, axis_ratio_b)`:
# the reference was generated with the positional shape constructors that came before the
# keyword ones.

using BiophysicalBehaviour, BiophysicalGeometry, HeatExchange, Unitful
using Setfield: @set

const BB = BiophysicalBehaviour

reference_fur() = (p = example_insulation_pars().dorsal; FibrousLayer(p.depth, p.diameter, p.density))
reference_body(shape) = Body(shape, CompositeInsulation(reference_fur(), FatLayer(0.1, 901.0u"kg/m^3")))
reference_physiology(shape) = example_heat_exchange_traits(;
    shape_pars = shape,
    conduction_pars_internal = example_conduction_pars_internal(; fat_fraction = 0.1),
)
const REFERENCE_BEHAVIOUR = BehavioralTraits(; thermoregulation = example_thermoregulation_limits())

function reference_organism(name)
    if name == :body
        shape = cylinder(1.0u"kg", 3.0)
        Organism(reference_body(shape), OrganismTraits(Endotherm(), reference_physiology(shape), REFERENCE_BEHAVIOUR))
    elseif name == :halves
        shape = half_cylinder(0.5u"kg", 3.0)
        flat = Attachment(Flat(), FullCover())
        body = CompositeBody(; parts = (; dorsal = reference_body(shape), ventral = reference_body(shape)),
            joins = (Join(dorsal = flat, ventral = flat),))
        physiology = (; dorsal = reference_physiology(shape), ventral = reference_physiology(shape))
        Organism(body, OrganismTraits(Endotherm(), physiology, REFERENCE_BEHAVIOUR; lung_part = :dorsal))
    else
        torso, head = cylinder(0.6u"kg", 3.0), cylinder(0.4u"kg", 3.0)
        body = CompositeBody(; parts = (; torso = reference_body(torso), head = reference_body(head)),
            joins = (Join(torso = Attachment(EndA(0.0u"m", 0.0), Disc(0.02u"m")),
                head = Attachment(EndB(0.0u"m", 0.0), Disc(0.02u"m"))),))
        physiology = (; torso = reference_physiology(torso), head = reference_physiology(head))
        couplings = name == :conductive ? (ConductiveCoupling(),) : ()
        Organism(body, OrganismTraits(Endotherm(), physiology, REFERENCE_BEHAVIOUR; lung_part = :torso, couplings))
    end
end

const REFERENCE_CASES = [
    (; organism, control, air, sun, wind)
    for organism in (:body, :halves, :shared, :conductive), control in (:rulebased, :ipopt),
        air in (5.0, 25.0, 40.0), sun in (false, true), wind in (0.1, 3.0)
]

function reference_thermoregulate(case)
    organism = reference_organism(case.organism)
    environment_vars = example_environment_vars(;
        air_temperature = u"K"(case.air * u"°C"),
        wind_speed = case.wind * u"m/s",
        global_radiation = case.sun ? 600.0u"W/m^2" : 0.0u"W/m^2",
    )
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    init = BB.initial_physiological_state(organism, environment_vars)
    control = case.control == :ipopt ? IPOPTControl() : BB.control_strategy(organism)
    thermoregulate(Endotherm(), control, organism, environment, init)
end

# Every number in a result, in SI units, keyed by its path.
function reference_values(out)
    values = Pair{String,Float64}[]
    _collect!(values, "", out, 0)
    values
end
_collect!(values, path, x::Unitful.Quantity, depth) = push!(values, path => Float64(ustrip(upreferred(x))))
_collect!(values, path, x::Real, depth) = push!(values, path => Float64(x))
_collect!(values, path, x::Union{Nothing,Symbol,AbstractString}, depth) = values
function _collect!(values, path, x, depth)
    depth > 4 && return values
    for (i, k) in enumerate(propertynames(x))
        key = x isa Tuple ? string(i) : string(k)
        _collect!(values, isempty(path) ? key : string(path, ".", key), getproperty(x, k), depth + 1)
    end
    values
end
