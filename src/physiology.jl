# Per-part physiology, the lung-part wrapper, and the lung-part index (plan §3.9, §4).
#
# An organism's physiology is a `NamedTuple` keyed by part name — one physiology
# object per part. A single `Body` broadcasts its one lumped `HeatExchangeTraits`
# into a one-entry NamedTuple keyed by `:body`, so every consumer writes the
# multi-part form and the single-body case falls out for free.
#
# Respiratory heat/mass exchange happens at the lung surface, which sits in one
# part (the torso, not the head). That part's physiology is wrapped in `LungPart`,
# which alone carries `panting_capacity`; non-lung parts have plain physiology and
# no panting field. `OrganismTraits.lung_part::Symbol` is the O(1) index into the
# physiology NamedTuple identifying it.

"""
    LungPart(physiology; panting_capacity = nothing)
    LungPart{Physiology,PantingCapacity}

Wraps the lung-hosting part's physiology, adding a `panting_capacity` field that
lives *only* on this wrapper — non-lung parts carry plain physiology. Dispatching
on `::LungPart` routes respiration, pulmonary evaporation, and O₂/CO₂ exchange to
the correct part. All `HeatExchange` physiology accessors forward transparently to
the wrapped physiology, so a `LungPart` is a drop-in for a plain physiology
everywhere except where the panting capacity or lung role matters.
"""
struct LungPart{Physiology,PantingCapacity}
    physiology::Physiology
    panting_capacity::PantingCapacity
end
LungPart(physiology; panting_capacity=nothing) = LungPart(physiology, panting_capacity)

# Peel the wrapper off to reach the underlying physiology; a no-op for plain parts.
unwrap_physiology(part::LungPart) = part.physiology
unwrap_physiology(part) = part

"""
    panting_capacity(part) -> capacity | nothing

The panting capacity of a part: the value stored on a `LungPart`, or `nothing` for
any part that does not host the lung.
"""
panting_capacity(part::LungPart) = part.panting_capacity
panting_capacity(::Any) = nothing

"""
    is_lung_part(part) -> Bool

Whether `part` is the lung-hosting part (i.e. a `LungPart`).
"""
is_lung_part(::LungPart) = true
is_lung_part(::Any) = false

# Forward every HeatExchange physiology accessor through the LungPart wrapper.
for accessor in (
    :shape_pars, :insulation_pars,
    :conduction_pars_external, :conduction_pars_internal,
    :convection_pars, :radiation_pars, :evaporation_pars,
    :hydraulic_pars, :respiration_pars, :metabolism_pars, :options,
)
    @eval HeatExchange.$accessor(part::LungPart) = HeatExchange.$accessor(part.physiology)
end

# ── Per-part physiology container ───────────────────────────────────────────

# A single lumped physiology is already the whole-organism physiology; treat it as
# the per-part physiology of the one part. A NamedTuple is already per-part.
_is_per_part_physiology(::NamedTuple) = true
_is_per_part_physiology(::Any) = false

"""
    broadcast_physiology(body, physiology; lung_part = SINGLE_PART_NAME,
                         panting_capacity = nothing) -> NamedTuple

Build the per-part physiology `NamedTuple` keyed by `body`'s part names.

- A single lumped `physiology` (a `HeatExchangeTraits`) is broadcast to every part;
  a plain `Body` yields `(; body = LungPart(physiology))`.
- A `NamedTuple` `physiology` (already per-part) is used as given.

Either way the `lung_part` entry is wrapped in `LungPart` (carrying
`panting_capacity`) and every other entry is left as plain physiology.
"""
function broadcast_physiology(body, physiology;
                              lung_part::Symbol=SINGLE_PART_NAME,
                              panting_capacity=nothing)
    names = part_names(body)
    per_part = _per_part_values(names, physiology)
    return _wrap_lung(names, per_part, lung_part, panting_capacity)
end

# Broadcast a lumped physiology across the part names, or pass a per-part NamedTuple through.
@inline _per_part_values(names::Tuple, physiology::NamedTuple) =
    NamedTuple{names}(map(n -> unwrap_physiology(physiology[n]), names))
@inline _per_part_values(names::Tuple, physiology) =
    NamedTuple{names}(ntuple(_ -> physiology, Val(length(names))))

# Wrap only the lung part in LungPart, leaving the rest as plain physiology.
@inline function _wrap_lung(names::Tuple, per_part::NamedTuple, lung_part::Symbol, panting_capacity)
    wrapped = ntuple(Val(length(names))) do i
        name = names[i]
        part = per_part[i]
        name === lung_part ? LungPart(part, panting_capacity) : part
    end
    return NamedTuple{names}(wrapped)
end

# ── Accessors on an organism / OrganismTraits ───────────────────────────────

"""
    physiology(organism) -> NamedTuple

The per-part physiology, keyed by part name, with the lung part wrapped in
`LungPart`. Broadcast on demand from a single lumped `heat_exchange`, or returned
directly when per-part physiology was supplied at construction.
"""
physiology(o::Organism) = physiology(HeatExchange.traits(o), HeatExchange.body(o))
function physiology(t::OrganismTraits, body)
    return _organism_physiology(t.heat_exchange, body, lung_part(t))
end
# Per-part physiology already stored (multi-part construction): use as-is.
_organism_physiology(physiology::NamedTuple, body, lung) =
    broadcast_physiology(body, physiology; lung_part=lung)
# Single lumped physiology (single-body construction): broadcast across parts.
_organism_physiology(physiology, body, lung) =
    broadcast_physiology(body, physiology; lung_part=lung)

"""
    part_physiology(organism, name) -> physiology

The physiology of the part named `name` (with the wrapper intact if it is the lung
part). `Val(name)` const-propagates for type stability.
"""
part_physiology(o::Organism, name::Symbol) = physiology(o)[name]

"""
    lung_part(organism) -> Symbol
    lung_part(traits::OrganismTraits) -> Symbol

The name of the lung-hosting part.
"""
lung_part(t::OrganismTraits) = t.lung_part
lung_part(o::Organism) = lung_part(HeatExchange.traits(o))

"""
    lung_physiology(organism) -> physiology

The unwrapped physiology of the lung part — the source for whole-organism
respiration and metabolism queries and for `lung_temperature`.
"""
lung_physiology(o::Organism) = unwrap_physiology(part_physiology(o, lung_part(o)))

"""
    pant_selector(organism) -> ByName

The default selector for `Pant`: the lung part only (§3.9). For a single `Body`
this is `ByName{(:body,)}`, equivalent to `WholeBody()`.
"""
pant_selector(o::Organism) = ByName{(lung_part(o),)}()
