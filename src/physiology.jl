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
    # The public keyword API takes a runtime `Symbol`; `Val(lung_part)` is the one
    # barrier that lifts it into the type for external callers. The organism path
    # instead calls `_broadcast_physiology` with the traits' stored `Val` directly.
    return _broadcast_physiology(body, physiology, Val(lung_part), panting_capacity)
end

# Lung part as a `Val` (name in the type), threaded straight to `_wrap_lung` with no
# runtime Symbol — so when the organism path passes the traits' stored `Val`, the
# whole physiology NamedTuple infers concretely.
@inline function _broadcast_physiology(body, physiology, lung::Val, panting_capacity)
    names = part_names(body)
    per_part = _per_part_values(names, physiology)
    return _wrap_lung(per_part, lung, panting_capacity)
end

# Broadcast a lumped physiology across the part names, or pass a per-part NamedTuple through.
@inline _per_part_values(names::Tuple, physiology::NamedTuple) =
    NamedTuple{names}(map(n -> unwrap_physiology(physiology[n]), names))
@inline _per_part_values(names::Tuple, physiology) =
    NamedTuple{names}(ntuple(_ -> physiology, Val(length(names))))

# Wrap only the lung part in LungPart, leaving the rest as plain physiology. The lung
# name is a `Val` type parameter (not a runtime Symbol), so `getfield` and `merge`
# specialise to the known field and the result is a concretely-typed NamedTuple. The
# earlier per-element `name === lung ? ...` compare branched on a runtime Symbol and
# so inferred `Union{LungPart{P},P}` in every slot.
@inline function _wrap_lung(per_part::NamedTuple, ::Val{lung}, panting_capacity) where {lung}
    lung_physiology = LungPart(getfield(per_part, lung), panting_capacity)
    return merge(per_part, NamedTuple{(lung,)}((lung_physiology,)))
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
    # Pass the stored `Val` lung part straight through — no Symbol round-trip — so the
    # result infers. `_broadcast_physiology` dispatches on the stored physiology: a
    # per-part `NamedTuple` (multi-part) is used as-is, a single lumped
    # `HeatExchangeTraits` (single-body) is broadcast across the parts.
    return _broadcast_physiology(body, t.heat_exchange, t.lung_part, nothing)
end

"""
    part_physiology(organism, name) -> physiology

The physiology of the part named `name` (with the wrapper intact if it is the lung
part). `Val(name)` const-propagates for type stability.
"""
part_physiology(o::Organism, name::Symbol) = physiology(o)[name]

"""
    lung_part(organism) -> Symbol
    lung_part(traits::OrganismTraits) -> Symbol

The name of the lung-hosting part. Stored as a `Val` in the traits (the name is a
type parameter), so this returns a compile-time-constant `Symbol` and the routing
built on it stays type-stable.
"""
@inline _lung_symbol(::Val{N}) where {N} = N
lung_part(t::OrganismTraits) = _lung_symbol(t.lung_part)
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
