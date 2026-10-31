# Part selectors and traversal primitives for multi-part organisms.
#
# A `PartSelector` names which anatomical parts a behavior targets — the
# whole body, a named subset, or a compartment. The traversal primitives
# (`map_parts`, `foldl_parts`, `set_part`) walk a body's parts, applying a
# transform only to the selected ones and leaving the rest untouched.
#
# A single `Body` is treated as a one-part composite whose sole part is
# named `:body`, so every downstream consumer writes the multi-part API and
# the single-body case falls out for free.
#
# Traversal is type-stable: part names live in the parts `NamedTuple`'s type
# and selector names in the selector's type parameter, so membership is a
# compile-time boolean and per-part values keep their concrete types.

# ── Selectors ─────────────────────────────────────────────────────────────

"""
    PartSelector

Supertype for part selectors. A selector, given a body, resolves to the
tuple of part names it targets (`select_names`).
"""
abstract type PartSelector end

"""
    WholeBody() <: PartSelector

Select every part. Degenerates to the single-body API on a plain `Body`.
"""
struct WholeBody <: PartSelector end

"""
    ByName{Names}() <: PartSelector
    ByName(:head, :torso)

Select the parts named in `Names` (a tuple of `Symbol`s).
"""
struct ByName{Names} <: PartSelector end
ByName(names::Symbol...) = ByName{names}()

"""
    Compartment{Name}() <: PartSelector
    Compartment(:torso_core)

Select the parts belonging to a thermal compartment. Until the compartment
graph lands (Phase 7) a compartment is identified with the like-named part;
`Vasodilate` / `Hyperthermia` use this selector.
"""
struct Compartment{Name} <: PartSelector end
Compartment(name::Symbol) = Compartment{name}()

# ── Effectors ─────────────────────────────────────────────────────────────
#
# Thermoregulatory effectors are singleton tags. A behavior is applied via
# `effect(op, selector, organism, limits)`, which resolves the selector to
# the targeted parts and folds the per-target update back into a new
# organism and new limits. The old scalar entry points (`piloerect`,
# `uncurl`, …) forward to `effect(op, WholeBody(), …)`.
#
# Scope (§3.10): `Piloerect` / `Uncurl` / `Sweat` are per-part;
# `Vasodilate` / `Hyperthermia` are per-compartment; `Pant` is routed to the
# lung part by default. Until per-part physiology lands (Phase 4) every
# effector operates on the whole single body.

"""
    Effector

Supertype for thermoregulatory effector tags. Dispatch target for
`effect` / `refresh!`.
"""
abstract type Effector end

struct Piloerect    <: Effector end
struct Uncurl       <: Effector end
struct Vasodilate   <: Effector end
struct Hyperthermia <: Effector end
struct Pant         <: Effector end
struct Sweat        <: Effector end

"""
    effect(op::Effector, selector::PartSelector, organism, limits) -> (new_limits, new_organism)

Apply effector `op` to the parts picked by `selector`, returning updated
limits and organism. Methods are defined per control strategy (rule-based in
`rulebased.jl`). The default selector for the scalar entry points is
`WholeBody()`.
"""
function effect end

# ── Body → parts normalisation ────────────────────────────────────────────

# The canonical name for the sole part of a plain `Body`. `lung_part` and
# single-part organism construction (§3.9) use the same convention.
const SINGLE_PART_NAME = :body

"""
    part_names(body) -> Tuple{Vararg{Symbol}}

The part names of a body, in order. A plain `Body` reports `(:body,)`.
"""
part_names(b::CompositeBody) = keys(b.parts)
part_names(::Body) = (SINGLE_PART_NAME,)

# Parts as a NamedTuple, keyed by name. A plain `Body` becomes a one-entry
# NamedTuple so the traversal primitives share one code path.
_parts(b::CompositeBody) = b.parts
_parts(b::Body) = NamedTuple{(SINGLE_PART_NAME,)}((b,))

# ── Selection resolution ──────────────────────────────────────────────────

"""
    select_names(selector, body) -> Tuple{Vararg{Symbol}}

Resolve a selector to the concrete tuple of targeted part names for `body`.
"""
select_names(::WholeBody, body) = part_names(body)
select_names(::ByName{Names}, body) where {Names} = Names
select_names(::Compartment{Name}, body) where {Name} = (Name,)

# Compile-time membership test, unrolled over the small selector tuple.
@inline _member(name::Symbol, ::Tuple{}) = false
@inline _member(name::Symbol, sel::Tuple) = name === sel[1] || _member(name, Base.tail(sel))

# ── map_parts ─────────────────────────────────────────────────────────────

"""
    map_parts(f, selector, body) -> body

Return a body of the same kind with `f(name, part)` substituted for each
selected part and every other part left unchanged. For a `CompositeBody`
the result is rebuilt through the constructor (poses re-solved); for a
`Body` the sole part is returned directly.
"""
function map_parts(f, selector::PartSelector, b::CompositeBody)
    sel = select_names(selector, b)
    newparts = _map_over(f, sel, b.parts)
    return CompositeBody(; parts = newparts, joins = b.joins)
end
function map_parts(f, selector::PartSelector, b::Body)
    sel = select_names(selector, b)
    return _member(SINGLE_PART_NAME, sel) ? f(SINGLE_PART_NAME, b) : b
end

# Build a new parts NamedTuple, applying f only to selected entries. Uses
# `ntuple(_, Val(N))` so each index is a compile-time literal — names and
# values keep their concrete types and no NamedTuple field is looked up by
# a runtime Symbol.
@inline function _map_over(f, sel::Tuple, parts::NamedTuple)
    names = keys(parts)
    vals = values(parts)
    newvals = ntuple(Val(length(names))) do i
        name = names[i]
        part = vals[i]
        _member(name, sel) ? f(name, part) : part
    end
    return NamedTuple{names}(newvals)
end

# ── foldl_parts ───────────────────────────────────────────────────────────

"""
    foldl_parts(f, init, selector, body) -> accumulator

Left-fold `f(accumulator, name, part)` over the selected parts, in order.
Unselected parts are skipped. Type-stable via tuple recursion.
"""
function foldl_parts(f, init, selector::PartSelector, body)
    sel = select_names(selector, body)
    parts = _parts(body)
    return _foldl_parts(f, init, sel, keys(parts), values(parts))
end

@inline _foldl_parts(f, acc, sel, ::Tuple{}, ::Tuple{}) = acc
@inline function _foldl_parts(f, acc, sel, names::Tuple, vals::Tuple)
    name = names[1]
    part = vals[1]
    acc = _member(name, sel) ? f(acc, name, part) : acc
    return _foldl_parts(f, acc, sel, Base.tail(names), Base.tail(vals))
end

# ── set_part ──────────────────────────────────────────────────────────────

"""
    set_part(body, name, new_part) -> body
    set_part(body, Val(name), new_part) -> body

Replace the part named `name` with `new_part`, returning a body of the same
kind. The `Val(name)` form is fully type-stable; the `Symbol` form is a
convenience that is stable when `name` const-propagates.
"""
set_part(b::AbstractBody, name::Symbol, new_part::AbstractBody) =
    set_part(b, Val(name), new_part)

@inline function set_part(b::CompositeBody, ::Val{N}, new_part::AbstractBody) where {N}
    newparts = merge(b.parts, NamedTuple{(N,)}((new_part,)))
    return CompositeBody(; parts = newparts, joins = b.joins)
end
# A plain `Body` has one part (`:body`); replacing it yields the new body.
@inline set_part(::Body, ::Val{SINGLE_PART_NAME}, new_part::AbstractBody) = new_part

# ── Compartment graph ───────────────────────────────────────────────────────

"""
    organism_compartment_graph(o::Organism) -> HeatExchange.CompartmentGraph

Build the thermal-compartment partition for `o` from its body topology and
`couplings(o)`. Parts joined by a `SharedCore` coupling equilibrate into one
compartment (sharing a single `core_temperature`); every other part is its own
compartment.

A plain `Body`, or a `CompositeBody` whose joins are all `SharedCore` (the
default when `couplings(o)` is empty), yields a single compartment — the one
regulated core that reproduces the single-body / dorsal-ventral behaviour.
"""
organism_compartment_graph(o::Organism) = _compartment_graph(HeatExchange.body(o), couplings(o))

# A plain `Body` is one part, one compartment, regardless of couplings.
_compartment_graph(::Body, ::Tuple) =
    HeatExchange.compartment_graph((SINGLE_PART_NAME,), ())

function _compartment_graph(body::CompositeBody, couplings::Tuple)
    resolved = _resolve_couplings(body.joins, couplings)
    return HeatExchange.compartment_graph(part_names(body), _shared_core_pairs(body.joins, resolved))
end

# Resolve the stored couplings against a body's joins: an empty tuple (the
# default) means "every join is SharedCore"; otherwise it is parallel to `joins`.
_resolve_couplings(joins::Tuple, ::Tuple{}) = map(_ -> HeatExchange.SharedCore(), joins)
_resolve_couplings(::Tuple, couplings::Tuple) = couplings

# The (parent, child) part-name pairs of the joins whose coupling contracts the
# two parts into one compartment (`SharedCore`). Runs once at graph construction,
# so a small mutable scratch vector is fine.
function _shared_core_pairs(joins::Tuple, couplings::Tuple)
    pairs = Tuple{Symbol,Symbol}[]
    for (join, coupling) in zip(joins, couplings)
        coupling isa HeatExchange.SharedCore && push!(pairs, join_partners(join))
    end
    return pairs
end
