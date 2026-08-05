# =============================================================================
# Precomputed per-part geometry cache (plan §3.5, Phase 5).
#
# The multi-part solve recomputes each part's surface area, silhouette, and
# characteristic dimension every outer-loop iteration inside `part_surface_setups`
# — pure functions of shape + insulation (+ the constant solar orientation).
# Across a thermoregulation loop these are invariant unless an effector reshapes
# the body, so we precompute them once and refresh only when geometry changes.
#
# Invalidation is by direct effector dispatch (§3.5): `refresh` is a no-op for
# the physiological effectors (`Vasodilate`, `Hyperthermia`, `Pant`, `Sweat`) —
# the compiler elides it at the call site — and rebuilds for `Piloerect` /
# `Uncurl`, which change fur depth or body axis ratio. The Tier-2 radiation cache
# (view factors, direct/diffuse areas) is a future slot; the current per-part
# view-factor decomposition is cheap and recomputed in `_part_surface_setup`.
# =============================================================================

"""
    ShapeCache{Parts}

Tier-1 geometry cache: a `NamedTuple` keyed by part name whose entries are
`(; total_area, silhouette_area, characteristic_dim)` — the geometry a per-part
surface solve needs, computed once from shape + insulation. Immutable and
stack-allocatable (all-`Unitful` leaves, fixed shape).
"""
struct ShapeCache{Parts<:NamedTuple}
    parts::Parts
end

# Per-part geometry entry. `solar_orientation` is fixed per organism, so the
# silhouette (its only dependency beyond the body) is cacheable.
@inline _shape_cache_entry(body, solar_orientation) = (;
    total_area         = BiophysicalGeometry.total_area(body),
    silhouette_area    = BiophysicalGeometry.silhouette_area(body, solar_orientation),
    characteristic_dim = characteristic_dimension(VolumeCubeRoot(), body),
)

"""
    precompute_shape_cache(organism) -> ShapeCache

Build the Tier-1 geometry cache for every part of `organism`, keyed by part name.
Each part's silhouette uses that part's own `radiation_pars.solar_orientation`.
"""
function precompute_shape_cache(o::Organism)
    body = HeatExchange.body(o)
    part_bodies = _parts(body)
    part_phys = physiology(o)
    entries = map(values(part_bodies), values(part_phys)) do part_body, phys
        _shape_cache_entry(part_body, HeatExchange.radiation_pars(phys).solar_orientation)
    end
    return ShapeCache(NamedTuple{part_names(body)}(entries))
end

"""
    refresh(cache, organism, effector::Effector) -> ShapeCache

Return the geometry cache valid after `effector` has been applied to `organism`.
Physiological effectors leave geometry untouched and return `cache` unchanged (a
no-op the compiler inlines away); `Piloerect` and `Uncurl` reshape the body, so
the cache is rebuilt from the updated organism.
"""
@inline refresh(cache::ShapeCache, ::Organism, ::Effector) = cache
@inline refresh(::ShapeCache, o::Organism, ::Piloerect) = precompute_shape_cache(o)
@inline refresh(::ShapeCache, o::Organism, ::Uncurl)    = precompute_shape_cache(o)

# Per-part cache lookup for `part_surface_setups`. A `nothing` cache yields a
# tuple of `nothing`s (compute-on-the-fly path); a `ShapeCache` yields its
# per-part entries in body order.
@inline _cache_entries(::Nothing, names::Tuple) = ntuple(_ -> nothing, Val(length(names)))
@inline _cache_entries(cache::ShapeCache, ::Tuple) = values(cache.parts)
