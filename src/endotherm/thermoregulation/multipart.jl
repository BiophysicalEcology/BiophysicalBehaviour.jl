# =============================================================================
# BB-side multi-part metabolic-rate solve (Phase 4 / Phase 6-7 integration).
#
# Builds the per-part `part_surface_setups` that `HeatExchange.solve_coupled_
# metabolic_rate` consumes directly from a multi-part `Organism` — its parts
# (a `CompositeBody`, or a single `Body` treated as the one `:body` part) plus
# its per-part physiology NamedTuple — and routes respiration, lung mass, and
# metabolism through the organism's `lung_part` (§3.9).
#
# This is the client-side generalisation of HeatExchange's `_pack_sides`: where
# `_pack_sides` splits one lumped body into a dorsal and a ventral pseudo-side
# and view-factor-weights their means, this walks genuine anatomical parts, each
# seeing the full sky/ground/bush/vegetation hemisphere (un-doubled) through its
# own `radiation_pars`, insulation through its own `insulation_pars`, and ground
# contact through its own `conduction_pars_external`. A single `Body` collapses
# to one part that reproduces the dorsal/ventral result in the symmetric case
# (the Phase 7 equivalence gate, now driven from the organism API).
#
# All parts share one regulated compartment core here (the all-`SharedCore`
# case). Floating compartments (independent leg/head cores via
# `solve_core_temperatures`) layer on top of this and are not wired yet.
# =============================================================================

"""
    solve_multipart_metabolic_rate(organism, environment,
                                   skin_temperature, insulation_temperature;
                                   smoothing = HeatExchange.HardBound())
        -> NamedTuple

Solve the metabolic heat flow of a multi-part `organism` whose parts share one
regulated core at the metabolism setpoint temperature.

Builds one `solve_part_surface` setup per part from that part's own physiology
and geometry (see `part_surface_setups`), then closes the whole-organism
respiration balance through the `lung_part`. `environment` carries
`environment_pars` and `environment_vars` exactly as `solve_metabolic_rate`
expects. Returns the `solve_coupled_metabolic_rate` NamedTuple
(`metabolic_heat_flow`, per-`parts` surface results, `net_metabolic_total`,
mean `skin_temperature` / `insulation_temperature`, `lung_temperature`,
`respiration_out`).
"""
function solve_multipart_metabolic_rate(
    organism::Organism,
    environment,
    skin_temperature,
    insulation_temperature;
    smoothing::HeatExchange.SmoothingStrategy=HeatExchange.HardBound(),
    cache=nothing,
)
    body = HeatExchange.body(organism)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars

    core_temperature = HeatExchange.metabolism_pars(organism).core_temperature

    setups = part_surface_setups(
        organism, environment_pars, environment_vars,
        core_temperature, skin_temperature, insulation_temperature; smoothing, cache,
    )

    # Respiration is a whole-organism balance closed at the lung surface: the
    # forwarding accessors resolve respiration/metabolism/options through the
    # lung part (§3.9), and the lung mass is that part's body mass.
    opts = HeatExchange.options(organism)
    metab_pars = HeatExchange.metabolism_pars(organism)
    resp_pars = HeatExchange.respiration_pars(organism)
    lung_mass = _lung_mass(body, lung_part(organism))

    return solve_coupled_metabolic_rate(;
        part_surface_setups=setups,
        core_temperature,
        skin_temperature,
        insulation_temperature,
        temperature_tolerance=opts.temperature_error_tolerance,
        respire=opts.respire,
        respiration_pars=resp_pars,
        lung_mass,
        air_temperature=environment_vars.air_temperature,
        atmos=AtmosphericConditions(environment_vars),
        gas_fractions=environment_pars.gas_fractions,
        metabolic_heat_flow_setpoint=metab_pars.metabolic_heat_flow,
        resp_tolerance=opts.resp_tolerance,
        smoothing,
    )
end

"""
    part_surface_setups(organism, environment_pars, environment_vars,
                        core_temperature, skin_temperature, insulation_temperature;
                        smoothing) -> Tuple

Build the tuple of per-part `solve_part_surface` setup NamedTuples for
`organism`, one per part in `body` order. Each setup pairs a part's `Body` with
its own physiology (`insulation_pars`, `conduction_pars_external` /
`_internal`, `radiation_pars`, `evaporation_pars`) and the covered join-patch
area subtracted from that part's exposed surface.
"""
function part_surface_setups(
    organism::Organism,
    environment_pars,
    environment_vars,
    core_temperature,
    skin_temperature,
    insulation_temperature;
    smoothing::HeatExchange.SmoothingStrategy=HeatExchange.HardBound(),
    cache=nothing,
)
    body = HeatExchange.body(organism)
    part_bodies = _parts(body)
    part_phys = physiology(organism)
    covered = _part_covered_areas(body)
    cache_entries = _cache_entries(cache, part_names(body))

    # `map` over the same-ordered NamedTuples' values returns a Tuple —
    # type-stable, no runtime Symbol indexing.
    return map(
        values(part_bodies), values(part_phys), values(covered), cache_entries,
    ) do part_body, phys, covered_area, cache_entry
        _part_surface_setup(
            part_body, phys, environment_pars, environment_vars,
            core_temperature, covered_area, skin_temperature, insulation_temperature, cache_entry; smoothing,
        )
    end
end

# Build one part's `solve_part_surface` setup. Mirrors the per-side setup in
# HeatExchange's `_pack_sides`, but for a genuine part: full (un-doubled)
# hemisphere view factors, the part's own insulation/conduction/radiation
# physiology, and the part's actual insulated `Body`. `cache_entry` is a Tier-1
# geometry entry (§3.5) or `nothing` to compute geometry on the fly.
function _part_surface_setup(
    part_body, phys, environment_pars, environment_vars,
    core_temperature, covered_area, skin_temperature, insulation_temperature, cache_entry; smoothing,
)
    ins_pars = HeatExchange.insulation_pars(phys)
    external = HeatExchange.conduction_pars_external(phys)
    internal = HeatExchange.conduction_pars_internal(phys)
    rad_pars = HeatExchange.radiation_pars(phys)
    evap_pars = HeatExchange.evaporation_pars(phys)

    # A bare part radiates/evaporates from skin: force full bare-skin exposure
    # when it carries no effective insulation (matches `_pack_sides`).
    avg_insulation_temperature = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation = insulation_properties(ins_pars, avg_insulation_temperature, rad_pars.ventral_fraction; smoothing)
    if insulation.insulation_test <= 0.0u"m" && evap_pars.bare_skin_fraction < 1.0
        evap_pars = AnimalEvaporationParameters(;
            skin_wetness=evap_pars.skin_wetness,
            insulation_wetness=evap_pars.insulation_wetness,
            eye_fraction=evap_pars.eye_fraction,
            bare_skin_fraction=1.0,
            insulation_fraction=evap_pars.insulation_fraction,
        )
    end

    # Hemisphere decomposition seen by this part (un-doubled — a genuine part
    # sees sky and ground at once, unlike a dorsal/ventral pseudo-side).
    vegetation_factor = rad_pars.sky_view_factor * environment_vars.shade
    sky_factor = rad_pars.sky_view_factor - vegetation_factor
    ground_factor = 1 - sky_factor - vegetation_factor
    bush_factor = rad_pars.bush_view_factor

    # Solar absorbed by this part over its own silhouette. Geometry comes from the
    # Tier-1 cache when supplied, else is computed on the fly (identical values).
    total_area, silhouette, characteristic_dim = _part_geometry(part_body, rad_pars, cache_entry)
    conduction_area = total_area * external.conduction_fraction
    absorptivities = Absorptivities(rad_pars, environment_pars)
    solar_view_factors = ViewFactors(sky_factor, ground_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    solar_flow = solar(
        part_body, absorptivities, solar_view_factors, solar_conditions, silhouette, conduction_area,
    ).solar_flow

    # Ground contact: a part touching the substrate loses heat by conduction;
    # a part clear of the ground has zero conductance coefficient.
    conductance_coefficient = if external.conduction_fraction > 0
        (conduction_area * environment_vars.substrate_conductivity) / environment_pars.conduction_depth
    else
        0.0u"W/K"
    end

    temperature = EnvironmentTemperatures(
        environment_vars.air_temperature,
        environment_vars.sky_temperature,
        environment_vars.ground_temperature,
        environment_vars.reference_air_temperature,   # shade-casting vegetation ≈ reference air
        environment_vars.bush_temperature,
        environment_vars.substrate_temperature,
    )
    view_factors = ViewFactors(sky_factor, ground_factor, bush_factor, vegetation_factor)
    packed_environment = (;
        temperature,
        view_factors,
        atmos=AtmosphericConditions(environment_vars),
        fluid=environment_pars.fluid,
        solar_flow,
        gas_fractions=environment_pars.gas_fractions,
        convection_enhancement=environment_pars.convection_enhancement,
    )
    part_traits = (;
        core_temperature,
        flesh_conductivity=internal.flesh_conductivity,
        fat_conductivity=internal.fat_conductivity,
        ϵ_body=rad_pars.body_emissivity_dorsal,
        skin_wetness=evap_pars.skin_wetness,
        insulation_wetness=evap_pars.insulation_wetness,
        bare_skin_fraction=evap_pars.bare_skin_fraction,
        eye_fraction=evap_pars.eye_fraction,
    )

    return (;
        body=part_body,
        insulation_pars=ins_pars,
        traits=part_traits,
        environment_vars=packed_environment,
        conduction_fraction=external.conduction_fraction,
        conductance_coefficient,
        ventral_fraction=rad_pars.ventral_fraction,
        longwave_depth_fraction=ins_pars.longwave_depth_fraction,
        covered_area,
        characteristic_dim,
    )
end

# Per-part geometry (total surface area, silhouette, characteristic dimension):
# read straight from the Tier-1 cache entry, or compute from the body when no
# cache is threaded through. Both paths return identical values.
@inline _part_geometry(part_body, rad_pars, ::Nothing) = (
    BiophysicalGeometry.total_area(part_body),
    BiophysicalGeometry.silhouette_area(part_body, rad_pars.solar_orientation),
    characteristic_dimension(VolumeCubeRoot(), part_body),
)
@inline _part_geometry(_, _rad_pars, cache_entry::NamedTuple) =
    (cache_entry.total_area, cache_entry.silhouette_area, cache_entry.characteristic_dim)

# Per-part covered join-patch area. A `CompositeBody` folds its joins into a
# per-part area NamedTuple; a plain `Body` has no joins, hence zero covered area.
_part_covered_areas(b::CompositeBody) = BiophysicalGeometry.covered_areas(b.parts, b.joins)
_part_covered_areas(b::Body) =
    NamedTuple{(SINGLE_PART_NAME,)}((zero(BiophysicalGeometry.total_area(b)),))

# The lung-hosting part's body mass, the effective air-exchanging mass for the
# respiration balance. `Val` keeps the NamedTuple lookup type-stable.
_lung_mass(body, name::Symbol) = _lung_mass(_parts(body), Val(name))
_lung_mass(part_bodies::NamedTuple, ::Val{N}) where {N} = getfield(part_bodies, N).shape.mass

# =============================================================================
# Per-part physiology mutation + multi-part rule-based effectors.
#
# Effectors modify the organism's *stored* per-part physiology (the raw
# `heat_exchange` NamedTuple — plain `HeatExchangeTraits`, before the on-demand
# `LungPart` wrapping in `physiology`). Per-part quantities (flesh_conductivity,
# skin_wetness) are applied across the selected parts; whole-organism quantities
# routed to the lung (core temperature, metabolic heat flow, pant rate) are
# applied to the `lung_part` only, via `pant_selector`. Piloerection and uncurl
# rebuild body geometry (fur depth / axis ratio) and so couple to the Tier-1
# geometry cache — deferred to Phase 5; the multi-part ladder below starts at
# vasodilation.
# =============================================================================

"""
    map_part_physiology(f, selector, organism) -> organism

Return `organism` with `f(physiology)` substituted for the stored physiology of
each part picked by `selector`, and every other part left unchanged. Operates on
the raw `heat_exchange` NamedTuple; type-stable via `ntuple(_, Val(N))`.
"""
function map_part_physiology(f, selector::PartSelector, o::Organism)
    sel = select_names(selector, HeatExchange.body(o))
    heat_exchange = HeatExchange.traits(o).heat_exchange
    return @set o.traits.heat_exchange = _map_physiology(f, sel, heat_exchange)
end

@inline function _map_physiology(f, sel::Tuple, heat_exchange::NamedTuple)
    names = keys(heat_exchange)
    vals = values(heat_exchange)
    newvals = ntuple(Val(length(names))) do i
        _member(names[i], sel) ? f(vals[i]) : vals[i]
    end
    return NamedTuple{names}(newvals)
end

# Whole-organism metabolic heat flow is produced at the lung; store it there.
_set_lung_metabolic_heat_flow(o::Organism, flow) =
    map_part_physiology(pant_selector(o), o) do physiology
        @set physiology.metabolism_pars.metabolic_heat_flow = flow
    end

# Vasodilation — raise tissue conductivity across every part (per-part effector).
function _vasodilate_multipart(o::Organism, flesh_conductivity_limits::SteppedParameter)
    flesh_conductivity = min(
        flesh_conductivity_limits.current + flesh_conductivity_limits.step,
        flesh_conductivity_limits.max,
    )
    flesh_conductivity_limits = @set flesh_conductivity_limits.current = flesh_conductivity
    o = map_part_physiology(WholeBody(), o) do physiology
        @set physiology.conduction_pars_internal.flesh_conductivity = flesh_conductivity
    end
    return flesh_conductivity_limits, o
end

# Hyperthermia — let the regulated core rise (whole-organism, routed to lung).
function _hyperthermia_multipart(o::Organism, core_temperature_limits::SteppedParameter, pant_cost)
    minimum_heat_flow = thermoregulation(o).minimum_heat_flow
    core_temperature = min(
        core_temperature_limits.current + core_temperature_limits.step,
        core_temperature_limits.max,
    )
    core_temperature_limits = @set core_temperature_limits.current = core_temperature
    metabolism = HeatExchange.metabolism_pars(o)
    minimum_heat_flow = (minimum_heat_flow + pant_cost) *
                        q10_scale(metabolism.q10, core_temperature, core_temperature_limits.reference)
    o = map_part_physiology(pant_selector(o), o) do physiology
        physiology = @set physiology.metabolism_pars.core_temperature = core_temperature
        @set physiology.metabolism_pars.metabolic_heat_flow = minimum_heat_flow
    end
    return core_temperature_limits, minimum_heat_flow, o
end

# Panting — evaporative cooling at the lung only (§3.9). The pant rate and the
# metabolic cost land on the lung part; no other part's flux is touched.
function _pant_multipart(o::Organism, panting_limits::PantingLimits)
    minimum_heat_flow = thermoregulation(o).minimum_heat_flow
    pant_rate_limits = panting_limits.pant
    pant_rate = min(pant_rate_limits.current + pant_rate_limits.step, pant_rate_limits.max)
    pant_cost = ((pant_rate - 1) / (pant_rate_limits.max + 1e-6 - 1)) *
                (panting_limits.multiplier - 1) * minimum_heat_flow
    panting_limits = @set panting_limits.pant.current = pant_rate
    panting_limits = @set panting_limits.cost = pant_cost
    metabolism = HeatExchange.metabolism_pars(o)
    minimum_heat_flow = (minimum_heat_flow + pant_cost) *
                        q10_scale(metabolism.q10, metabolism.core_temperature, panting_limits.core_temperature_ref)
    o = map_part_physiology(pant_selector(o), o) do physiology
        physiology = @set physiology.respiration_pars.pant = pant_rate
        @set physiology.metabolism_pars.metabolic_heat_flow = minimum_heat_flow
    end
    return panting_limits, minimum_heat_flow, o
end

# Sweating — raise skin wetness across every part (per-part effector).
function _sweat_multipart(o::Organism, skin_wetness_limits::SteppedParameter)
    skin_wetness = min(skin_wetness_limits.current + skin_wetness_limits.step, skin_wetness_limits.max)
    skin_wetness_limits = @set skin_wetness_limits.current = skin_wetness
    o = map_part_physiology(WholeBody(), o) do physiology
        @set physiology.evaporation_pars.skin_wetness = skin_wetness
    end
    return skin_wetness_limits, o
end

"""
    thermoregulate(::Endotherm, ::RuleBasedSequentialControl, organism, environment, init)

Rule-based sequential control for a multi-part (`CompositeBody`) organism.

Drives `solve_multipart_metabolic_rate` and applies the effector ladder —
vasodilation → hyperthermia → panting (lung only) → sweating — until the
metabolic heat flow reaches the minimum, mirroring the single-body loop but with
per-part physiology and lung-routed respiration. Piloerection and uncurl (which
reshape body geometry) are deferred to Phase 5; a multi-part organism whose fur
must flatten before other responses is not yet supported by this ladder.

Returns the final `solve_multipart_metabolic_rate` NamedTuple.
"""
function thermoregulate(
    ::Endotherm,
    control::RuleBasedSequentialControl,
    o::Organism{<:CompositeBody},
    environment::NamedTuple,
    init::NamedTuple,
)
    (; skin_temperature, insulation_temperature) = init
    limits = thermoregulation(o)
    (; mode, tolerance, max_iterations) = control

    flesh_conductivity_limits = limits.flesh_conductivity
    core_temperature_limits   = limits.core_temperature
    panting_limits            = limits.panting
    skin_wetness_limits       = limits.skin_wetness

    # Tier-1 geometry cache (§3.5): computed once, refreshed only when an effector
    # reshapes the body. The ladder below is all physiological (Vasodilate,
    # Hyperthermia, Pant, Sweat), so every `refresh` is a no-op and the cache
    # survives the whole loop.
    cache = precompute_shape_cache(o)

    out = solve_multipart_metabolic_rate(o, environment, skin_temperature, insulation_temperature; cache)
    skin_temperature       = out.skin_temperature
    insulation_temperature = out.insulation_temperature
    metabolic_heat_flow    = out.metabolic_heat_flow
    minimum_heat_flow      = limits.minimum_heat_flow

    iteration = 0
    while metabolic_heat_flow < minimum_heat_flow * (1 - tolerance)
        iteration += 1
        if iteration > max_iterations
            @warn "max_iterations exceeded"
            return out
        end

        if flesh_conductivity_limits.current < flesh_conductivity_limits.max
            flesh_conductivity_limits, o = _vasodilate_multipart(o, flesh_conductivity_limits)
            cache = refresh(cache, o, Vasodilate())

        elseif core_temperature_limits.current < core_temperature_limits.max
            core_temperature_limits, minimum_heat_flow, o =
                _hyperthermia_multipart(o, core_temperature_limits, panting_limits.cost)
            cache = refresh(cache, o, Hyperthermia())
            if simultaneous_pant(mode) && panting_limits.pant.current < panting_limits.pant.max
                panting_limits, minimum_heat_flow, o = _pant_multipart(o, panting_limits)
                cache = refresh(cache, o, Pant())
            end
            if simultaneous_sweat(mode)
                if (skin_wetness_limits.current > skin_wetness_limits.max) || (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return out
                end
                skin_wetness_limits, o = _sweat_multipart(o, skin_wetness_limits)
                cache = refresh(cache, o, Sweat())
            end

        elseif panting_limits.pant.current < panting_limits.pant.max
            panting_limits, minimum_heat_flow, o = _pant_multipart(o, panting_limits)
            cache = refresh(cache, o, Pant())
            if simultaneous_sweat(mode)
                if (skin_wetness_limits.current > skin_wetness_limits.max) || (skin_wetness_limits.step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return out
                end
                skin_wetness_limits, o = _sweat_multipart(o, skin_wetness_limits)
                cache = refresh(cache, o, Sweat())
            end

        else
            if (skin_wetness_limits.current > skin_wetness_limits.max) || (skin_wetness_limits.step <= 0)
                return out
            end
            skin_wetness_limits, o = _sweat_multipart(o, skin_wetness_limits)
            cache = refresh(cache, o, Sweat())
        end

        out = solve_multipart_metabolic_rate(o, environment, skin_temperature, insulation_temperature; cache)
        skin_temperature       = out.skin_temperature
        insulation_temperature = out.insulation_temperature
        metabolic_heat_flow    = out.metabolic_heat_flow
    end

    return out
end
