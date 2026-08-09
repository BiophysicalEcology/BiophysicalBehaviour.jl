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
    # Dispatch on the couplings *type*: the default empty tuple (every organism that
    # doesn't opt in) is all-`SharedCore` — one regulated core — and goes straight to
    # the single-compartment path with its concrete return type, no graph built. Only
    # organisms with explicit couplings pay for the graph and the multi-compartment
    # branch.
    return _solve_multipart(couplings(organism), organism, environment,
        skin_temperature, insulation_temperature; smoothing, cache)
end

# Default (empty couplings) → single regulated compartment shared by all parts.
_solve_multipart(::Tuple{}, organism, environment, skin_temperature, insulation_temperature; smoothing, cache) =
    _solve_single_compartment(organism, environment, skin_temperature, insulation_temperature; smoothing, cache)

# Explicit couplings → build the compartment graph and branch on its compartment
# count. All-`SharedCore` still collapses to the single-compartment path.
function _solve_multipart(::Tuple, organism, environment, skin_temperature, insulation_temperature; smoothing, cache)
    graph = organism_compartment_graph(organism)
    if HeatExchange.num_compartments(graph) == 1
        return _solve_single_compartment(organism, environment, skin_temperature, insulation_temperature; smoothing, cache)
    else
        return _solve_multicompartment(graph, organism, environment,
            skin_temperature, insulation_temperature; smoothing, cache)
    end
end

# One regulated compartment: every part shares the setpoint core (single `Body`, or a
# composite whose joins are all `SharedCore`). The dorsal/ventral-equivalent path.
function _solve_single_compartment(
    organism::Organism, environment, skin_temperature, insulation_temperature; smoothing, cache,
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

# =============================================================================
# Multi-compartment solve (floating ConductiveCoupling compartments).
#
# The regulated compartment (hosting the lung / metabolism setpoint) is pinned at
# the setpoint core; its metabolic heat is the free output closed by respiration.
# Every other compartment "floats": it produces its parts' basal metabolic heat and
# its core settles wherever conduction to its own skins and to its neighbours
# (across `ConductiveCoupling` joins, including the pinned regulated core) balances.
#
# Outer fixed point: solve each part's surface at its compartment's current core →
# aggregate per compartment → `solve_regulated_core_temperatures` for the floating
# cores → repeat until the cores stop moving. Then the regulated compartment's
# metabolic heat is closed by `solve_coupled_metabolic_rate`, with the conductive
# heat it sheds to floating neighbours passed as `extra_net_metabolic`.
# =============================================================================
function _solve_multicompartment(
    graph::HeatExchange.CompartmentGraph{P,N,K}, organism::Organism, environment,
    skin_temperature, insulation_temperature; smoothing, cache,
) where {P,N,K}
    body = HeatExchange.body(organism)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    setpoint = HeatExchange.metabolism_pars(organism).core_temperature
    opts = HeatExchange.options(organism)
    metab_pars = HeatExchange.metabolism_pars(organism)
    resp_pars = HeatExchange.respiration_pars(organism)
    lung_mass = _lung_mass(body, lung_part(organism))
    tol = opts.temperature_error_tolerance

    setups = part_surface_setups(
        organism, environment_pars, environment_vars,
        setpoint, skin_temperature, insulation_temperature; smoothing, cache,
    )

    part_phys = physiology(organism)
    part_compartment = graph.part_compartment                                   # NTuple{N,Int}
    basal = map(phys -> HeatExchange.metabolism_pars(phys).metabolic_heat_flow, values(part_phys))
    k_by_name = map(phys -> HeatExchange.conduction_pars_internal(phys).flesh_conductivity, part_phys)
    regulated = HeatExchange.compartment_of(graph, lung_part(organism))
    entries = _conductance_entries(body, organism, graph, k_by_name)

    # Fixed-point iteration over the floating compartment cores (regulated pinned).
    cores = ntuple(_ -> setpoint, Val(K))
    parts = _solve_parts_at_cores(setups, part_compartment, cores, skin_temperature, insulation_temperature, tol, smoothing)
    for _ in 1:100
        agg = _compartment_aggregates(Val(K), part_compartment, parts, basal)
        new_vec = HeatExchange.solve_regulated_core_temperatures(
            graph, regulated, entries, agg.net_generation,
            agg.flesh_conductance, agg.flesh_weighted_skin, setpoint)
        new_cores = ntuple(i -> new_vec[i], Val(K))
        Δ = maximum(ntuple(i -> abs(new_cores[i] - cores[i]), Val(K)))
        cores = new_cores
        parts = _solve_parts_at_cores(setups, part_compartment, cores, skin_temperature, insulation_temperature, tol, smoothing)
        Δ < tol && break
    end

    # Close the regulated compartment's metabolic heat via respiration, charging it
    # the conductive heat it sheds to the floating compartments.
    coupling_loss = _regulated_coupling_loss(entries, cores, regulated)
    reg_setups = _select_compartment_setups(setups, part_compartment, regulated)
    reg_out = solve_coupled_metabolic_rate(;
        part_surface_setups=reg_setups,
        core_temperature=setpoint,
        skin_temperature,
        insulation_temperature,
        temperature_tolerance=tol,
        respire=opts.respire,
        respiration_pars=resp_pars,
        lung_mass,
        air_temperature=environment_vars.air_temperature,
        atmos=AtmosphericConditions(environment_vars),
        gas_fractions=environment_pars.gas_fractions,
        metabolic_heat_flow_setpoint=metab_pars.metabolic_heat_flow,
        resp_tolerance=opts.resp_tolerance,
        extra_net_metabolic=coupling_loss,
        smoothing,
    )

    number_of_parts = length(parts)
    skin_mean = sum(p -> p.skin_temperature, parts) / number_of_parts
    insulation_mean = sum(p -> p.insulation_temperature, parts) / number_of_parts

    return (;
        metabolic_heat_flow=reg_out.metabolic_heat_flow,
        parts,
        net_metabolic_total=reg_out.net_metabolic_total,
        skin_temperature=skin_mean,
        insulation_temperature=insulation_mean,
        lung_temperature=reg_out.lung_temperature,
        respiration_out=reg_out.respiration_out,
        compartment_core_temperatures=cores,
    )
end

# Solve every part's surface at its compartment's current core temperature.
function _solve_parts_at_cores(setups, part_compartment, cores, skin, insul, tol, smoothing)
    return map(setups, part_compartment) do setup, c
        _solve_part_at_core(setup, cores[c], skin, insul, tol, smoothing)
    end
end

# One part's surface solve at an imposed core (overriding the setpoint baked into
# the setup's traits).
function _solve_part_at_core(setup, core, skin, insul, tol, smoothing)
    traits = merge(setup.traits, (; core_temperature=core))
    return solve_part_surface(;
        merge(setup, (; traits))...,
        skin_temperature=skin,
        insulation_temperature=insul,
        temperature_tolerance=tol,
        smoothing,
    )
end

# Per-compartment aggregates the compartment core solve needs: total core→skin flesh
# conductance, flesh-conductance-weighted skin temperature, and basal generation.
function _compartment_aggregates(::Val{K}, part_compartment, parts, basal) where {K}
    fc0 = zero(first(parts).flesh_conductance)
    fws0 = zero(first(parts).flesh_conductance * first(parts).skin_temperature)
    ng0 = zero(first(basal))
    flesh_conductance = ntuple(Val(K)) do k
        s = fc0
        for i in eachindex(parts)
            part_compartment[i] == k && (s += parts[i].flesh_conductance)
        end
        s
    end
    flesh_weighted_skin = ntuple(Val(K)) do k
        s = fws0
        for i in eachindex(parts)
            part_compartment[i] == k && (s += parts[i].flesh_conductance * parts[i].skin_temperature)
        end
        s
    end
    net_generation = ntuple(Val(K)) do k
        s = ng0
        for i in eachindex(basal)
            part_compartment[i] == k && (s += basal[i])
        end
        s
    end
    return (; flesh_conductance, flesh_weighted_skin, net_generation)
end

# Inter-compartment conductance triples `(compartment_i, compartment_j, G)`, one per
# join whose coupling conducts across a compartment boundary. `SharedCore` joins (both
# parts in one compartment) contribute nothing and are skipped.
function _conductance_entries(body::CompositeBody, organism, graph, k_by_name)
    resolved = _resolve_couplings(body.joins, couplings(organism))
    entries = Tuple{Int,Int,typeof(1.0u"W/K")}[]
    for (join, coupling) in zip(body.joins, resolved)
        parent, child = join_partners(join)
        ci = HeatExchange.compartment_of(graph, parent)
        cj = HeatExchange.compartment_of(graph, child)
        ci == cj && continue
        area = BiophysicalGeometry.join_area(join, body)
        dist = BiophysicalGeometry.internal_distance(join, body)
        G = HeatExchange.contribution_to_conductance(
            coupling, area, dist.parent, dist.child, k_by_name[parent], k_by_name[child])
        G === nothing && continue
        push!(entries, (ci, cj, G))
    end
    return entries
end

# The setups of the parts in `target` compartment (a runtime-length tuple).
_select_compartment_setups(setups, part_compartment, target) =
    Tuple(setups[i] for i in eachindex(setups) if part_compartment[i] == target)

# Conductive heat the regulated compartment sheds to its floating neighbours:
# Σ G·(core_regulated − core_neighbour) over the joins touching `regulated`.
function _regulated_coupling_loss(entries, cores, regulated)
    loss = 0.0u"W"
    for (i, j, G) in entries
        if i == regulated
            loss += G * (cores[i] - cores[j])
        elseif j == regulated
            loss += G * (cores[j] - cores[i])
        end
    end
    return loss
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

# A single `Body` stores one lumped physiology (not a per-part NamedTuple); its sole
# part is `SINGLE_PART_NAME`, so apply `f` when that part is selected. Mirrors
# `broadcast_physiology`'s single-body handling so `map_part_physiology` works on a
# plain `Body` organism, not only a `CompositeBody`.
@inline _map_physiology(f, sel::Tuple, heat_exchange) =
    _member(SINGLE_PART_NAME, sel) ? f(heat_exchange) : heat_exchange

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

# =============================================================================
# Single-part endotherm output assembler (Phase 8).
#
# Produces the full `ThermoregulationOutput` (thermoregulation / morphology /
# energy_flows / mass_flows) that downstream consumers read, from a *single-part*
# multipart solve (a single `Body` → one part = the whole animal). Morphology is a
# pure function of the body, computed with the same expressions the legacy
# dorsal/ventral assembler used, so it matches to numerical tolerance. Energy and
# mass flows come from the one part's `flows` + the whole-organism respiration; the
# genuine per-part flows already integrate the part's actual view factors, so
# longwave in/out read straight off them (longwave_in = Σ incoming radiation,
# longwave_out = longwave_in + net longwave) with no dorsal/ventral reweighting.
# For a single body there are no sides, so the `dorsal`/`ventral` sub-tuples both
# report the one part (kept for output-shape compatibility until callers migrate).
# =============================================================================
function _assemble_endotherm_output(organism::Organism, environment, out;
                                    smoothing::HeatExchange.SmoothingStrategy=HeatExchange.HardBound())
    body = HeatExchange.body(organism)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    ins_pars = HeatExchange.insulation_pars(organism)
    rad_pars = HeatExchange.radiation_pars(organism)
    evap_pars = HeatExchange.evaporation_pars(organism)
    resp_pars = HeatExchange.respiration_pars(organism)
    metab_pars = HeatExchange.metabolism_pars(organism)
    external = HeatExchange.conduction_pars_external(organism)
    internal = HeatExchange.conduction_pars_internal(organism)

    p = out.parts[1]
    resp = out.respiration_out
    core_temperature = metab_pars.core_temperature
    skin_temperature = out.skin_temperature
    insulation_temperature = out.insulation_temperature

    # --- morphology (body-only; same expressions as the legacy assembler) ---
    cond_frac = external.conduction_fraction
    total_area = BiophysicalGeometry.total_area(body)
    area_convection = total_area * (1 - cond_frac)
    vegetation_factor = rad_pars.sky_view_factor * environment_vars.shade
    sky_factor = rad_pars.sky_view_factor - vegetation_factor
    ground_factor = 1 - sky_factor - vegetation_factor
    fat = FatLayer(internal.fat_fraction, internal.fat_density)
    morphology = (;
        total_area,
        area_skin = BiophysicalGeometry.skin_area(body),
        area_evaporation = BiophysicalGeometry.evaporation_area(body),
        area_convection,
        area_conduction = total_area * cond_frac,
        area_silhouette = HeatExchange.silhouette_area(body, rad_pars.solar_orientation),
        sky_view_factor = sky_factor,
        ground_view_factor = ground_factor,
        volume = body.geometry.volume,
        volume_flesh = BiophysicalGeometry.flesh_volume(body),
        characteristic_dimension = HeatExchange.characteristic_dimension(HeatExchange.VolumeCubeRoot(), body),
        fat_mass = body.shape.mass * fat.fraction,
        body.geometry.length...,
    )

    # --- energy flows (one part's flows + respiration) ---
    # Gross outgoing longwave = ε·σ·A·T⁴ over the radiating surface (the part sees the
    # full hemisphere, so its view factors sum to one — no dorsal/ventral doubling).
    # Gross incoming = gross outgoing − net longwave (`flows.longwave` is the net loss).
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    insulation_test = insulation_properties(
        ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction; smoothing
    ).insulation_test
    has_insulation = insulation_test > 0.0u"m"
    radiating_surface_temperature = has_insulation ? p.insulation_temperature : p.skin_temperature
    # Gross longwave radiates from the same surface the internal physics uses:
    # `area_convection = total_area·(1 − conduction_fraction)` (heat_balance.jl). The
    # ground-contact patch's exchange with the substrate (conduction and near-field
    # radiation across the air layer, which can't be separated in practice) is bundled
    # into the conduction term, so it is excluded from the *environmental* longwave here
    # to avoid double-counting — using `total_area` would over-state gross in/out by
    # 1/(1 − conduction_fraction) (the net, longwave_flow_in = out − p.flows.longwave,
    # is unaffected either way).
    radiating_area = area_convection
    longwave_flow_out = rad_pars.body_emissivity_dorsal * σ * radiating_area *
                        radiating_surface_temperature^4
    longwave_flow_in = longwave_flow_out - p.flows.longwave
    respiration_heat_flow = resp === nothing ? zero(out.metabolic_heat_flow) : resp.respiration_heat_flow
    evaporation_heat_flow = p.flows.skin_evaporation + p.flows.insulation_evaporation + respiration_heat_flow
    metabolic_heat_flow = out.metabolic_heat_flow
    heat_balance_val = p.flows.solar + longwave_flow_in + metabolic_heat_flow -
                       longwave_flow_out - p.flows.convection - evaporation_heat_flow - p.flows.conduction
    side_flows = (;
        solar = p.flows.solar, convection = p.flows.convection, conduction = p.flows.conduction,
        net_metabolic = p.flows.net_metabolic,
        skin_evaporation = p.flows.skin_evaporation,
        insulation_evaporation = p.flows.insulation_evaporation,
        longwave = p.flows.longwave, sky_radiation = p.flows.sky_radiation,
        bush_radiation = p.flows.bush_radiation,
        vegetation_radiation = p.flows.vegetation_radiation,
        ground_radiation = p.flows.ground_radiation,
    )
    energy_flows = (;
        solar_flow = p.flows.solar, longwave_flow_in, metabolic_heat_flow,
        evaporation_heat_flow, respiration_heat_flow, longwave_flow_out,
        convection_heat_flow = p.flows.convection, conduction_flow = p.flows.conduction,
        heat_balance = heat_balance_val, balance = heat_balance_val,
        ntry = p.ntry, success = p.success,
        dorsal = side_flows, ventral = side_flows,
    )

    # --- mass flows ---
    latent_heat = HeatExchange.enthalpy_of_vaporisation(environment_vars.air_temperature)
    m_sweat = u"g/hr"(p.flows.skin_evaporation / latent_heat)
    respiration_mass_flow = resp === nothing ? nothing : resp.respiration_mass_flow
    m_evap = respiration_mass_flow !== nothing ? u"g/hr"(respiration_mass_flow + m_sweat) : u"g/hr"(m_sweat)
    mass_flows = (;
        air_flow = resp === nothing ? missing : resp.air_flow,
        oxygen_flow_standard = resp === nothing ? missing : resp.oxygen_flow_standard,
        m_evap, respiration_mass_flow, m_sweat,
        molar_fluxes_in = resp === nothing ? missing : resp.molar_fluxes_in,
        molar_fluxes_out = resp === nothing ? missing : resp.molar_fluxes_out,
    )

    # --- thermoregulation ---
    insulation_final = insulation_properties(
        ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction; smoothing)
    insulation_depth = ins_pars.dorsal.depth * (1 - rad_pars.ventral_fraction) +
                       ins_pars.ventral.depth * rad_pars.ventral_fraction
    axis_ratio_b = body.shape isa Sphere ? 1.0 : body.shape.axis_ratio_b
    side_treg = (;
        skin_temperature = p.skin_temperature,
        insulation_temperature = p.insulation_temperature,
        insulation_conductivity = p.insulation_conductivity,
        insulation_depth,
    )
    thermoregulation = (;
        core_temperature, skin_temperature, insulation_temperature,
        lung_temperature = out.lung_temperature, insulation_depth, axis_ratio_b,
        pant = resp_pars.pant, skin_wetness = evap_pars.skin_wetness,
        flesh_conductivity = internal.flesh_conductivity,
        insulation_conductivity_effective = insulation_final.conductivities.average,
        insulation_conductivity_compressed = insulation_final.conductivity_compressed,
        q10 = metab_pars.q10,
        dorsal = side_treg, ventral = side_treg,
    )

    return ThermoregulationOutput(
        ThermoregulationState(thermoregulation),
        MorphologyState(morphology),
        EnergyFlowState(energy_flows),
        MassFlowState(mass_flows),
    )
end
