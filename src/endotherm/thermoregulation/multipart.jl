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
)
    body = HeatExchange.body(organism)
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars

    core_temperature = HeatExchange.metabolism_pars(organism).core_temperature

    setups = part_surface_setups(
        organism, environment_pars, environment_vars,
        core_temperature, skin_temperature, insulation_temperature; smoothing,
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
)
    body = HeatExchange.body(organism)
    part_bodies = _parts(body)
    part_phys = physiology(organism)
    covered = _part_covered_areas(body)

    # `map` over the three same-ordered NamedTuples' values returns a Tuple —
    # type-stable, no runtime Symbol indexing.
    return map(
        values(part_bodies), values(part_phys), values(covered),
    ) do part_body, phys, covered_area
        _part_surface_setup(
            part_body, phys, environment_pars, environment_vars,
            core_temperature, covered_area, skin_temperature, insulation_temperature; smoothing,
        )
    end
end

# Build one part's `solve_part_surface` setup. Mirrors the per-side setup in
# HeatExchange's `_pack_sides`, but for a genuine part: full (un-doubled)
# hemisphere view factors, the part's own insulation/conduction/radiation
# physiology, and the part's actual insulated `Body`.
function _part_surface_setup(
    part_body, phys, environment_pars, environment_vars,
    core_temperature, covered_area, skin_temperature, insulation_temperature; smoothing,
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

    # Solar absorbed by this part over its own silhouette.
    total_area = BiophysicalGeometry.total_area(part_body)
    conduction_area = total_area * external.conduction_fraction
    silhouette = BiophysicalGeometry.silhouette_area(part_body, rad_pars.solar_orientation)
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
    )
end

# Per-part covered join-patch area. A `CompositeBody` folds its joins into a
# per-part area NamedTuple; a plain `Body` has no joins, hence zero covered area.
_part_covered_areas(b::CompositeBody) = BiophysicalGeometry.covered_areas(b.parts, b.joins)
_part_covered_areas(b::Body) =
    NamedTuple{(SINGLE_PART_NAME,)}((zero(BiophysicalGeometry.total_area(b)),))

# The lung-hosting part's body mass, the effective air-exchanging mass for the
# respiration balance. `Val` keeps the NamedTuple lookup type-stable.
_lung_mass(body, name::Symbol) = _lung_mass(_parts(body), Val(name))
_lung_mass(part_bodies::NamedTuple, ::Val{N}) where {N} = getfield(part_bodies, N).shape.mass
