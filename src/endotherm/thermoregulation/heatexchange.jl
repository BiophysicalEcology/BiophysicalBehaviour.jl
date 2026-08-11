# =============================================================================
# TODO: This file is an interim copy of helpers that should live in
# `HeatExchange.jl`. The function signatures here are written exactly as we
# intend to call them once HeatExchange exports them — when that happens,
# delete this file and add the names to a `using HeatExchange:` import in
# `BiophysicalBehaviour.jl`. No call sites should need to change.
#
# These helpers exist because `solve_metabolic_rate` in HeatExchange currently
# inlines all of this setup, and IPOPT-based thermoregulation (`ipopt.jl`)
# needs the same pieces but composed differently (single mean-weighted body
# rather than separate dorsal/ventral solves).
# =============================================================================

"""
    view_factor_decomposition(rad_pars, environment_vars) -> NamedTuple

Decompose the sky view factor into sky and vegetation components based on
overhead shade, then derive the ground component as the residual. Bush comes
straight from `rad_pars`.

Returns a NamedTuple with fields `sky`, `vegetation`, `ground`, `bush`.
"""
function view_factor_decomposition(rad_pars, environment_vars)
    vegetation = rad_pars.sky_view_factor * environment_vars.shade
    sky        = rad_pars.sky_view_factor - vegetation
    ground     = 1 - sky - vegetation
    bush       = rad_pars.bush_view_factor
    return (; sky, vegetation, ground, bush)
end

"""
    solar_heat_flows(body, rad_pars, environment_pars, environment_vars,
                     view_factors, conduction_fraction) -> NamedTuple

Compute solar absorption and split it into dorsal and ventral components.

Wraps `HeatExchange.solar` and applies the per-side adjustments used inside
`solve_metabolic_rate`: the dorsal flow assumes the sun is visible from both
directions (2× direct + 2× sky-scattered), and the ventral flow recovers the
ground-scattered radiation that was attenuated by the conduction-fraction
factor inside `solar`.

`view_factors` is the decomposition returned by `view_factor_decomposition`.

Returns a NamedTuple with fields `total`, `direct`, `sky`, `substrate`,
`dorsal`, `ventral`.
"""
function solar_heat_flows(body, rad_pars, environment_pars, environment_vars,
                          view_factors, conduction_fraction)
    total_area      = BiophysicalGeometry.total_area(body)
    conduction_area = total_area * conduction_fraction
    area_silhouette = silhouette(body, rad_pars.solar_orientation)
    absorptivities  = Absorptivities(rad_pars, environment_pars)
    solar_view_fac  = ViewFactors(view_factors.sky, view_factors.ground, 0.0, 0.0)
    solar_conds     = SolarConditions(environment_vars)
    result          = solar(body, absorptivities, solar_view_fac, solar_conds,
                            area_silhouette, conduction_area)
    if result.solar_flow > 0.0u"W"
        dorsal  = 2.0 * result.direct_flow + result.solar_sky_flow * 2.0
        ventral = (result.solar_substrate_flow / (1.0 - view_factors.sky - view_factors.vegetation)) *
                  (1.0 - 2.0 * conduction_fraction)
    else
        dorsal  = 0.0u"W"
        ventral = 0.0u"W"
    end
    return (; total = result.solar_flow,
              direct = result.direct_flow,
              sky = result.solar_sky_flow,
              substrate = result.solar_substrate_flow,
              dorsal, ventral)
end

"""
    pack_heat_balance_env(environment_vars, environment_pars, view_factors,
                          vegetation_temperature, solar_flow) -> NamedTuple

Pack the `environment_vars` NamedTuple consumed by `HeatExchange.heat_balance`.

`vegetation_temperature` is supplied separately because callers commonly
override `environment_vars.vegetation_temperature` with
`environment_vars.reference_air_temperature` (the standard assumption for
shade-casting vegetation).
"""
function pack_heat_balance_env(environment_vars, environment_pars, view_factors,
                               vegetation_temperature, solar_flow)
    temperature = EnvironmentTemperatures(
        environment_vars.air_temperature,
        environment_vars.sky_temperature,
        environment_vars.ground_temperature,
        vegetation_temperature,
        environment_vars.bush_temperature,
        environment_vars.substrate_temperature,
    )
    atmos = AtmosphericConditions(environment_vars)
    return (;
        temperature,
        view_factors,
        atmos,
        fluid                  = environment_pars.fluid,
        solar_flow,
        gas_fractions          = environment_pars.gas_fractions,
        convection_enhancement = environment_pars.convection_enhancement,
    )
end

"""
    pack_heat_balance_traits(internal_conduction, evaporation_pars, body_emissivity)
        -> NamedTuple

Pack the `traits` NamedTuple consumed by `HeatExchange.heat_balance`.
`body_emissivity` is supplied separately because callers may use a single
mean-weighted emissivity instead of choosing dorsal or ventral.
"""
function pack_heat_balance_traits(internal_conduction, evaporation_pars, body_emissivity)
    return (;
        fat_conductivity   = internal_conduction.fat_conductivity,
        flesh_conductivity = internal_conduction.flesh_conductivity,
        ϵ_body             = body_emissivity,
        skin_wetness       = evaporation_pars.skin_wetness,
        insulation_wetness = evaporation_pars.insulation_wetness,
        bare_skin_fraction = evaporation_pars.bare_skin_fraction,
        eye_fraction       = evaporation_pars.eye_fraction,
    )
end

"""
    longwave_output_flows(body, rad_pars, conduction_fraction, view_factors,
                          dorsal_surface_temperature, ventral_surface_temperature)
        -> NamedTuple

Compute outgoing longwave radiation for the dorsal and ventral surfaces.

Dorsal radiates over the full body surface area; ventral radiates over the
fraction of surface area not in conductive contact with the substrate.
`view_factors` is the decomposition returned by `view_factor_decomposition`.

Callers combine the two components via dorsal/ventral weights appropriate to
their model (e.g. `sky + vegetation` for dorsal, `1 - dorsal_weight` for
ventral).

Returns a NamedTuple with fields `dorsal`, `ventral`.
"""
function longwave_output_flows(body, rad_pars, conduction_fraction, view_factors,
                               dorsal_surface_temperature, ventral_surface_temperature)
    σ            = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    total_area   = BiophysicalGeometry.total_area(body)
    ventral_area = total_area * (1 - conduction_fraction)
    dorsal  = 2 * view_factors.sky * σ * rad_pars.body_emissivity_dorsal *
              total_area * dorsal_surface_temperature^4
    ventral = 2 * view_factors.ground * σ * rad_pars.body_emissivity_ventral *
              ventral_area * ventral_surface_temperature^4
    return (; dorsal, ventral)
end
