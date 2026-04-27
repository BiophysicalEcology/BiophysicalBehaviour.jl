# =============================================================================
# Shared helpers used across endotherm thermoregulation control strategies.
# Top-level dispatchers, mode dispatch helpers, body-reconstruction helper
# (used by rule-based effectors and IPOPT residuals), and initial-state
# helper for the inner physiological loop.
# =============================================================================

# ----------------------------------------------------------------------------
# Inner-loop output type
# ----------------------------------------------------------------------------

"""
    ThermoregulationOutput(thermoregulation, morphology, energy_flows, mass_flows)

Return type of the inner physiological thermoregulation loop, produced by both
the rule-based (`solve_metabolic_rate`) and IPOPT control paths.

Each field is a NamedTuple:

- `thermoregulation`: temperature states and effector values at the solution.
- `morphology`: body surface areas and geometric measures.
- `energy_flows`: heat fluxes (solar, longwave, generated, evaporation, etc.).
- `mass_flows`: respiration and sweat mass flows.

The behavioural outer loop embeds one of these as the `endotherm_out` field
of its own (wider) NamedTuple result.
"""
struct ThermoregulationOutput{T, M, E, F}
    thermoregulation::T
    morphology::M
    energy_flows::E
    mass_flows::F
end

# ----------------------------------------------------------------------------
# Mode dispatch helpers
# ----------------------------------------------------------------------------

"""
    simultaneous_pant(mode::AbstractThermoregulationMode) -> Bool

Return true if the mode allows panting as an effector running in parallel with
hyperthermia.
"""
simultaneous_pant(::CoreFirst) = false
simultaneous_pant(::CoreAndPantingFirst) = true
simultaneous_pant(::CorePantingSweatingFirst) = true

"""
    simultaneous_sweat(mode::AbstractThermoregulationMode) -> Bool

Return true if the mode allows sweating as an effector running in parallel
with hyperthermia and panting.
"""
simultaneous_sweat(::CoreFirst) = false
simultaneous_sweat(::CoreAndPantingFirst) = false
simultaneous_sweat(::CorePantingSweatingFirst) = true

# ----------------------------------------------------------------------------
# Top-level thermoregulate dispatchers
# ----------------------------------------------------------------------------

"""
    thermoregulate(organism, environment, generated_heat_flow, skin_temperature,
                   insulation_temperature)

Run the thermoregulation loop to find heat balance.

Dispatches on the organism's thermal strategy (`Endotherm`, `Ectotherm`,
`Heterotherm`). For endotherms, further dispatches on the control strategy
(`RuleBasedSequentialControl`, `IPOPTControl`).
"""
function thermoregulate(
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    thermoregulate(
        thermal_strategy(organism),
        organism,
        environment,
        generated_heat_flow,
        skin_temperature,
        insulation_temperature,
    )
end

function thermoregulate(
    ::Endotherm,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    thermoregulate(
        Endotherm(),
        control_strategy(organism),
        organism,
        environment,
        generated_heat_flow,
        skin_temperature,
        insulation_temperature,
    )
end

function thermoregulate(
    ::Heterotherm,
    organism::Organism,
    environment::NamedTuple,
    generated_heat_flow,
    skin_temperature,
    insulation_temperature,
)
    error("Heterotherm thermoregulation not yet implemented")
end

# ----------------------------------------------------------------------------
# Q10 metabolic scaling
# ----------------------------------------------------------------------------

"""
    q10_scale(q10, T, T_reference) -> Real

Q10 metabolic-rate multiplier: `q10 ^ ((T − T_reference) / 10)`.

Inputs must be Unitful temperatures. Scientific code in this package never
operates on stripped temperatures; if a caller has Float64 K it has crossed
a unit boundary it shouldn't have crossed.
"""
q10_scale(q10, T::Unitful.Quantity, T_reference::Unitful.Quantity) =
    q10 ^ (ustrip(u"K", T - T_reference) / 10)

# ----------------------------------------------------------------------------
# Body reconstruction helper
# ----------------------------------------------------------------------------

"""
    rebuild_body(shape, fur::Fur, fat) -> Body
    rebuild_body(shape, fur_depth, fur_diameter, fur_density, fat) -> Body

Construct a `Body` from a shape, a fur insulation layer, and a fat layer.
The 5-argument form constructs the `Fur` from its primitive fields. Used by
rule-based effectors that mutate the organism (`piloerect`, `uncurl`) and by
IPOPT to build trial bodies inside the residual function.
"""
rebuild_body(shape, fur::Fur, fat) = Body(shape, CompositeInsulation(fur, fat))
rebuild_body(shape, fur_depth, fur_diameter, fur_density, fat) =
    rebuild_body(shape, Fur(fur_depth, fur_diameter, fur_density), fat)

# ----------------------------------------------------------------------------
# Initial state for the inner physiological thermoregulate
# ----------------------------------------------------------------------------

"""
    initial_physiological_state(organism, environment_vars) ->
        (; generated_heat_flow, skin_temperature, insulation_temperature)

Standard initial guesses for the inner physiological thermoregulate loop:
no pre-assumed metabolic heat, skin 3 K below the setpoint core temperature,
and insulation surface at air temperature.
"""
function initial_physiological_state(organism::Organism, environment_vars)
    (;
        generated_heat_flow    = zero(thermoregulation(organism).Q_minimum_ref),
        skin_temperature       = HeatExchange.metabolism_pars(organism).core_temperature - 3u"K",
        insulation_temperature = environment_vars.air_temperature,
    )
end
