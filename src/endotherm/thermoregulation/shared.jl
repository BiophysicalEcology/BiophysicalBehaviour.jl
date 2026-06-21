# =============================================================================
# Shared helpers used across endotherm thermoregulation control strategies.
# Top-level dispatchers, mode dispatch helpers, body-reconstruction helper
# (used by rule-based effectors and IPOPT residuals), and initial-state
# helper for the inner physiological loop.
# =============================================================================

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
    thermoregulate(organism, environment, init)

Run the thermoregulation loop to find heat balance.

`init` is a NamedTuple of initial guesses with fields `metabolic_heat_flow`,
`skin_temperature`, and `insulation_temperature` — typically produced by
`initial_physiological_state(organism, environment_vars)`.

Dispatches on the organism's thermal strategy (`Endotherm`, `Ectotherm`,
`Heterotherm`). For endotherms, further dispatches on the control strategy
(`RuleBasedSequentialControl`, `IPOPTControl`).
"""
thermoregulate(organism::Organism, environment::NamedTuple, init::NamedTuple) =
    thermoregulate(thermal_strategy(organism), organism, environment, init)
thermoregulate(strategy::Endotherm, organism::Organism, environment::NamedTuple, init::NamedTuple) =
    thermoregulate(strategy, control_strategy(organism), organism, environment, init)
thermoregulate(::Heterotherm, organism::Organism, environment::NamedTuple, init::NamedTuple) =
    error("Heterotherm thermoregulation not yet implemented")

# ----------------------------------------------------------------------------
# Q10 metabolic scaling
# ----------------------------------------------------------------------------

"""
    q10_scale(q10, temperature, reference_temperature) -> Real

Q10 metabolic-rate multiplier: `q10 ^ ((temperature − reference_temperature) / 10)`.

Inputs must be Unitful temperatures. Scientific code in this package never
operates on stripped temperatures; if a caller has Float64 K it has crossed
a unit boundary it shouldn't have crossed.
"""
q10_scale(q10, temperature::Unitful.Quantity, reference_temperature::Unitful.Quantity) =
    q10 ^ (ustrip(u"K", temperature - reference_temperature) / 10)

# ----------------------------------------------------------------------------
# Body reconstruction helper
# ----------------------------------------------------------------------------

"""
    rebuild_body(shape, fur::FibrousLayer, fat) -> Body
    rebuild_body(shape, fur_depth, fur_diameter, fur_density, fat) -> Body

Construct a `Body` from a shape, a fur insulation layer, and a fat layer.
The 5-argument form constructs the `FibrousLayer` from its primitive fields. Used by
rule-based effectors that mutate the organism (`piloerect`, `uncurl`) and by
IPOPT to build trial bodies inside the residual function.
"""
rebuild_body(shape, fur::FibrousLayer, fat) = Body(shape, CompositeInsulation(fur, fat))
rebuild_body(shape, fur_depth, fur_diameter, fur_density, fat) =
    rebuild_body(shape, FibrousLayer(fur_depth, fur_diameter, fur_density), fat)

# ----------------------------------------------------------------------------
# Initial state for the inner physiological thermoregulate
# ----------------------------------------------------------------------------

"""
    initial_physiological_state(organism, environment_vars) ->
        (; metabolic_heat_flow, skin_temperature, insulation_temperature)

Standard initial guesses for the inner physiological thermoregulate loop:
no pre-assumed metabolic heat, skin 3 K below the setpoint core temperature,
and insulation surface at air temperature.
"""
function initial_physiological_state(organism::Organism, environment_vars)
    (;
        metabolic_heat_flow    = zero(thermoregulation(organism).minimum_heat_flow),
        skin_temperature       = HeatExchange.metabolism_pars(organism).core_temperature - 3u"K",
        insulation_temperature = environment_vars.air_temperature,
    )
end
