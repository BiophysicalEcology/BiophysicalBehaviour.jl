"""
    example_ectotherm_conduction_pars_external(; conduction_fraction=0.1) → ExternalConductionParameters

Create example `ExternalConductionParameters` for an ectotherm with NicheMapR-style defaults.
`conduction_fraction=0.1` corresponds to `pct_cond=10` in NicheMapR's ectotherm model.
"""
function example_ectotherm_conduction_pars_external(;
    conduction_fraction=0.1,
)
    ExternalConductionParameters(; conduction_fraction)
end

"""
    example_ectotherm_metabolism_pars(; kwargs...) → MetabolismParameters

Create example `MetabolismParameters` for an ectotherm.

Defaults to `AndrewsPough2()` (Andrews & Pough 1985 squamate metabolic rate equation).
`Q_metabolism` is not used directly — metabolic heat is computed from the model
given body mass and temperature.
"""
function example_ectotherm_metabolism_pars(;
    Q_metabolism=0.0u"W",
    model=AndrewsPough2(),
)
    MetabolismParameters(; Q_metabolism, model)
end

"""
    example_ectotherm_behavioral_limits(; kwargs...) → EctothermBehavioralLimits

Create example `EctothermBehavioralLimits` with default values matching the
NicheMapR R `ectotherm()` function defaults (ectotherm.R).

# Keyword arguments

| Parameter                    | Default   | NicheMapR equivalent          |
|------------------------------|-----------|-------------------------------|
| `shade_min`                  | 0.0       | minshade / 100                |
| `shade_max`                  | 0.9       | maxshade / 100                |
| `shade_step`                 | 0.03      | delta_shade / 100             |
| `depth_min_underground`      | 2         | mindepth (shallowest underground node)|
| `depth_max`                  | all nodes | maxdepth (default: use all)   |
| `height_max`                 | all nodes | (default: use all heights)    |
| `alpha_min`                  | 0.9       | alpha_min                     |
| `alpha_max`                  | 0.9       | alpha_max                     |
| `alpha_step`                 | 0.0       | (derived: alpha_max-alpha_min)|
| `pant_max`                   | 1.0       | (panting multiplier max)      |
| `pant_step`                  | 0.1       | (panting step)                |
| `T_target`                   | 30.0 °C   | T_pref (TPREF, rises to T_F_max)|
| `T_target_step`              | 1.0 K     | TBIG (TPREF increment, 0.5 in NicheMapR THERMOREG)|
| `T_active_min`               | 24.0 °C   | T_F_min                       |
| `T_active_max`               | 34.0 °C   | T_F_max                       |
| `T_bask`                     | 17.5 °C   | T_B_min                       |
| `T_critical_min`             | 6.0 °C    | CT_min                        |
| `T_critical_max`             | 40.0 °C   | CT_max                        |
| `T_emerge`                   | 17.5 °C   | T_RB_min                      |
| `can_retreat_underground`    | `true`    | burrow                        |
| `can_climb`                  | `false`   | climb                         |
| `can_seek_shade`             | `true`    | shade_seek                    |
| `can_change_absorptivity`    | `false`   | alpha_max ≠ alpha_min         |
| `can_solar_orient`           | `true`    | postur                        |
| `can_press_to_ground`        | `true`    | pct_cond > 0                  |
| `can_pant`                   | `false`   | panting                       |
| `underground_tb_equals_soil` | `true`    | (NicheMapR BELOWGROUND.f)     |
| `burrow_shade_mode`          | `MaxShadeOnly()` | shdburrow=2 (NicheMapR) |
"""
function example_ectotherm_behavioral_limits(;
    # Shade
    shade_min        = 0.0,
    shade_max        = 0.9,
    shade_step       = 0.03,
    # Solar absorptivity (NicheMapR alpha_min / alpha_max)
    alpha_min        = 0.9,
    alpha_max        = 0.9,
    alpha_step       = alpha_max - alpha_min,   # jumps to max in one step by default
    # Panting
    pant_max         = 1.0,
    pant_step        = 0.1,
    # Depth (soil node indices)
    depth_min_underground = 2,      # mindepth: shallowest accessible underground node
    depth_max        = typemax(Int), # maxdepth: default = use all available depth nodes
    # Height (atmospheric profile node indices)
    height_max       = typemax(Int), # default = use all available height nodes
    # Thermal thresholds
    T_target      = u"K"(30.0u"°C"),   # TPREF: starting target temp (rises to T_active_max)
    T_target_step = 0.5u"K",           # TBIG: step size for TPREF increment (0.5 in NicheMapR THERMOREG)
    T_active_min  = u"K"(24.0u"°C"),
    T_active_max  = u"K"(34.0u"°C"),
    T_bask           = u"K"(17.5u"°C"),
    T_critical_min   = u"K"(6.0u"°C"),
    T_critical_max   = u"K"(40.0u"°C"),
    T_emerge         = u"K"(15.0u"°C"),
    # Capability flags
    can_retreat_underground = true,
    can_climb               = false,
    can_seek_shade          = true,
    can_change_absorptivity = false,
    can_solar_orient        = true,
    can_press_to_ground     = true,
    can_pant                = false,
    underground_tb_equals_soil = true,
    burrow_shade_mode          = MaxShadeOnly(),
    Δsoil_signal                = 0.0u"K/hr",
    # Control
    thermoregulation_mode = CoreFirst(),
    tolerance             = 0.005,
    max_iterations        = 1000,
)
    control = RuleBasedSequentialControl(;
        mode=thermoregulation_mode,
        tolerance,
        max_iterations,
    )
    shade_param = SteppedParameter(;
        current   = shade_min,
        reference = shade_min,
        max       = shade_max,
        step      = shade_step,
    )
    depth_param = SteppedParameter(;
        current   = 1,
        reference = 1,
        max       = depth_max,
        step      = 1,
    )
    height_param = SteppedParameter(;
        current   = 1,
        reference = 1,
        max       = height_max,
        step      = 1,
    )
    absorptivity_param = SteppedParameter(;
        current   = alpha_min,
        reference = alpha_min,
        max       = alpha_max,
        step      = alpha_step,
    )
    pant_param = SteppedParameter(;
        current   = 1.0,
        reference = 1.0,
        max       = pant_max,
        step      = pant_step,
    )
    T_target_param = SteppedParameter(;
        current   = T_target,
        reference = T_target,
        max       = T_active_max,
        step      = T_target_step,
    )
    EctothermBehavioralLimits(;
        control,
        shade            = shade_param,
        depth            = depth_param,
        depth_min_underground,
        height           = height_param,
        absorptivity = absorptivity_param,
        pant_rate    = pant_param,
        T_target     = T_target_param,
        T_active_min,
        T_active_max,
        T_bask,
        T_critical_min,
        T_critical_max,
        T_emerge,
        can_retreat_underground,
        can_climb,
        can_seek_shade,
        can_change_absorptivity,
        can_solar_orient,
        can_press_to_ground,
        can_pant,
        underground_tb_equals_soil,
        burrow_shade_mode,
        Δsoil_signal,
    )
end

"""
    example_ectotherm_behavioral_traits(; kwargs...) → BehavioralTraits

Create example `BehavioralTraits` for an ectotherm.

# Keyword arguments
- `activity_period`: Activity period type. Default `Diurnal()` (diurn=1, nocturn=0, crepus=0
  in NicheMapR defaults). Use `CombinedActivity(Diurnal(), Crepuscular())` etc.
  for combined periods.
- All kwargs from `example_ectotherm_behavioral_limits` are forwarded.
"""
function example_ectotherm_behavioral_traits(;
    activity_period = Diurnal(),
    kwargs...,
)
    limits = example_ectotherm_behavioral_limits(; kwargs...)
    BehavioralTraits(; thermoregulation=limits, activity_period)
end

"""
    example_ectotherm_radiation_pars(; kwargs...) → RadiationParameters

Create example `RadiationParameters` for an ectotherm with NicheMapR defaults.

| Parameter         | Default        | NicheMapR equivalent |
|-------------------|----------------|----------------------|
| `α_body_dorsal`   | 0.85           | alpha_max = 0.85     |
| `α_body_ventral`  | 0.85           | alpha_min = 0.85     |
| `ϵ_body_dorsal`   | 0.95           | epsilon = 0.95       |
| `ϵ_body_ventral`  | 0.95           | epsilon = 0.95       |
| `F_sky`           | 0.4            | fatosk = 0.4         |
| `F_ground`        | 0.4            | fatosb = 0.4         |
| `F_bush`          | 0.0            | (none)               |
| `ventral_fraction`| 0.5            | (none)               |
| `solar_orientation`| `Intermediate()` | postur = 0         |
"""
function example_ectotherm_radiation_pars(;
    α_body_dorsal    = 0.85,
    α_body_ventral   = 0.85,
    ϵ_body_dorsal    = 0.95,
    ϵ_body_ventral   = 0.95,
    F_sky            = 0.4,
    F_ground         = 0.4,
    F_bush           = 0.0,
    ventral_fraction = 0.5,
    solar_orientation = Intermediate(),
)
    RadiationParameters(;
        α_body_dorsal,
        α_body_ventral,
        ϵ_body_dorsal,
        ϵ_body_ventral,
        F_sky,
        F_ground,
        F_bush,
        ventral_fraction,
        solar_orientation,
    )
end

"""
    example_ectotherm_evaporation_pars(; kwargs...) → EvaporationParameters

Create example `EvaporationParameters` for an ectotherm with NicheMapR defaults.

| Parameter            | Default  | NicheMapR equivalent         |
|----------------------|----------|------------------------------|
| `skin_wetness`       | 0.001    | pct_wet = 0.1 (%)            |
| `eye_fraction`       | 0.0003   | pct_eyes = 0.03 (%)          |
| `bare_skin_fraction` | 1.0      | naked skin (no insulation)   |
| `insulation_wetness` | 0.0      | (no insulation)              |
| `insulation_fraction`| 0.0      | (no insulation)              |
"""
function example_ectotherm_evaporation_pars(;
    skin_wetness        = 0.001,
    eye_fraction        = 0.0003,
    bare_skin_fraction  = 1.0,
    insulation_wetness  = 0.0,
    insulation_fraction = 0.0,
)
    EvaporationParameters(;
        skin_wetness,
        insulation_wetness,
        eye_fraction,
        bare_skin_fraction,
        insulation_fraction,
    )
end

"""
    example_ectotherm_respiration_pars(; kwargs...) → RespirationParameters

Create example `RespirationParameters` for an ectotherm with NicheMapR defaults.

| Parameter        | Default  | NicheMapR equivalent                             |
|------------------|----------|--------------------------------------------------|
| `fO2_extract`    | 0.2      | F_O2 = 20 (%)                                    |
| `rq`             | 0.8      | RQ = 0.8                                         |
| `pant`           | 1.0      | pantmax = 1 (no panting)                         |
| `Δ_breath`       | 0.1 K    | delta_air = 0.1 (°C above Tair)                  |
| `rh_exit`        | 1.0      | (saturated exhaled air)                          |
| `mouth_fraction` | 0.05     | pct_mouth/100: added to skin_wetness when panting|
"""
function example_ectotherm_respiration_pars(;
    fO2_extract    = 0.2,
    rq             = 0.8,
    pant           = 1.0,
    Δ_breath       = 0.1u"K",
    rh_exit        = 1.0,
    mouth_fraction = 0.05,
)
    RespirationParameters(; fO2_extract, rq, pant, Δ_breath, rh_exit, mouth_fraction)
end

"""
    example_ectotherm_hydraulic_pars(; kwargs...) → HydraulicParameters

Create example `HydraulicParameters` for an ectotherm with NicheMapR defaults.

| Parameter              | Default        | NicheMapR equivalent          |
|------------------------|----------------|-------------------------------|
| `water_potential`      | -707.0 J/kg    | psi_body = -707               |
| `hydraulic_conductance`| 0.0 …          | (no liquid exchange by default)|
| `specific_hydration`   | 0.000304 …     | (standard value)              |
"""
function example_ectotherm_hydraulic_pars(;
    water_potential      = -707.0u"J/kg",
    hydraulic_conductance = 0.0u"kg / (m^2 * s * (J/kg))",
    specific_hydration   = 0.000304u"m^3 / (m^3 * (J/kg))",
)
    HydraulicParameters(; water_potential, hydraulic_conductance, specific_hydration)
end

"""
    example_ectotherm_conduction_pars_internal(; kwargs...) → InternalConductionParameters

Create example `InternalConductionParameters` for an ectotherm with NicheMapR defaults.

| Parameter      | Default       | NicheMapR equivalent |
|----------------|---------------|----------------------|
| `k_flesh`      | 0.5 W/m/K     | k_flesh = 0.5        |
| `fat_fraction` | 0.0           | rinsul = 0           |
"""
function example_ectotherm_conduction_pars_internal(;
    k_flesh      = 0.5u"W/m/K",
    fat_fraction = 0.0,
    k_fat        = 0.23u"W/m/K",
    ρ_fat        = 901.0u"kg/m^3",
)
    InternalConductionParameters(; fat_fraction, k_flesh, k_fat, ρ_fat)
end

"""
    example_ectotherm_heat_exchange_traits(; kwargs...) → HeatExchangeTraits

Create example `HeatExchangeTraits` for an ectotherm using NicheMapR defaults throughout.
All component parameter functions can be overridden individually.
"""
function example_ectotherm_heat_exchange_traits(;
    shape_pars               = DesertIguana(40.0u"g", 1000.0u"kg/m^3"),
    conduction_pars_external = example_ectotherm_conduction_pars_external(),
    conduction_pars_internal = example_ectotherm_conduction_pars_internal(),
    convection_pars          = ConvectionParameters(),
    radiation_pars           = example_ectotherm_radiation_pars(),
    evaporation_pars         = example_ectotherm_evaporation_pars(),
    hydraulic_pars           = example_ectotherm_hydraulic_pars(),
    respiration_pars         = example_ectotherm_respiration_pars(),
    metabolism_pars          = example_ectotherm_metabolism_pars(),
    options                  = example_metabolic_rate_options(),
    insulation_pars          = example_insulation_pars(),
)
    HeatExchange.HeatExchangeTraits(
        shape_pars,
        insulation_pars,
        conduction_pars_external,
        conduction_pars_internal,
        radiation_pars,
        convection_pars,
        evaporation_pars,
        hydraulic_pars,
        respiration_pars,
        metabolism_pars,
        options,
    )
end

"""
    example_ectotherm_organism_traits(; kwargs...) → OrganismTraits

Create example `OrganismTraits` for an ectotherm organism.

Uses `Ectotherm()` as the thermal strategy. Physical traits default to a naked
40 g lizard using `DesertIguana` as the shape, matching NicheMapR's `shape=3`
(desert iguana) default. Surface areas are computed from allometric functions
(Porter and Tracy 1984); body dimensions for convection and conduction use the
cylinder approximation from the original NicheMapR Fortran (see `DesertIguana`).

Pass `shape_pars` explicitly to use a different shape (e.g. `Ellipsoid` or
`LeopardFrog`).

Keyword arguments are forwarded to `example_ectotherm_behavioral_traits`
(e.g. `can_climb`, `can_retreat_underground`, thermal thresholds). To customise heat exchange traits,
pass `heat_exchange` or `shape_pars` explicitly.
"""
function example_ectotherm_organism_traits(;
    thermal_strategy = Ectotherm(),
    heat_exchange    = nothing,
    behavior         = nothing,
    shape_pars       = DesertIguana(40.0u"g", 1000.0u"kg/m^3"),
    kwargs...,
)
    he    = isnothing(heat_exchange) ? example_ectotherm_heat_exchange_traits(; shape_pars) : heat_exchange
    behav = isnothing(behavior)      ? example_ectotherm_behavioral_traits(; kwargs...) : behavior
    OrganismTraits(thermal_strategy, he, behav)
end
