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
| `absorptivity_min`           | 0.9       | alpha_min                     |
| `absorptivity_max`           | 0.9       | alpha_max                     |
| `absorptivity_step`          | 0.0       | (derived: alpha_max-alpha_min)|
| `pant_max`                   | 1.0       | (panting multiplier max)      |
| `pant_step`                  | 0.1       | (panting step)                |
| `target_temperature`              | 30.0 °C   | T_pref (TPREF, initial preferred temp; rises to active_temperature_max)|
| `target_temperature_step`         | 0.5 K     | TBIG (TPREF increment, 0.5 in NicheMapR THERMOREG)|
| `active_temperature_min`          | 24.0 °C   | T_F_min                       |
| `active_temperature_max`          | 34.0 °C   | T_F_max                       |
| `basking_temperature_min`         | 17.5 °C   | T_B_min                       |
| `escape_temperature_min`        | 6.0 °C    | CT_min                        |
| `escape_temperature_max`          | 40.0 °C   | (escape-behaviour trigger)    |
| `emerge_temperature_min`          | 15.0 °C   | T_RB_min                      |
| `can_retreat_underground`         | `true`    | burrow                        |
| `can_climb`                       | `false`   | climb                         |
| `can_seek_shade`                  | `true`    | shade_seek                    |
| `can_change_absorptivity`         | `false`   | alpha_max ≠ alpha_min         |
| `can_solar_orient`                | `true`    | postur                        |
| `can_press_to_ground`             | `true`    | pct_cond > 0                  |
| `can_pant`                        | `false`   | panting                       |
| `solve_underground`               | `false`   | (NicheMapR BELOWGROUND.f)     |
| `burrow_shade_mode`               | `MaxShadeOnly()` | shdburrow=2 (NicheMapR) |
"""
function example_ectotherm_behavioral_limits(;
    # Shade
    shade_min        = 0.0,
    shade_max        = 0.9,
    shade_step       = 0.03,
    # Solar absorptivity (NicheMapR alpha_min / alpha_max)
    absorptivity_min  = 0.9,
    absorptivity_max  = 0.9,
    absorptivity_step = absorptivity_max - absorptivity_min,   # jumps to max in one step by default
    # Panting
    pant_max         = 1.0,
    pant_step        = 0.1,
    # Depth (soil node indices)
    depth_min_underground = 2,       # mindepth: shallowest accessible underground node
    depth_max             = typemax(Int), # maxdepth: default = use all available depth nodes
    # Height (atmospheric profile node indices)
    height_max            = typemax(Int), # default = use all available height nodes
    # Thermal thresholds
    target_temperature      = u"K"(30.0u"°C"),   # TPREF: starting preferred temp (rises to active_temperature_max)
    target_temperature_step = 0.5u"K",           # TBIG: step size for TPREF increment (0.5 in NicheMapR THERMOREG)
    active_temperature_min  = u"K"(24.0u"°C"),
    active_temperature_max  = u"K"(34.0u"°C"),
    basking_temperature_min    = u"K"(17.5u"°C"),
    escape_temperature_min   = u"K"(6.0u"°C"),
    escape_temperature_max     = u"K"(40.0u"°C"),
    emerge_temperature_min     = u"K"(15.0u"°C"),
    # Capability flags
    can_retreat_underground      = true,
    can_climb                    = false,
    can_seek_shade               = true,
    can_change_absorptivity      = false,
    can_solar_orient             = true,
    can_press_to_ground          = true,
    can_pant                     = false,
    solve_underground                = false,
    burrow_shade_mode                = MaxShadeOnly(),
    emerge_signal                = 0.0u"K/hr",
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
    shade = SteppedParameter(;
        current   = shade_min,
        reference = shade_min,
        max       = shade_max,
        step      = shade_step,
    )
    depth = SteppedParameter(;
        current   = 1,
        reference = 1,
        max       = depth_max,
        step      = 1,
    )
    height = SteppedParameter(;
        current   = 1,
        reference = 1,
        max       = height_max,
        step      = 1,
    )
    absorptivity = SteppedParameter(;
        current   = absorptivity_min,
        reference = absorptivity_min,
        max       = absorptivity_max,
        step      = absorptivity_step,
    )
    pant_rate = SteppedParameter(;
        current   = 1.0,
        reference = 1.0,
        max       = pant_max,
        step      = pant_step,
    )
    target_temperature = SteppedParameter(;
        current   = target_temperature,
        reference = target_temperature,
        max       = active_temperature_max,
        step      = target_temperature_step,
    )
    limits = EctothermBehavioralLimits(;
        control,
        shade,
        depth,
        depth_min_underground,
        height,
        absorptivity,
        pant_rate,
        target_temperature,
        active_temperature_min,
        active_temperature_max,
        basking_temperature_min,
        escape_temperature_min,
        escape_temperature_max,
        emerge_temperature_min,
        can_retreat_underground,
        can_climb,
        can_seek_shade,
        can_change_absorptivity,
        can_solar_orient,
        can_press_to_ground,
        can_pant,
        solve_underground,
        burrow_shade_mode,
        emerge_signal,
    )
    _validate_ectotherm_thresholds(limits)
    return limits
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
