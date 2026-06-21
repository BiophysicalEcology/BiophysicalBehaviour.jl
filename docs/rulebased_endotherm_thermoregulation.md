# Rule-Based Endotherm Thermoregulation

**Source files:**
- `src/endotherm/thermoregulation/rulebased.jl` — effectors and the sequential control loop
- `src/endotherm/thermoregulation/shared.jl` — dispatch chain, Q10 scaling, body reconstruction, initial-state helper
- `src/endotherm/endotherm_traits.jl` — `ThermoregulationLimits`, `InsulationLimits`, `PantingLimits`, `SteppedParameter`

---

## Background

When an endotherm is thermally stressed it recruits effectors roughly in order of cost: it
first adjusts its insulation (cheap, reversible), then its posture, then its peripheral blood
flow, and only resorts to evaporative cooling (expensive in water and energy) as a last resort.
This biological priority hierarchy is directly encoded here as a sequential loop: the model
applies one effector step at a time — in the order piloerection → posture → vasodilation →
hyperthermia → panting → sweating — solving the full heat balance after each step and stopping
as soon as the balance is restored.

In reality some effectors are deployed in parallel rather than strictly in series: many species
pant and allow hyperthermia simultaneously, and some also sweat at the same time. The `mode`
parameter controls this: `CoreFirst` applies effectors strictly one at a time, while
`CoreAndPantingFirst` activates panting alongside hyperthermia, and `CorePantingSweatingFirst`
activates all three together. The choice of mode is therefore a biological hypothesis about
co-activation, not just a performance tuning option.

The approach is transparent and deterministic, making it easy to trace exactly which effectors
were recruited and why. There are, however, some practical limitations to be aware of. The
result can be sensitive to the choice of step size for each effector: large steps may overshoot
the balance point and cause the loop to terminate with a coarser effector combination than the
animal would actually use, while very small steps require many more heat-balance solves to reach
the same endpoint. Step sizes should therefore be chosen to reflect the biological resolution at
which the animal can actually adjust each effector. Speed is also asymmetric: in cold conditions the loop does not run at all (the animal simply
uses maximum insulation and the solver computes the required thermogenesis in one call). In hot
conditions the loop sweeps from whatever the starting effector state is, working through
combinations one step at a time toward increasing heat loss. Insulation depth is always reset
to its maximum at the start of each call (the code explicitly sets it to the fully erected
state), so the loop always starts from the most heat-conserving insulation configuration. Other
effectors — posture, vasodilation, panting, sweating — begin from the `current` value stored
in each `SteppedParameter`, which the caller controls. For a resting animal this is typically
the most heat-conserving posture with no vasodilation, giving the longest possible sweep. For
an active animal already in an extended posture with elevated tissue conductivity, the loop
starts partway along the sequence and reaches the balance point in fewer steps. The number of
iterations — and thus computation time — therefore depends on both the severity of the heat
load and how far the starting state already is from the heat-conserving extreme. Very hot
environments with a resting starting state and small step sizes can require hundreds of
heat-balance solves per time step.

The alternative, `IPOPTControl`, replaces the fixed sequence with a mathematical optimiser that
finds the best combination of all effectors simultaneously (see
[ipopt_endotherm_thermoregulation.md](ipopt_endotherm_thermoregulation.md)).

---

## Overview

`RuleBasedSequentialControl` applies physiological thermoregulation effectors in a fixed
priority order until the animal's heat balance closes above its minimum metabolic rate.
This approach directly encodes the biological priority hierarchy (insulation first, then
postural changes, then circulatory adjustments, then evaporative cooling) in a transparent
step-by-step loop.

It is the default control strategy and complements `IPOPTControl`, which finds the globally
optimal effector combination via nonlinear programming. Both use the same biophysical models
from `HeatExchange.jl`.

---

## Entry Points

```julia
# Short form — dispatches on the control strategy stored in ThermoregulationLimits
thermoregulate(organism, environment, init)

# Explicit form
thermoregulate(Endotherm(), RuleBasedSequentialControl(), organism, environment, init)
```

`thermoregulate(organism, ...)` dispatches via:

```julia
thermoregulate(organism, env, init) →
    thermoregulate(thermal_strategy(organism), organism, env, init) →
        thermoregulate(::Endotherm, control_strategy(organism), organism, env, init)
```

The control strategy is read from `organism.traits.behavioral_traits.thermoregulation.control`.

`init` is a NamedTuple with fields `metabolic_heat_flow`, `skin_temperature`, and
`insulation_temperature`. Produce standard defaults with:

```julia
init = initial_physiological_state(organism, environment_vars)
# metabolic_heat_flow = 0 W, skin_temperature = setpoint − 3 K, insulation_temperature = air_temperature
```

---

## Physiological Logic

The rule-based loop is driven by a single criterion:

```julia
while metabolic_heat_flow < minimum_heat_flow * (1 - tolerance)
    # apply next effector ...
end
```

`metabolic_heat_flow` is the value returned by `solve_metabolic_rate` — the heat generation
required to close the energy balance given the current organism state. `minimum_heat_flow` is
the basal metabolic rate (the minimum the organism can generate).

**Cold conditions** (`metabolic_heat_flow ≥ minimum_heat_flow`): the loop does not run. The
organism starts with fully erected insulation (maximum piloerection depth) and `solve_metabolic_rate`
computes whatever elevated thermogenesis is required to close the balance. No step-by-step
effector adjustment is needed — the solver finds the consistent metabolic rate directly.

**Hot conditions** (`metabolic_heat_flow < minimum_heat_flow`): even at minimum metabolic rate,
the organism gains more heat from the environment than it can dissipate. The loop applies cooling
effectors one step at a time until `metabolic_heat_flow` rises to meet `minimum_heat_flow`,
meaning the heat balance can now close at or above basal rate.

---

## Piloerection Starting Condition

Before the first `solve_metabolic_rate` call, the insulation depth is set to its maximum
(fully erected) if the step size is positive:

```julia
insulation_limits.dorsal.current  = insulation_limits.dorsal.max
insulation_limits.ventral.current = insulation_limits.ventral.max
```

This ensures the animal starts in its best-insulated state regardless of the `init` values.
The loop then flattens insulation as needed for the hot case.

---

## Effectors and Application Order

Six effectors are applied in priority order during the hot loop:

| Priority | Effector | Function | Effect |
|---|---|---|---|
| 1 | Piloerection | `piloerect` | Flatten fur/feathers; reduce insulation depth toward reference |
| 2 | Posture | `uncurl` | Increase body aspect ratio; expose more surface area |
| 3 | Vasodilation | `vasodilate` | Increase flesh conductivity; push more heat to skin |
| 4 | Hyperthermia | `hyperthermia` | Allow core temperature to rise; increases thermal gradient |
| 5 | Panting | `pant` | Increase ventilation rate for respiratory evaporative cooling |
| 6 | Sweating | `sweat` | Increase skin wetness for cutaneous evaporative cooling |

Each effector adjusts its parameter by one `step` per loop iteration and returns the updated
`organism`. The loop calls `solve_metabolic_rate` after each step to check whether the balance
is restored.

### Piloerect

```julia
piloerect(organism, insulation_limits) → (insulation_limits, organism)
```

Decrements dorsal and ventral insulation depth by `step × fibre_length` toward their
`reference` (minimum) values. Rebuilds the body geometry (`rebuild_body`) with the new
mean insulation depth. Applied first — reducing insulation is the least costly cooling strategy.

### Uncurl

```julia
uncurl(organism, aspect_ratio_limits) → (aspect_ratio_limits, organism)
```

Increases `axis_ratio_b` by `step` toward `max`, elongating the body to increase surface-area-to-
volume ratio. No-ops for a `Sphere` (setting `current = max` immediately). Rebuilds body geometry.

### Vasodilate

```julia
vasodilate(organism, flesh_conductivity_limits) → (flesh_conductivity_limits, organism)
```

Increases `flesh_conductivity` by `step` toward `max`. Higher conductivity routes more
metabolic heat from core to skin, improving heat loss through convection and radiation.

### Hyperthermia

```julia
hyperthermia(organism, core_temperature_limits, pant_cost) → (core_temperature_limits, minimum_heat_flow, organism)
```

Increments the core temperature setpoint by `step` toward `max`. A higher core temperature
increases the core-to-skin gradient, driving more heat outward. Also Q10-scales
`minimum_heat_flow` to account for the temperature dependence of metabolic reactions:

```julia
minimum_heat_flow = (minimum_heat_flow + pant_cost) × Q10^((core_temp − core_temp_ref) / 10)
```

Depending on `mode`, panting and/or sweating may run simultaneously with hyperthermia
(see Control Modes below).

### Pant

```julia
pant(organism, panting_limits) → (panting_limits, minimum_heat_flow, organism)
```

Increments `pant_rate` by `step` toward `max`. Panting increases ventilation rate for
respiratory evaporative cooling. The panting cost (additional metabolic heat from breathing
muscles) is tracked in `panting_limits.cost` and added to `minimum_heat_flow`. Q10-scales
`minimum_heat_flow` using the current core temperature.

### Sweat

```julia
sweat(organism, skin_wetness_limits) → (skin_wetness_limits, organism)
```

Increments `skin_wetness` by `step` toward `max`. Higher skin wetness increases cutaneous
evaporative heat loss. Applied last — sweating is the most costly strategy in terms of water
loss and is reserved for when other effectors are exhausted.

---

## Control Modes

The `mode` field of `RuleBasedSequentialControl` controls whether panting and sweating activate
simultaneously with hyperthermia or only sequentially:

| Mode | Panting with hyperthermia? | Sweating with hyperthermia? |
|---|---|---|
| `CoreFirst()` | No | No |
| `CoreAndPantingFirst()` | Yes | No |
| `CorePantingSweatingFirst()` | Yes | Yes |

`CorePantingSweatingFirst` is the typical choice for birds and many mammals. With this mode,
as soon as hyperthermia is triggered, panting and sweating also step forward each iteration —
all three activate in parallel.

Set `mode` in `RuleBasedSequentialControl`:

```julia
ThermoregulationLimits(;
    control = RuleBasedSequentialControl(;
        mode           = CorePantingSweatingFirst(),
        tolerance      = 0.005,
        max_iterations = 1000,
    ),
    ...
)
```

---

## Q10 Metabolic Scaling

```julia
q10_scale(q10, temperature, reference_temperature) = q10 ^ ((temperature − reference_temperature) / 10)
```

Q10 scaling is applied in both `hyperthermia` and `pant` to keep `minimum_heat_flow` consistent
with the temperature dependence of biochemical reaction rates. As core temperature rises,
minimum metabolic rate rises proportionally — preventing the solver from finding spurious
solutions where high temperature alone closes the heat balance at artificially low metabolic rates.

---

## Key Structs

### `SteppedParameter`

```julia
SteppedParameter(; current, reference=current, max, step)
```

Tracks a parameter that adjusts incrementally. `current` is the active value; `reference` is
the baseline (minimum for the hot case); `max` is the ceiling; `step` is the per-iteration
adjustment size.

### `InsulationLimits`

```julia
InsulationLimits(; dorsal::SteppedParameter, ventral::SteppedParameter)
```

Holds separate `SteppedParameter` instances for dorsal and ventral insulation depth. The
`step` field is multiplied by fibre length inside `piloerect` so it represents a fractional
change in fibre length per iteration.

### `PantingLimits`

```julia
PantingLimits(;
    pant::SteppedParameter,    # panting rate (dimensionless; 1 = resting)
    cost = 0.0u"W",            # current panting metabolic cost
    multiplier = 1.05,         # cost at max panting as fraction of minimum_heat_flow
    core_temperature_ref,      # reference core temperature for Q10
)
```

### `ThermoregulationLimits`

The top-level struct aggregating all effector limits plus the control strategy and IPOPT
penalty weights. Stored at `organism.traits.behavioral_traits.thermoregulation`.

```julia
ThermoregulationLimits(;
    control            = RuleBasedSequentialControl(),
    minimum_heat_flow,          # basal metabolic rate (W)
    insulation,                 # InsulationLimits
    aspect_ratio_factor,        # SteppedParameter (body shape)
    flesh_conductivity,         # SteppedParameter (vasodilation)
    core_temperature,           # SteppedParameter (hyperthermia)
    panting,                    # PantingLimits
    skin_wetness,               # SteppedParameter (sweating)
    # IPOPT penalty weights (ignored by RuleBasedSequentialControl):
    core_temperature_penalty  = 1.0,
    metabolic_heat_penalty    = 0.1,
    panting_penalty           = 1.0,
    skin_wetness_penalty      = 1.0,
    gradient_penalty          = 0.0,
    target_core_skin_gradient = 2.0,
)
```

---

## Output

`thermoregulate` returns the `endotherm_out` NamedTuple from the final `solve_metabolic_rate`
call. Key sub-fields:

```julia
endotherm_out.thermoregulation   # temperatures, conductivities, insulation depths, pant, skin_wetness
endotherm_out.energy_flows       # metabolic_heat_flow, solar_flow, convection_heat_flow, balance, ...
endotherm_out.mass_flows         # air_flow, m_evap, m_sweat, respiration_mass_flow, ...
endotherm_out.morphology         # total_area, sky_view_factor, ground_view_factor, ...
```

The format is identical to the IPOPT solver output, so both solvers can be swapped without
changing downstream code.

---

## Example Usage

```julia
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful

# Build organism with thermoregulation limits
thermoregulation_limits = ThermoregulationLimits(;
    control = RuleBasedSequentialControl(;
        mode           = CorePantingSweatingFirst(),
        tolerance      = 0.005,
        max_iterations = 1000,
    ),
    minimum_heat_flow = metabolic_rate(McKechnieWolf(), shape_pars.mass),
    insulation = InsulationLimits(;
        dorsal  = SteppedParameter(; current = 16.0u"mm", reference = 5.9u"mm", max = 16.0u"mm", step = 0.1),
        ventral = SteppedParameter(; current = 15.9u"mm", reference = 5.7u"mm", max = 15.9u"mm", step = 0.1),
    ),
    aspect_ratio_factor = SteppedParameter(; current = 1.1, max = 5.0, step = 0.1),
    flesh_conductivity  = SteppedParameter(; current = 0.9u"W/m/K", max = 2.8u"W/m/K", step = 0.1u"W/m/K"),
    core_temperature    = SteppedParameter(; current = 311.15u"K", reference = 311.15u"K", max = 316.15u"K", step = 0.1u"K"),
    panting = PantingLimits(;
        pant                = SteppedParameter(; current = 1.0, max = 15.0, step = 0.01),
        multiplier          = 1.05,
        core_temperature_ref = 311.15u"K",
    ),
    skin_wetness = SteppedParameter(; current = 0.005, max = 0.05, step = 0.0025),
)

# Run a sweep across air temperatures
results = NamedTuple[]

for air_temperature in air_temperatures
    environment_vars = example_environment_vars(; air_temperature, ...)
    environment      = (; environment_pars, environment_vars)
    init             = initial_physiological_state(organism, environment_vars)

    # Reset init every call — no carry-forward for rule-based (effectors reset via organism rebuild)
    out = thermoregulate(organism, environment, init)

    push!(results, (
        core_temperature    = out.thermoregulation.core_temperature,
        metabolic_heat_flow = out.energy_flows.metabolic_heat_flow,
        pant                = out.thermoregulation.pant,
        skin_wetness        = out.thermoregulation.skin_wetness,
    ))
end
```

See `examples/budgerigar.jl` for a complete worked example including observed data comparison.

---

## Comparison with IPOPT Control

| Aspect | `RuleBasedSequentialControl` | `IPOPTControl` |
|---|---|---|
| **Philosophy** | Fixed priority order; one effector at a time | Simultaneous optimisation across all effectors |
| **Cold response** | Implicit: start at max insulation; solver computes metabolic rate | Explicit: all effectors jointly adjusted to minimise objective |
| **Hot response** | Sequential: piloerect → uncurl → vasodilate → hyperthermia → pant → sweat | Order emerges from penalty weights |
| **Mode control** | `CoreFirst`, `CoreAndPantingFirst`, `CorePantingSweatingFirst` | Relative penalty weights (`panting_penalty` vs `skin_wetness_penalty`) |
| **Tuning** | Step sizes and iteration order | Five scalar penalty weights |
| **Speed** | Fast (deterministic iterations) | Comparable; exact Hessian via Enzyme reduces IPOPT iterations |
| **Guarantee** | Always converges; may not be globally optimal | NLP optimum within bounds; may not converge in extreme cases |
| **Warm-starting** | Partial: caller sets `current` on each `SteppedParameter` and passes previous `skin_temperature` / `insulation_temperature` via `init` | `IPOPTSolverCache` enables automatic primal+dual warm-start across a sweep |
