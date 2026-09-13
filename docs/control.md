# Control: Metric/Bound/Condition Primitives

**Source files:**
- `src/control/conditions.jl` — `AbstractCondition` interface, `BelowBound`/`AboveBound`, `AnyDirection`/`RisingDirection`/`FallingDirection`, `signal_value`/`signal_rate`, `NeverController`
- `src/control/metrics.jl` — `AbstractMetric`: `RawSignal`, `RawProgress`, `Accumulate`
- `src/control/bounds.jl` — `AbstractBound`: `FixedBound`
- `src/control/threshold.jl` — `ThresholdController` (metric-vs-bound comparison)
- `src/control/function_controller.jl` — `FunctionController` (escape hatch for custom logic)
- `src/control/proportional.jl` — `ProportionalController` (raw fraction toward a target)
- `src/control/composition.jl` — `AnyController`/`AllController` (fuzzy OR/AND over sub-conditions)
- `src/control/step_state.jl` — `step_state` (generic NamedTuple-rate integration)
- `src/control/structure.jl` — leaf-level `describe_*`/`node_label` (shared by both structure printers below)

This is the package's central abstraction: a dynamical system whose regime changes at a threshold
on some measured or accumulated quantity. Two domains are built on it — `src/arrest/`
(dormancy/diapause/quiescence, see `docs/arrest.md`) and `src/stages/` (ordered life-cycle staging,
see `docs/stages.md`) — as siblings that both depend on this layer and not on each other.

The running example throughout this document is *Chortoicetes terminifera* (the Australian
plague locust) egg dormancy, which has two independent pauses:

- **Diapause** — induced only while developmental progress is within a narrow window, and only
  while accumulated chilling and diapause duration both stay under their own limits.
- **Quiescence** — a desiccation-triggered pause during two separate developmental windows,
  active whenever hydration drops below a tolerance threshold.

## Control-theory vocabulary, in biological terms

| Control-theory term | In this module | Chortoicetes example |
|---|---|---|
| Process variable | `Metric` — what a condition reads (`RawSignal`, `RawProgress`, `Accumulate`) | Developmental progress, accumulated chill-hours, hydration |
| Setpoint / reference | `Bound` (`FixedBound`) | A diapause window edge (progress `0.45`), a chill-hour limit (`720.0`) |
| Error signal | The signed gap between metric and bound | How far past (or short of) a threshold cue the current state is |
| Comparator | `BelowBound()`/`AboveBound()` | "Hydration below the desiccation tolerance", "progress above the window's low edge" |
| Relay / bang-bang control | `HardBound()` | A sharp switch: crossed a cue, or hasn't |
| Proportional / smoothed control | `SmoothBound(ε)` | A graded response near a threshold instead of an instant switch |
| Derivative (rate) term | `RisingDirection()`/`FallingDirection()` | A cue keyed to the rate of change of a signal, not just its value |
| Controller output | `controller_level ∈ [0,1]` | How strongly a pathway is currently engaged |
| Integrator | `Accumulate` | The chill-hour counter and the diapause-duration counter |
| AND gate (fuzzy min) | `AllController` | Diapause requires window position AND low chill AND low duration, all at once |
| OR gate (fuzzy max) | `AnyController` | Two separate quiescence windows are alternative triggers |
| Event / zero-crossing detection | `trigger_conditions`, `register_callback`, a host's `ContinuousCallback` | The precise hour a threshold is crossed |

## Core idea: level as a continuous multiplier

`controller_level(condition, own_state, progress, signals, model, container_state)` returns a
value in `[0, 1]`, meant to be read directly as a multiplier on a process rate — not passed
through a boolean gate. `1` means fully engaged, `0` unrestricted, values between are partial
engagement. A host derives whatever boolean it needs from the level itself; there is no built-in
boolean convenience wrapper, since a single fixed cutoff can't represent both a graded response
(`SmoothBound`) and an exactly-binary one (`HardBound`) with the same rule.

## Metrics (`AbstractMetric`)

- `RawSignal(:x)` — reads `signals.x` directly, no state of its own. Accepts either a bare value
  or a `(; value, rate)` NamedTuple; the latter is required only for `RisingDirection`/`FallingDirection`.
- `RawProgress()` — reads `progress` directly (not a signal). Has no notion of rate (progress is
  driven externally by the host, not tracked as own-state); pairing it with
  `RisingDirection`/`FallingDirection` raises a clear error.
- `Accumulate(rate, init=0.0)` — integrates a rate functor
  `(progress, signals, model, container_state) -> rate` into its own accumulator via `step_state`.
  `container_state` is the *whole* enclosing state tree (an `arrest_state` or `stage_state`), not
  just this metric's own slice — the hook that lets one condition's rate read another's running
  total elsewhere in the composition.

Extending: implement `metric_value`/`metric_rate_value` (and `metric_state`/`metric_rate` if the
metric needs its own state). Missing methods raise a clear `"no metric_value method for ..."`
error rather than a bare `MethodError`.

## Bounds (`AbstractBound`)

- `FixedBound(value)` — a constant threshold, `Unitful` or plain `Float64`.

Extending: implement `bound_value` (and `bound_state` if the bound needs its own state).

## `ThresholdController`

```julia
ThresholdController(; metric, bound, direction=AnyDirection(), comparison=BelowBound(),
                     smoothing=HardBound(), scale=1.0)
```

Fully parametric `Base.@kwdef` struct. `scale` matters for `SmoothBound`: a `Unitful` metric/bound
combined with `SmoothBound` needs an explicit `scale` in matching units (e.g. `scale=1.0u"hr"`) —
the default `scale=1.0` is dimensionless. `register_callback(c)` is `true` only for `HardBound`
(the only real discontinuity worth root-finding); `trigger_conditions(c)` returns the zero-crossing
condition (the signed gap) for a host's `ContinuousCallback`, or `()` otherwise.

**`safe_step`'s numeric type is not input-type-preserving.** `HeatExchange.safe_step` (both
`HardBound` and `SmoothBound`) currently always returns `Float64`, regardless of the metric/bound's
own numeric type. Code that must preserve a non-`Float64` controller level end-to-end (e.g.
`Float32`, a dual/autodiff type) should use `FunctionController`'s `level_function` instead, which
returns whatever type its body produces.

## `FunctionController`

```julia
FunctionController(; condition, level_function=nothing, needs_callback=true)
```

Escape hatch for logic that isn't expressible as metric-vs-bound. `condition(own_state, progress,
signals, model, container_state)` returns a signed distance (same convention as
`ThresholdController`'s gap: `>= 0` means the condition holds), giving a hard `1.0`/`0.0` level by
default; pass `level_function` for a continuous level instead, of any numeric type.

## `ProportionalController`

```julia
ProportionalController(; metric, target)
```

Reads a metric's raw fraction toward `target`, clamped to `[0,1]`, read directly as the level —
unlike `ThresholdController` (even under `SmoothBound`, which only smooths a narrow window right
at the crossing), this ramps across the metric's whole range from zero to target. Always
`register_callback == false` (a smooth, monotone ramp has no crossing event to root-find).

## Composite conditions: `AnyController` / `AllController`

```julia
AnyController(; a=..., b=..., ...)   # OR: max over sub-condition levels
AllController(; a=..., b=..., ...)   # AND: min over sub-condition levels
```

A named `NamedTuple` of sub-conditions (not a plain `Tuple`) — state threads through
`ComponentArrays` by name. Level combination is Zadeh fuzzy-logic `max`/`min`, reducing to boolean
any/all at the `{0,1}` extremes — always hard-edged even if every child uses `SmoothBound` (see
Known limitations).

## `step_state`

`step_state(state, rate, dt)` recursively advances a nested `NamedTuple` state by `rate * dt`; an
empty-rate branch (`NamedTuple()`) is left unchanged rather than erroring. Shared by `arrest/`'s
`advance_arrest` and `stages/`'s `advance_stage` — both produce rate NamedTuples with the same
shape as their state, and both apply them the same way.

## Extending the module

Every abstract type here (`AbstractMetric`, `AbstractBound`, `AbstractComparison`,
`AbstractDirection`, `AbstractCondition`) is extended by implementing plain functions dispatched on
the new concrete type, from outside the package. None of these functions are exported by default;
reach them qualified (`BiophysicalBehaviour.metric_value`) or via an explicit import list.

| Abstract type | Methods to implement | Mandatory? |
|---|---|---|
| `AbstractMetric` | `metric_value` | always |
| | `metric_rate_value` | only if used with `RisingDirection`/`FallingDirection` |
| | `metric_state`, `metric_rate` | only if the metric needs its own state |
| `AbstractBound` | `bound_value` | always |
| | `bound_state` | only if the bound needs its own state |
| `AbstractComparison` | `signed_gap(::YourType, value, bound)` | always |
| `AbstractDirection` | `direction_gate(::YourType, rate)` | always |
| `AbstractCondition` | `controller_level` | always |
| | `initial_controller_state`, `controller_rate` | only if the condition needs its own state |
| | `trigger_conditions`, `register_callback` | only for continuous-ODE hosts needing root-finding |

A custom `AbstractCondition` must return `controller_level` in `[0, 1]` — a documented
precondition, not runtime-checked. Each type also has a `describe_*` (or `node_label` for a new
leaf `AbstractCondition`) method used only by the structure printers; skipping it raises a clear
error only for diagrams containing your type, leaving simulation unaffected.

## Known limitations

- **Composition doesn't propagate smoothing.** `AnyController`/`AllController` (and `arrest/`'s
  `AnyArrestModel`/`AllArrestModel`) combine child levels with a hard `max`/`min` — a kink at the
  crossover regardless of whether every child uses `SmoothBound`. A fully AD-safe composed model
  would need a `SmoothingStrategy` field on these types, routed through
  `HeatExchange.safe_max`/`safe_min` instead of `Base.maximum`/`minimum`.
- **No stochastic pathway.** The package is deterministic end-to-end. If randomness is needed, the
  recommended approach is a Monte Carlo wrapper *outside* the model (draw a threshold once per
  individual, feed it in as a `FixedBound`, run many deterministic sims) rather than threading
  randomness through the condition/bound abstraction itself.
