# Arrest: Generalized Dormancy/Diapause/Quiescence Controllers

**Source files:**
- `src/arrest/controllers.jl` — `AbstractArrestController` interface, `BelowBound`/`AboveBound`, `AnyDirection`/`RisingDirection`/`FallingDirection`, `signal_value`/`signal_rate`, `NeverController`
- `src/arrest/metrics.jl` — `AbstractMetric`: `RawSignal`, `RawProgress`, `Accumulate`
- `src/arrest/bounds.jl` — `AbstractBound`: `FixedBound`
- `src/arrest/threshold.jl` — `ThresholdController` (metric-vs-bound comparison)
- `src/arrest/function_controller.jl` — `FunctionController` (escape hatch for custom logic)
- `src/arrest/composition.jl` — `AnyController`/`AllController` (fuzzy OR/AND over sub-controllers)
- `src/arrest/arrest_model.jl` — `AbstractArrestModel`: `ComposedArrest`, `AnyArrestModel`, `AllArrestModel`, `step_state`
- `src/arrest/component.jl` — `arrest_component` (`ComponentArray` wrapper for named ODE-state access)
- `src/arrest/structure.jl` — `print_arrest_structure` (composition-tree diagram)

This module expresses dormancy — diapause, quiescence, or any other developmental pause — as a
small set of composable building blocks: a **metric** (what's measured), a **bound** (what it's
compared against), and a **level** (how strongly the resulting response is engaged, `0`–`1`).
Those blocks compose into controllers, and controllers compose into whole arrest models.

The running example throughout this document is *Chortoicetes terminifera* (the Australian
plague locust) egg dormancy, which has two independent pauses:

- **Diapause** — induced only while developmental progress is within a narrow window, and only
  while accumulated chilling and diapause duration both stay under their own limits; either limit
  being exceeded ends diapause.
- **Quiescence** — a desiccation-triggered pause during two separate developmental windows,
  active whenever hydration drops below a tolerance threshold.

An egg is arrested if either pathway is active — diapause and quiescence are alternative routes
to the same outcome, not a single combined rule.

## Control-theory vocabulary, in biological terms

The building blocks borrow control-theory language because that's what they structurally are —
comparators, gates, and integrators wired together. The table below maps each term to its
biological reading, using Chortoicetes throughout.

| Control-theory term | In this module | Chortoicetes example |
|---|---|---|
| Process variable | `Metric` — what a controller reads (`RawSignal`, `RawProgress`, `Accumulate`) | Developmental progress, accumulated chill-hours, hydration |
| Setpoint / reference | `Bound` (`FixedBound`) | A diapause window edge (progress `0.45`), a chill-hour limit (`720.0`), a desiccation tolerance (`0.6`) |
| Error signal | The signed gap between metric and bound | How far past (or short of) a threshold cue the egg currently is |
| Comparator | `BelowBound()`/`AboveBound()` | "Hydration below the desiccation tolerance" (`BelowBound`), "progress above the window's low edge" (`AboveBound`) |
| Relay / bang-bang control | `HardBound()` | A sharp switch: the egg either has, or hasn't, crossed a cue |
| Proportional / smoothed control | `SmoothBound(ε)` | A graded response near a threshold instead of an instant switch — not used by Chortoicetes, whose cues are treated as sharp |
| Derivative (rate) term | `RisingDirection()`/`FallingDirection()` direction gating | A cue keyed to the rate of change of a signal, not just its value — not used by Chortoicetes, but common in real phenology models (e.g. photoperiod-change-rate cues) |
| Controller output | `controller_level ∈ [0,1]` | How strongly a given pathway is currently engaged |
| Integrator | `Accumulate` | The chill-hour counter and the diapause-duration counter — degree-day-style accumulators |
| AND gate (fuzzy min) | `AllController` | Diapause requires window position AND low accumulated chill AND low accumulated duration, all at once |
| OR gate (fuzzy max) | `AnyController` / `AnyArrestModel` | Two separate quiescence windows are alternative triggers; diapause OR quiescence are alternative dormancy pathways |
| Set/reset latch | `ComposedArrest`'s `induction`/`breakage` | An induction cue that engages a pathway vs. a termination cue that releases it |
| Event / zero-crossing detection | `trigger_conditions`, `register_callback`, a host's `ContinuousCallback` | The precise hour a threshold is crossed — timing a fixed time-step loop can only approximate |

## Core idea: level as a continuous multiplier

`arrest_level(model, arrest_state, progress, signals)` returns a `Float64` in `[0, 1]`, meant to
be read directly as a multiplier on a process rate — `development_rate * (1 - arrest_level)` —
rather than passed through a boolean gate. `1` means fully arrested (rate driven to zero), `0`
means unrestricted, and values in between represent partial arrest. There is no boolean
convenience wrapper anywhere in the module: a single fixed cutoff (e.g. `controller_level > 0.5`) can't
represent both a graded response (`SmoothBound`) and an exactly-binary one (`HardBound`) with the
same rule, so a host reads the multiplier directly and derives whatever boolean it needs itself.
For Chortoicetes under `HardBound` throughout, `arrest_level` only ever takes the values `0.0` or
`1.0` — but the interface is the same regardless.

## Three composition layers

1. **Metric vs. bound** → a signed gap (`ThresholdController`'s internal `_gap`) — e.g. accumulated
   chill-hours vs. the `720.0`-hour limit.
2. **Smoothing** turns the gap into a `controller_level` — `HardBound()` (exact step) or `SmoothBound(ε)`
   (differentiable transition).
3. **Controllers compose** via `AllController`/`AnyController` (fuzzy AND/OR); **arrest models
   compose** an induction/breakage pair via `ComposedArrest`, or whole models via
   `AnyArrestModel`/`AllArrestModel`.

## Metrics (`AbstractMetric`)

- `RawSignal(:x)` — reads `signals.x` directly, no state of its own. Chortoicetes uses
  `RawSignal(:hydration)` for the quiescence trigger. Accepts either a bare value or a
  `(; value, rate)` NamedTuple; the latter is required only if a controller uses `RisingDirection`/
  `FallingDirection` direction gating.
- `RawProgress()` — reads `progress` directly (not a signal): both the diapause window
  (`0.45 < progress < 0.50`) and the two quiescence windows (`0.25`–`0.30`, `0.45`–`0.50`) are
  expressed this way. `RawProgress` has no notion of rate (progress is driven externally by the
  host, not tracked as own-state); pairing it with `RisingDirection`/`FallingDirection` raises a clear error.
- `Accumulate(rate, init=0.0)` — integrates a rate functor
  `(progress, signals, model, arrest_state) -> rate` into its own accumulator (see below).
  Chortoicetes uses two: a chill-hour counter and a diapause-duration counter.

Extending: implement `metric_value`/`metric_rate_value` (and `metric_state`/`metric_rate` if the
metric needs its own state) for a new `AbstractMetric` subtype — see "Extending the module"
below for a worked example. Missing methods raise a clear `"no metric_value method for ..."`
error rather than a bare `MethodError`.

### `Accumulate`: a functor-driven integrator

A **functor** is just a callable value — anything that can be invoked with parentheses,
regardless of what kind of object it is. In Julia this covers ordinary functions, closures
(anonymous functions that have captured surrounding variables), and even plain `struct`s made
callable by defining a method on them (`(x::MyType)(args...) = ...`). `Accumulate{R,I}` is
generic over `R` precisely so `rate` can be any of these — the struct only requires that it be
callable with four arguments, never that it be a `Function`. In practice, in this module `rate`
is always a closure.

`Accumulate` (`src/arrest/metrics.jl`) holds two fields: `rate` (the functor) and `init` (the
starting value). Three methods wire it into the controller machinery:

```julia
metric_state(m::Accumulate) = (; accumulator=m.init)
metric_rate(m::Accumulate, own_state, progress, signals, model, arrest_state) =
    (; accumulator=m.rate(progress, signals, model, arrest_state))
metric_value(m::Accumulate, own_state, progress, signals) = own_state.accumulator
```

`metric_state` sets the running total to `init`. `metric_rate` calls the functor every step
with the current `(progress, signals, model, arrest_state)` and reports the result as the
accumulator's instantaneous rate of change. `metric_value` — what actually gets compared against
the `Bound` — just reads the running total. All three are exported, public functions — see
"Extending the module" below.

The integration itself isn't done here at all: it happens in the generic `step_state`, which
performs the Euler step `accumulator_new = accumulator_old + rate * dt`. `Accumulate` only
supplies the derivative each tick (via the functor); `step_state` accumulates it. That's exactly
an integrator block in control-theory terms — an input signal in, a running total out — and it's
why the doc comment calls `rate` a "rate functor": it's not the accumulated value itself, just
the thing that produces its rate of change on demand.

For the two Chortoicetes counters:

- **Chill**: `Accumulate((p, sig, m, as) -> ustrip(u"K", sig.temperature) < 285.65 ? 1.0 : 0.0, 0.0)`.
  The closure ignores everything except `sig` (the signals argument) and returns `1.0` — add one
  hour — whenever the current temperature is below the cold threshold, else `0.0`.
- **Duration**: `Accumulate(duration_rate, 0.0)`, where

  ```julia
  duration_rate(progress, signals, model, arrest_state) =
      (diapause_lo < progress < diapause_hi &&
       arrest_state.induction.chill.metric.accumulator <= cold_hour_threshold &&
       arrest_state.induction.duration.metric.accumulator <= diapause_hour_threshold) ? 1.0 : 0.0
  ```

  This is the interesting case: `arrest_state` here isn't the duration controller's own little
  state slice, it's the *whole model's* state tree, passed down to every rate/level call for
  exactly this reason. `duration_rate` reaches sideways into a sibling controller's accumulator
  (`...induction.chill.metric.accumulator`) to decide whether it should keep ticking — duration
  only accrues while progress is inside the window *and* chill hasn't yet exceeded its own limit.
  This is why the functor signature includes the full `arrest_state` and not just the metric's
  own state: it's the hook that lets one `Accumulate` (or any controller) read another's running
  total elsewhere in the composition tree, which is what makes coupled rules like this one
  expressible without a special case in the framework itself.

## Bounds (`AbstractBound`)

- `FixedBound(value)` — a constant threshold, `Unitful` or plain `Float64`. Every bound in the
  Chortoicetes model is a `FixedBound`: window edges (`0.25`, `0.30`, `0.45`, `0.50`), the chill
  and duration limits (`720.0`, `1240.0`), and the desiccation tolerance (`0.6`).

Extending: implement `bound_value` (and `bound_state` if the bound needs its own state). Missing
methods raise a clear `"no bound_value method for ..."` error — see "Extending the module" below.

## `ThresholdController`

```julia
ThresholdController(; metric, bound, direction=AnyDirection(), comparison=BelowBound(),
                     smoothing=HardBound(), scale=1.0)
```

The chill-accumulation leaf of the diapause controller:

```julia
ThresholdController(
    metric=Accumulate((p, sig, m, as) -> ustrip(u"K", sig.temperature) < 285.65 ? 1.0 : 0.0, 0.0),
    bound=FixedBound(720.0), comparison=BelowBound(),
)
```

Fully parametric `Base.@kwdef` struct — every field, including `scale`, is a type parameter.
`scale` matters for `SmoothBound`: `safe_abs(SmoothBound(ε), gap; scale)` needs `ε * scale` to
carry the same units as `gap`, so a `Unitful` metric/bound combined with `SmoothBound` needs an
explicit `scale` in matching units (e.g. `scale=1.0u"hr"` for an hour-accumulator bound) — the
default `scale=1.0` is dimensionless and only correct for dimensionless metrics.

`register_callback(c)` is `true` only for `HardBound` (the only real discontinuity worth
root-finding); `trigger_conditions(c)` returns the zero-crossing condition (the same signed gap)
for a host's `ContinuousCallback`, or `()` otherwise.

## `FunctionController`

```julia
FunctionController(; condition, level_function=nothing, needs_callback=true)
```

Escape hatch for logic that isn't expressible as metric-vs-bound — genuinely custom or
order-dependent conditions. `condition(own_state, progress, signals, model, arrest_state)`
returns a signed distance (same convention as `ThresholdController`'s gap: `>= 0` means the
condition holds), giving a hard `1.0`/`0.0` level by default; pass `level_function` for a continuous
level instead. Not needed by Chortoicetes — every one of its rules fits the metric-vs-bound shape.

## Composite controllers: `AnyController` / `AllController`

```julia
AnyController(; a=..., b=..., ...)   # OR: max over sub-controller levels
AllController(; a=..., b=..., ...)   # AND: min over sub-controller levels
```

Diapause induction is a four-way `AllController` — window position (two edges), chill, and
duration must *all* hold at once:

```julia
AllController(
    lo=ThresholdController(metric=RawProgress(), bound=FixedBound(0.45), comparison=AboveBound()),
    hi=ThresholdController(metric=RawProgress(), bound=FixedBound(0.50), comparison=BelowBound()),
    chill=ThresholdController(metric=Accumulate(...), bound=FixedBound(720.0), comparison=BelowBound()),
    duration=ThresholdController(metric=Accumulate(...), bound=FixedBound(1240.0), comparison=BelowBound()),
)
```

Quiescence induction is an `AnyController` over two `AllController` windows — either window is
sufficient:

```julia
AnyController(
    w1=window_controller(0.25, 0.30),
    w2=window_controller(0.45, 0.50),
)
```

A named `NamedTuple` of sub-controllers (not a plain `Tuple`) — state threads through
`ComponentArrays` by name (`u.arrest.induction.chill...`). Level combination is Zadeh fuzzy-logic
`max`/`min`, reducing to boolean any/all at the `{0,1}` extremes. This reduction is always
hard-edged — see "Known limitations" below.

## Arrest models (`AbstractArrestModel`)

- `ComposedArrest(; induction, breakage=NeverController())` —
  `arrest_level = induction_level * (1 - breakage_level)`, a fuzzy AND-NOT. Chortoicetes uses the
  default `breakage=NeverController()` for *both* pathways: rather than a dedicated termination
  controller, diapause's own chill/duration limits are folded into the induction `AllController`
  itself — once either accumulator exceeds its limit, that leg of the AND drops to `0` and
  induction (hence `arrest_level`) drops with it. A dedicated `breakage` controller is the
  alternative shape for the same biological idea (a release cue tracked independently of the
  entry cue) — either is valid; Chortoicetes happens to use the first.
- `AnyArrestModel(; a=..., b=...)` / `AllArrestModel(; a=..., b=...)` — whole-model OR/AND
  composition. Chortoicetes is `AnyArrestModel(diapause=diapause_model,
  quiescence=quiescence_model)`: diapause and quiescence are independent pathways, not a single
  induction/breakage pair, so they compose one level up from `ComposedArrest`.

## State and stepping

- `initial_arrest_state(model)` builds a NamedTuple state tree matching the model's composition.
  `arrest_component(model)` wraps it in a `ComponentArray` for named ODE-state access.
- `advance_arrest(model, arrest_state, progress, signals)` returns the matching rate NamedTuple —
  for Chortoicetes, this is where the chill and duration counters' `+1`/hour rates are computed.
- `step_state(state, rate, dt)` recursively advances `state` by `rate * dt`; an empty-rate branch
  (`NamedTuple()`, e.g. a `FixedBound`'s state) is left unchanged rather than erroring.

The whole tree is built with pure functional NamedTuple construction — no in-place mutation
anywhere — so it's trivially Enzyme-differentiable without needing `Duplicated` buffers, unlike
the mutable-cache style used in the IPOPT thermoregulation path (`src/endotherm/thermoregulation/ipopt.jl`).

## ODE integration hooks

Not exported by default (see "Extending the module") — no host wires arrest into a
`ContinuousCallback` yet, so this path is untested; reach it via
`BiophysicalBehaviour.trigger_conditions`/`arrest_conditions` until it is.

`trigger_conditions(model)` / `arrest_conditions(model)` flatten every `HardBound` controller's
zero-crossing condition into a fixed-length `Tuple` (via `_flatten_tuples`, fully unrolled and
concretely typed at compile time), ready to splice into a host's own `ContinuousCallback`
condition tuple for type-stable root-finding — the exact hour diapause's chill accumulator
crosses `720.0`, for instance. `SmoothBound` controllers never register (no discontinuity to
find); `register_callback` reports whether any leaf in a composition needs one at all.

## `print_arrest_structure`

```julia
print_arrest_structure(model, label="model")
```

Not a simulation result — walks the composition tree and prints an indented text diagram (AND/OR
gate at each composite, a human-readable `metric comparison bound` condition at each leaf).
Running it on the Chortoicetes model:

```
└─ chortoicetes (OR)
   ├─ diapause (induction AND NOT breakage)
   │  └─ induction (AND)
   │     ├─ lo: progress > 0.45
   │     ├─ hi: progress < 0.5
   │     ├─ chill: accumulator < 720.0
   │     └─ duration: accumulator < 1240.0
   └─ quiescence (induction AND NOT breakage)
      └─ induction (OR)
         ├─ w1 (AND)
         │  ├─ lo: progress > 0.25
         │  ├─ hi: progress < 0.3
         │  └─ wet: signal :hydration < 0.6
         └─ w2 (AND)
            ├─ lo: progress > 0.45
            ├─ hi: progress < 0.5
            └─ wet: signal :hydration < 0.6
```

Both `diapause` and `quiescence` show `(induction AND NOT breakage)` even though neither declares
an explicit `breakage` controller — that's `ComposedArrest`'s default `NeverController()` printing
literally, since its `controller_level` is always `0`, `1 - breakage` is always `1`, and the AND-NOT
collapses to plain induction. The diagram is otherwise a direct read of the model definition
above: one four-way AND for diapause, an OR of two three-way ANDs for quiescence, and an OR
across both pathways at the top.

## Cross-check against an independent reference

`test/arrest.jl`'s "Chortoicetes equivalence" testset drives the composable model above and a
second, independently written implementation of the same rules (plain functions over
`progress`/`chill`/`duration`/`hydration`, with no controller abstraction at all) through the same
2000-hour synthetic forcing series, asserting `arrest_level(...) > 0.5` matches the independent
implementation's boolean result at every step with zero mismatches.

## Extending the module

Every abstract type in this module (`AbstractMetric`, `AbstractBound`, `AbstractComparison`,
`AbstractDirection`) is extended the same way: implement one or more plain functions dispatched
on the new concrete type, from outside the package. None of these functions are underscore-
prefixed or otherwise internal — but they're deliberately **not exported** by default, since
nothing outside the package extends them yet (adding an export later is non-breaking; removing
one isn't). Reach them either qualified (`BiophysicalBehaviour.metric_value`) or with an explicit
import list: `using BiophysicalBehaviour: AbstractMetric, metric_state, metric_rate, metric_value`.

| Abstract type | Methods to implement | Mandatory? |
|---|---|---|
| `AbstractMetric` | `metric_value` | always |
| | `metric_rate_value` | only if used with `RisingDirection`/`FallingDirection` |
| | `metric_state`, `metric_rate` | only if the metric needs its own state (default: none, no rate) |
| `AbstractBound` | `bound_value` | always |
| | `bound_state` | only if the bound needs its own state (default: none) |
| `AbstractComparison` | `signed_gap(::YourType, value, bound)` | always |
| `AbstractDirection` | `direction_gate(::YourType, rate)` | always |

Each also has a `describe_*` (or `node_label` for a new leaf `AbstractArrestController`) method
used only by `print_arrest_structure`; skipping it means the diagnostic raises a clear
`"no describe_metric method for ..."`-style error for models containing your type, while
simulation itself is unaffected.

As a worked example, here's a metric with its own state — an exponentially-smoothed signal that
averages out noisy hourly readings before thresholding, something none of `RawSignal`/
`RawProgress`/`Accumulate` can express, since it needs to remember a running average between
calls:

```julia
using BiophysicalBehaviour
using BiophysicalBehaviour: AbstractMetric, metric_state, metric_rate, metric_value, signal_value, describe_metric

struct SmoothedSignal{S<:Symbol,T,U} <: AbstractMetric
    signal::S
    halflife::T
    init::U      # seed value, same units as the signal -- mirrors Accumulate's own `init`
end
SmoothedSignal(signal, halflife) = SmoothedSignal(signal, halflife, 0.0)

metric_state(m::SmoothedSignal) = (; average=m.init)

function metric_rate(m::SmoothedSignal, own_state, progress, signals, model, arrest_state)
    x = signal_value(getfield(signals, m.signal))
    return (; average=(x - own_state.average) / m.halflife)
end

metric_value(m::SmoothedSignal, own_state, progress, signals) = own_state.average

# optional: keeps `print_arrest_structure` working for this metric too
describe_metric(m::SmoothedSignal) = "smoothed signal :$(m.signal)"
```

Used exactly like a built-in metric — `init` is seeded in the same units as the signal
(`290.0u"K"`) so that `metric_rate`'s `(x - average) / halflife` and `step_state`'s
`average + rate * dt` both stay dimensionally consistent:

```julia
ThresholdController(metric=SmoothedSignal(:temperature, 24.0u"hr", 290.0u"K"), bound=FixedBound(285.0u"K"), comparison=BelowBound())
```

`AbstractBound`, `AbstractComparison`, and `AbstractDirection` follow the identical shape — e.g. a
custom `AbstractBound` whose threshold is looked up from a species-parameter table at
construction time (via `bound_state`) rather than supplied as a literal `FixedBound(value)`.

## Known limitations / not (yet) built here

- **Composition doesn't propagate smoothing.** `AnyController`/`AllController` and
  `AnyArrestModel`/`AllArrestModel` combine child levels with a hard `max`/`min` — a kink at the
  crossover regardless of whether every child controller uses `SmoothBound`. A fully AD-safe
  composed model would need a `SmoothingStrategy` field on these four types, routed through
  `HeatExchange.safe_max`/`safe_min` instead of `Base.maximum`/`minimum`.
- **Not yet wired into any simulation driver.** None of `src/transient/`, `src/ectotherm/`,
  `src/endotherm/` currently read `arrest_level`.
- **No stochastic pathway.** The package is deterministic end-to-end; there's no RNG-backed bound
  or metric (e.g. population-variation thresholds, hazard-rate triggers). If that's needed later,
  the recommended approach is a Monte Carlo wrapper *outside* the model — draw a threshold once
  per individual, feed it in as a `FixedBound`, run many deterministic sims — rather than
  threading randomness through the controller/bound abstraction itself.
