# Arrest: Generalized Dormancy/Diapause/Quiescence Controllers

**Source files:**
- `src/arrest/arrest_model.jl` — `AbstractArrestModel`: `ComposedArrest`, `AnyArrestModel`, `AllArrestModel`
- `src/arrest/component.jl` — `arrest_component` (`ComponentArray` wrapper for named ODE-state access)
- `src/arrest/structure.jl` — `print_arrest_structure` (composition-tree diagram)

The underlying metric/bound/condition primitives (`ThresholdController`, `FunctionController`,
`AnyController`/`AllController`, `step_state`, etc.) live in `src/control/` — see `docs/control.md`
for those, including the control-theory vocabulary table and the `[0,1]` level contract. This
document covers only what's specific to *arrest*: composing conditions into an induction/breakage
pair, and composing whole arrest models.

The running example throughout this document is *Chortoicetes terminifera* (the Australian
plague locust) egg dormancy, which has two independent pauses:

- **Diapause** — induced only while developmental progress is within a narrow window, and only
  while accumulated chilling and diapause duration both stay under their own limits; either limit
  being exceeded ends diapause.
- **Quiescence** — a desiccation-triggered pause during two separate developmental windows,
  active whenever hydration drops below a tolerance threshold.

An egg is arrested if either pathway is active — diapause and quiescence are alternative routes
to the same outcome, not a single combined rule. The full composed model, built entirely from
`docs/control.md`'s primitives:

```julia
diapause_induction = AllController(
    lo=ThresholdController(metric=RawProgress(), bound=FixedBound(0.45), comparison=AboveBound()),
    hi=ThresholdController(metric=RawProgress(), bound=FixedBound(0.50), comparison=BelowBound()),
    chill=ThresholdController(metric=Accumulate(...), bound=FixedBound(720.0), comparison=BelowBound()),
    duration=ThresholdController(metric=Accumulate(...), bound=FixedBound(1240.0), comparison=BelowBound()),
)
diapause_model = ComposedArrest(induction=diapause_induction)

quiescence_model = ComposedArrest(induction=AnyController(
    w1=window_controller(0.25, 0.30),
    w2=window_controller(0.45, 0.50),
))

chortoicetes_model = AnyArrestModel(diapause=diapause_model, quiescence=quiescence_model)
```

## Core idea: level as a continuous multiplier

`arrest_level(model, arrest_state, progress, signals)` returns a value in `[0, 1]`, read directly
as a multiplier on a process rate — `development_rate * (1 - arrest_level)`. `1` means fully
arrested, `0` unrestricted. For Chortoicetes under `HardBound` throughout, `arrest_level` only
ever takes the values `0.0` or `1.0` — but the interface is the same regardless of smoothing.

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
  composition, same Zadeh fuzzy `max`/`min` idiom as `AnyController`/`AllController`. Chortoicetes
  is `AnyArrestModel(diapause=..., quiescence=...)`: diapause and quiescence are independent
  pathways, not a single induction/breakage pair, so they compose one level up from `ComposedArrest`.

## State and stepping

- `initial_arrest_state(model)` builds a NamedTuple state tree matching the model's composition.
  `arrest_component(model)` wraps it in a `ComponentArray` for named ODE-state access.
- `advance_arrest(model, arrest_state, progress, signals)` returns the matching rate NamedTuple —
  for Chortoicetes, this is where the chill and duration counters' `+1`/hour rates are computed.
- `step_state(state, rate, dt)` (from `control/`) recursively advances `state` by `rate * dt`.

The whole tree is built with pure functional NamedTuple construction — no in-place mutation
anywhere — so it's trivially Enzyme-differentiable without needing `Duplicated` buffers, unlike
the mutable-cache style used in the IPOPT thermoregulation path (`src/endotherm/thermoregulation/ipopt.jl`).

## ODE integration hooks

Not exported by default; reach via `BiophysicalBehaviour.trigger_conditions`/`arrest_conditions`.
`arrest_conditions(model)` flattens every `HardBound` leaf condition's zero-crossing condition into
a fixed-length `Tuple` (via `_flatten_tuples`, fully unrolled and concretely typed at compile
time), ready to splice into a host's own `ContinuousCallback` condition tuple for type-stable
root-finding. `SmoothBound` conditions never register (no discontinuity to find);
`register_callback` reports whether any leaf in a composition needs one at all. This package's own
test suite exercises the underlying `trigger_conditions`/`register_callback` protocol against a
real `VectorContinuousCallback` in `test/stages.jl` (`stage_conditions` is built from the same
`control/` primitives as `arrest_conditions`); a downstream package (Chortoicetes) exercises
`arrest_conditions` itself the same way.

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
collapses to plain induction.

## Cross-check against an independent reference

`test/arrest.jl`'s "Chortoicetes equivalence" testset drives the composable model above and a
second, independently written implementation of the same rules (plain functions over
`progress`/`chill`/`duration`/`hydration`, with no controller abstraction at all) through the same
2000-hour synthetic forcing series, asserting `arrest_level(...) > 0.5` matches the independent
implementation's boolean result at every step with zero mismatches.

## Known limitations

See `docs/control.md`'s "Known limitations" for the composition-smoothing and stochasticity
points, which apply here unchanged. Arrest-specific: none of `src/transient/`, `src/ectotherm/`,
`src/endotherm/` currently read `arrest_level` — no host in this package wires arrest into a
simulation driver yet (Chortoicetes, downstream, does).
