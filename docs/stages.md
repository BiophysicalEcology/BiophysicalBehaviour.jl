# Stages: Ordered Life-Cycle Topology

**Source files:**
- `src/stages/stage_sequence.jl` — `AbstractStage`, `AbstractTransition`, `StageSequence`, the
  `stage_traits`/`transition_controller`/`stage_key`/`transition_key` accessor protocol
- `src/stages/lookup.jl` — type-keyed `getindex`/`haskey`/`get`
- `src/stages/weights.jl` — `stage_weights`, `current_stage`, `stage_value`
- `src/stages/state.jl` — `initial_stage_state`, `advance_stage`
- `src/stages/conditions.jl` — `stage_conditions`
- `src/stages/component.jl` — `stage_component`
- `src/stages/structure.jl` — `print_stage_structure`

This module expresses an ordered life cycle (egg → instar → instar → ..., or any other
domain's staged progression) as a sequence of stages separated by transitions, where each
transition is a `control/` condition (see `docs/control.md`) — the same threshold/level machinery
`arrest/` uses for dormancy. `stages/` ships no domain vocabulary: `AbstractStage`/
`AbstractTransition` are abstract only, and a downstream package (or, for this package's own
tests, `test/stages.jl`) supplies concrete subtypes such as `Egg`/`Instar{N}`/`Hatch`/`Molt{N}`.

## Topology

A `StageSequence` holds exactly one more stage than transitions — `stage_1 --[transition_1]-->
stage_2 --[transition_2]--> ... --> stage_n`. Both construction paths go through the same
validated inner constructor:

```julia
StageSequence(initial_stage, transition_1 => stage_2, transition_2 => stage_3, ...)  # chain form
StageSequence((stage_1, stage_2, ..., stage_n), (transition_1, transition_2, ..., transition_{n-1}))  # raw tuple-pair form
```

Construction throws `ArgumentError` if the cardinality is wrong, or if any two stages (or any two
transitions) share a `stage_key`/`transition_key`. A singleton sequence (one stage, zero
transitions) is valid.

## The accessor protocol

```julia
abstract type AbstractStage{T} end
abstract type AbstractTransition{T<:AbstractCondition} end

stage_traits(stage::AbstractStage) = stage.traits
transition_controller(transition::AbstractTransition) = transition.controller
```

The accessor functions, not the field names, are the extensibility contract. A domain type is not
required to be literally named `traits`/`controller` internally — override `stage_traits`/
`transition_controller` instead (`test/stages.jl`'s `ExternalStage`/`ExternalTransition` fixtures,
with fields `data`/`gate`, exist specifically to prove this). Every generic function in `stages/`
calls the accessor, never `.traits`/`.controller` directly.

`stage_key`/`transition_key` follow the same pattern, with one addition: the primary method
dispatches on `Type{S}`, not the instance —

```julia
stage_key(::Type{S}) where {S<:AbstractStage} = nameof(S)
stage_key(stage::AbstractStage) = stage_key(typeof(stage))
```

— because `StageSequence` resolves and caches these keys once, at construction, as type
parameters (`StageSequence{StageKeys,TransitionKeys,...}`), so every hot-path `NamedTuple{...}`
construction downstream (`stage_weights`, `initial_stage_state`, `advance_stage`) reads a literal
key tuple directly off the type, with no per-call key re-derivation. Overriding only the instance
method would silently defeat this. The default `nameof(S)` collides across different parametric
instantiations of one type (`Instar{1}` and `Instar{2}` are both `:Instar`) — override the
`Type{S}` method to fold the parameter in, e.g. `stage_key(::Type{<:Instar{N}}) where N =
Symbol(:Instar, N)`, whenever more than one instantiation of the same parametric stage type
appears in one sequence. (`Symbol(a, b)`-style construction itself allocates, but this only runs
once per sequence, at construction — never in a hot path.)

## Type-keyed lookup

```julia
sequence[SomeStageType]          # unwraps via stage_traits; ArgumentError if absent
sequence[SomeTransitionType]     # unwraps via transition_controller
get(sequence, SomeType, default)
haskey(sequence, SomeType)
```

Matched via plain `isa`, which already does the right thing for partial parametric types:
`Instar{3}` is `Instar{3,T} where T`, so `instance isa Instar{3}` matches any concrete
`Instar{3,...}` and correctly rejects `Instar{5,...}`. `haskey` tests presence of the matched
stage/transition object itself, never its unwrapped (possibly-`nothing`) data, so a
present-but-`nothing` value is never confused with absence.

## Weighting: sequential gating

For transition levels `l_1, ..., l_{n-1}`:

```
w_1 = 1 - l_1
w_i = l_1 * l_2 * ... * l_{i-1} * (1 - l_i)   for 1 < i < n
w_n = l_1 * l_2 * ... * l_{n-1}
```

Nonnegative for levels in `[0,1]`, sums to 1 algebraically. Numeric-type-preserving: no literal
`Float64` anywhere in `stage_weights`, so `Float32`/`Unitful`/dual-number controller levels pass
through unchanged (aside from `HeatExchange.safe_step`'s own `Float64`-only behavior — see
`docs/control.md`).

- `current_stage(sequence, stage_state, progress, signals)` — the single stage with the highest
  weight (exact ties resolve to the earlier stage). Its return type is an inherent
  `Union` over the sequence's stage types (indexing a heterogeneous `Tuple` by a runtime index),
  not a type-stability defect.
- `stage_value(sequence, stage_state, progress, signals, accessor)` — the weighted sum of
  `accessor(stage_traits(stage), progress, signals)` over every stage. `accessor` is called for
  *every* stage unconditionally (not just the nominal "current" one, and not skipped for
  near-zero weight), trading a few extra calls for a guaranteed concrete return type — so
  `accessor` must be valid for every stage in the sequence, not only the one currently dominant.

Use `current_stage` for a hard pick (dispatch on "what stage is this"); use `stage_value` for a
blend of stage-specific numeric/array/`Unitful` data across a graded transition.

## No latching: stateless recomputation

`StageSequence` computes stage identity fresh from current controller levels on every call — it
holds no latch/history state. `stage_weights`/`current_stage` therefore already recompute
correctly if a driving level decreases and later recovers (DEB theory's "rejuvenation": maturity
can decrease under starvation). What this design does **not** support: a genuinely
history-dependent query like "has this individual ever reached puberty, even if maturity has
since regressed" — that needs additional latch state, added externally by the host; it is not
part of `StageSequence`'s contract.

Interpreting the result as *ordered life-cycle milestones* (stage `N` implies stages `1..N-1` were
already passed) requires the host to ensure monotonic, biologically-ordered controller levels —
this is the caller's responsibility, not something `StageSequence` checks or enforces.

## `StageSequence` vs. survival/death

`StageSequence` transitions move an individual into a *continuing* stage — development keeps
going, under different traits. A survival/death criterion (deterministic, like a hard temperature
or desiccation limit, or probabilistic, like a cumulative-hazard senescence model) *terminates*
the individual — nothing continues afterward. These are different things regardless of whether the
termination criterion itself is a hard threshold or a hazard process; `StageSequence` does not
attempt to represent termination as a stage. A complete simulation runs both in parallel:
`StageSequence` picks which traits/governing equations apply; an independent survival model
decides whether the individual is still alive to have them applied at all.

## Conditions driven by a different state variable

A transition's condition need not be driven by the sequence's own `progress` argument — a
`RawSignal`-based condition reads an independently-tracked state variable instead (DEB theory's
separate aging/damage accumulator driving senescence, as opposed to the maturity variable driving
ordinary stage transitions, is the motivating case). `progress` and `signals` are both threaded
through to every `controller_level` call; which one a given transition's condition actually reads
is a per-condition choice, not fixed by `stages/` itself.

## `stage_conditions` and ODE integration

```julia
stage_conditions(sequence)  # -> flattened Tuple of (stage_state, progress, signals) -> Real closures
```

Returns bare condition closures with **no attached transition identity** — mapping a fired
callback index back to a transition is the host's responsibility, exactly like
`arrest_conditions`. A minimal host-side pattern (see `test/stages.jl`, test 14, for a complete
worked example against a real `VectorContinuousCallback`):

```julia
conditions = stage_conditions(sequence)
function condition!(out, u, t, integrator)
    for (i, cond) in enumerate(conditions)
        out[i] = cond(stage_state, u[1], signals)
    end
end
function affect!(integrator, i)
    # `i` is the position of the fired condition in `conditions`, in the same
    # order stage_conditions built them (transition order, then each
    # transition's own trigger_conditions order) -- the host owns this mapping.
end
```
