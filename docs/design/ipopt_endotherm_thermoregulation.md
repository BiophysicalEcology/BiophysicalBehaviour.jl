# IPOPT-Based Endotherm Thermoregulation

**Source files:**
- `src/endotherm/thermoregulation/ipopt.jl` — shared direct-Ipopt + Enzyme solver machinery (cache, callbacks, solve loop, `thermoregulate` dispatch)
- `src/endotherm/thermoregulation/multipart_nlp.jl` — the `MultipartNLP` strategy: variable/residual templates, objective, bounds, scaling, output assembly
- `src/endotherm/endotherm_traits.jl` — `ThermoregulationLimits` struct and objective-weight fields
- `HeatExchange.jl/src/nlp_interface.jl` — `NLPStrategy` abstraction and the `nlp_pack` entry point; physics residuals come from `part_surface.jl`

> **Recommended reading:** readers unfamiliar with the effectors, heat balance loop, and
> biological priority hierarchy should read
> [rulebased_endotherm_thermoregulation.md](rulebased_endotherm_thermoregulation.md) first.
> The rule-based approach is the simpler of the two and provides the biological context that
> makes the IPOPT formulation meaningful.

---

## Background

### The optimisation problem

An endotherm facing thermal stress has several physiological and behavioural options: erect or
flatten its fur, change its posture, dilate or constrict peripheral blood vessels, pant, or
sweat. Each option has a cost (metabolic energy, water loss, time) and each shifts the heat
balance by a different amount. The question of which combination is best — or even feasible —
is not obvious when several effectors interact non-linearly through the same heat balance
equations.

One answer is to apply effectors in a fixed biological priority order, one at a time, until
balance is restored. That is what `RuleBasedSequentialControl` does and it is intuitive and
transparent. The alternative explored here is to let a mathematical optimiser find the best
combination simultaneously, without prescribing the order. In this framing, objective weights
take the place of hard-wired priority rules: setting `panting_weight` lower than
`skin_wetness_weight`, for example, makes the solver prefer panting before sweating — the
same effect as placing panting earlier in the sequential rule-based sequence, but with a smooth
trade-off rather than a step-wise one. This makes the approach well suited to sensitivity
analyses and to species where the biological priority order is uncertain.

Formally, we treat the effectors (insulation depth, body shape, flesh conductivity, panting
rate, skin wetness, metabolic rate) as the **primal variables** — the quantities the solver
actually adjusts. We define an **objective function** that measures how bad the current state
is: primarily how far core temperature is from the setpoint, with smaller penalties for active
cooling (panting, sweating). The solver searches for the primal variable values that minimise
this objective. The primal variables and the objective are sometimes called the *primal problem*.

### Constraints and feasibility

The solver cannot simply minimise the objective without restriction — it must respect physics.
The heat balance equations impose hard **equality constraints**: net heat flow at every surface
must be zero (energy is conserved). An additional **inequality constraint** enforces that
metabolic heat generation cannot fall below its Q10-scaled minimum even when the animal is hot.

These constraints define the **feasible region** — the set of effector combinations that are
physically possible. Solutions outside this region are not just suboptimal; they are physically
meaningless. Naive gradient descent ignores constraints and often wanders outside the feasible
region. IPOPT (Interior Point OPTimizer) avoids this by adding a **barrier function**: a term
that grows toward infinity as any constraint boundary is approached. This transforms each hard
constraint into a smooth hill in the objective landscape that the solver is strongly repelled
from. The feasible region is not enforced as a wall the solver bounces off, but as a geometry
embedded in the search space — the barrier reshapes the hills and valleys of the objective so
that the unconstrained minimum always lies strictly inside the feasible interior. As the barrier
is gradually weakened across iterations, the solution migrates toward the true constrained
optimum. This is why the method is called *interior point*: the solver remains inside feasible
space throughout, approaching the boundary only in the limit.

### Derivatives, Jacobians, and Hessians

To descend toward the minimum efficiently, the solver needs to know the **slope** of the
objective and the slope of each constraint residual, with respect to each primal variable. A
single slope in one dimension becomes a vector of partial derivatives in many dimensions. When
there are multiple outputs (several constraint residuals) and multiple inputs (several
effectors), the full table of partial derivatives — one row per output, one column per input —
is called the **Jacobian matrix**. It answers the question: "if I nudge flesh conductivity by a
tiny amount, how does each heat-balance residual change?"

Knowing the slope tells you which direction is downhill, but not how steeply the landscape
curves. The **Hessian matrix** captures that curvature: it is the table of second partial
derivatives of the objective (how the slope itself changes as you move). With curvature
information, the solver can take much larger, more confident steps — analogous to knowing not
just whether a hillside tilts but how quickly it steepens. Computing the full Hessian naively
requires evaluating the objective many times with small perturbations to every variable; here it
is computed analytically via automatic differentiation (see below).

### Lagrange multipliers and the dual problem

A convenient way to handle constraints is to fold them into a single combined function called
the **Lagrangian**: the objective plus a weighted sum of the constraint residuals, where each
weight is called a **Lagrange multiplier** (or dual variable). At the true solution, the
gradient of the Lagrangian with respect to the primal variables is zero — meaning the objective
gradient is exactly balanced by the constraint gradients, scaled by the multipliers. The
multipliers therefore measure the *shadow value* of each constraint: how much better the
objective could be if that physical constraint were slightly relaxed. A large multiplier on the
energy-balance constraint, for instance, means the solution is tightly limited by energy
conservation at that point. The set of Lagrange multipliers is called the *dual* of the primal
problem, and IPOPT tracks both primal and dual variables simultaneously across iterations.

Storing the dual variables (multipliers) at the end of a solve is what makes **warm-starting**
possible: when the next environmental condition is only slightly different, both the primal
effector values and the dual multipliers from the previous solve are excellent starting
guesses, and IPOPT converges in far fewer iterations than it would from a cold start.

### Automatic differentiation

Computing the Jacobian and Hessian analytically by hand for a complex heat balance model is
impractical. Symbolic differentiation (like a computer algebra system) is feasible but
produces unwieldy expressions. Finite differences (perturbing each variable by a small amount
and measuring the change in output) are simple but slow and accumulate numerical error.

**Automatic differentiation** (autodiff) is a third approach that exploits the fact that every
computer function, no matter how complex, is ultimately a composition of elementary operations
(addition, multiplication, exponentiation, etc.) for which exact derivatives are known. By
applying the chain rule mechanically at each operation as the function executes, the computer
produces exact derivatives (to floating-point precision) at the cost of roughly the same
computation as the function itself. Here, the Enzyme.jl library performs this on the Julia
heat-balance code directly — no algebraic simplification and no finite-difference approximation.

A particularly important property of Julia is that Enzyme operates on compiled Julia code
regardless of which package it came from. The heat balance functions live in `HeatExchange.jl`,
the geometry functions in `BiophysicalGeometry.jl`, and the thermoregulation control logic here
in `BiophysicalBehaviour.jl` — yet Enzyme can differentiate through all of them in a single
pass, following the computation across package boundaries without any special instrumentation.
In most other scientific computing languages, automatic differentiation is limited to code
within a single library or requires wrapping every external call by hand; Julia's design makes
multi-package AD a natural consequence of how the language compiles.


---

## Overview

The IPOPT thermoregulation solver finds the physiologically optimal combination of behavioural
and physiological effectors (piloerection, vasodilation, panting, sweating, body shape) that
minimises thermal discomfort, given the physical constraints imposed by the heat balance equations.

It is an alternative to `RuleBasedSequentialControl`, which applies effectors in a fixed priority
order. Rather than following a rigid sequence of rules, the IPOPT solver finds the globally
optimal combination in a single nonlinear program (NLP). Both approaches use the same underlying
biophysical models from `HeatExchange.jl`.

The entry point is:

```julia
thermoregulate(Endotherm(), IPOPTControl(), organism, environment, init)
```

where `init` is a NamedTuple with fields `metabolic_heat_flow`, `skin_temperature`, and
`insulation_temperature` (see `initial_physiological_state` for the standard defaults).

For repeated solves over a sweep (e.g., across an air temperature range), use
`IPOPTSolverCache` to enable primal+dual warm-starting and reuse callback allocations:

```julia
cache = IPOPTSolverCache(IPOPTControl(), organism, environment, init)
out   = thermoregulate(Endotherm(), IPOPTControl(), organism, environment, init; cache)
```

---

## Control-Theory Framing

Following standard control theory, the NLP decision variables split into two categories.
Some are shared across the whole organism; others exist once *per body part*, so the count
scales with the organism (`3 + 4·N` variables for `N` parts).

**Control variables (effectors, u)** — things the animal actively adjusts:
| Variable | Scope | Symbol | Description |
|---|---|---|---|
| Metabolic heat generation | whole-organism | `log(metabolic_heat_flow)` | Log-space for numerical scaling; exponentiated inside the solver |
| Panting rate | whole-organism | `panting_rate` | Multiplier on baseline ventilation (1 = resting) |
| Flesh conductivity | per-part | `flesh_conductivity` | Tissue conductivity (vasoconstriction/vasodilation) |
| Skin wetness | per-part | `skin_wetness` | Fractional wetted surface area (sweating/cutaneous evaporation) |

Geometry effectors (piloerection `insulation_depth`, uncurling `axis_ratio`) are **deferred** —
they would rebuild each part's body inside the differentiated residual, which the multi-part
rule-based ladder also avoids. The per-part `setups` are therefore fixed across the solve.

**State variables (x)** — outcomes determined by the effectors and the heat balance:
| Variable | Scope | Description |
|---|---|---|
| Core temperature | whole-organism | `core_temperature` (K) — one regulated core shared across all parts |
| Skin temperature | per-part | `skin_temperature` (K) |
| Insulation surface temperature | per-part | `insulation_temperature` (K) |

In principle, state variables are fully determined by the effectors via the heat balance
constraints. Including them explicitly as NLP variables (with equality constraints to enforce
consistency) improves numerical stability and allows IPOPT to find feasible points more reliably.

---

## NLP Strategy

There is one NLP formulation, `MultipartNLP()` (the default `IPOPTControl.nlp_strategy`). It
solves a genuine multi-part organism: one regulated core shared across all parts, per-part
surface temperatures and effectors, and a single whole-organism respiration balance closed at
the lung. A single `Body` is the one-part case; a `CompositeBody` (torso/head/legs) adds two
equality constraints per part.

| Variables | Constraints |
|---|---|
| `3 + 4·N` (`N` = number of parts) | `2·N + 2` (`2·N` per-part equalities + 1 whole-organism equality + 1 Q10 inequality) |

The decision-variable layout is **not** hand-indexed: it is a nested `NamedTuple` template with
`Unitful` leaves, and `Flatten.flatten` / `Flatten.reconstruct` bridge it to the flat
`Vector{Float64}` Ipopt optimises. Variables, bounds, scaling and residuals are all templates of
the same shape flattened in the same order, so a bound can never drift out of alignment with the
variable it bounds, and adding a part or a per-part effector is a change to the *structure*, never
to a pile of `x[3 + 4i]` offsets.

---

## Constraints

For an `N`-part organism the residual vector is `2·N + 2` long: two equalities per part, then the
whole-organism respiration balance, then the Q10 inequality (last, so it is the single `≥ 0`
constraint).

**Per part (2·N equalities, = 0)** — from `HeatExchange.part_surface_residuals`:

1. **Surface balance** (W): `residual_energy_balance − residual_internal_conduction = 0` — pure
   surface heat exchange (solar, longwave, convection, conduction, net internal heat flow);
   metabolic and respiration terms cancel algebraically, so respiration is counted once
   whole-organism (below) rather than per part.
2. **Skin temperature consistency** (K): skin temperature is consistent with the core–skin
   gradient given flesh conductivity and body geometry.

**Whole-organism energy balance (= 0)**

```
metabolic_heat_flow − respiration_heat_flow = Σ_parts net_metabolic_heat_internal
```

Respiration is closed once at the lung, using the mean skin temperature across parts for the lung
temperature. This mirrors `solve_coupled_metabolic_rate` at its converged point.

**Q10 metabolic scaling inequality (≥ 0)**

```
generated_heat_flow ≥ minimum_heat_flow × Q10^((core_temperature − setpoint) / 10)
```

This enforces that minimum metabolic rate rises with core temperature during hyperthermia,
consistent with the Q10 temperature coefficient of biochemical reactions. In cold conditions it
has no effect (generated_heat_flow rises naturally to close the energy deficit). In hot conditions
it ensures the solver does not suppress generated_heat_flow below its temperature-corrected minimum.

The Q10 value is taken from `MetabolismParameters.q10` in `HeatExchange.jl`.

---

## Objective Function

The objective minimises thermal discomfort by penalising deviations from setpoint and use of
active cooling effectors. All terms are normalised to [0, 1] so the objective weights are directly
comparable:

```
J = w_core × ((core_temperature − setpoint) / Δcore_temperature)²
  + w_met  × ((generated_heat_flow − heat_flow_min) / Δheat_flow)²   ← regularisation only
  + w_grad × ((core_temperature − skin_temperature − Δ*) / Δcore_temperature)²  ← optional; default disabled
  + w_pant × ((pant − 1) / Δpant)²
  + w_wet  × ((skin_wetness − skin_wetness_min) / Δskin_wetness)²
```

where `Δ` denotes the range of each variable over its physiological limits.

### Why generated_heat_flow is not the primary target

`generated_heat_flow` is a state variable. In cold conditions, keeping `core_temperature` near
setpoint already forces `generated_heat_flow` upward (less metabolic heat → larger energy deficit
→ colder core → objective penalty). No explicit penalty on `generated_heat_flow` is needed
for thermogenesis.

In hot conditions, the Q10 inequality constraint (constraint 4) prevents `generated_heat_flow`
from dropping below its temperature-scaled minimum, so the objective cannot drive it to zero.

The small `metabolic_heat_weight` (default 0.1) acts as a **regularisation** term only: it
breaks an otherwise degenerate manifold in cold conditions where combinations of high panting
and high `generated_heat_flow` satisfy the energy balance equally well. A small weight is
enough; the Q10 constraint takes over in hot conditions.

### Objective weights and their effects

The weights are stored as fields of `ThermoregulationLimits`:

| Field | Default | Effect |
|---|---|---|
| `core_temperature_weight` | 1.0 | Lower → core temperature allowed to deviate more before effectors are exhausted |
| `metabolic_heat_weight` | 0.1 | Small regularisation to prevent high-panting/high-generated_heat_flow degeneracy in cold |
| `panting_weight` | 1.0 | Lower → panting activates sooner |
| `skin_wetness_weight` | 1.0 | Higher than `panting_weight` → panting before sweating (birds/rabbits); lower → sweating first (humans) |
| `gradient_weight` | 0.0 | Non-zero → penalises deviation from `target_core_skin_gradient` (K); disabled by default |
| `target_core_skin_gradient` | 2.0 | Target `core_temperature − skin_temperature` (K); only used when `gradient_weight > 0` |

---

## Comparison with Rule-Based Sequential Control

| Aspect | `RuleBasedSequentialControl` | `IPOPTControl` |
|---|---|---|
| **Philosophy** | Fixed priority order; apply one effector at a time until balance is achieved | Simultaneous optimisation across all effectors |
| **Cold response** | Piloerect → curl → vasoconstrict → thermogenesis (strict order) | All effectors adjusted together; `metabolic_heat_weight` regularisation biases toward the sequential ordering |
| **Hot response** | Vasodilate → allow hyperthermia → pant → sweat (strict order) | Order emerges from relative objective weights |
| **Panting vs sweating** | Set by `mode` (e.g. `CorePantingSweatingFirst`) | Set by `panting_weight` vs `skin_wetness_weight` relative magnitudes |
| **Q10 scaling** | Applied explicitly: `generated_heat_flow_min` rises with Q10 at each step | Enforced via inequality constraint 4 |
| **Core temperature** | Set explicitly by `hyperthermia()` when energy balance cannot close | Free to rise within bounds; penalised by `core_temperature_weight` |
| **Speed** | Fast (iterative, closed-form steps) | Comparable in practice; both approaches take similar wall time per temperature point |
| **Tuning** | Logic rules and step sizes | Five scalar objective weights |
| **Guarantee** | Always produces a result; may not be globally optimal | Returns the NLP optimum within bounds; may not converge in extreme cases |

---

## Solver Configuration

IPOPT is invoked directly via `Ipopt.jl` (not through `Optimization.jl` / `OptimizationIpopt`).
Going direct exposes `mult_g`, `mult_x_L`, `mult_x_U` on `IpoptProblem`, which is what enables
dual warm-starting across solves.

Key solver options applied per call:

```julia
Ipopt.AddIpoptNumOption(prob, "acceptable_tol",        1e-3)
Ipopt.AddIpoptIntOption(prob, "acceptable_iter",       5)
Ipopt.AddIpoptNumOption(prob, "tol",                   1e-4)
Ipopt.AddIpoptIntOption(prob, "max_iter",              300)
Ipopt.AddIpoptStrOption(prob, "mu_strategy",           "adaptive")
Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method",    "user-scaling")
```

**Exact Hessian (forward-over-reverse Enzyme).** The solver registers a real `evaluate_hessian`
callback computed via Enzyme forward-over-reverse automatic differentiation on the Lagrangian. IPOPT uses this exact
second-order information (dense lower-triangle) rather than an L-BFGS approximation. This is
why `hessian_approximation` is left at its default `"exact"` — no override is applied.

**Adaptive barrier parameter** (`mu_strategy = "adaptive"`): IPOPT picks the barrier parameter
`mu` per-iteration based on the optimality-error / complementarity ratio instead of the default
`monotone` schedule. On smooth NLPs (heat balance with continuous physical relationships) this
typically converges in fewer iterations.

**User-supplied variable scaling** (`nlp_scaling_method = "user-scaling"`): per-variable scale
factors are set via `Ipopt.SetIpoptProblemScaling`. Temperatures scale by 1/300, small
quantities like skin_wetness and insulation_depth by 50, and log(metabolic_heat_flow) by 1.
This gives the solver a well-conditioned internal problem from iteration zero and skips IPOPT's
own gradient-based scaling heuristic.

**Differentiation strategy.** Gradients and constraint Jacobians are computed via
[`Enzyme.jl`](https://github.com/EnzymeAD/Enzyme.jl): reverse-mode for the scalar objective
gradient and the Lagrangian inner pass; forward-mode for the Jacobian rows and outer Hessian
columns. All five Ipopt callbacks (`evaluate_objective`, `evaluate_constraints`,
`evaluate_objective_gradient`, `evaluate_constraint_jacobian`, `evaluate_hessian`) are
implemented as named functor structs (`EvaluateObjective`, `EvaluateConstraints`, etc.)
subtyping `Function`, which keeps type stability explicit and makes stack traces readable.

---

## Warm-Start Cache

`IPOPTSolverCache` reuses callback functors and scratch buffers across a sweep and enables
primal+dual warm-starting from the previous solve:

```julia
cache = IPOPTSolverCache(control, organism, environment, init)

for (air_temperature, ...) in ...
    # update organism and environment ...
    out = thermoregulate(Endotherm(), IPOPTControl(), organism, environment, init; cache)
    init = (;
        metabolic_heat_flow    = out.energy_flows.metabolic_heat_flow,
        skin_temperature       = out.thermoregulation.skin_temperature,
        insulation_temperature = out.thermoregulation.insulation_temperature,
    )
end
```

The cache warm-starts **both** primal and dual:

- **Primal**: the previous solve's full effector vector (all 9 or 11 variables, clamped into
  the new per-call bounds) is fed in as the next solve's initial `x`. The fresh (no-cache) path
  seeds only the three `init` fields; the cache seeds all decision variables from the prior
  solution — a strictly better starting point when the trajectory changes smoothly.
- **Dual**: Lagrange multipliers from the previous solve (`mult_g` on constraints, `mult_x_L` /
  `mult_x_U` on bounds) are restored into the next `IpoptProblem` and IPOPT is instructed to
  use them via `warm_start_init_point = "yes"`.

Because IPOPT bakes bounds into the C-side `IpoptProblem` at construction time, a new
`IpoptProblem` is created each call (cheap: one malloc + struct init), while Julia-side
primal/dual state is carried in the cache across rebuilds.

**Invalidate the cache** (build a new one) if the organism's body shape, NLP strategy, or
smoothing strategy changes — the functors bake in those types. Routine per-iteration changes
(environment temperatures, q10, relative humidity) reuse the cache without rebuilding.

Use `reset_warm_start!(cache)` to discard the stored solution and force a cold start on the
next solve — useful after a large discontinuous step in environmental conditions.

The cache must not be shared between threads.

---

## Initial Conditions

The solver requires initial values for the decision variables. The standard helper is:

```julia
init = initial_physiological_state(organism, environment_vars)
# returns (; metabolic_heat_flow, skin_temperature, insulation_temperature)
```

which seeds `metabolic_heat_flow = 0 W` (clamped to `heat_flow_min` inside the solver),
`skin_temperature = setpoint − 3 K`, and `insulation_temperature = air_temperature`.

When running a sweep, **carry-forward initialisation** is recommended: the solution from each
step becomes the starting point for the next, which is the natural choice for a smooth sequence:

```julia
init = (;
    metabolic_heat_flow    = minimum_heat_flow,
    skin_temperature       = core_temperature - 3.0u"K",
    insulation_temperature = air_temperatures[1],
)

for (air_temperature, ...) in zip(air_temperatures, ...)
    # ... build organism and environment ...
    out = thermoregulate(Endotherm(), IPOPTControl(), organism, environment, init)

    init = (;
        metabolic_heat_flow    = out.energy_flows.metabolic_heat_flow,
        skin_temperature       = out.thermoregulation.skin_temperature,
        insulation_temperature = out.thermoregulation.insulation_temperature,
    )
end
```

With `IPOPTSolverCache`, warm-starting is handled internally: the full prior solution vector
is restored automatically, so manually threading `init` through the loop is not needed for
warm-starting (though passing updated `init` still sets the fallback cold-start point).

---

## HeatExchange.jl NLP Interface

The IPOPT solver does **not** call `solve_metabolic_rate` (or `solve_coupled_metabolic_rate`).
Those internally iterate to find consistent temperatures given fixed physiological parameters; the
optimiser cannot use them as constraints because the temperatures *are* decision variables —
passing a trial `core_temperature` in would trigger a nested solve that ignores the optimiser's
current guess. Instead it evaluates the physics *non-iteratively* at the optimiser's trial point.

- **`nlp_pack(strategy, organism, environment, skin_temperature_init, insulation_temperature_init; smoothing)`**
  (HeatExchange, extended in `multipart_nlp.jl`) — packs the fixed per-part physics once, returning
  a `MultipartNLPPacked` with the per-part `part_surface_residuals` setups, the whole-organism
  respiration inputs, and the `Unitful` variable template. `NLPStrategy` and `nlp_pack` are the
  only names HeatExchange exposes for this; the strategy and all solver policy live in
  BiophysicalBehaviour.

- **`HeatExchange.part_surface_residuals(setup, core_temperature, skin_temperature, insulation_temperature,
  metabolic_heat_flow; ...)`** — evaluates one part's surface-balance and skin-temperature residuals
  for a trial set of decision variables, without any iterative solve. `_heat_balance_residuals!`
  walks the parts calling this, then closes the whole-organism respiration balance once. This is the
  same primitive the iterative coupled solver uses, so no physics is duplicated.

Output assembly is `_assemble(::MultipartNLPPacked, …)` (`multipart_nlp.jl`): it reconstructs the
converged `x` into named per-part fields and evaluates `part_surface_residuals` once more at the
solution for the final per-part heat flows.

---

## Example Usage

See `examples/budgerigar.jl` for a complete example comparing both solvers across 0–50°C.
`RuleBasedSequentialControl` runs first (independently), and the IPOPT solver runs in a separate
pass using carry-forward initialisation. The two solvers are entirely independent — neither is
used as input to the other.

```julia
# RuleBasedSequentialControl pass (resets init every call — no carry-forward)
for (air_temperature, relative_humidity, q10) in zip(air_temperatures, ...)
    init = (;
        metabolic_heat_flow    = 0.0u"W",
        skin_temperature       = core_temperature - 3.0u"K",
        insulation_temperature = air_temperature,
    )
    out = thermoregulate(organism_loop, environment_loop, init)
    push!(results, ...)
end

# IPOPT pass (carry-forward initialisation across the sweep)
metabolic_heat_flow_ipopt    = min_metabolic_heat_flow
skin_temperature_ipopt       = core_temperature - 3.0u"K"
insulation_temperature_ipopt = air_temperatures[1]

for (air_temperature, relative_humidity, q10) in zip(air_temperatures, ...)
    out = thermoregulate(
        Endotherm(), IPOPTControl(), organism_loop, environment_loop,
        (;
            metabolic_heat_flow    = metabolic_heat_flow_ipopt,
            skin_temperature       = skin_temperature_ipopt,
            insulation_temperature = insulation_temperature_ipopt,
        ),
    )
    global metabolic_heat_flow_ipopt    = out.energy_flows.metabolic_heat_flow
    global skin_temperature_ipopt       = out.thermoregulation.skin_temperature
    global insulation_temperature_ipopt = out.thermoregulation.insulation_temperature
    push!(ipopt_results, ...)
end
```

Typical tuning guidance for objective weights:

- **Birds (panting-first, no sweating):** `panting_weight = 0.1`, `skin_wetness_weight = 0.1`
  (skin_wetness_max near zero, so the wetness penalty matters little)
- **Sweating mammals (humans):** `skin_wetness_weight < panting_weight`
- **Panting mammals (rabbits, dogs):** `skin_wetness_weight > panting_weight`
- **Cold-temperature fit too high:** increase `metabolic_heat_weight` (stronger bias against
  thermogenesis before behavioral effectors are exhausted)

---

## Limitations and Future Work

- **Variable scaling:** `_scaling` authors the per-variable Ipopt scale factors in the variable
  *structure* and flattens it, so they stay aligned with the variables by construction (no index
  drift). The scale *magnitudes* are still hand-chosen constants; a future refinement could derive
  them from each variable's bounds.
- **Per-part sub-side symmetry:** each part is a lumped body with a single skin and insulation
  temperature, so asymmetric loading *within* a part (sun above, ground below) is averaged rather
  than resolved. A multi-part organism can capture asymmetry *between* parts, but not within one.
- **Deferred geometry effectors:** piloerection (`insulation_depth`) and uncurling (`axis_ratio`)
  are not decision variables — they would rebuild each part's body inside the differentiated
  residual. Adding them requires a per-part body-rebuild path that Enzyme can differentiate.
- **Global optimality:** IPOPT finds a local optimum of the NLP. For well-posed problems the
  objective is approximately convex and the local optimum is unique, but unusual initial
  conditions or extreme parameter combinations may yield suboptimal solutions.
- **Unitful and automatic differentiation:** `Unitful.jl` attaches physical units to every
  quantity as part of its type. Enzyme propagates derivatives by pairing each number with a
  companion derivative value (a "dual number"), and Unitful's type machinery does not compose
  cleanly with these dual numbers — a derivative of a temperature is not itself a temperature,
  so the units become ambiguous. The solution is a hard boundary at the interface between the
  optimizer and the physics: the optimizer works entirely in plain `Float64` (no units), hands
  stripped numbers to the heat-balance functions, which reattach units internally, then strip
  units from the residuals before returning plain `Float64` back to Enzyme. This boundary is
  enforced in `_heat_balance_residuals!`, the only function that crosses between the two worlds.
  The result is that published scientific units are preserved everywhere they matter while the
  autodiff machinery operates cleanly on bare numbers.
