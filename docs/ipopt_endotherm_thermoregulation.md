# IPOPT-Based Endotherm Thermoregulation

**Source files:**
- `src/endotherm/thermoregulation/ipopt.jl` — NLP formulation, solver setup, output assembly
- `src/endotherm/endotherm_traits.jl` — `ThermoregulationLimits` struct and penalty fields
- `HeatExchange.jl/src/nlp_interface.jl` — `nlp_pack`, `nlp_residuals`, `nlp_assemble_output` functions used by the NLP solver

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

---

## Control-Theory Framing

Following standard control theory, the nine NLP decision variables are split into two categories:

**Control variables (effectors, u)** — things the animal actively adjusts:
| Variable | Symbol | Description |
|---|---|---|
| Metabolic heat generation | `log(generated_heat_flow)` | Log-space for numerical scaling; exponentiated inside the solver |
| Flesh conductivity | `flesh_conductivity` | Tissue conductivity (vasoconstriction/vasodilation) |
| Panting rate | `pant` | Multiplier on baseline ventilation (1 = resting) |
| Skin wetness | `skin_wetness` | Fractional wetted surface area (sweating/cutaneous evaporation) |
| Insulation depth | `insulation_depth` | Fur/feather erection depth (piloerection) |
| Body shape | `aspect_ratio_factor` | Ellipsoid aspect ratio (curling/elongating posture) |

**State variables (x)** — outcomes determined by the effectors and the heat balance:
| Variable | Description |
|---|---|
| Core temperature | `core_temperature` (K) |
| Skin temperature | `skin_temperature` (K) |
| Insulation surface temperature | `insulation_temperature` (K) |

In principle, state variables are fully determined by the effectors via the heat balance
constraints. Including them explicitly as NLP variables (with equality constraints to enforce
consistency) improves numerical stability and allows IPOPT to find feasible points more reliably.

---

## NLP Strategies

Two NLP formulations are available via `IPOPTControl(nlp_strategy=...)`:

| Strategy | Decision variables | Constraints | Description |
|---|---|---|---|
| `WeightedMeanNLP()` | 9 | 4 (3 equality + 1 Q10) | Single mean-weighted body; dorsal/ventral merged by view-factor weights |
| `MultiSidedNLP()` | 11 | 6 (5 equality + 1 Q10) | Explicit dorsal and ventral sides; two skin temperatures, two insulation temperatures |

`WeightedMeanNLP` is the default. `MultiSidedNLP` gives more accurate results under strongly
asymmetric conditions (high solar loading) at the cost of a slightly larger NLP.

---

## Constraints

### `WeightedMeanNLP`: 4 constraints

**Constraints 1–3: Heat balance equalities (= 0)**

Provided by `HeatExchange.nlp_residuals(::WeightedMeanNLPPacked, ...)`, which returns three residuals:

1. **Energy balance residual** (W): net heat flow at the insulation surface must be zero.
2. **Internal conduction residual** (W): heat flow from core to skin through flesh and fat must
   match the heat flow from skin through insulation to the environment.
3. **Skin temperature consistency** (K): skin temperature must be consistent with the
   core–skin gradient given flesh conductivity and body geometry.

**Constraint 4: Q10 metabolic scaling inequality (≥ 0)**

```
generated_heat_flow ≥ minimum_heat_flow × Q10^((core_temperature − setpoint) / 10)
```

### `MultiSidedNLP`: 6 constraints

**Constraints 1–2: Dorsal surface physics (= 0)**

1. **Dorsal surface balance** (W): `residual_energy_balance_d − residual_internal_conduction_d = 0`
   — pure surface heat exchange (solar, longwave, convection, conduction, net internal heat flow);
   metabolic and respiration terms cancel algebraically.
2. **Dorsal skin temperature** (K): `residual_skin_temperature_d = 0`

**Constraints 3–4: Ventral surface physics (= 0)**

3. **Ventral surface balance** (W): `residual_energy_balance_v − residual_internal_conduction_v = 0`
4. **Ventral skin temperature** (K): `residual_skin_temperature_v = 0`

**Constraint 5: Whole-organism energy balance (= 0)**

```
metabolic_heat_flow − respiration_heat_flow = dmult × net_metabolic_heat_internal_d + vmult × net_metabolic_heat_internal_v
```

where `dmult = sky_view_factor + vegetation_view_factor` and `vmult = 1 − dmult`. This mirrors
the validated rule-based multi-sided solver criterion exactly.

**Constraint 6: Q10 metabolic scaling inequality (≥ 0)**

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
active cooling effectors. All terms are normalised to [0, 1] so the penalty weights are directly
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

The small `metabolic_heat_penalty` (default 0.1) acts as a **regularisation** term only: it
breaks an otherwise degenerate manifold in cold conditions where combinations of high panting
and high `generated_heat_flow` satisfy the energy balance equally well. A small weight is
enough; the Q10 constraint takes over in hot conditions.

### Penalty weights and their effects

The weights are stored as fields of `ThermoregulationLimits`:

| Field | Default | Effect |
|---|---|---|
| `core_temperature_penalty` | 1.0 | Lower → core temperature allowed to deviate more before effectors are exhausted |
| `metabolic_heat_penalty` | 0.1 | Small regularisation to prevent high-panting/high-generated_heat_flow degeneracy in cold |
| `panting_penalty` | 1.0 | Lower → panting activates sooner |
| `skin_wetness_penalty` | 1.0 | Higher than `panting_penalty` → panting before sweating (birds/rabbits); lower → sweating first (humans) |
| `gradient_penalty` | 0.0 | Non-zero → penalises deviation from `target_core_skin_gradient` (K); disabled by default |
| `target_core_skin_gradient` | 2.0 | Target `core_temperature − skin_temperature` (K); only used when `gradient_penalty > 0` |

---

## Comparison with Rule-Based Sequential Control

| Aspect | `RuleBasedSequentialControl` | `IPOPTControl` |
|---|---|---|
| **Philosophy** | Fixed priority order; apply one effector at a time until balance is achieved | Simultaneous optimisation across all effectors |
| **Cold response** | Piloerect → curl → vasoconstrict → thermogenesis (strict order) | All effectors adjusted together; `metabolic_heat_penalty` regularisation biases toward the rule-based ordering |
| **Hot response** | Vasodilate → allow hyperthermia → pant → sweat (strict order) | Order emerges from relative penalty weights |
| **Panting vs sweating** | Set by `mode` (e.g. `CorePantingSweatingFirst`) | Set by `panting_penalty` vs `skin_wetness_penalty` relative magnitudes |
| **Q10 scaling** | Applied explicitly: `generated_heat_flow_min` rises with Q10 at each step | Enforced via inequality constraint 4 |
| **Core temperature** | Set explicitly by `hyperthermia()` when energy balance cannot close | Free to rise within bounds; penalised by `core_temperature_penalty` |
| **Speed** | Fast (iterative, closed-form steps) | Comparable in practice; both approaches take similar wall time per temperature point |
| **Tuning** | Logic rules and step sizes | Five scalar penalty weights |
| **Guarantee** | Always produces a result; may not be globally optimal | Returns the NLP optimum within bounds; may not converge in extreme cases |

---

## Solver Configuration

IPOPT is invoked via `OptimizationIpopt.jl` (part of the `Optimization.jl` ecosystem):

```julia
solve(optimization_prob,
    IpoptOptimizer(;
        hessian_approximation = "limited-memory",   # L-BFGS; avoids O(n²) Hessian evaluations
        acceptable_tol        = 1e-3,
        acceptable_iter       = 5,
    );
    reltol   = 1e-4,
    maxiters = 300,
)
```

**L-BFGS Hessian approximation** (`hessian_approximation="limited-memory"`) is critical for
performance. The exact Hessian requires O(n²) = 81 finite-difference evaluations per IPOPT
iteration; L-BFGS constructs a low-rank approximation from gradient differences instead,
reducing the cost substantially.

**Differentiation strategy.** Gradients and constraint Jacobians are computed via
[`Enzyme.jl`](https://github.com/EnzymeAD/Enzyme.jl): reverse-mode for the scalar objective
gradient and forward-mode for the constraint Jacobian.

---

## Initial Conditions

The IPOPT solver requires initial values for the nine decision variables. When running a sequence
of temperature points (as in the example scripts), **carry-forward initialisation** is used: the
solution from each step becomes the starting point for the next. This avoids dependence on the
rule-based solver and is the natural choice for a sequence that changes smoothly across
environmental conditions:

```julia
init = (;
    metabolic_heat_flow    = minimum_heat_flow,                # start at minimum metabolic rate
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

For a single solve (e.g., during testing), any physiologically plausible initial point works —
for example `metabolic_heat_flow = 0.0u"W"` (clamped to `heat_flow_min` inside the solver),
`skin_temperature` near `setpoint_temperature - 3 K`, and `insulation_temperature` near ambient temperature.

---

## Mean-Weighted Body Approximation

The IPOPT solver uses a single mean-weighted body (dorsal/ventral average) rather than the
full dorsal/ventral split used by `solve_metabolic_rate`. Weights are derived from view factors:

```
dorsal_weight  = sky_view_factor + vegetation_view_factor
ventral_weight = 1 − dorsal_weight
```

Insulation depth, fibre properties, view factors, and conduction coefficients are all
weighted averages of the dorsal and ventral values. This is the same weighting used internally
by `solve_metabolic_rate` and keeps the NLP problem dimensionality manageable.

---

## HeatExchange.jl NLP Interface

The IPOPT solver does **not** call `solve_metabolic_rate`. Instead it uses three functions from
`HeatExchange.jl/src/nlp_interface.jl`:

- **`nlp_pack(strategy, organism, environment, skin_temperature_init, insulation_temperature_init)`** —
  pre-computes all geometry and environment quantities that are fixed for a given hour, returning a
  packed parameter struct (`WeightedMeanNLPPacked` or `MultiSidedNLPPacked`).

- **`nlp_residuals(packed, core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow, ...)`** —
  evaluates the heat balance residuals for a trial set of decision variables without any iterative
  solve. Returns a `NamedTuple` with `residuals` (tuple of physics residuals in W/K) plus
  individual heat flow components.

- **`nlp_assemble_output(packed, organism, environment, ...)`** —
  reconstructs the full thermoregulation output `NamedTuple` from the solver solution, including
  respiration, mass flows, and energy flows in the same format as `solve_metabolic_rate`.

This distinction is critical: `solve_metabolic_rate` internally iterates to find consistent
temperatures given a fixed set of physiological parameters. The IPOPT solver cannot use it as a
constraint because the temperatures themselves are decision variables — passing a trial
`core_temperature` to `solve_metabolic_rate` would trigger a nested solve that ignores the
optimizer's current guess. `nlp_residuals` instead accepts all temperature and effector values
as explicit arguments and simply returns the residuals, leaving the root-finding entirely to IPOPT.

---

## Example Usage

See `examples/budgerigar.jl` for a complete example comparing both solvers across 0–50°C.
The rule-based solver runs first (independently), and the IPOPT solver runs in a separate pass
using carry-forward initialisation. The two solvers are entirely independent — neither is used
as input to the other.

```julia
# Rule-based pass
for (air_temperature, rh, q10) in zip(air_temperatures, ...)
    init = (; metabolic_heat_flow = 0.0u"W", skin_temperature, insulation_temperature)
    out = thermoregulate(organism, environment, init)
    push!(results, ...)
end

# IPOPT pass (carry-forward initialisation)
init = (;
    metabolic_heat_flow    = minimum_heat_flow,
    skin_temperature       = core_temperature - 3.0u"K",
    insulation_temperature = air_temperatures[1],
)

for (air_temperature, rh, q10) in zip(air_temperatures, ...)
    out = thermoregulate(Endotherm(), IPOPTControl(), organism, environment, init)
    init = (;
        metabolic_heat_flow    = out.energy_flows.metabolic_heat_flow,
        skin_temperature       = out.thermoregulation.skin_temperature,
        insulation_temperature = out.thermoregulation.insulation_temperature,
    )
    push!(ipopt_results, ...)
end
```

Typical tuning guidance for penalty weights:

- **Birds (panting-first, no sweating):** `panting_penalty = 0.1`, `skin_wetness_penalty = 0.1`
  (skin_wetness_max near zero, so the wetness penalty matters little)
- **Sweating mammals (humans):** `skin_wetness_penalty < panting_penalty`
- **Panting mammals (rabbits, dogs):** `skin_wetness_penalty > panting_penalty`
- **Cold-temperature fit too high:** increase `metabolic_heat_penalty` (stronger bias against
  thermogenesis before behavioral effectors are exhausted)

---

## Limitations and Future Work

- **Automatic differentiation:** derivatives are computed via `Enzyme.jl` (reverse-mode for the objective gradient, forward-mode for the constraint Jacobian). The Hessian and constraint Hessian callbacks are stubs because IPOPT runs with the L-BFGS Hessian approximation; supplying exact second-order information may improve convergence further.
- **Dorsal/ventral symmetry:** the mean-weighted body approximation merges dorsal and ventral sides. The full `solve_metabolic_rate` computes them separately. In strongly asymmetric conditions (e.g. high solar loading on dorsal surface) this may introduce small errors.
- **Global optimality:** IPOPT finds a local optimum of the NLP. For well-posed problems the objective is approximately convex and the local optimum is unique, but unusual initial conditions or extreme parameter combinations may yield suboptimal solutions.
