# IPOPT-Based Endotherm Thermoregulation

**Source files:**
- `src/endotherm/ipopt_thermoregulation.jl` — NLP formulation, solver setup, output assembly (branch: `IPOPT-implementation`)
- `src/endotherm/endotherm_traits.jl` — `ThermoregulationLimits` struct and penalty fields
- `HeatExchange.jl/src/endotherm/heat_balance.jl` — `heat_balance` function used as NLP constraints (branch: `IPOPT-preparation`)

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
thermoregulate(Endotherm(), IPOPTControl(), organism, environment, generated_heat_flow_init, skin_temperature_init, insulation_temperature_init)
```

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
| Core temperature | `T_core` (K) |
| Skin temperature | `T_skin` (K) |
| Insulation surface temperature | `T_insulation` (K) |

In principle, state variables are fully determined by the effectors via the heat balance
constraints. Including them explicitly as NLP variables (with equality constraints to enforce
consistency) improves numerical stability and allows IPOPT to find feasible points more reliably.

---

## Constraints

Four constraints are imposed on the NLP:

### Constraints 1–3: Heat balance equalities (= 0)

These are provided by `HeatExchange.heat_balance(T_core, T_skin, T_insulation, Q_gen; ...)`, which
returns three residuals:

1. **Energy balance residual** (W): net heat flow at the insulation surface must be zero.
2. **Internal conduction residual** (W): heat flow from core to skin through flesh and fat must
   match the heat flow from skin through insulation to the environment.
3. **Skin temperature consistency** (K): skin temperature must be consistent with the
   core–skin gradient given flesh conductivity and body geometry.

These three equality constraints make the three temperature state variables outcomes of the
six effectors, not free choices.

### Constraint 4: Q10 metabolic scaling inequality (≥ 0)

```
generated_heat_flow ≥ heat_flow_min × Q10^((core_temperature − setpoint) / 10)
```

This enforces that minimum metabolic rate rises with core temperature during hyperthermia,
consistent with the Q10 temperature coefficient of biochemical reactions. In cold conditions it
has no effect (Q_gen rises naturally to close the energy deficit). In hot conditions it ensures
the solver does not suppress Q_gen below its temperature-corrected minimum.

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

### Why Q_gen is not the primary target

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
| `core_temperature_penalty` | 1.0 | Lower → T_core allowed to deviate more before effectors are exhausted |
| `metabolic_heat_penalty` | 0.1 | Small regularisation to prevent high-panting/high-Q_gen degeneracy in cold |
| `panting_penalty` | 1.0 | Lower → panting activates sooner |
| `skin_wetness_penalty` | 1.0 | Higher than `panting_penalty` → panting before sweating (birds/rabbits); lower → sweating first (humans) |
| `gradient_penalty` | 0.0 | Non-zero → penalises deviation from `target_core_skin_gradient` (K); disabled by default |
| `target_core_skin_gradient` | 2.0 | Target T_core − T_skin (K); only used when `gradient_penalty > 0` |

---

## Comparison with Rule-Based Sequential Control

| Aspect | `RuleBasedSequentialControl` | `IPOPTControl` |
|---|---|---|
| **Philosophy** | Fixed priority order; apply one effector at a time until balance is achieved | Simultaneous optimisation across all effectors |
| **Cold response** | Piloerect → curl → vasoconstrict → thermogenesis (strict order) | All effectors adjusted together; `metabolic_heat_penalty` regularisation biases toward the rule-based ordering |
| **Hot response** | Vasodilate → allow hyperthermia → pant → sweat (strict order) | Order emerges from relative penalty weights |
| **Panting vs sweating** | Set by `mode` (e.g. `CorePantingSweatingFirst`) | Set by `panting_penalty` vs `skin_wetness_penalty` relative magnitudes |
| **Q10 scaling** | Applied explicitly: `Q_gen_min` rises with Q10 at each step | Enforced via inequality constraint 4 |
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

**Why finite differences?** `Unitful.jl` units propagate through `HeatExchange.heat_balance`,
making the residuals and objective incompatible with Julia automatic differentiation packages
(Zygote, ForwardDiff). All gradients, Jacobians, and Hessians are computed via `FiniteDiff.jl`.
The `hess` and `cons_h` callbacks are registered (required by `IpoptOptimizer`) but never
called at runtime when L-BFGS is active.

---

## Initial Conditions

The IPOPT solver requires initial values for the nine decision variables. When running a sequence
of temperature points (as in the example scripts), **carry-forward initialisation** is used: the
solution from each step becomes the starting point for the next. This avoids dependence on the
rule-based solver and is the natural choice for a sequence that changes smoothly across
environmental conditions:

```julia
generated_heat_flow_ipopt    = Q_minimum                          # start at minimum metabolic rate
skin_temperature_ipopt       = core_temperature - 3.0u"K"
insulation_temperature_ipopt = air_temperatures[1]

for (T_air, ...) in zip(air_temperatures, ...)
    # ... build organism and environment ...
    out = thermoregulate(Endotherm(), IPOPTControl(), organism, environment,
                         generated_heat_flow_ipopt, skin_temperature_ipopt, insulation_temperature_ipopt)

    generated_heat_flow_ipopt    = out.energy_flows.generated_heat_flow
    skin_temperature_ipopt       = out.thermoregulation.skin_temperature
    insulation_temperature_ipopt = out.thermoregulation.insulation_temperature
end
```

For a single solve (e.g., during testing), any physiologically plausible initial point works —
for example `generated_heat_flow_init = 0.0u"W"` (clamped to `heat_flow_min` inside the solver),
`skin_temperature_init` near `T_setpoint - 3 K`, and `insulation_temperature_init` near ambient temperature.

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

## HeatExchange.jl: `heat_balance` Function

The IPOPT solver does **not** call `solve_metabolic_rate`. Instead it calls a dedicated lower-level
function `HeatExchange.heat_balance` (added on the `IPOPT-preparation` branch), which evaluates
the heat balance residuals for an arbitrary set of trial temperatures and effectors without
performing any iterative solve of its own.

This distinction is critical: `solve_metabolic_rate` internally iterates to find consistent
temperatures given a fixed set of physiological parameters. The IPOPT solver cannot use it as a
constraint because the temperatures themselves are decision variables — passing a trial `T_core`
to `solve_metabolic_rate` would trigger a nested solve that ignores the optimizer's current
guess. `heat_balance` instead accepts all four temperature/heat values as explicit arguments and
simply returns the residuals, leaving the root-finding entirely to IPOPT.

The function signature is:

```julia
HeatExchange.heat_balance(
    core_temperature, skin_temperature, insulation_temperature, generated_heat_flow;
    body, insulation_pars, insulation, geometry_vars, environment_vars, traits,
    resp_pars, k_flesh, pant, skin_wetness,
)
```

It returns a `NamedTuple` with:
- `residual_energy_balance` (W) — net heat flow at the body surface (must equal zero at equilibrium)
- `residual_internal_conduction` (W) — internal core-to-skin conduction imbalance (must equal zero)
- `residual_skin_temperature` (K) — skin temperature consistency residual (must equal zero)
- `radiation_heat_flow`, `convection_heat_flow`, `conduction_heat_flow` — individual fluxes
- `skin_evaporation_heat_flow`, `insulation_evaporation_heat_flow`, `respiration_heat_flow`

---

## Example Usage

See `examples/budgerigar.jl` for a complete example comparing both solvers across 0–50°C.
The rule-based solver runs first (independently), and the IPOPT solver runs in a separate pass
using carry-forward initialisation. The two solvers are entirely independent — neither is used
as input to the other.

```julia
# Rule-based pass
for (T_air, rh, q10) in zip(air_temperatures, ...)
    out = thermoregulate(organism, environment, 0.0u"W", skin_temperature, T_insulation)
    push!(results, ...)
end

# IPOPT pass (carry-forward initialisation)
generated_heat_flow_ipopt    = Q_minimum
skin_temperature_ipopt       = core_temperature - 3.0u"K"
insulation_temperature_ipopt = air_temperatures[1]

for (T_air, rh, q10) in zip(air_temperatures, ...)
    out = thermoregulate(Endotherm(), IPOPTControl(), organism, environment,
                         generated_heat_flow_ipopt, skin_temperature_ipopt, insulation_temperature_ipopt)
    generated_heat_flow_ipopt    = out.energy_flows.generated_heat_flow
    skin_temperature_ipopt       = out.thermoregulation.skin_temperature
    insulation_temperature_ipopt = out.thermoregulation.insulation_temperature
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

- **Automatic differentiation:** `Unitful.jl` units propagate through `HeatExchange.heat_balance`, making it incompatible with `ForwardDiff` or `Zygote`. All derivatives are computed by finite differences via `FiniteDiff.jl`. Stripping units from the inner-loop heat balance computation would allow exact derivatives and potentially faster convergence.
- **Dorsal/ventral symmetry:** the mean-weighted body approximation merges dorsal and ventral sides. The full `solve_metabolic_rate` computes them separately. In strongly asymmetric conditions (e.g. high solar loading on dorsal surface) this may introduce small errors.
- **Global optimality:** IPOPT finds a local optimum of the NLP. For well-posed problems the objective is approximately convex and the local optimum is unique, but unusual initial conditions or extreme parameter combinations may yield suboptimal solutions.
