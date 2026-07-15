# Transient (Lumped-Capacitance) Body Temperature

**Source files:**
- `HeatExchange.jl/src/transient_lumped.jl` — `onelump`/`twolump` physics (`dT/dt`)
- `HeatExchange.jl/src/internal_temperature.jl` — `internal_gradient_shape_factor`
- `src/transient/forcing.jl` — `EnvironmentForcing` (time-varying environment)
- `src/transient/simulate.jl` — `simulate_onelump`/`simulate_twolump`
- `src/transient/ectotherm/behavioral_driver.jl` — `simulate_diurnal_behavior`

Ported from NicheMapR's `onelump.R`/`onelump_var.R`/`twolump.R`/`trans_behav.R`. Unlike the
rest of the package, which solves steady-state heat balance (`HeatExchange.heat_balance`,
zero thermal mass), this gives the organism real thermal inertia — `dT_b/dt` driven by a
lumped heat-capacitance model, integrated with `OrdinaryDiffEq.jl`.

## Physics / solver split

`HeatExchange.jl` exposes `onelump`/`twolump` as pure, solver-agnostic functions — no ODE
type appears in their signature. `BiophysicalBehaviour.jl` is the only place that knows about
`OrdinaryDiffEqTsit5`, mirroring the existing `nlp_interface.jl` (physics) / `ipopt.jl`
(solver) split for the steady-state endotherm path. This keeps the derivative reusable later
as a dynamics constraint in an optimal-control formulation without rewriting it.

`onelump` has two methods on the same name (`core_temperature, t, ...` for the ODE derivative;
`t::AbstractVector, core_temperature_init, ...` for the closed-form trajectory under a
constant environment) rather than a `_rhs`-suffixed sibling — the closed form's
`final_core_temperature`/`initial_rate` outputs are the cross-check invariants used in
`test/onelump_derivative.jl`, since no R reference data exists for these yet.

## Deliberate deviations from the R reference

- **`twolump`'s outer surface temperature `Ts` is algebraic, not a third ODE state.** The R
  reference integrates it via a `deSolve` relaxation trick (`dTs = Ts_algebraic - Ts_current`)
  purely so the solver has something to converge — it has no real thermal mass. Here it's
  solved fresh every call from the current shell temperature `Tsk` (linearisation point for
  convection/radiation), giving a 2-state ODE (`Tc`, `Tsk`) instead of 3. `test/transient.jl`
  and the interactive checks in this port confirm the resulting trajectories and steady states
  are sane; no R-generated reference data exists yet to confirm bit-for-bit parity.
- **No solar-constant cap.** `onelump.R` caps direct beam at 1367 W/m² and zeros it above 89°
  zenith; `HeatExchange.solar()` doesn't have this cap and it was deliberately not added, to
  keep the transient port consistent with the rest of the package rather than replicate an R
  quirk. Divergence is only material at extreme (near-horizon) sun angles.
- **Convection correlations** match the R reference exactly for `Cylinder`, `Ellipsoid`,
  `DesertIguana`; `Plate` and `LeopardFrog` use `HeatExchange`'s existing (differently
  calibrated, already-tested) correlations rather than R's. `twolump` only supports
  `Cylinder`/`Ellipsoid` anyway, matching the R source.

## Behavioral driver (`simulate_diurnal_behavior`)

Event-driven port of `trans_behav.R`'s sleep → bask → forage ⇄ cool state machine, using
`OrdinaryDiffEqTsit5.ContinuousCallback` for each transition (root-finding events in the R
reference) instead of R's `deSolve` rootfun/events mechanism. Each phase is its own
`ODEProblem`, solved to its terminating event; state and time carry forward into the next
phase's problem (bout-by-bout, following the same pattern `c:/git/DEBtool_J.jl` uses for
its lifecycle-transition simulations).

Reuses existing types rather than inventing parallel ones: `EctothermBehavioralLimits`'
`active_temperature_min/max` and `basking_temperature_min` fields are the phase thresholds
(R's `T_F_min`/`T_F_max`/`T_B_min`), and `NormalToSun`/`Intermediate` are the postures.

`simulate_diurnal_behavior` dispatches on the initial state, matching R's `lump` parameter:
a plain temperature runs the one-lump model (`HeatExchange.onelump`); a
`(; core_temperature, shell_temperature)` NamedTuple runs two-lump (`HeatExchange.twolump`,
needs a `shell_thickness` keyword). Behavioral thresholds always act on core temperature.

**Simplified relative to `trans_behav.R`:**
- Four phases (`SleepPhase`, `BaskPhase`, `ForagePhase`, `CoolPhase`) cover the core cycle;
  `OrganismState` only distinguishes `Resting`/`Basking`/`Active`, so `Sleep` and `Cool` both
  report as `Resting` — a diagnostic-detail simplification noted in the original plan, not
  revisited here.
- Burrow/depth/height/absorptivity/panting behaviours (`EctothermBehavioralLimits`' other
  fields) aren't wired into this driver — `trans_behav.R` itself doesn't use them either; they
  belong to the separate, existing steady-state `ectothermy.jl` loop.
- Diurnal activity only. `is_active`'s `Nocturnal`/`Crepuscular`/`CombinedActivity` dispatch
  isn't threaded through the phase transitions yet.
- `test/trans_behav_r.jl` compares against a real `trans_behav.R` run (`test/R/trans_behav.R`,
  real `micro_global` forcing saved to `test/data/trans_behav/`) — `max_Tb` within ~1°C of R's.
  `test/transient.jl` additionally checks a synthetic-forcing scenario for sanity.

### Event-detection notes

- Phases are resolved for zero-duration handoffs before integrating (a callback starting
  already on its own root has no sign change to find).
- `ContinuousCallback` uses `affect_neg!=nothing`: only the upward crossing ends a phase.
- Each phase integrates in bounded chunks (`bout_chunk`, default 1 hour) so a chunk that
  starts on a flat span (e.g. zenith pinned at 90° overnight) gets re-checked at the next
  boundary rather than running unchecked.
- `CoolPhase`'s resume threshold follows `trans_behav.R`'s `forage` event: `T_F_min`, or
  shade air temperature + 1K if that's higher (unreachable otherwise in warm shade),
  falling back to `T_F_min` if shade is already near `T_F_max`.
- `SOLR` is unattenuated in both `metout`/`shadmet` — apply canopy shading via
  `EnvironmentalVarsVec`'s `shade` fraction, not the CSV values.
- `_next_phase(::ForagePhase, ...)` compares `_active_max_signal`/`_active_min_signal` to each
  other rather than each to zero independently — a fast body landing near `active_max` from
  root-finding noise could otherwise flip to `BaskPhase`, which immediately bounces back to
  `ForagePhase` with zero elapsed time, stalling the simulation until `max_bouts` is exhausted.

## Not built here (deferred)

- Unifying this with the existing steady-state `ectothermy.jl` loop via a future
  `BodyTemperatureModel` trait (`SteadyState` vs `LumpedCapacitance`).
- A JuMP/IPOPT dynamic-optimization layer reusing `onelump`/`twolump` as collocation
  constraints — the physics/solver split above is what makes this possible later, not
  something built now.
- NicheMapR's custom-shape option (`geom=5`) — see the main plan for why.
