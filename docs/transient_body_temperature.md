# Transient (Lumped-Capacitance) Body Temperature

**Source files:**
- `HeatExchange.jl/src/transient_lumped.jl` — `onelump`/`twolump` physics (`dT/dt`)
- `HeatExchange.jl/src/internal_temperature.jl` — `internal_gradient_shape_factor`
- `HeatExchange.jl/src/transient_lumped_endotherm.jl` — `endotherm_onelump` physics (`dT_core/dt`)
- `src/transient/forcing.jl` — `EnvironmentForcing` (time-varying environment)
- `src/transient/simulate.jl` — `simulate_onelump`/`simulate_twolump`
- `src/transient/ectotherm/behavioral_driver.jl` — `simulate_diurnal_behavior`
- `src/transient/endotherm/simulate.jl` — `simulate_endotherm_onelump`

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

## Endotherm transient (fixed effectors, one-lump)

`endotherm_onelump` adds core-temperature thermal mass to the existing steady-state endotherm
physics, without porting new physics from scratch: `HeatExchange.heat_balance`/`_pack_sides`
already compute exactly what a transient RHS needs, just framed as residuals driven to zero
rather than as a rate.

- **`Naked` bodies** reuse `heat_balance(core_temperature, organism, e)` verbatim — its
  `metabolic_heat_flow` already comes from `organism`'s own `metabolism_pars.model` (a genuine
  forward model, e.g. `Kleiber`/`McKechnieWolf`/`AndrewsPough2`), so `.energy_balance.heat_balance`
  (heat in − heat out at the given, possibly off-equilibrium, core temperature) is already the
  correct ODE numerator.
- **Insulated bodies** (`FibrousLayer`/`CompositeInsulation`) need `metabolic_heat_flow` as an
  exogenous keyword (`Quantity` or `Function(core_temperature)`) — `solve_metabolic_rate`
  zbrent-solves it as an unknown at equilibrium; calling that inside the ODE loop would force
  instantaneous equilibrium every step, defeating the point of adding thermal mass. Skin/
  insulation temperature are still solved algebraically each call via `_pack_sides` (not a
  second ODE state, the same treatment `twolump`'s outer surface temperature gets). The
  numerator is `metabolic_heat_flow − respiration_heat_flow − net_metabolic_heat_internal`
  (`HeatExchange.heat_balance`'s `residual_internal_conduction`, matching `respiration()`'s own
  `.balance` residual — `MetabolicRates.sum` must be `net_metabolic_heat_internal`, not
  `metabolic_heat_flow`, to match this).
- **Thermal capacitance**: `flesh_volume(body) * density * flesh_specific_heat` (flesh only,
  fat/insulation excluded) — mirrors `twolump`'s core/shell split, not `onelump`'s whole-body
  mass (which never had a fat layer to exclude).
- **`minimum_metabolic_heat`** defaults to `metabolism_pars(organism).metabolic_heat_flow`
  (the organism's configured basal rate), matching `solve_metabolic_rate`'s own convention —
  it's a physiological floor, not zero.
- **Fixed effectors only**: insulation depth, panting, skin wetness, flesh conductivity are
  read from `organism`'s traits and held constant for the whole simulation — no periodic
  re-solving via `RuleBasedSequentialControl`/`IPOPTControl`. Recommended initialization:
  call `solve_metabolic_rate` once, externally, to seed a physiologically self-consistent
  `core_temperature_init` before simulating (see `examples/endotherm_transient.jl`).
- **Validation**: no R reference exists for this (new capability, not a NicheMapR port) — the
  primary check (`test/endotherm_onelump.jl`) is that feeding `endotherm_onelump` the exact
  equilibrium `metabolic_heat_flow` `solve_metabolic_rate` independently finds gives
  `core_temperature_rate ≈ 0`.

`simulate_endotherm_onelump` mirrors `simulate_onelump`'s structure and returns
`core_temperature`/`core_temperature_rate` (both first-class, not buried in diagnostics) plus
`skin_temperature` and per-timestep `diagnostics` (energy/mass flows), reconstructed post-hoc
rather than via a solver callback — same simple approach `twolump` already uses. There's no
automatic critical-temperature stopping: "how long can this be sustained" is answered by
inspecting the returned trajectory (e.g. `findfirst(>=(critical_temperature),
core_temperature)`), as in `examples/endotherm_transient.jl`.

## Not built here (deferred)

- Unifying this with the existing steady-state `ectothermy.jl` loop via a future
  `BodyTemperatureModel` trait (`SteadyState` vs `LumpedCapacitance`).
- A JuMP/IPOPT dynamic-optimization layer reusing `onelump`/`twolump` as collocation
  constraints — the physics/solver split above is what makes this possible later, not
  something built now.
- NicheMapR's custom-shape option (`geom=5`) — see the main plan for why.
- Periodic effector re-solving mid-simulation for endotherms, `endotherm_twolump`,
  torpor/heterothermy phases, and an endotherm analogue of `simulate_diurnal_behavior`'s
  event-driven effector switching.
