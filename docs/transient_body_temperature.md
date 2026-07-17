# Transient (Lumped-Capacitance) Body Temperature

**Source files:**
- `HeatExchange.jl/src/transient.jl` — `onelump`/`twolump` physics (`dT/dt`), shared with
  steady state via `heat_balance`/`_radiative_convective_flows`
- `HeatExchange.jl/src/internal_temperature.jl` — `internal_gradient_shape_factor`
- `src/transient/forcing.jl` — `EnvironmentForcing` (time-varying environment)
- `src/transient/simulate.jl` — `simulate_onelump`/`simulate_twolump`
- `src/transient/ectotherm/behavioral_driver.jl` — `simulate_diurnal_behavior`
- `src/transient/endotherm/behavioral_driver.jl` — `simulate_endotherm_activity_cycle`

Ported from NicheMapR's `onelump.R`/`onelump_var.R`/`twolump.R`/`trans_behav.R`. Unlike the
rest of the package, which solves steady-state heat balance (`HeatExchange.heat_balance`,
zero thermal mass), this gives the organism real thermal inertia — `dT_b/dt` driven by a
lumped heat-capacitance model, integrated with `OrdinaryDiffEq.jl`.

## Physics dispatch

`onelump` and `twolump` are outer functions that extract `insulation(body(organism))` (or
shape, for `twolump`) and forward to an inner method dispatching on that value — the same
pattern `heat_balance` uses. `onelump`'s `Naked` branch reuses `heat_balance` verbatim,
including evaporative and respiratory heat loss; its insulated branch
(`FibrousLayer`/`CompositeInsulation`) needs an exogenous `metabolic_heat_flow`. `twolump`
(`Cylinder`/`Ellipsoid`-only) reuses
`_radiative_convective_flows` for its surface exchange via a caller-selectable
`SurfaceSolveStrategy`:
- `LinearizedSurface()` (default): linearizes about the current shell temperature, fast,
  has a closed-form `final_core_temperature`.
- `RootFindSurface()`: exact root-find each call, no closed form.

`simulate_diurnal_behavior` and the underlying `_body_temperature_rate` dispatch (one-lump
vs. two-lump) pass `surface_solve` straight through to `twolump`.

`simulate_onelump` handles both `Naked` and insulated bodies (one wrapper); `simulate_twolump`
is `Cylinder`/`Ellipsoid`-only.

## Physics / solver split

`HeatExchange.jl` exposes `onelump`/`twolump` as pure, solver-agnostic functions — no ODE
type appears in their signature. `BiophysicalBehaviour.jl` is the only place that knows about
`OrdinaryDiffEqTsit5`, mirroring the existing `nlp_interface.jl` (physics) / `ipopt.jl`
(solver) split for the steady-state endotherm path. This keeps the derivative reusable later
as a dynamics constraint in an optimal-control formulation without rewriting it.

## Deliberate deviations from the R reference

- **`twolump`'s outer surface temperature `Ts` is algebraic, not a third ODE state.** The R
  reference integrates it via a `deSolve` relaxation trick (`dTs = Ts_algebraic - Ts_current`)
  purely so the solver has something to converge — it has no real thermal mass. Here it's
  solved fresh every call from the current shell temperature `Tsk`, giving a 2-state ODE
  (`Tc`, `Tsk`) instead of 3.
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
`active_temperature_min/max`, `basking_temperature_min`, `emerge_temperature_min`,
`can_climb`, `can_retreat_underground` fields are the phase thresholds/gates (R's
`T_F_min`/`T_F_max`/`T_B_min`/`T_RB_min`/`CLIMB`/`BURROW`), and `NormalToSun`/`Intermediate`
are the postures — the same fields the steady-state `ectothermy.jl` loop already reads.

`simulate_diurnal_behavior` dispatches on the initial state, matching R's `lump` parameter:
a plain temperature runs the one-lump model (`HeatExchange.onelump`); a `(; core_temperature,
shell_temperature)` NamedTuple runs two-lump (`HeatExchange.twolump`, needs a
`shell_thickness` keyword). Behavioral thresholds always act on core temperature.

**Seven phases**: `SleepPhase`/`BurrowPhase` (inactive site), `BaskPhase`, `ForagePhase`,
`CoolPhase`/`ClimbPhase` (too-hot escape site), `RefugePhase` (deeper too-hot escalation).
`BurrowPhase` substitutes for `SleepPhase` whenever `limits.can_retreat_underground` and an
`underground_forcing` are both given; `ClimbPhase` substitutes for `CoolPhase` whenever
`limits.can_climb` and a `climb_forcing` are both given. Both substitutions are unconditional
(decided once from capability + forcing presence, not re-evaluated as a competing choice each
bout) and apply for the whole run. `sleep_forcing` (default `shade_forcing`) sets the
inactive-period site independently — pass `sun_forcing` for an unshaded sleeper. `initial_phase`
overrides the default day/night starting logic to start in any reachable phase. `activity_period`
(default `Diurnal()`) generalizes which part of the day is the active window: `_activity_signal`
is a continuous analogue of `is_active` (`thermoregulation.jl`), reused for `Nocturnal`/
`Crepuscular`/`CombinedActivity` — `Bask`/`Forage` run during the active window, `night_phase`
during the inactive one, regardless of which part of the day that maps to. `ResponsiveActivity`
isn't supported (no continuous signal can be derived from an arbitrary boolean function).

When `has_burrow` (same gate as `BurrowPhase`), `CoolPhase`/`ClimbPhase` can additionally
escalate to `RefugePhase` (same `underground_forcing` site) if core temperature keeps rising
past `limits.escape_temperature_max` despite shade/climbing — the same threshold
`select_depth` uses in the steady-state loop to pick an underground node. `RefugePhase` resumes
straight to `ForagePhase` once core temperature drops back to `active_temperature_max -
active_min_hysteresis`.

Unlike the steady-state loop, which searches many height/depth *nodes* each hour
(`AvailableEnvironments`, `select_depth`, `climb`), `ClimbPhase`/`BurrowPhase` are each locked
to one caller-chosen `EnvironmentForcing` for the whole run — the same simplification already
accepted for `sun_forcing`/`shade_forcing` (one shade fraction, not a continuous search), not a
port of the multi-node search (which needs discrete hourly indexing the continuous ODE driver
doesn't have mid-integration).

**Reused vs. deliberately unused `EctothermBehavioralLimits` fields**: `can_climb`,
`can_retreat_underground`, `emerge_temperature_min` are now read by this driver.
`emerge_signal` (K/hr soil-trend gate), `burrow_shade_mode`, `depth`/`height`/
`depth_min_underground` remain unused — see "Not built here" below.

**Simplified relative to `trans_behav.R`:**
- `OrganismState` only distinguishes `Resting`/`Basking`/`Active`; `Sleep`/`Cool`/`Climb`/
  `Burrow` all report as `Resting` (location doesn't change activity classification, matching
  the steady-state loop) — fine-grained site identity is in the `phase` trajectory field
  instead.
- `depth`/`height`/`absorptivity`/`panting` `SteppedParameter`s aren't wired into this driver —
  `climb_forcing`/`underground_forcing` are each a single fixed site, not a searched node.
- `is_active`'s `Diurnal`/`Nocturnal`/`Crepuscular`/`CombinedActivity` dispatch is reused via
  `activity_period` (see above); `ResponsiveActivity` isn't supported (needs a continuous
  signal, not an arbitrary boolean function). `sunlight`/global radiation isn't part of the
  signal (only zenith), matching the existing "SOLR unattenuated" simplification below.
- `test/trans_behav_r.jl` compares against a real `trans_behav.R` run (`test/R/trans_behav.R`,
  real `micro_global` forcing saved to `test/data/trans_behav/`) — `max_Tb` within ~2°C of R's,
  a loose tolerance that absorbs `onelump`'s evaporative loss and dorsal/ventral
  absorptivity/emissivity split. `test/transient.jl` additionally checks a synthetic-forcing
  scenario for sanity.

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
- `active_min_hysteresis` (default `0.15u"K"`, matching `trans_behav.R`'s `T_F_min ± 0.15°C`)
  separates `BaskPhase`'s exit threshold from `ForagePhase`'s cold-exit threshold, both otherwise
  anchored at the same `active_temperature_min`. Without it, a body settling exactly on that
  value oscillates `BaskPhase`⇄`ForagePhase` forever at zero elapsed time, exhausting
  `max_bouts` — unlike the `ForagePhase` fix above, a relative-signal comparison can't fix this
  since the two thresholds are legitimately the same number.
- `_validate_thresholds` checks `escape_temperature_min < emerge_temperature_min <
  basking_temperature_min < active_temperature_min < active_temperature_max <
  escape_temperature_max` (with the `active_min_hysteresis` margin below
  `active_temperature_min`) once per call. The base ordering (`_validate_ectotherm_thresholds`,
  `src/ectotherm/ectotherm_traits.jl`) is shared with the steady-state loop (`thermoregulate`
  and `example_ectotherm_behavioral_limits` both call it too), so a misconfigured
  `EctothermBehavioralLimits` raises an `ArgumentError` up front instead of silently stalling
  or behaving nonsensically in either driver.
- `BurrowPhase`'s emergence check compares the animal's own core temperature (the ODE state)
  against `emerge_temperature_min`, not the underground forcing's temperature — the transient
  driver never pins core temperature to the environment underground (unlike the steady-state
  loop's `solve_underground=false` shortcut), so only the animal's own state answers "is it
  warm enough to move."
- The day/night boundary (`activity_signal=0`) needs two independent fixes, each protecting
  against a different failure mode:
  - **`activity_hysteresis`** (default `0.1u"°"`) separates `SleepPhase`/`BurrowPhase`'s exit
    (`activity_signal <= -activity_hysteresis`) from `BaskPhase`/`ForagePhase`/`CoolPhase`/
    `ClimbPhase`'s day-arrived check (`activity_signal` compared at the raw boundary) — two
    *different* phases anchored at the same `activity_signal=0`, the same failure class
    `active_min_hysteresis` fixes for temperature. Without it, an animal already warm going
    into the boundary can hand off and bounce straight back with zero elapsed time.
  - **Relative routing** in `_next_phase(::BaskPhase/ForagePhase/CoolPhase/ClimbPhase, ...)`
    compares the day-arrived signal against that *same* phase's own temperature-based exit
    reason, not against zero independently — this is `activity_hysteresis` alone can't fix,
    because it's one phase's own routing ambiguity (which of its several exit reasons fired),
    not two phases sharing a threshold. Needed for `Crepuscular`/`CombinedActivity` in
    particular: their active window is narrow enough (10° for `Crepuscular`) that a root can
    land within floating-point noise of the boundary on almost every crossing, and an
    independent `>=0` check reads that noise inconsistently.

## Endotherm transient (fixed effectors, one-lump)

`onelump`'s insulated branch adds core-temperature thermal mass to the existing steady-state
endotherm physics, without porting new physics from scratch: `HeatExchange.heat_balance`/
`_pack_sides` already compute exactly what a transient RHS needs, just framed as residuals
driven to zero rather than as a rate.

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
  fat/insulation excluded) for both branches, mirroring `twolump`'s core/shell split.
- **`minimum_metabolic_heat`** defaults to `metabolism_pars(organism).metabolic_heat_flow`
  (the organism's configured basal rate), matching `solve_metabolic_rate`'s own convention —
  it's a physiological floor, not zero.
- **Fixed effectors only**: insulation depth, panting, skin wetness, flesh conductivity are
  read from `organism`'s traits and held constant for the whole simulation — no periodic
  re-solving via `RuleBasedSequentialControl`/`IPOPTControl`. Recommended initialization:
  call `solve_metabolic_rate` once, externally, to seed a physiologically self-consistent
  `core_temperature_init` before simulating (see `examples/endotherm_transient.jl`).
- **Validation**: no R reference exists for this (new capability, not a NicheMapR port) — the
  primary check (`HeatExchange.jl/test/onelump.jl`) is that feeding `onelump` the exact
  equilibrium `metabolic_heat_flow` `solve_metabolic_rate` independently finds gives
  `core_temperature_rate ≈ 0`.

`simulate_onelump` returns `core_temperature`/`core_temperature_rate` (both first-class, not
buried in diagnostics) plus `skin_temperature` and per-timestep `diagnostics` (energy/mass
flows), reconstructed post-hoc rather than via a solver callback — same simple approach
`twolump` already uses. There's no automatic critical-temperature stopping: "how long can
this be sustained" is answered by inspecting the returned trajectory (e.g.
`findfirst(>=(critical_temperature), core_temperature)`), as in `examples/endotherm_transient.jl`.

### Activity/rest cycling (`simulate_endotherm_activity_cycle`)

Event-driven analogue of `simulate_diurnal_behavior` for endotherms: alternates an
`EndothermActivePhase` (elevated metabolic rate, e.g. flight or running) and
`EndothermRestingPhase` (basal rate) purely as a function of core temperature, using
`ContinuousCallback` the same way. No day/night cycle — this is an activity/rest thermal
cycle, not a full diel budget (deliberately out of scope; see below). Insulated bodies only,
since `onelump`'s `metabolic_heat_flow` override only exists on that branch.

Two independent `EnvironmentForcing` streams (`active_forcing`/`resting_forcing`) — same
two-forcing pattern as the ectotherm driver's `sun_forcing`/`shade_forcing` — is what gives
this "height choice" for free: build `active_forcing` from reference-height microclimate
columns (e.g. NicheMapR's `TAREF`/`VREF`) for a flying animal and `resting_forcing` from
local-height columns (`TALOC`/`VLOC`) for a perched/walking one, with no change to the driver
itself (see `examples/endotherm_transient.jl`, which does exactly this from real
`test/R/trans_behav.R` output). The example also applies effective wind speed as
`max(movement_speed, ambient_wind_at_that_height)`, not movement speed alone — NicheMapR's own
(untested) butterfly/flight model in `ectotherm.R` has the same two levers under different
names (`flymetab`, `flyspeed`).

Because there's no third phase to disambiguate (unlike the ectotherm driver's `ForagePhase`,
which had three possible exits), `_next_phase` here is unconditional alternation — the
floating-point boundary-noise bug fixed in the ectotherm driver doesn't have an analogue to
fix.

## Not built here (deferred)

- Unifying this with the existing steady-state `ectothermy.jl` loop via a future
  `BodyTemperatureModel` trait (`SteadyState` vs `LumpedCapacitance`).
- Adopting Microclimate.jl's `CommonSolve.solve`/`init`/`solve!` pattern for `onelump`/
  `twolump` — would need the same redesign applied to `heat_balance`/`solve_temperature`/
  `solve_metabolic_rate` to be coherent; tracked as a future whole-package follow-on, not a
  partial application here.
- A JuMP/IPOPT dynamic-optimization layer reusing `onelump`/`twolump` as collocation
  constraints — the physics/solver split above is what makes this possible later, not
  something built now. Forward-mode AD (`ForwardDiff.jl`) through a single `simulate_onelump`
  call should already work today (no events, `SmoothBound` smoothing) as a cheap way to
  validate the chain before investing in this.
- **Optimal movement height/speed** — `simulate_endotherm_activity_cycle` takes height and
  movement speed as fixed inputs (baked into `active_forcing`); finding the height/speed that
  optimizes some objective (e.g. maximise distance covered without exceeding a thermal limit)
  is a natural follow-on once the AD/optimal-control layer above exists.
- Periodic effector re-solving mid-simulation, torpor/heterothermy phases, and a full diurnal
  cycle (Sleep at night) for the activity/rest driver — deliberately out of scope for the same
  "keep effectors fixed, no day/night" reasons as above.
- An ectotherm analogue of the activity/rest driver (`onelump` + a flight/movement metabolic
  term + movement-speed wind exposure, mirroring `ectotherm.R`'s `flyer`/`flymetab`/`flyspeed`/
  `flyhigh`) — the trigger logic in NicheMapR lives in compiled Fortran, not inspectable from
  the R wrapper, so this needs its own design pass rather than porting undocumented behaviour.
- `emerge_signal` (K/hr soil-warming/cooling trend gate on emergence) — a discrete
  finite-difference construct in the R source that doesn't map cleanly onto a
  `ContinuousCallback` root; most configs use the disabled default anyway.
- `burrow_shade_mode`-driven blended underground forcing, and graded/multi-node
  `ClimbPhase`/`BurrowPhase` site search (matching the steady-state loop's `select_depth`/
  `climb` over many height/depth nodes) — `climb_forcing`/`underground_forcing` are each a
  single fixed site for now.
