# BiophysicalBehaviour — Multi-part refactor plan

Working design + ordered execution roadmap for updating `BiophysicalBehaviour` to
handle arbitrary multi-part organism geometries defined in `BiophysicalGeometry`,
alongside coordinated refactors to `HeatExchange` and `BiophysicalGeometry`.

We control all four packages (`BiophysicalGeometry`, `HeatExchange`,
`BiophysicalBehaviour`, and, tangentially, `DEBTOOL_J` for state-management
patterns and `ModelParameters` for parameter tooling). Nothing is registered.
Backwards compatibility is not a constraint. Simplicity, type stability, and
zero-allocation hot paths are.

## Status (updated 2026-08-05)

- **Phase 0-3** — DONE (family hierarchies, join accessors, part selectors,
  effector refactor).
- **Phase 4** — DONE. Per-part physiology + `LungPart` + `lung_part` index;
  co-located per-part **limits container** deferred to 7.5 (its only consumer is
  the NLP builder); `ThermoregulationLimits` gained the four magic-number fields
  + `weight_for`; per-part effectors (`map_part_physiology`) with panting/core
  routed to the lung part; multi-part `RuleBasedSequentialControl` loop driving
  the coupled solve. `ThermoregulationLimits` *slimming* (removing
  insulation/flesh_conductivity/etc.) stays with the Phase 8 deletion pass since
  the single-body loop + 328-case regression still consume them.
- **Phase 6-7** — Core physics landed in `HeatExchange` (`solve_part_heat_balance`,
  `HeatCoupling`, `CompartmentGraph`, `solve_part_surface`,
  `solve_coupled_metabolic_rate`, shell-fraction fix). `CommonSolve`
  `HeatBalanceProblem`/`HeatBalanceSolver` interface (§3.6) and `Vasodilate`
  conductance-matrix caching not yet wired. **BB-side integration DONE**:
  `solve_multipart_metabolic_rate` / `part_surface_setups` build per-part setups
  from a multi-part `Organism` and route through `lung_part`. Floating
  compartments (independent leg/head cores via `solve_core_temperatures`) not
  yet wired — current solve assumes one all-`SharedCore` compartment.
- **Phase 5** — DONE (Tier-1). `ShapeCache` per-part geometry cache
  (`total_area`, `silhouette_area`, `characteristic_dim`), built once by
  `precompute_shape_cache` and invalidated by direct effector dispatch
  (`refresh` — no-op for physiological effectors, rebuild for `Piloerect` /
  `Uncurl`). Threaded through `solve_multipart_metabolic_rate` /
  `part_surface_setups` (cached == uncached bitwise). Tier-2 radiation cache and
  the formal JET/zero-allocation gate on the whole `thermoregulate` call remain
  a later optimisation pass (the outer `physiology` broadcast still allocates).
- **Phase 9** — Example (`examples/dog_thermoregulation.jl`) + multi-part tests
  (`test/multipart_solve.jl`) landed for the single-compartment path.
- **Phase 7.5** (multi-part NLP + `IPOPTControl`), **Phase 8** (dorsal/ventral
  deletion), **Phase 10** (cleanup) — REMAINING. Multi-part `Piloerect` /
  `Uncurl` (which reshape geometry and drive the Phase 5 cache rebuild) land with
  the per-part insulation-limits container in 7.5. Phase 8 is gated on migrating
  the 328-case endotherm regression off the dorsal/ventral reference model;
  Phase 7.5 is a new Flatten-templated NLP subsystem across both packages plus
  the Ipopt extension.

## Naming conventions used throughout

**Fully spelled identifiers.** No `T_` / `Q_` / `k_` prefixes, no `Ra` /
`Pr` / `Nu` for dimensionless numbers, no `evap` / `resp` / `pars` / `env`
/ `geom` abbreviations — e.g. `T_skin` → `skin_temperature`, `Q_gen` →
`metabolic_heat`, `k_flesh` → `flesh_conductivity`, `Nu` →
`nusselt_number`, `IR` → `infrared`. Existing abbreviated identifiers are
renamed by the phase that touches them. When this document quotes
existing code (current field names, NicheMapR reference material), the
original identifiers are kept as-is for accuracy.

---

## 1. Goal

Support organisms composed of an arbitrary number of anatomical parts — head,
torso, legs, etc. — each with its own shape, insulation, physiology, and
thermal state. Parts are joined at explicit interfaces (already modelled by
`BiophysicalGeometry.CompositeBody` + `Join`). Heat exchanges through those
interfaces according to per-join coupling policy. External heat exchange
(radiation, convection, evaporation) happens per part.

The existing single-body, dorsal/ventral-split model becomes a special case:
one part with `Naked`, or two `HalfCylinder` parts with a `SharedCore` join.

---

## 2. Current architecture (diagnosis)

Three layers assume a single lumped body:

- **Effectors** (`src/endotherm/thermoregulation.jl`) rebuild
  `organism.body = Body(...)` — hard-fails on `CompositeBody`.
- **State** (`src/endotherm/homeothermy.jl:98-229`) threads skin /
  insulation / core temperatures and metabolic heat as scalars;
  `initialise_state` is a `NullBehavior` stub (`src/organism.jl:126`).
- **Physiology traits** (`src/organism.jl:243-257`) — a single lumped
  `HeatExchangeTraits` with `insulation_pars.dorsal/ventral`; no
  per-part concept.

`HeatExchange` is tightly coupled to `BiophysicalGeometry`: direct scalar
accessor calls (`total_area`, `silhouette_area`, `evaporation_area`,
`flesh_volume`, `characteristic_dimension`) sprinkled through
`biophysics.jl` / `heat_balance.jl` / `nlp_interface.jl` /
`insulated/utils.jl`; shape dispatch for closed-form conduction and
Nusselt correlations; and a deep dorsal/ventral split (`_pack_sides` in
`insulated/utils.jl:60-230`, `MultiSided`, `Dorsal`/`Ventral` singletons,
dorsal/ventral fields on `Radiation`/`Insulation`/`Absorptivity`
parameters, `ventral_fraction`, `mean_skin_temperature.jl`) —
roughly half the `insulated/` directory exists to service the split.

---

## 3. Target architecture

### 3.1 Shared shape and insulation hierarchies (in `BiophysicalGeometry`)

Introduce physics-relevant abstract intermediates between `AbstractShape` and
concrete shapes. Same for insulation. Single hierarchy, zero mapping tables,
every consumer dispatches on the same types.

```julia
# BiophysicalGeometry

abstract type AbstractShape end
abstract type AbstractCylindrical <: AbstractShape end
abstract type AbstractSpherical   <: AbstractShape end
abstract type AbstractEllipsoidal <: AbstractShape end
abstract type AbstractSlab        <: AbstractShape end

struct Cylinder      <: AbstractCylindrical ... end
struct HalfCylinder  <: AbstractCylindrical ... end
struct Cone          <: AbstractCylindrical ... end
struct Sphere        <: AbstractSpherical   ... end
struct Ellipsoid     <: AbstractEllipsoidal ... end
struct HalfEllipsoid <: AbstractEllipsoidal ... end
struct Plate         <: AbstractSlab        ... end
struct TriMesh       <: AbstractShape       ... end    # no thermal dispatch (yet)

abstract type AbstractInsulation end
abstract type AbstractFibrous <: AbstractInsulation end    # fur, feathers, hair
abstract type AbstractLayer   <: AbstractInsulation end    # fat, muscle, skin

struct Naked                                          <: AbstractInsulation ... end
struct Fur                                            <: AbstractFibrous    ... end
struct Fat                                            <: AbstractLayer      ... end
struct CompositeInsulation{Layers<:Tuple}             <: AbstractInsulation ... end
```

`TriMesh` deliberately has no family membership — attempting a thermal solve on
a mesh part will error with a missing-method dispatch. Silent fallback to the
wrong physics correlation is worse than an explicit signal.

Free wins from the family types: `Feathers <: AbstractFibrous`,
`Muscle <: AbstractLayer`, `SkinLayer <: AbstractLayer` become drop-in structs
with no new methods. `HalfEllipsoid` inherits ellipsoid thermal physics.

### 3.2 Multi-part composition and joins

`BiophysicalGeometry.CompositeBody{Root,Parts<:NamedTuple,Joins<:Tuple,RootPose<:NamedTuple}`
already gives us the right primitive: parts as a `NamedTuple` keyed by
anatomical name, joins as a `Tuple` with parent/child part names lifted into
type parameters, poses computed at construction. Adopt as-is.

Small additions needed in `BiophysicalGeometry`:

- `join_area(join, body) -> Area` — `π · radius^2` for `Disc(radius)`,
  delegates to `surface_area` for `FullCover`.
- `join_position(join, body) -> NTuple{3,Length}` — world-frame center of the
  patch, from existing `Pose`.
- `join_partners(join) -> (Symbol, Symbol)` — expose part names currently only
  in type parameters.
- `internal_distance(shape, surface_location) -> Length` — centroid-to-surface
  path length for lumped-resistance calculations. Closed-form for
  cylinder/sphere/ellipsoid; a bit of work for half-shapes; skip mesh.

### 3.3 Heat coupling per join

The `HeatCoupling` hierarchy describes *how heat exchanges between parts* —
that is physics, not geometry and not behavior. It lives in **`HeatExchange`**.
`BiophysicalGeometry` stays pure geometry (usable by aerodynamics /
mechanical consumers with no thermal opinion). `BiophysicalBehaviour` merely
*stores* a coupling per join in `OrganismTraits.couplings::NamedTuple`,
parallel to `body.joins` and keyed the same way, and hands it to
`HeatExchange` at problem-construction time.

```julia
abstract type HeatCoupling end
struct SharedCore <: HeatCoupling end          # parts equilibrate — contract into one node

struct ConductiveCoupling{InterfaceConductivity} <: HeatCoupling
    interface_conductivity::InterfaceConductivity
end
```

Joins without an entry in `couplings` are insulated by default (zero
contribution to conductance matrix and heat load) — no explicit
`InsulatedJoin` sentinel is needed.

Example dog:
```julia
couplings = (
    dorsal_ventral = SharedCore(),          # torso halves, one blood pool
    dorsal_head    = ConductiveCoupling(),  # derive from flesh_conductivity + geometry
    ventral_leg_fl = ConductiveCoupling(),
    ventral_leg_fr = ConductiveCoupling(),
    ventral_leg_bl = ConductiveCoupling(),
    ventral_leg_br = ConductiveCoupling(),
)
```

The two derivation modes are distinguished at the type-parameter level,
not by a runtime `nothing` check. `ConductiveCoupling()` constructs
`ConductiveCoupling{Nothing}`; `contribution_to_conductance_matrix(::ConductiveCoupling{Nothing}, ...)`
derives series resistance
`total_resistance = distance_parent / (conductivity_parent * area) + distance_child / (conductivity_child * area)`
from each part's own `flesh_conductivity` and `internal_distance` to the
join. `ConductiveCoupling(k)` constructs `ConductiveCoupling{typeof(k)}`;
the `<:Real` method uses the explicit override. No branch in the hot loop.

### 3.4 Compartment solve

The organism is a **graph of thermal nodes**. Today each part contributes one
node (its core); once `LayeredShape` (section 8.1) lands, a layered part
contributes one node per shell. Design the graph in terms of nodes, not parts,
from the start.

`CompartmentGraph` and the coupling-contribution methods live in
**`HeatExchange`** (they encode physics topology, not behavior). Union-find
over `SharedCore` edges partitions nodes into thermal compartments; each
compartment has one temperature. Topology is compile-time known (part names,
coupling types, layer counts all in type parameters), so build once at
`HeatBalanceProblem` construction (section 3.6) and freeze into a
`CompartmentGraph{Groups,EdgeIndices,NumCompartments}`.

Couplings contribute to the solve via a two-method interface — new physics
(e.g., `BloodPerfusionCoupling` from section 8.4) is added as new methods, not
by rewriting the outer loop:

```julia
contribution_to_conductance_matrix(::SharedCore, ...)                    # zero; edge already contracted
contribution_to_conductance_matrix(::ConductiveCoupling{Nothing}, ...)   # derived series resistance
contribution_to_conductance_matrix(::ConductiveCoupling{<:Real}, ...)    # explicit override
contribution_to_heat_load(::HeatCoupling, ...)                           # right-hand-side shifts
                                                                         # (nonzero for perfusion)
```

Per outer-loop iteration:
1. Given current per-part `skin_temperature` and `insulation_temperature`,
   compute per-part `environmental_heat` by calling `HeatExchange` once per
   part.
2. Build the compartment heat-load vector:
   `heat_load_vector[compartment] = sum over parts in compartment of
    (metabolic_heat[part] − environmental_heat[part])`.
3. Solve the linear system
   `conductance_matrix * core_temperatures = heat_load_vector`. The
   `conductance_matrix` contains contributions from `ConductiveCoupling` joins
   only. An `SMatrix{NumCompartments,NumCompartments}` solve — stack-allocated,
   no heap allocation.
4. Propagate compartment `core_temperatures` back through flesh to update
   per-part `skin_temperature` for the next iteration.

The `conductance_matrix` is a per-organism constant, recomputed only when
`Vasodilate` fires.

Degenerate cases fall out for free:
- 1-part `Body` → 1 compartment → 1-by-1 solve → same numerics as today.
- All-`SharedCore` composite → 1 compartment → external solve per part,
  single core.

### 3.5 Two-tier cache on the organism

Caches live on the organism (mutable-style, via `@set`). Fixed-shape
`NamedTuple`s keyed by part name, all `Float64` (or `Unitful`) — type-stable
and stack-allocatable.

**Tier 1 — Shape cache:**
`(; skin_area, insulation_area, evaporation_area, flesh_volume, characteristic_dimension)`
per part. Function of shape + insulation only. Invalidated by `Piloerect` (fur
depth changes insulation surface) and `Uncurl`.

**Tier 2 — Radiation cache:**
`(; area_direct, area_diffuse_up, area_diffuse_down, view_factor_ground, view_factor_sky)`
per part. Function of shape + part pose + sun direction + sky geometry.
Invalidated by shape changes and by environment ticks (sun moved). Not
invalidated by `Piloerect`, `Vasodilate`, `Hyperthermia`, `Pant`, `Sweat` —
the vast majority of thermoregulation-loop iterations.

Posture-change behaviors (`Roll`, `Curl`) are out of scope for this refactor.

**Cache invalidation via direct effector dispatch.** The effector type
is concrete at the call site, so `refresh!` dispatches on it directly —
no `Val` wrapper, no intermediate `invalidates` trait function:

```julia
@inline refresh!(cache, body, environment, ::Effector)    = cache   # default no-op
@inline refresh!(cache, body, environment, ::Piloerect)   = ...     # refresh shape
@inline refresh!(cache, body, environment, ::Uncurl)      = ...     # refresh shape + radiation
@inline refresh!(cache, body, environment, ::Vasodilate)  = ...     # refresh conductance matrix
```

`Vasodilate`, `Hyperthermia`, `Pant`, `Sweat` fall through to the
`::Effector` no-op and the compiler inlines them away at the call site.
Only `Piloerect`, `Uncurl`, and `Vasodilate` do real work.

### 3.6 HeatExchange as pure physics + `CommonSolve` interface

`HeatExchange` loses all geometry function calls and all dorsal/ventral
machinery, and exposes the standard **`CommonSolve.jl`** interface: `init`,
`solve!`, `reinit!`, `solve`. This matches how the rest of the SciML
ecosystem (`DifferentialEquations`, `LinearSolve`, `NonlinearSolve`,
`Optimization`) works, and gives us two usage modes with no code duplication.

**Problem type** — a self-contained multi-part heat-balance problem
(topology, physics, environment only — no solver state):

```julia
struct HeatBalanceProblem{Geometry, Physiology, Environment, Couplings, CompartmentGraph}
    geometry::Geometry                    # per-part shape + insulation + Tier 1/2 cache scalars
    physiology::Physiology                # per-part physiological parameters
    environment::Environment              # air / ground / sky temperatures, wind,
                                          # solar_direct, solar_diffuse, ...
    couplings::Couplings                  # NamedTuple keyed by join name, values
                                          # <: HeatCoupling
    compartment_graph::CompartmentGraph   # precomputed from couplings (section 3.4)
end
```

The problem holds *everything* required to define the physics — no hidden
state, no callbacks into `BiophysicalGeometry`. Initial state flows into
the solver via `init` / `reinit!`, not the problem (SciML convention).

**Solver type** — mutable workspace + current state, reused across successive
solves. Per-solve output metadata (`converged`, `iterations`, `retcode`,
`stats`) belongs on the returned `Solution`, not the solver:

```julia
mutable struct HeatBalanceSolver{Problem, Workspace, State}
    problem::Problem
    workspace::Workspace   # Jacobian / root-finding buffers / conductance matrix cache
    state::State           # current best estimate of skin, insulation, core temperatures
end
```

**Interface** (implemented against `CommonSolve.jl`):

```julia
CommonSolve.init(problem::HeatBalanceProblem; initial_state, kwargs...) -> HeatBalanceSolver
CommonSolve.solve!(solver::HeatBalanceSolver) -> Solution   # run to convergence in place
CommonSolve.reinit!(solver::HeatBalanceSolver, new_problem::HeatBalanceProblem;
                    warm_start::Bool = true)                # reset for a new problem,
                                                            # optionally reusing state
CommonSolve.solve(problem::HeatBalanceProblem; initial_state, kwargs...) =
    solve!(init(problem; initial_state, kwargs...))         # one-shot convenience
```

**Two usage modes fall out:**

1. **One-shot from downstream code** (e.g., someone modelling a static
   organism / environment pair, or a NicheMapR-style time-stepping driver):
   ```julia
   solution = solve(HeatBalanceProblem(geometry, physiology, environment,
                                       couplings); initial_state)
   ```
   Full multi-part solve, no thermoregulation control loop, no
   `BiophysicalBehaviour` needed.

2. **Stepwise from `BiophysicalBehaviour`'s control loop**:
   ```julia
   solver = init(HeatBalanceProblem(from_organism(organism), ...))
   solve!(solver)                                   # first solve
   while solver.state.metabolic_heat < minimum_metabolic_heat
       organism = apply_effector(effector, organism)          # immutable
       reinit!(solver, HeatBalanceProblem(from_organism(organism), ...))
       solve!(solver)                               # warm-started from previous state
   end
   ```
   Workspace and state persist across effector iterations; warm-starting is a
   significant performance win because small effector changes → small state
   changes.

**Where the per-part physics primitive lives:** the family-dispatched
functions below (`nusselt_free`, `surface_temperature`, `insulation_conductance`)
are the atoms used inside `solve!`. They stay scalar-in, scalar-out — the
solver calls them once per part per iteration and threads results through the
`CompartmentGraph` (section 3.4). External consumers wanting a *single-part*
solve can still call them directly, but the `CommonSolve` API is the primary
surface.

**Signatures** — the per-part primitive:

```julia
solve_part_heat_balance(
    physiology,        # insulation parameters, flesh_conductivity,
                       # metabolism parameters, options
    environment,       # air_temperature, ground_temperature, sky_temperature,
                       # wind, solar_direct, solar_diffuse, ...
    geometry,          # (; shape, insulation,
                       #    skin_area, insulation_area, evaporation_area,
                       #    contact_area, flesh_volume, characteristic_dimension,
                       #    area_direct, area_diffuse_up, area_diffuse_down)
    skin_temperature, insulation_temperature, core_temperature_boundary,
) -> (; metabolic_heat, skin_temperature, insulation_temperature,
        convective_heat, radiative_heat, evaporative_heat, respiratory_heat, ...)
```

Called **once per `Body` part per solver iteration**. Dispatches on
`geometry.shape::AbstractCylindrical` etc. — reads shape fields directly
(`shape.radius`, `shape.semi_axes`, ...), no function calls into
`BiophysicalGeometry`. `HeatExchange` still `using BiophysicalGeometry` for
the abstract types, but never calls heavy geometric functions in hot paths.
(This function replaces today's `solve_metabolic_rate`, which fused per-part
physics with dorsal/ventral bookkeeping; renaming makes clear that it is
strictly a *single part* primitive.)

Family dispatch in `HeatExchange` — `nusselt_free`,
`surface_temperature`, `insulation_conductance` each get one method per
family abstract (`AbstractCylindrical`, `AbstractSpherical`,
`AbstractEllipsoidal`, `AbstractSlab`; `Naked`, `AbstractFibrous`,
`AbstractLayer`, `CompositeInsulation`). Adding a new concrete shape or
insulation is a one-liner in `BiophysicalGeometry`; physics works
automatically via `<:` dispatch.

### 3.7 NLP interface + strategy dispatch

**Split across two packages.** `HeatExchange/src/nlp_interface.jl` is the
**physics side** — three pure functions (`nlp_pack`, `nlp_residuals`,
`nlp_assemble_output`) that any nonlinear solver can drive, with two
strategies today (`WeightedMeanNLP` — 3 residuals, mean-weighted;
`MultiSidedNLP` — 5 residuals, per-side).
`BiophysicalBehaviour/src/endotherm/thermoregulation/ipopt.jl` (on the
`update-IPOPT-implementation-to-new-HeatExchange-names-and-NLP-hook`
branch) is the **policy side** — objective, bounds, initial values, Q10
metabolic-scaling inequality constraint. The tech stack is
**`Optimization.jl` + `OptimizationIpopt.jl`** (SciML): the driver builds an
`OptimizationFunction` and `OptimizationProblem`, then calls
`SciMLBase.solve` with `IpoptOptimizer(; hessian_approximation =
"limited-memory", ...)`. Optimization.jl already implements the
`CommonSolve` interface, so it composes with the design in section 3.6
directly — no separate JuMP model to persist.

The intellectual content: effectors (piloerection, uncurl, vasodilate,
panting, sweating) become continuous **decision variables** with box
bounds from `ThermoregulationLimits`, plus per-effector quadratic penalties
in the objective. Ipopt walks the effector space in one solve — an
alternative to the rule-based sequential ladder in
`RuleBasedSequentialControl`.

**Numerical invariants worth preserving from the current branch:**

- **Hard `Float64` ↔ `Unitful` boundary at `_heat_balance_residuals_*!`.**
  Solver world is unitless (`Float64`); physics world is `Unitful`. Units
  are attached at the top of the residual function and stripped at the
  bottom. Zero `Unitful` inside the optimizer's hot loop. Preserve this
  layer under multi-part; extend the strip/attach machinery to unpack the
  flat decision vector into per-part / per-compartment `NamedTuple`s of
  `Unitful` quantities.
- **Log-transform on `metabolic_heat_flow`** (`decision[k] =
  log(metabolic_heat_flow)`). Keeps the variable in a numerically stable
  range and enforces positivity implicitly. Preserve.
- **Objective is a weighted sum of quadratic penalties on** core-temperature
  deviation from setpoint, metabolic-heat magnitude, core-skin gradient,
  panting rate, skin wetness. Each normalised by a `_range` to keep terms
  O(1). Extend to per-part / per-compartment where semantically per-part
  (e.g. penalise `insulation_depth` deviation from reference per part).
- **Q10 inequality constraint**
  `metabolic_heat_flow ≥ minimum_heat_flow · q10 ^ ((core_temperature −
  setpoint_temperature) / 10)`. This is the physically-required lower bound
  on metabolic heat production; without it the optimiser drives metabolic
  heat to zero. Stays as a single whole-organism inequality (metabolism is
  routed to `lung_part`; the constraint is on total metabolic heat, not
  per-part).

**Under multi-part composition, `WeightedMeanNLP` and `MultiSidedNLP`
collapse into a single strategy.** Dorsal/ventral become just two parts; the
mean-weighting hack disappears. Residuals and decision variables become
natural functions of the compile-time topology (`NumParts`,
`NumCompartments`):

- **Residuals:**
  `2 · NumParts` (per-part energy balance + skin-temperature closure)
  `+ NumCompartments` (per-compartment core balance)
  `+ 1` (whole-organism metabolic − respiration = net internal flow).
- **Decision variables:**
  `2 · NumParts` (per-part skin & insulation temperatures)
  `+ 2 · NumCompartments` (per-compartment `core_temperature` and
  `flesh_conductivity`)
  `+ 3 · NumParts` (per-part `insulation_depth`, `skin_wetness`)
  `+ 2` (`metabolic_heat_flow` and `pant`, whole-organism, routed to
  `lung_part`).

For a 6-part dog with `SharedCore` on torso halves (5 compartments): **18
residuals, 42 decision variables**. Well-posed underdetermined system with
an objective to minimize (e.g., total metabolic effort or a weighted sum of
effector "costs").

**Variable-layout mechanism: `NamedTuple` template + `Flatten.jl`.**
The layout is *not* a hand-maintained compile-time table of names, sizes,
and index ranges — it's a nested `NamedTuple` template whose shape is
derived from the `HeatBalanceProblem`'s topology (`body.parts`, `couplings`,
`compartment_graph`) and whose leaves carry the `Unitful` value type of each
decision variable. `Flatten.jl` is the pack / unpack engine, exactly as
`DEBTOOL_J`'s `StateReconstructor`
(`DEBTOOL_J/src/animals/ode/utils.jl:39`) uses it to bridge nested state and
DifferentialEquations' flat `u` vector.

```julia
nlp_template(problem) =
    (; parts = (; head    = (; skin_temperature       = 300.0u"K",
                                insulation_temperature = 295.0u"K",
                                insulation_depth       = 0.01u"m",
                                skin_wetness           = 0.1),
                  torso_dorsal   = (; ...), torso_ventral = (; ...),
                  leg_fl = (; ...), leg_fr = (; ...), ...),
       compartments = (; torso_core = (; core_temperature   = 310.0u"K",
                                          flesh_conductivity = 0.5u"W/m/K"),
                         head_core  = (; ...), leg_fl_core = (; ...), ...),
       whole_organism = (; log_metabolic_heat_flow = log(0.5),
                           pant                    = 1.0))
```

Pack/unpack is three `Flatten` calls mirroring `DEBTOOL_J`'s `to_vec` /
`to_obj`: `unpack(template, x) = Flatten.reconstruct(template, x,
Number)` (Float64 → Unitful, per-field convert via
`typeof(template_leaf)`); `pack(obj) = SVector(map(ustrip,
Flatten.flatten(obj, Number)))` (Unitful → Float64); `initial(template)
= pack(template)`.

**The NLP API in `HeatExchange`** (no JuMP / optimizer / `Optimization.jl`
dep; one new dep: `Flatten`):

```julia
nlp_template(problem::HeatBalanceProblem)  -> NamedTuple
nlp_pack(problem::HeatBalanceProblem)      -> PackedNLP{Problem, Template, Workspace}
nlp_residuals!(residuals, packed, x)       -> residuals
nlp_variable_bounds(problem, limits)       -> (lower_template, upper_template)
nlp_assemble_output(packed, x)             -> Solution
```

`nlp_residuals!` is the only hot function. It calls
`Flatten.reconstruct` once at the top (Float64 → Unitful), runs
`solve_part_heat_balance` per part, and `ustrip`s each residual as it's
written. Same physics primitives (`solve_part_heat_balance`,
`contribution_to_conductance_matrix`, `contribution_to_heat_load`) as
the iterative solver — **zero physics duplication** between paths.

**Strategy dispatch — layered on the SciML stack.** Because
Optimization.jl already implements `CommonSolve`, an NLP
thermoregulation solve is *already* a `solve(::OptimizationProblem)`
call. `BiophysicalBehaviour` owns the driver:

```julia
struct IPOPTControl{Objective, Options} <: AbstractControlStrategy
    objective::Objective   # penalty weights, targets, ranges
    options::Options       # Ipopt options (tol, max_iter, hessian_approximation, ...)
end
```

`thermoregulate(::IPOPTControl, organism, environment, ...)` builds a
`HeatBalanceProblem`, calls `HeatExchange.nlp_pack` /
`nlp_variable_bounds`, flattens the templates to `Vector{Float64}` (same
`pack` path used everywhere), wraps into an `OptimizationFunction` +
`OptimizationProblem`, and calls `SciMLBase.solve(_, IpoptOptimizer(;
control.options...))`. `RuleBasedSequentialControl` continues to drive
the `HeatBalanceSolver{IterativeStrategy}` from 3.6.

**Package layout:** `HeatExchange` has no `Optimization.jl` /
`OptimizationIpopt` dep. `BiophysicalBehaviour` depends on
`Optimization.jl` + `SciMLBase`; `OptimizationIpopt` + `Ipopt` sit in
`ext/BiophysicalBehaviourIpoptExt.jl` (loaded when `OptimizationIpopt`
is present) so users of the rule-based path don't have to build Ipopt.
Alternative optimisers (KNITRO, MadNLP, Percival) plug in as sibling
extensions.

**`AbstractControlStrategy` maps to solver machinery:**
`RuleBasedSequentialControl` → `reinit!` / `solve!` on the iterative
`HeatBalanceSolver`. `IPOPTControl` → `OptimizationProblem` via the NLP
API + `SciMLBase.solve` with `IpoptOptimizer`. Both share the
`HeatBalanceProblem` and per-part `solve_part_heat_balance` primitive —
no forked effector code path.

**Coefficient plumbing — one home per fact, walked at solve time.**
Today's `_run_ipopt` reads `limits.insulation.dorsal.reference` — a
piloerection range stored on `ThermoregulationLimits`, parallel to
`organism.body.insulation.dorsal`. Two homes for one fact; tolerable
for one shape, unmaintainable for N parts. The refactor moves every
per-part fact to a single home — the part itself — and template
builders walk the organism at solve time.

*Where each fact lives:*

| Fact | Home |
|---|---|
| Insulation properties, piloerection range | `part.physiology.insulation`, `.insulation_limits` |
| Skin wetness range, flesh conductivity range | `part.physiology.skin_wetness_range`, `.flesh_conductivity_limits` |
| Panting capacity | field on the `LungPart` wrapper around the lung part's physiology (§3.9); not present on other parts |
| Setpoint, Q10, minimum heat flow, penalty weights, magic-number replacements | `ThermoregulationLimits` (scalars — same as today) |

For the common case (uniform tuning across parts) the user surface is
exactly what it is today — scalar fields on `ThermoregulationLimits`
plus the per-part physiology the user already writes.

Per-part variation is opt-in: any scalar penalty field on
`ThermoregulationLimits` accepts either `Float64` (broadcast to every
leaf) or a `NamedTuple{PartNames}` (per-part override). A `weight_for(w,
part)` helper handles both forms.

`nlp_template`, `nlp_variable_bounds`, and `nlp_objective_templates`
iterate `problem.body.parts` and derive their `NamedTuple` outputs from
per-part physiology plus whole-organism scalars — the user writes none
of them. The objective is one generic reduction: `sum over leaves of
w · ((x − ref) / range)²`, plus a possibly-empty list of pairwise
`cross_penalty(kind, obj)` terms (today's `gradient_penalty` reduces to
a `CoreSkinGradient{part, compartment}` entry). Ranges auto-derive from
bounds as `upper − lower`, guarded by `minimum_normalisation_range`.

*Magic numbers move to named `ThermoregulationLimits` fields:*

| Current inline value | New field |
|---|---|
| `air_temperature - 5.0` | `skin_temperature_undershoot::Temperature = 5.0u"K"` (subtracted from the *coldest* radiant boundary — air, sky, or ground) |
| `core_temperature_max + 5.0` | `skin_temperature_core_overshoot::Temperature = 5.0u"K"` |
| `heat_flow_min_W * 20.0` | `metabolic_heat_flow_max_multiplier::Float64 = 20.0` |
| `max(_, 1e-6)` guards | `minimum_normalisation_range::Float64 = 1e-6` |

The `insulation_depth_max` / `aspect_ratio_min` initial heuristics move
to `default_initial_template(problem, environment, limits)`, which
consults each part's own physiology. Warm-starting is trivial: the
previous solve's output *is* a `NamedTuple` of the same shape, so
`initial_template = previous_solution`.

The Q10 constraint stays a single whole-organism inequality reading two
named leaves off the reconstructed template — no index arithmetic — and
tracks whichever part hosts the lung via `lung_part` from 3.9.

### 3.8 State schema (per-part, DEBTOOL_J-style)

State follows `DEBTOOL_J`'s pattern — nested `NamedTuple` template fixed at
the organism type level, produced by a generated `initialise_state`. Per-part
state is `(; skin_temperature, insulation_temperature)`; compartment state
is `(; core_temperature)`. State flows through the outer loop as an
immutable `NamedTuple`. For future DiffEq / ForwardDiff paths, a
`StateReconstructor`-style wrapper (`DEBTOOL_J/src/animals/ode/utils.jl`)
flattens to `SVector` and rebuilds on demand — not needed immediately.

### 3.9 Lung-part identification

Respiratory heat and mass exchange happens primarily at the **lung
surface**, which sits in the torso (not the head — the head hosts the
airway opening but not alveolar gas exchange). `HomoTherm.R` gets this
right: it computes `TLUNG` from `parts[[2]]` and sets `RESPIRE = 0` on
all other parts.

Wrap the lung-hosting part's physiology in a `LungPart{Physiology,
PantingCapacity}` type that adds a `panting_capacity` field. Only that
one part carries the wrapper — non-lung parts have plain physiology and
no zero-panting field. Dispatch on `::LungPart` routes respiration,
metabolic O₂/CO₂ exchange, and pulmonary evaporation to the correct
part; there is no per-part `panting_capacity` field on other parts.

Add `OrganismTraits.lung_part::Symbol` as an O(1) index into
`body.parts` (e.g. `:torso`, or `:body` for single-part organisms).
Semantics: `lung_temperature` derives from `lung_part`'s
`core_temperature`; `Pant`'s default selector is
`ByName{(lung_part,)}`; respiratory evaporative heat loss is subtracted
from `lung_part`'s energy balance. Single-`Body` organisms wrap the
sole part in `LungPart`. Nasal passages are a separate site (deferred
with `nasal_part` to §8.4). Add both in Phase 4.

### 3.10 Effector + selector interface

Selectors let a behavior target the whole body, a subset, or a single
part. `WholeBody()` degenerates to the current single-body API.

```julia
abstract type PartSelector end
struct WholeBody         <: PartSelector end
struct ByName{Names}     <: PartSelector end
struct Compartment{Name} <: PartSelector end     # Vasodilate / Hyperthermia

abstract type Effector end
struct Piloerect <: Effector end
struct Uncurl <: Effector end
struct Vasodilate <: Effector end
struct Hyperthermia <: Effector end
struct Pant <: Effector end
struct Sweat <: Effector end

effect(op, selector, organism, limits) -> (new_limits, new_organism)
apply(op, part, part_limits)           -> (new_part, new_part_limits)

piloerect(organism, limits) = effect(Piloerect(), WholeBody(), organism, limits)
```

Scope: `Piloerect` / `Uncurl` / `Sweat` are per-part; `Vasodilate` /
`Hyperthermia` are per-compartment; `Pant` is routed to `lung_part` by
default.

---

## 4. Dependency graph

```
BiophysicalGeometry
  └── owns AbstractShape hierarchy + family abstracts,
      AbstractInsulation hierarchy + family abstracts,
      Body, CompositeBody, Join, Pose, silhouette, meshes

HeatExchange
  └── depends on BiophysicalGeometry (types only — no function calls in hot path)
  └── depends on CommonSolve
  └── depends on Flatten (NLP variable-layout pack / unpack — same pattern
      DEBTOOL_J uses for state reconstruction; see 3.7)
  └── owns thermal physics (dispatched on family abstracts),
      HeatCoupling hierarchy, CompartmentGraph,
      HeatBalanceProblem + HeatBalanceSolver,
      CommonSolve interface (init / solve! / reinit! / solve),
      multi-part NLP interface (nlp_template / nlp_pack / nlp_residuals! /
      nlp_variable_bounds / nlp_assemble_output)
  └── does NOT depend on Optimization.jl / OptimizationIpopt / JuMP.
      NLP interface is `Optimization.jl`-agnostic (pure Julia function
      contracts).

BiophysicalBehaviour additional deps for the IPOPT control strategy:
  └── depends on Optimization + SciMLBase (control-strategy driver machinery)
  └── OptimizationIpopt + Ipopt live either as direct deps or in
      ext/BiophysicalBehaviourIpoptExt.jl (loaded when OptimizationIpopt is
      present). Extension only defines `thermoregulate(::IPOPTControl, ...)`.

BiophysicalBehaviour
  └── depends on both, plus CommonSolve for the client-side interface
  └── owns Body cache (Tier 1 + Tier 2), effectors, selectors,
      thermoregulation control loop, radiation geometry (silhouette caching),
      state schema
  └── drives HeatExchange via init / reinit! / solve! in the control loop
```

Ecosystem reference (not modified by this refactor):
- `DEBTOOL_J` — source of state-management patterns (`StateReconstructor` /
  `Flatten.reconstruct` for pack / unpack, `initialise_state`,
  fold-over-transitions). Directly borrowed by the NLP variable-layout
  mechanism in 3.7 to generalise IPOPT over arbitrary multi-part organisms
  without a bespoke index table.
- `ModelParameters` — source of `Flatten` / `Setfield` patterns (already used).

---

## 5. Ordered execution roadmap

Each phase has entry criteria (what must be true to start), work (the changes),
and exit criteria (how we know it's done). Phases can be committed
independently; the numeric regression test after Phase 6 gates the deletion
phase.

### Phase 0 — Shape and insulation family hierarchy
**Package:** `BiophysicalGeometry`
**Work:**
- Add `AbstractCylindrical`, `AbstractSpherical`, `AbstractEllipsoidal`,
  `AbstractSlab` between `AbstractShape` and concrete shapes.
- Add `AbstractFibrous`, `AbstractLayer` between `AbstractInsulation` and
  concrete insulations.
- Update concrete `struct` declarations to use the family supertype.

**Exit:** `BiophysicalGeometry` tests pass unchanged (dispatch on
`AbstractShape` still resolves; family types are additive).

### Phase 1 — Join geometry accessors
**Package:** `BiophysicalGeometry`
**Work:**
- `join_area(join, body)`, `join_position(join, body)`, `join_partners(join)`,
  `internal_distance(shape, surface_location)`.
- Unit tests for each on the existing `examples/dog.jl` composite body.

**Exit:** New accessors covered by tests. No changes to existing API.

### Phase 2 — Part selectors and traversal primitives
**Package:** `BiophysicalBehaviour`
**Work:** New `src/parts.jl`:
- `PartSelector` hierarchy: `WholeBody`, `ByName{Names}`, `Compartment{Name}`.
- `map_parts(f, selector, body)`, `foldl_parts(f, init, selector, body)`.
- `set_part(body, name, new_part)`.
- Degeneracy: `Body` treated as a 1-part `CompositeBody` for traversal.

**Exit:** Traversal tests round-trip on both `Body` and `CompositeBody`.

### Phase 3 — Effector refactor (single-body semantics preserved)
**Package:** `BiophysicalBehaviour`
**Work:** Rewrite `src/endotherm/thermoregulation.jl` around
`effect(op, selector, organism, limits)` + `apply(op, part, part_limits)`:
- Effector tags: `Piloerect`, `Uncurl`, `Vasodilate`, `Hyperthermia`, `Pant`,
  `Sweat`.
- `piloerect(organism, limits) = effect(Piloerect(), WholeBody(), organism, limits)`
  keeps the old scalar API.
- Per-part `apply` methods for `Piloerect`, `Uncurl`, `Sweat`.
- Compartment-scoped `apply` for `Vasodilate`, `Hyperthermia`.
- Rename existing abbreviated identifiers (`T_skin`, `Q_gen`, etc.) inside
  the touched files to their spelled-out forms (see naming conventions above).

**Exit:** Existing `test/endotherm.jl` passes unchanged. No numerical drift.

### Phase 4 — Per-part physiology and limits
**Package:** `BiophysicalBehaviour`
**Work:**
- `OrganismTraits.physiology` becomes `NamedTuple` keyed by part name; falls
  back to a scalar broadcaster for single-`Body` construction.
- **Co-locate per-part limits with per-part physiology** (see 3.7
  "Coefficient plumbing"). Each entry in `OrganismTraits.physiology`
  carries its own `insulation_limits`, `skin_wetness_range`,
  `flesh_conductivity_limits` — the same struct that owns the fur owns
  the piloerection range for that fur. Delete the parallel per-part
  limits containers this replaces (`InsulationLimits.dorsal` /
  `.ventral`, etc.) as part of Phase 8.
- `ThermoregulationLimits` keeps only whole-organism scalars: control
  strategy, setpoint / Q10 / minimum heat flow, penalty weights, and
  the four magic-number-replacement fields
  (`skin_temperature_undershoot`, `skin_temperature_core_overshoot`,
  `metabolic_heat_flow_max_multiplier`, `minimum_normalisation_range`).
  Each penalty-weight field accepts either a `Float64` (broadcast to
  every relevant leaf) or a `NamedTuple{PartNames}` (opt-in per-part
  override); the template builder honours whichever form was given.
- Wrap the lung-hosting part's physiology in `LungPart(physiology,
  panting_capacity)` (§3.9). `panting_capacity` lives only on the
  wrapper — nowhere else. Add `OrganismTraits.lung_part::Symbol` as an
  O(1) index into `body.parts`. Route `Pant`'s default selector, plus
  respiratory evaporative heat loss and O₂/CO₂ exchange attribution,
  through `lung_part`. Derive `lung_temperature` from `lung_part`'s
  `core_temperature`. Single-`Body` organisms wrap the sole part.

**Exit:** Existing tests pass. New tests construct a 2-part organism
with distinct per-part physiology *and* per-part insulation limits set
from the parts themselves — no `ThermoregulationLimits` per-part table,
no `panting_capacity` on non-lung parts. Panting on a multi-part
organism affects only the lung part's evaporative flux and energy
balance.

### Phase 5 — Cache layer + `invalidates` traits
**Package:** `BiophysicalBehaviour`
**Work:** New `src/cache.jl`:
- `PrecomputedCache{Parts}` — Tier 1 (shape) + Tier 2 (radiation) `NamedTuple`s
  keyed by part name.
- `refresh!(cache, body, environment, ::Effector)` dispatched directly on
  effector type. Default method returns `cache` unchanged (compiler
  elides it at the call site); `::Piloerect`, `::Uncurl`, `::Vasodilate`
  do actual work.
- Cache is a field on `Organism`, updated via `@set` after each effector.

**Exit:**
- `JET.@report_opt` on `thermoregulate(organism, environment, ...)` returns no
  dynamic dispatches, runtime type-instability warnings, or method errors
  across the whole call graph (not just the outer frame).
- `@allocated` on a single effector-loop iteration is 0.
- Spot-check with `Cthulhu.@descend` on the effector loop and any hot inner
  method — confirm inferred types are concrete throughout, no `Any` /
  `Union{...}` in the descended frames.

### Phase 6 — HeatExchange decoupling + `CommonSolve` interface
**Package:** `HeatExchange` (major), `BiophysicalBehaviour` (integration)
**Work:**
- `HeatExchange` dispatch changes from concrete shapes to family abstracts.
- Introduce the per-part physics primitive `solve_part_heat_balance` —
  scalar-in, scalar-out, taking a `geometry` `NamedTuple` prepared by
  `BiophysicalBehaviour`. Replaces the geometry-coupled `solve_metabolic_rate`.
- Add `CommonSolve.jl` as a direct dependency of `HeatExchange`. Implement
  `HeatBalanceProblem`, `HeatBalanceSolver`, `CommonSolve.init`,
  `CommonSolve.solve!`, `CommonSolve.reinit!`, `CommonSolve.solve`.
  The `solve!` implementation iterates: per-part `solve_part_heat_balance`
  → assemble `heat_load_vector` → compartment solve (`SMatrix`) → propagate
  back to per-part `skin_temperature` → check convergence.
- Remove all `total_area(body)`, `silhouette_area(body, ...)`,
  `evaporation_area(body)`, `flesh_volume(body)`,
  `characteristic_dimension(body)` calls from `biophysics.jl`,
  `heat_balance.jl`, `nlp_interface.jl`, `insulated/*`.
- Rename touched identifiers to spelled-out forms.
- New `src/radiation_geometry.jl` in `BiophysicalBehaviour`: consumes
  `BiophysicalGeometry.silhouette_rasterized` (and cheaper closed-form
  silhouettes where applicable), produces Tier 2 cache entries.
- `BiophysicalBehaviour` wires the effector loop around `init` / `reinit!` /
  `solve!` on a `HeatBalanceSolver`. Warm-starts each iteration by default.

**Exit:**
- `grep -r "BiophysicalGeometry\." HeatExchange/src/` returns only import
  statements and abstract-type references — no function calls.
- Standalone `solve(HeatBalanceProblem(...))` runs a full multi-part solve
  with no `BiophysicalBehaviour` present in the environment.
- Existing single-body tests still pass (dorsal/ventral machinery still
  present at this point but bypassed by the new path).

### Phase 7 — Compartment graph and coupled solve
**Package:** `HeatExchange` (major), `BiophysicalBehaviour` (integration)
**Work:**
- In `HeatExchange`: define `HeatCoupling` hierarchy (`SharedCore`,
  `ConductiveCoupling{InterfaceConductivity}`). Joins without an entry
  in `couplings` default to zero contribution — no `InsulatedJoin`
  sentinel type.
- In `HeatExchange`: define `CompartmentGraph{Groups,EdgeIndices,NumCompartments}`,
  built from a `NamedTuple` of couplings via union-find over `SharedCore`
  edges. Add `contribution_to_conductance_matrix` and
  `contribution_to_heat_load` methods dispatched on `HeatCoupling` subtypes.
- In `HeatExchange`: extend `HeatBalanceProblem` (Phase 6) with `couplings`
  and `compartment_graph` fields; extend `solve!` to build the
  `conductance_matrix` (`SMatrix`, cached on the solver's workspace, only
  rebuilt when a `Vasodilate` reinit! passes a new `flesh_conductivity`) and
  solve
  `core_temperatures = conductance_matrix \ heat_load_vector` each iteration.
- In `BiophysicalBehaviour`: `OrganismTraits.couplings::NamedTuple` parallel
  to `body.joins`, populated at organism construction. Passed straight
  through to `HeatBalanceProblem`.
- Convergence criterion becomes per-compartment inside `solve!`.

**Exit:**
- 1-part organism produces identical numerics to Phase 6.
- 2-part `CompositeBody` of two `HalfCylinder`s with matched fur and
  `SharedCore` reproduces today's dorsal/ventral results within tolerance
  (equivalence test — this is the gate for Phase 8).

### Phase 7.5 — Multi-part NLP interface + IPOPTControl
**Package:** `HeatExchange` (major), `BiophysicalBehaviour` (major)
**Entry criterion:** Phase 7 iterative solve landed and passes equivalence
tests. The IPOPT branch
(`update-IPOPT-implementation-to-new-HeatExchange-names-and-NLP-hook`) is
the starting reference — its `src/endotherm/thermoregulation/ipopt.jl`
already carries the numerical invariants worth preserving (Float64 /
Unitful boundary, log-transform on `metabolic_heat_flow`, quadratic
penalty objective, Q10 inequality). This phase generalises that driver to
multi-part.
**Work:**
- In `HeatExchange`: rewrite `src/nlp_interface.jl` around the multi-part
  API — `nlp_template`, `nlp_pack`, `nlp_residuals!`,
  `nlp_variable_bounds`, `nlp_assemble_output`. Delete
  `WeightedMeanNLP{Packed}`, `MultiSidedNLP{Packed}`, and their
  strategy-specific `nlp_pack` / `nlp_residuals` / `nlp_assemble_output`
  methods — the multi-part `NamedTuple`-template API supersedes both.
  Variable layout is a nested `NamedTuple` template with `Unitful` leaf
  values, built from `problem.body.parts`, `problem.couplings`, and
  `problem.compartment_graph`. Pack / unpack is `Flatten.flatten` /
  `Flatten.reconstruct` — the same mechanism `DEBTOOL_J`'s
  `StateReconstructor` (`DEBTOOL_J/src/animals/ode/utils.jl:39`) uses to
  bridge nested state and DifferentialEquations' flat `u` vector. Add
  `Flatten` as a direct `HeatExchange` dep. All physics reuses
  `solve_part_heat_balance`, `contribution_to_conductance_matrix`,
  `contribution_to_heat_load` from Phases 6-7. Zero physics duplication.
  `HeatExchange` gains no `Optimization.jl` dependency.
- In `BiophysicalBehaviour`: rewrite the branch's
  `src/endotherm/thermoregulation/ipopt.jl` for multi-part. Delete the
  strategy-specific `_objective_value_weighted` /
  `_objective_value_multisided` and `_heat_balance_residuals_weighted!` /
  `_heat_balance_residuals_multisided!` pairs (with their hardcoded
  `effectors[1]`..`effectors[N]` indices and per-formulation bounds
  vectors). Replace with a single template-driven objective and residual
  closure that reads named fields off the reconstructed template — no
  index arithmetic. Introduce
  `IPOPTControl{Objective, Options} <: AbstractControlStrategy`. Its
  `thermoregulate` method builds an `OptimizationFunction` +
  `OptimizationProblem` from the `HeatBalanceProblem` and `IPOPTControl`,
  then calls `SciMLBase.solve` with `IpoptOptimizer`.
- Preserve numerical invariants from section 3.7: Float64 / Unitful
  boundary at the single `Flatten.reconstruct` call inside
  `nlp_residuals!` (attach units) and the matching `ustrip` on residual
  writes (strip units); log-transform on `metabolic_heat_flow`
  (represented as a `log_metabolic_heat_flow` leaf on the template);
  per-part quadratic penalties in the objective (extending the current
  5-term objective to per-part `insulation_depth` deviation, per-part
  `skin_wetness` penalty, etc.); single whole-organism Q10 inequality.
- Bounds, weights, references, and ranges are all machine-derived —
  the user authors none of them. `nlp_variable_bounds(problem, limits)`
  and `nlp_objective_templates(problem, limits)` walk
  `problem.body.parts` and read per-part limits from each part's own
  physiology (Phase 4); whole-organism scalars come from
  `ThermoregulationLimits`. The three coefficient templates
  (`penalty_weights_template`, `reference_values_template`,
  `normalisation_ranges_template`) match the shape of `nlp_template`
  and are produced by the builder — user writes them nowhere.
  `normalisation_ranges_template` auto-derives as
  `upper_leaf - lower_leaf`, guarded by `minimum_normalisation_range`.
- Replace the five-term hand-written objective with the generic
  `_leaf_objective` reduction plus a possibly-empty `cross_penalty`
  list on `IPOPTControl.objective`.
- `ThermoregulationLimits` gains named fields for the current inline
  magic numbers (`skin_temperature_undershoot`,
  `skin_temperature_core_overshoot`,
  `metabolic_heat_flow_max_multiplier`, `minimum_normalisation_range`);
  penalty-weight fields accept `Float64` (broadcast) or
  `NamedTuple{PartNames}` (opt-in per-part override) via a
  `weight_for(w, part)` helper.
- Replace inline `insulation_depth_max` / `aspect_ratio_min` initial
  heuristics with a `default_initial_template(problem, environment,
  limits)` helper. Wire warm-starting to reuse the previous solve's
  assembled `NamedTuple` output directly as the next initial template.
- `RuleBasedSequentialControl` continues to use the `IterativeStrategy`
  `HeatBalanceSolver` from section 3.6 unchanged.
- Move `OptimizationIpopt` and `Ipopt` into
  `ext/BiophysicalBehaviourIpoptExt.jl` (loaded when `OptimizationIpopt`
  is present), so users of the rule-based path don't have to build Ipopt.
  Only `thermoregulate(::IPOPTControl, ...)` lives in the extension.

**Exit:**
- On a 1-part organism, `IPOPTControl(...)` produces the same converged
  state as `RuleBasedSequentialControl` (numerical equivalence within
  Ipopt's tolerance).
- On the 6-part dog, `IPOPTControl` converges with per-topology residual
  and variable counts (18 residuals, 42 decision variables for the
  `SharedCore`-torso dog); solution matches iterative-strategy solution
  on the same problem within tolerance.
- No `Optimization.jl` / `OptimizationIpopt` / `Ipopt` imports appear in
  `HeatExchange/src/`. `HeatExchange` remains loadable without them.
- `BiophysicalBehaviour` loads without `OptimizationIpopt` installed;
  `IPOPTControl` is available only when the extension has loaded.

### Phase 8 — Deletion pass
**Package:** `HeatExchange` (major), `BiophysicalBehaviour` (minor)
**Entry criterion:** Phase 7 equivalence test passes.
**Work:**
- Delete `traits.jl:21-28` (`MultiSided`, `SingleBody`, `EvaluationStrategy`).
- Delete `insulated/types.jl` (`Dorsal`, `Ventral`, `BodySide`).
- Delete `insulated/utils.jl:60-230` (`_pack_sides`).
- Delete `insulated/mean_skin_temperature.jl`.
- Delete `RadiationParameters.ventral_fraction` and all uses.
- Delete `dorsal`/`ventral` fields on `AbsorptivityParameters`,
  `EmissivityParameters`, `InsulationParameters` — becomes single-valued.
- Delete `solve_temperature.jl` `MultiSided` branch.
- Delete `WeightedMeanNLP`, `MultiSidedNLP`, `WeightedMeanNLPPacked`,
  `MultiSidedNLPPacked` and their `nlp_pack` / `nlp_residuals` /
  `nlp_assemble_output` methods from `nlp_interface.jl`. The multi-part
  API from Phase 7.5 supersedes both. Delete the mean-weighting block
  (`nlp_interface.jl:74-97`) and the multi-side assembly
  (`nlp_interface.jl:102-133`) that depend on `_pack_sides`.
- Migrate example organism parameterisations to the composite form.
- Migrate test assertions from `T_skin_dorsal`/`T_skin_ventral` to
  `state.parts.<name>.skin_temperature`.

**Exit:**
- `grep -ri "dorsal\|ventral\|ventral_fraction" HeatExchange/src/` returns
  zero results.
- All tests pass with the composite representation.

### Phase 9 — Multi-part example + tests
**Package:** `BiophysicalBehaviour`
**Work:**
- Port `BiophysicalGeometry/examples/dog.jl` as
  `examples/dog_thermoregulation.jl`.
- `test/multipart.jl`:
  (i) single-`Body` regression identical to Phase 6 results;
  (ii) 2-part composite with `SharedCore` matches today's dorsal/ventral
       (equivalence test from Phase 7 gate);
  (iii) selector round-trip (`ByName{(:head,)}` changes only the head).

**Exit:** All tests pass. Dog example runs and produces sensible output.

### Phase 10 — Cleanup
**Package:** `BiophysicalBehaviour`
**Work:**
- Delete `ectotherm.jl` at the repo root (stray 456-line file, not `include`d
  anywhere).
- Remove now-unused imports in `src/BiophysicalBehaviour.jl:12`.
- Update `README.md`.

---

## 6. Global acceptance criteria

Beyond per-phase exits:

- **Full-call-graph type stability**, not just outer frame.
  `@code_warntype` only inspects the top-level return type — it will not
  catch a boxed variable or dynamic dispatch inside a called function. Use
  instead:
  - **`JET.@report_opt thermoregulate(organism, environment, ...)`** — static
    analysis across the whole call graph. Should report no dynamic dispatches,
    no runtime `Union{...}` splits, no method errors.
  - **`Cthulhu.@descend`** on the effector loop and on any hot inner method
    (per-part `solve_metabolic_rate`, `refresh!`, compartment solve) — walk
    the call tree and confirm inferred types are concrete at each level.
  - **Sampled profiling** — `Profile.@profile thermoregulate(...)` plus
    `ProfileView.view()` (red bars indicate type-unstable / dynamically
    dispatched frames). Also usable from `Profile.print()` text output when
    a GUI isn't available.
- `BenchmarkTools.@benchmark thermoregulate($organism, $environment)` on a
  1-part body reports zero allocations per iteration and stays within a few
  percent of the pre-refactor timing.
- On a 6-compartment dog composite, the compartment
  `SMatrix{NumCompartments,NumCompartments} \ SVector{NumCompartments}` solve
  contributes negligible overhead compared to the per-part
  `solve_metabolic_rate` calls.
- No `Ref`, no `Dict`, no growing `Vector`, no runtime closures in the
  effector loop.
- **No new abbreviated identifiers** in code added or edited during any phase.
  Identifiers in touched files are renamed to spelled-out forms as part of the
  same phase.

---

## 7. Explicit non-goals (for this refactor)

- `PDEControl` strategy — remains stubbed. The state schema and cache design
  are compatible with a future PDE path but no PDE implementation now.
- `Ectotherm` and `Heterotherm` `thermoregulate` — remain stubs.
- Posture change behaviors (`Roll`, `Curl`, `Stand`, `LieDown`) — out of scope.
  Tier 2 cache design accommodates them for the future.
- `TriMesh` thermal physics — no family membership, will error usefully at
  dispatch time.
- Backwards compatibility for `ventral_fraction`-calibrated organisms —
  document the translation recipe (`ventral_fraction = x` ≡ two
  `HalfCylinder` parts with mass fractions `x` / `1 − x` joined by
  `SharedCore`) but do not provide shims.
- The four extensions in section 8 — see below.

## 8. Future extensions (design must not preclude)

Four extensions we know we'll need eventually. Not built here; each one lists
the aspect of today's design that must stay open so a future add-on doesn't
force a rewrite.

- **8.1 Concentric layered shapes** (`LayeredShape{Shape}`, e.g. brain → CSF
  → skull → skin). Requires: per-part state schema extensible from
  `(skin_temperature, insulation_temperature)` to include
  `layer_temperatures::NTuple{NumLayers}`; compartment graph nodes at
  layer-granularity, not part-granularity (union-find over nodes, not parts —
  Phase 7).
- **8.2 Part-to-part view factors** (leg occludes belly's view of sky).
  Requires: Tier 2 radiation cache with room for a per-part incoming-IR
  vector; compartment solve extensible from block-diagonal (independent
  parts) to a full linear system so cross-part radiation off-diagonals can
  be added.
- **8.3 Local convection per part** (leeward parts see less wind).
  Requires: Tier 2 cache accepts an extra per-part `flow_factor` scalar
  (default 1.0). No structural change.
- **8.4 Selective brain cooling** (nasal countercurrent + carotid rete).
  Requires: `HeatCoupling` hierarchy open to new subtypes via
  `contribution_to_conductance_matrix` /
  `contribution_to_heat_load` methods (so `BloodPerfusionCoupling` is a
  new method, not a rewrite); `lung_part` from 3.9 in place so nasal /
  pulmonary flux can be distinguished when a future `nasal_part` lands.

- **8.5 Animal-in-shelter coupling** (e.g. possum in a nestbox, animal in
  a burrow). The shelter is a *second heat-exchanging body* — its own
  geometry, wall insulation, and heat balance — that both modifies the
  microclimate the animal sees (attenuates wind, replaces sky/ground
  radiant boundaries with interior-wall temperatures, traps humidity) and
  is itself warmed by the animal. Their temperatures are mutually
  dependent, so it is a coupled solve, not a one-way environment tweak.
  Old monolith code for this exists but was never re-integrated post-split;
  not built here. Requires two things to stay open, both already implied by
  the "nodes, not parts" principle (§3.4) plus 8.2/8.4:
  (a) a `HeatBalanceProblem`'s thermal-node graph is **not hardwired to a
  single organism** — it can hold nodes belonging to more than one body
  (organism parts + shelter walls), coupled by non-flesh `HeatCoupling`
  subtypes (`RadiativeEnclosureCoupling`, `CavityAirCoupling`) added as new
  `contribution_to_*` methods, not a rewrite;
  (b) the per-node **environment/boundary is resolvable from other nodes'
  state** — an animal node's radiant/convective boundary can be a shelter
  wall node's temperature (a solve variable), not only a fixed global
  `environment`. Keep `environment` access in `solve_part_heat_balance`
  behind a per-part boundary lookup so a shelter can later inject
  node-derived boundaries. The Tier 2 radiation cache's per-part
  incoming-IR vector (8.2) already carries enclosure view factors.

**Cross-cutting principle:** state schema, compartment graph, radiation
cache, and coupling interface are all designed additively — new node types,
coupling types, and cache fields plug in without rewriting the outer loop.
Every phase-8 exit criterion is checked against "would this survive one of
these extensions?"
