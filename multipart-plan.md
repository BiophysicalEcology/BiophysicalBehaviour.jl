# BiophysicalBehaviour — Multi-part refactor plan

Working design + ordered execution roadmap for updating `BiophysicalBehaviour` to
handle arbitrary multi-part organism geometries defined in `BiophysicalGeometry`,
alongside coordinated refactors to `HeatExchange` and `BiophysicalGeometry`.

We control all four packages (`BiophysicalGeometry`, `HeatExchange`,
`BiophysicalBehaviour`, and, tangentially, `DEBTOOL_J` for state-management
patterns and `ModelParameters` for parameter tooling). Nothing is registered.
Backwards compatibility is not a constraint. Simplicity, type stability, and
zero-allocation hot paths are.

## Naming conventions used throughout

**Fully spelled identifiers.** No `T_` for temperature, no `Q_` for heat, no
`k_` for conductivity, no `Ra`/`Pr`/`Nu` for dimensionless numbers, no `G`/`b`
for matrix/vector, no `evap`/`resp`/`pars`/`env`/`geom` abbreviations.
Existing code that uses abbreviated identifiers will be renamed as part of the
phase that touches it. When this document quotes existing code or names
(e.g. current field names, NicheMapR reference material), the original
identifiers are kept as-is for accuracy.

Standard target substitutions used throughout this plan:

| Abbreviated | Full |
|---|---|
| `T_skin`, `T_insulation`, `T_core`, `T_air`, `T_ground`, `T_sky` | `skin_temperature`, `insulation_temperature`, `core_temperature`, `air_temperature`, `ground_temperature`, `sky_temperature` |
| `T_lung`, `TLUNG`, `T_brain`, `T_layers` | `lung_temperature`, `brain_temperature`, `layer_temperatures` |
| `Q_gen`, `Q_metab`, `Q_env`, `Q_conv`, `Q_rad`, `Q_evap`, `Q_resp` | `metabolic_heat`, `environmental_heat`, `convective_heat`, `radiative_heat`, `evaporative_heat`, `respiratory_heat` |
| `Q_minimum` | `minimum_metabolic_heat` |
| `k_flesh`, `k_interface`, `k` | `flesh_conductivity`, `interface_conductivity`, `conductivity` |
| `Ra`, `Pr`, `Nu`, `Re` | `rayleigh_number`, `prandtl_number`, `nusselt_number`, `reynolds_number` |
| `G` (matrix), `b` (vector) | `conductance_matrix`, `heat_load_vector` |
| `IR`, `SBC` | `infrared`, `selective brain cooling` |
| `evap_area`, `evap_source` | `evaporation_area`, `evaporation_source` |
| `geom`, `env`, `pars` (new code) | `geometry`, `environment`, `parameters` |
| `A_direct`, `A_diffuse_up`, `A_diffuse_down` | `area_direct`, `area_diffuse_up`, `area_diffuse_down` |

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

- **Effectors** (`src/endotherm/thermoregulation.jl`) — `piloerect`, `uncurl`,
  `vasodilate`, `pant`, `sweat` each read `HeatExchange.<foo>pars(organism)`
  and rebuild `organism.body = Body(...)`. They will hard-fail on
  `CompositeBody` (no `.insulation` field) and would erase part structure even
  if they didn't.
- **State** (`src/endotherm/homeothermy.jl:98-229`) — skin, insulation, core
  temperatures and generated metabolic heat are threaded as scalars (existing
  code names them `T_skin`, `T_insulation`, `T_core`, `Q_gen`).
  `initialise_state` exists only as a `NullBehavior` stub
  (`src/organism.jl:126`).
- **Physiology traits** (`src/organism.jl:243-257`) — a single lumped
  `HeatExchangeTraits` with `insulation_pars.dorsal/ventral` fields. No concept
  of per-part insulation, metabolism, or radiation.

`HeatExchange` is tightly coupled to `BiophysicalGeometry`:

- Direct scalar accessor calls (`total_area`, `silhouette_area`,
  `evaporation_area`, `flesh_volume`, `characteristic_dimension`) sprinkled
  through `biophysics.jl`, `heat_balance.jl`, `nlp_interface.jl`,
  `insulated/utils.jl`.
- Shape dispatch for closed-form conduction (`internal_temperature.jl:28-65`,
  cylinder/plate/ellipsoid) and Nusselt correlations
  (`biophysics.jl:313-392`).
- Deep dorsal/ventral split — `_pack_sides` in `insulated/utils.jl:60-230`,
  `MultiSided` in `traits.jl:21-28`, `Dorsal`/`Ventral` singletons, dorsal/ventral
  fields on `RadiationParameters`, `InsulationParameters`, `AbsorptivityParameters`,
  `ventral_fraction`, `mean_skin_temperature.jl`. Roughly half the `insulated/`
  directory exists to service the split.

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
    interface_conductivity::InterfaceConductivity   # optional override;
                                                    # `nothing` → derive from parts
end

struct InsulatedJoin <: HeatCoupling end       # no heat crosses
```

Example dog:
```julia
couplings = (
    dorsal_ventral = SharedCore(),                # torso halves, one blood pool
    dorsal_head    = ConductiveCoupling(nothing), # derive from flesh_conductivity + geometry
    ventral_leg_fl = ConductiveCoupling(nothing),
    ventral_leg_fr = ConductiveCoupling(nothing),
    ventral_leg_bl = ConductiveCoupling(nothing),
    ventral_leg_br = ConductiveCoupling(nothing),
)
```

Coupling policy for `ConductiveCoupling(nothing)` defaults to a series
resistance
`total_resistance = distance_parent / (conductivity_parent * area) + distance_child / (conductivity_child * area)`
computed from each part's own `flesh_conductivity` and `internal_distance` to
the join. `ConductiveCoupling(interface_conductivity)` is an override.

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
contribution_to_conductance_matrix(::SharedCore, ...)          # zero; edge already contracted
contribution_to_conductance_matrix(::ConductiveCoupling, ...)  # symmetric off-diagonal
contribution_to_conductance_matrix(::InsulatedJoin, ...)       # zero
contribution_to_heat_load(::HeatCoupling, ...)                 # right-hand-side shifts
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

**Cache invalidation via `Val`-dispatch:**

```julia
invalidates(::Piloerect)    = (:shape,)
invalidates(::Uncurl)       = (:shape, :radiation)
invalidates(::Vasodilate)   = (:conductance,)          # conductance matrix
invalidates(::Hyperthermia) = ()
invalidates(::Pant)         = ()
invalidates(::Sweat)        = ()

@inline refresh!(cache, body, environment, ::Val{()}) = cache   # compiler elides
@inline refresh!(cache, body, environment, ::Val{(:shape,)}) = ...
```

In the hot loop the compiler sees `Val{()}` for `Vasodilate`, `Hyperthermia`,
`Pant`, `Sweat` and eliminates the call. Only `Piloerect` and `Uncurl` do work.

### 3.6 HeatExchange as pure physics + `CommonSolve` interface

`HeatExchange` loses all geometry function calls and all dorsal/ventral
machinery, and exposes the standard **`CommonSolve.jl`** interface: `init`,
`solve!`, `reinit!`, `solve`. This matches how the rest of the SciML
ecosystem (`DifferentialEquations`, `LinearSolve`, `NonlinearSolve`,
`Optimization`) works, and gives us two usage modes with no code duplication.

**Problem type** — a self-contained multi-part heat-balance problem:

```julia
struct HeatBalanceProblem{Geometry, Physiology, Environment, Couplings,
                         CompartmentGraph, InitialState}
    geometry::Geometry                    # per-part shape + insulation + Tier 1/2 cache scalars
    physiology::Physiology                # per-part physiological parameters
    environment::Environment              # air / ground / sky temperatures, wind,
                                          # solar_direct, solar_diffuse, ...
    couplings::Couplings                  # NamedTuple keyed by join name, values
                                          # <: HeatCoupling
    compartment_graph::CompartmentGraph   # precomputed from couplings (section 3.4)
    initial_state::InitialState           # per-part skin/insulation, per-compartment core
end
```

The problem holds *everything* required to compute a full multi-part coupled
heat balance. No hidden state; no callbacks into `BiophysicalGeometry`.

**Solver type** — mutable workspace + current state, reused across successive
solves:

```julia
mutable struct HeatBalanceSolver{Problem, Workspace, State}
    problem::Problem
    workspace::Workspace   # Jacobian / root-finding buffers / conductance matrix cache
    state::State           # current best estimate of skin, insulation, core temperatures
    converged::Bool
    iterations::Int
end
```

**Interface** (implemented against `CommonSolve.jl`):

```julia
CommonSolve.init(problem::HeatBalanceProblem; kwargs...) -> HeatBalanceSolver
CommonSolve.solve!(solver::HeatBalanceSolver) -> Solution   # run to convergence in place
CommonSolve.reinit!(solver::HeatBalanceSolver, new_problem::HeatBalanceProblem;
                    warm_start::Bool = true)                # reset for a new problem,
                                                            # optionally reusing state
CommonSolve.solve(problem::HeatBalanceProblem; kwargs...) =
    solve!(init(problem; kwargs...))                        # one-shot convenience
```

**Two usage modes fall out:**

1. **One-shot from downstream code** (e.g., someone modelling a static
   organism / environment pair, or a NicheMapR-style time-stepping driver):
   ```julia
   solution = solve(HeatBalanceProblem(geometry, physiology, environment,
                                       couplings, initial_state))
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

Family dispatch in `HeatExchange`:

```julia
nusselt_free(::AbstractCylindrical, rayleigh_number, prandtl_number) = ...  # McAdams
nusselt_free(::AbstractSpherical,   rayleigh_number, prandtl_number) = ...
nusselt_free(::AbstractEllipsoidal, rayleigh_number, prandtl_number) = ...
nusselt_free(::AbstractSlab,        rayleigh_number, prandtl_number) = ...  # Gates plate

surface_temperature(::AbstractCylindrical, radius,         heat_per_volume, conductivity) =
    heat_per_volume * radius^2 / (4 * conductivity)
surface_temperature(::AbstractSpherical,   radius,         heat_per_volume, conductivity) =
    heat_per_volume * radius^2 / (6 * conductivity)
surface_temperature(::AbstractSlab,        half_thickness, heat_per_volume, conductivity) =
    heat_per_volume * half_thickness^2 / (2 * conductivity)
surface_temperature(::AbstractEllipsoidal, semi_axes,      heat_per_volume, conductivity) =
    ...                                                    # exact three-axis formula

insulation_conductance(::Naked, args...)                     = zero(...)
insulation_conductance(::AbstractFibrous, args...)           = ...  # fur, feathers, hair
insulation_conductance(::AbstractLayer,   args...)           = ...  # fat, muscle
insulation_conductance(c::CompositeInsulation, args...)      = sum(insulation_conductance, c.layers)
```

Adding a new concrete shape or insulation is a one-liner in
`BiophysicalGeometry`. Physics works automatically via family dispatch.

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
  `+ 3 · NumParts` (per-part `insulation_depth`, `skin_wetness`,
  `aspect_ratio`)
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
                                skin_wetness           = 0.1,
                                aspect_ratio           = 1.0),
                  torso_dorsal   = (; ...), torso_ventral = (; ...),
                  leg_fl = (; ...), leg_fr = (; ...), ...),
       compartments = (; torso_core = (; core_temperature   = 310.0u"K",
                                          flesh_conductivity = 0.5u"W/m/K"),
                         head_core  = (; ...), leg_fl_core = (; ...), ...),
       whole_organism = (; log_metabolic_heat_flow = log(0.5),
                           pant                    = 1.0))
```

The whole pack / unpack toolkit is three `Flatten` calls, mirroring
`DEBTOOL_J`'s `to_vec` / `to_obj`:

```julia
# Float64 vector (optimizer world) → NamedTuple of Unitful quantities
unpack(template, x) = Flatten.reconstruct(template, x, Number)
                      # per-field convert via typeof(template_leaf),
                      # attaching units back on

# NamedTuple of Unitful quantities → Float64 vector
pack(obj)           = SVector(map(ustrip, Flatten.flatten(obj, Number)))

# template itself → initial Float64 decision vector
initial(template)   = pack(template)
```

**The NLP API in `HeatExchange`** (no JuMP dep, no optimizer dep, no
`Optimization.jl` dep; one new dep: `Flatten`):

```julia
nlp_template(problem::HeatBalanceProblem) -> NamedTuple

nlp_pack(problem::HeatBalanceProblem) -> PackedNLP{Problem, Template, Workspace}
    # bundles problem + template + any per-solve constant physics scalars

nlp_residuals!(residuals::AbstractVector,
               packed::PackedNLP,
               decision_variables::AbstractVector) -> residuals

nlp_variable_bounds(problem::HeatBalanceProblem, limits) ->
    (lower_template, upper_template)
    # two NamedTuples with the same shape as `nlp_template(problem)`;
    # flatten each to a Vector{Float64} via the same `pack` path

nlp_assemble_output(packed::PackedNLP, decision_variables) -> Solution
```

`nlp_residuals!` is the only hot function. Structure:

```julia
function nlp_residuals!(residuals, packed, x)
    obj = Flatten.reconstruct(packed.template, x, Number)   # Unitful, nested
    part_results = map(part -> solve_part_heat_balance(packed.problem, part,
                                                       obj.parts[part], obj),
                       keys(obj.parts))
    fill_compartment_residuals!(residuals, packed, obj, part_results)
    fill_whole_organism_residuals!(residuals, packed, obj, part_results)
    return residuals
end
```

One `Flatten.reconstruct` at the top (Float64 → Unitful), `ustrip` on each
residual as it's written to `residuals`. The optimizer never sees a Unitful
quantity. Same `solve_part_heat_balance` and
`contribution_to_conductance_matrix` / `contribution_to_heat_load` methods
as the iterative solver — **zero physics duplication** between paths.

**Adding a part, compartment, or per-part decision variable is a one-line
change to the template.** Residual counts and variable counts fall out of
the template shape. Bounds and objective-weight tables use the same
template pattern (parallel `NamedTuple`s of matching shape), so a new
decision variable adds one field in three places (template, lower-bound
template, upper-bound template) rather than N places across an indexed
layout table.

**Strategy dispatch — layered on the SciML stack.** Because Optimization.jl
already implements `CommonSolve`, an NLP thermoregulation solve is *already*
a `solve(::OptimizationProblem)` call. We stay on the same rails:

```julia
# HeatExchange stays free of Optimization.jl / Ipopt. It only exposes:
nlp_template(problem)                     -> NamedTuple    # Flatten-backed layout
nlp_pack(problem)                         -> PackedNLP
nlp_residuals!(residuals, packed, x)      -> residuals
nlp_variable_bounds(problem, limits)      -> (lower_template, upper_template)
nlp_assemble_output(packed, x)            -> Solution

# BiophysicalBehaviour owns policy and the driver:
struct IPOPTControl{Objective, Options} <: AbstractControlStrategy
    objective::Objective   # penalty weights, targets, ranges
    options::Options       # Ipopt options (tol, max_iter, hessian_approximation, ...)
end

# thermoregulate on IPOPTControl builds an OptimizationProblem from the
# HeatBalanceProblem + IPOPTControl, then calls SciMLBase.solve:
function thermoregulate(::IPOPTControl, organism, environment, ...)
    heat_balance_problem = HeatBalanceProblem(from_organism(organism), environment, ...)
    packed = HeatExchange.nlp_pack(heat_balance_problem)
    template = HeatExchange.nlp_template(heat_balance_problem)
    lower_template, upper_template =
        HeatExchange.nlp_variable_bounds(heat_balance_problem, limits)

    # Flatten the three parallel NamedTuples to Float64 vectors — same code path.
    lower = collect(map(ustrip, Flatten.flatten(lower_template, Number)))
    upper = collect(map(ustrip, Flatten.flatten(upper_template, Number)))
    x0    = collect(map(ustrip, Flatten.flatten(template,       Number)))

    optimization_function = OptimizationFunction(
        (x, p) -> _objective(x, control.objective, template),
        SciMLBase.NoAD();
        cons = (r, x, p) -> _residuals!(r, x, packed, control),
        grad = _grad!, hess = _hess!, cons_j = _cons_j!, cons_h = _cons_h!,
    )
    optimization_problem = OptimizationProblem(
        optimization_function, x0;
        lb = lower, ub = upper,
        lcons = equality_and_inequality_lower_bounds,
        ucons = equality_and_inequality_upper_bounds,
    )
    SciMLBase.solve(optimization_problem, IpoptOptimizer(; control.options...))
end
```

For the iterative path, `RuleBasedSequentialControl` continues to drive the
`HeatBalanceSolver{IterativeStrategy}` from section 3.6.

**Package layout:**

- `HeatExchange` — no `Optimization.jl` dep, no `OptimizationIpopt` dep.
  Just the physics primitives.
- `BiophysicalBehaviour` — depends directly on `Optimization.jl` and
  `SciMLBase` (control-strategy driver machinery). `OptimizationIpopt` and
  `Ipopt` can be direct deps or moved behind a package extension
  (`ext/BiophysicalBehaviourIpoptExt.jl` loaded when `OptimizationIpopt` is
  present) so that users who only want the rule-based path don't have to
  build Ipopt. The extension only defines
  `thermoregulate(::IPOPTControl, ...)`; every other type and function is
  in the main module.
- Alternative optimisers (KNITRO, MadNLP, Percival) plug in the same way —
  each is a separate extension that adds another `::AbstractControlStrategy`
  method or reuses `IPOPTControl` with a different optimiser factory.

**`AbstractControlStrategy` maps to solver machinery:**

- `RuleBasedSequentialControl` → rebuild `HeatBalanceProblem` after each
  effector via `reinit!`, `solve!` on the `IterativeStrategy`
  `HeatBalanceSolver` from section 3.6.
- `IPOPTControl{Objective, Options}` → build `HeatBalanceProblem` once,
  build an `OptimizationProblem` via `HeatExchange`'s NLP interface, call
  `SciMLBase.solve` with `IpoptOptimizer`.
- `PDEControl` (future) → another strategy, same rails.

All three share the `HeatBalanceProblem` and the per-part
`solve_part_heat_balance` primitive. No forked effector code path.

**Coefficient plumbing — one home per fact, walked at solve time.**

The scaling concern is real: today's `_run_ipopt` reads
`limits.insulation.dorsal.reference` — a piloerection range stored on
`ThermoregulationLimits`, *parallel to and disconnected from* the
actual insulation on the body at `organism.body.insulation.dorsal`.
That's two homes for one fact. For one shape it's tolerable
duplication; for N parts it would force the user to author N per-part
limits entries kept in lock-step with N per-part insulation entries.
That doesn't scale.

The refactor eliminates the parallelism by moving each per-part fact
to exactly one home — the part itself — and keeping only genuinely
whole-organism knobs on `ThermoregulationLimits`. Template builders
walk the organism at solve time; the user authors nothing that lives
in more than one place.

*Where each fact lives:*

| Fact | Home | Notes |
|---|---|---|
| Insulation properties (depth, conductivity, ...) | `part.physiology.insulation` | already there today |
| Insulation piloerection range (min / ref / max depth) | `part.physiology.insulation_limits` | co-located with the insulation |
| Skin wetness range | `part.physiology.skin_wetness_range` | per part; typically zero for parts without sweat glands |
| Panting capacity | `part.physiology.panting_capacity` | zero for non-lung parts (auto-derived from `lung_part` in 3.9) |
| Flesh conductivity range | `part.physiology.flesh_conductivity_limits` | co-located with conductivity |
| Aspect-ratio bounds | dispatched on `part.shape` | `aspect_ratio_bounds(::Sphere) = (1,1)`; no user input |
| Whole-organism setpoint, Q10, minimum heat flow | `ThermoregulationLimits` (scalar) | same as today |
| Whole-organism penalty weights | `ThermoregulationLimits` (scalar; opt-in NamedTuple override) | same as today for uniform case |
| Skin-temperature slack constants, heat-flow multiplier, range guard | `ThermoregulationLimits` (scalar) | replaces inline magic numbers (table below) |

Per-part physiology already becomes a `NamedTuple` keyed by part name in
Phase 4; this section extends each entry to carry its own limits
alongside its physiology. **For the common case (uniform tuning across
parts) the user surface is exactly what it is today** — a handful of
scalar fields on `ThermoregulationLimits`, plus the per-part physiology
structs the user already writes to define the organism.

*Per-part variation is opt-in.* Anywhere a scalar sits on
`ThermoregulationLimits`, swapping it for a `NamedTuple{PartNames}`
gives per-part control:

```julia
# uniform (default):
ThermoregulationLimits(skin_wetness_penalty = 1.0, ...)

# per-part (opt-in):
ThermoregulationLimits(skin_wetness_penalty =
    (head = 0.0, torso_dorsal = 1.0, torso_ventral = 5.0,
     leg_fl = 2.0, leg_fr = 2.0, leg_bl = 2.0, leg_br = 2.0), ...)
```

The template builder honours whichever form was given — same code path.

*Template builders walk the parts; the user writes none of them.*
`nlp_template`, `nlp_variable_bounds`, and `nlp_objective_templates`
each iterate `problem.body.parts` and derive their `NamedTuple` output
from the part's own physiology plus whole-organism limits:

```julia
function nlp_variable_bounds(problem, limits)
    parts_lower = map(problem.body.parts, problem.body.parts_physiology) do part, physiology
        skin_low, skin_high = skin_temperature_bounds(part, problem.environment, limits)
        (; skin_temperature       = skin_low,
           insulation_temperature = skin_low,
           insulation_depth       = physiology.insulation_limits.reference,  # <- from part
           skin_wetness           = physiology.skin_wetness_range.reference, # <- from part
           aspect_ratio           = aspect_ratio_bounds(part.shape, limits)[1])
    end
    ...
end

function nlp_objective_templates(problem, limits)
    parts_weights = map(problem.body.parts) do part
        (; skin_wetness = weight_for(limits.skin_wetness_penalty, nameof(part)),
           ...)
    end
    ...
end

weight_for(w::Real, _)             = w                          # broadcast scalar
weight_for(w::NamedTuple, part)    = w[part]                    # honour per-part override
```

One code path per template, all parts handled uniformly. Adding a part
adds no new user surface — its physiology carries everything the
builders need. Adding a per-part penalty override is a scalar →
`NamedTuple` swap on one field.

*The resulting objective is one generic reduction over leaves:*

```julia
function _leaf_objective(x, template, weights, references, ranges)
    xs  = Flatten.flatten(Flatten.reconstruct(template, x, Number), Number)
    ws  = Flatten.flatten(weights,    Number)
    rs  = Flatten.flatten(references, Number)
    ns  = Flatten.flatten(ranges,     Number)
    sum(w * ((ustrip(x) - ustrip(r)) / ustrip(n))^2
        for (w, x, r, n) in zip(ws, xs, rs, ns))
end
```

Every current objective term reduces to one leaf weight —
`core_temperature_penalty` → `weights.compartments[c].core_temperature`,
`panting_penalty` → `weights.whole_organism.pant`, etc. On a 1-part
organism the scalar broadcasting produces numerically identical
objectives to the current branch.

*Cross-part / cross-compartment penalties.* The current
`gradient_penalty` term reads two leaves (`core - skin`) against a
target — pairwise, not per-leaf. Under multi-part, "core – skin" for
which pair? A small explicit list on `IPOPTControl.objective`:

```julia
struct CoreSkinGradient{PartName, CompartmentName, Weight, Target, Range}
    part::PartName; compartment::CompartmentName
    weight::Weight; target::Target; range::Range
end
```

Objective sums leaf terms plus `cross_penalty(kind, obj) -> Float64`
per entry. Empty list → pure sum-over-leaves objective. New cross-term
physics (per-limb balance, brain-core gradient once selective brain
cooling lands) is a new struct + method, not a rewrite.

*Shape-derivable bounds go via per-shape dispatch — no `isa` branches
and no user input:*

```julia
aspect_ratio_bounds(::Sphere,              limits) = (1.0, 1.0)   # frozen
aspect_ratio_bounds(::AbstractCylindrical, limits) =
    (Float64(limits.aspect_ratio_factor.reference),
     Float64(limits.aspect_ratio_factor.max))

skin_temperature_bounds(part, environment, limits) =
    (min(environment.air_temperature,
         environment.sky_temperature,
         environment.ground_temperature) - limits.skin_temperature_undershoot,
     limits.core_temperature.max + limits.skin_temperature_core_overshoot)
```

A new shape family adds one `aspect_ratio_bounds` method; no other
touchpoint.

*Magic numbers move to named `ThermoregulationLimits` fields* —
tunable, discoverable, unitful, still whole-organism scalar:

| Current inline value | New named field on `ThermoregulationLimits` |
|---|---|
| `air_temperature - 5.0` | `skin_temperature_undershoot::Temperature = 5.0u"K"` (subtracted from the *coldest* radiant boundary — air, sky, or ground) |
| `core_temperature_max + 5.0` | `skin_temperature_core_overshoot::Temperature = 5.0u"K"` |
| `heat_flow_min_W * 20.0` | `metabolic_heat_flow_max_multiplier::Float64 = 20.0` (or a `maximum_heat_flow::Power` field on the lung part's physiology) |
| `panting_rate_min = 1.0` | already the reference on the lung part's `panting_capacity`; expose explicitly |
| `max(_, 1e-6)` guards | `minimum_normalisation_range::Float64 = 1e-6` (single field; applied inside the auto-derivation of `normalisation_ranges_template` from bounds) |

*Initial values via a helper; warm-start via template reuse.* The
inline `insulation_depth_max` / `aspect_ratio_min` heuristics move to
`default_initial_template(problem, environment, limits) -> NamedTuple`,
which fills each per-part leaf from a small set of principled rules
(cold environment → erected fur & curled posture; hot → the opposite)
by consulting each part's own physiology limits. Warm-starting across
effector iterations or across time steps is trivial: the previous
solve's assembled output *is* a `NamedTuple` of the same shape, so
`initial_template = previous_solution` — no repacking.

*Consolidation with Phase 8's deletion pass.* Today's
`limits.insulation.dorsal` and `body.insulation.dorsal` — two homes for
one fact — both fold into `part.physiology.insulation_limits` on a
single `torso_dorsal` part; the `.dorsal` / `.ventral` field indirection
disappears with the rest of the dorsal/ventral machinery in Phase 8.

*Q10 constraint* stays a single whole-organism inequality; the residual
reads two named leaves off the reconstructed template — no index
arithmetic — and tracks whichever part hosts the lung via `lung_part`
from 3.9:

```julia
q10_residual(obj, limits, lung_part) =
    exp(obj.whole_organism.log_metabolic_heat_flow) -
    limits.minimum_heat_flow *
        limits.q10 ^ ((obj.compartments[compartment_of(lung_part)].core_temperature
                       - limits.setpoint_temperature) / 10u"K")
```

### 3.8 State schema (per-part, DEBTOOL_J-style)

Follow `DEBTOOL_J`'s pattern: nested `NamedTuple` template fixed at the
organism type level, generated function that guarantees compile-time shape.

```julia
@generated function initialise_state(::Type{O}) where {O<:Organism}
    ...  # derives keys from body.parts NamedTuple type
    quote
        (; parts = NamedTuple{$partnames}(map(_ -> _part_state_template(), $partnames)),
           compartments = NamedTuple{$compartmentnames}(...))
    end
end
```

Per-part state is `(; skin_temperature, insulation_temperature)`. Compartment
state is `(; core_temperature)`. State flows through the outer loop as an
immutable `NamedTuple`. For future DiffEq / ForwardDiff paths, a
`StateReconstructor`-style wrapper (borrowed from `DEBTOOL_J`'s
`src/animals/ode/utils.jl`) flattens to `SVector` and rebuilds on demand — not
needed immediately but the schema makes it easy.

### 3.9 Lung-part identification

Respiratory heat and mass exchange happens primarily at the **lung surface**,
which sits in the torso — not in the head. The head hosts the airway opening
(nares/mouth) but is not where alveolar gas exchange, air warming, or the bulk
of respiratory water loss occur. `HomoTherm.R` gets this right — it computes
lung temperature from `parts[[2]]` (the trunk, in a `TLUNG` variable) and
disables respiration on all other parts (`RESPIRE = 0` in the per-part `endoR`
loop).

Add to `OrganismTraits`:

```julia
lung_part::Symbol   # e.g. :torso, or :body for single-part organisms
```

Semantics:
- `lung_temperature` is derived from `lung_part`'s `core_temperature` via the
  appropriate closed-form conduction (part shape + interior geometry).
- The `Pant` effector's default selector becomes `ByName{(lung_part,)}`.
  Panting increases respiratory rate, which increases evaporative flux at the
  lung surface — attributed to `lung_part`.
- Respiratory evaporative heat loss is subtracted from `lung_part`'s energy
  balance, not spread across the body.
- Metabolic O₂ / CO₂ exchange is attributed to `lung_part`.
- For single-`Body` (1-part) organisms, `lung_part = :body` matches the sole
  part name and everything degenerates cleanly.

**Separate from lung site:** the *nasal* passages (in the head) also
contribute a smaller share of heat/water exchange via nasal countercurrent
recovery, and are the anatomical basis for selective brain cooling (see
section 8.4 — selective brain cooling uses nasal-vein evaporative cooling to
chill arterial blood bound for the brain). A `nasal_part::Symbol` field is
deliberately deferred until 8.4 lands; nothing in the current refactor needs
it.

Add `lung_part` in Phase 4 alongside per-part physiology.

### 3.10 Effector + selector interface

Selectors let a behavior target the whole body, a subset, or a single part.
`WholeBody()` degenerates to the current single-body API.

```julia
abstract type PartSelector end
struct WholeBody           <: PartSelector end
struct ByName{Names}       <: PartSelector end
struct Compartment{Name}   <: PartSelector end     # for e.g. Vasodilate / Hyperthermia

# Effectors as small dispatch tags:
abstract type Effector end
struct Piloerect    <: Effector end
struct Uncurl       <: Effector end
struct Vasodilate   <: Effector end
struct Hyperthermia <: Effector end
struct Pant         <: Effector end
struct Sweat        <: Effector end

effect(op, selector, organism, limits) -> (new_limits, new_organism)
apply(op, part, part_limits)           -> (new_part, new_part_limits)   # per-part atom

# Existing scalar API kept:
piloerect(organism, limits) = effect(Piloerect(), WholeBody(), organism, limits)
```

Effector scope by anatomy:
- `Piloerect`, `Uncurl`, `Sweat` → per-part (each part has its own fur / shape /
  sweat glands).
- `Vasodilate`, `Hyperthermia` → per-compartment (perfusion and core
  temperature live at the compartment level, not the part level).
- `Pant` → per-part where anatomically real — routed to `lung_part` by default.

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
      HeatBalanceStrategy hierarchy (IterativeStrategy, NLPStrategy),
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
  `panting_capacity`, `flesh_conductivity_limits` — the same struct that
  owns the fur owns the piloerection range for that fur. Delete the
  parallel per-part limits containers this replaces
  (`InsulationLimits.dorsal` / `.ventral`, etc.) as part of Phase 8.
- `ThermoregulationLimits` keeps only whole-organism scalars: control
  strategy, setpoint / Q10 / minimum heat flow, penalty weights, and
  the four magic-number-replacement fields
  (`skin_temperature_undershoot`, `skin_temperature_core_overshoot`,
  `metabolic_heat_flow_max_multiplier`, `minimum_normalisation_range`).
  Each penalty-weight field accepts either a `Float64` (broadcast to
  every relevant leaf) or a `NamedTuple{PartNames}` (opt-in per-part
  override); the template builder honours whichever form was given.
- Add `OrganismTraits.lung_part::Symbol` (see 3.8). Route `Pant`'s default
  selector through it. Route respiratory evaporative heat loss and O₂/CO₂
  exchange attribution through it. Derive `lung_temperature` from `lung_part`'s
  `core_temperature`. Single-`Body` organisms default to the sole part name.
  Auto-derive each part's `panting_capacity` from `lung_part` (zero on
  non-lung parts).

**Exit:** Existing tests pass. New tests construct a 2-part organism with
distinct per-part physiology *and* per-part insulation limits set from
the parts themselves — no `ThermoregulationLimits` per-part table.
Panting on a multi-part organism affects only the lung part's evaporative
flux and energy balance.

### Phase 5 — Cache layer + `invalidates` traits
**Package:** `BiophysicalBehaviour`
**Work:** New `src/cache.jl`:
- `PrecomputedCache{Parts}` — Tier 1 (shape) + Tier 2 (radiation) `NamedTuple`s
  keyed by part name.
- `invalidates(::Effector)` trait table.
- `refresh!(cache, body, environment, ::Val{keys})` methods; `Val{()}` is a
  compiler no-op.
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
  `ConductiveCoupling{InterfaceConductivity}`, `InsulatedJoin`).
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
- Replace the `isa Sphere` special-case with per-shape
  `aspect_ratio_bounds` dispatch; the builder consults each part's
  own physiology for `insulation_depth_bounds` / `skin_wetness_bounds`
  (no top-level `limits.insulation[part]` indirection).
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
- `ModelParameters.Param` wrapping — not adopted. Use `Flatten` /
  `ConstructionBase` / `Setfield` directly (already in the dep tree).
- The four extensions in section 8 — see below.

## 8. Future extensions (design must not preclude)

Four physical extensions we know we will need eventually. Not built in this
refactor, but the architecture must accommodate them without a redesign. Each
subsection sketches the interface so today's decisions don't paint us into a
corner.

### 8.1 Concentric layered shapes

**Motivation:** Head as brain → cerebrospinal fluid → skull → skin. Torso as
lung → ribcage → muscle → fat → skin. `HomoTherm.R` already handles two
layers (fat + flesh) inside one part; we want arbitrary radial layering.

**Sketch:**

```julia
struct LayeredShape{Shape<:AbstractShape, NumLayers, Boundaries, Conductivities}
        <: AbstractShape
    outer::Shape                            # outermost shape, defines family for dispatch
    boundaries::Boundaries                  # NTuple{NumLayers-1} of inner boundary radii
                                            # (or fractions)
    conductivities::Conductivities          # NTuple{NumLayers} of per-layer conductivity values
end
```

`LayeredShape{Shape}` inherits `Shape`'s family supertype
(`LayeredShape{Cylinder} <: AbstractCylindrical` etc.) via explicit type
parameter constraint, so `HeatExchange`'s family dispatch still works — just
with additional per-shell conduction calculations for the same shape family.
Cylinder, sphere, and ellipsoid all have closed-form radial conduction.

**Impact on today's design:**
- Per-part state schema must be extensible from
  `(; skin_temperature, insulation_temperature)` to
  `(; skin_temperature, insulation_temperature,
      layer_temperatures::NTuple{NumLayers})`. Design `initialise_state` to
  derive shape from the part's shape type, so a `LayeredShape` naturally
  contributes extra state slots.
- Compartment graph nodes are today "one per part core"; must allow
  "`NumLayers` per part" (one per shell). Union-find groups nodes, not parts.
  Design accordingly from Phase 7.

### 8.2 Part-to-part view factors

**Motivation:** A leg partially occludes the belly's view of the sky; the
belly radiates infrared to the leg above the leg's `skin_temperature`.
Ignoring this mis-attributes ~10-20% of infrared exchange in animals with
limbs close to the body or in curled postures.

**Sketch:**

```julia
struct ViewFactors{NumParts, ElementType}
    parts::SMatrix{NumParts, NumParts, ElementType}   # entry [i, j]: fraction of
                                                      # part i's outgoing infrared
                                                      # reaching part j
    sky::SVector{NumParts, ElementType}               # view_factor_to_sky[i]
    ground::SVector{NumParts, ElementType}            # view_factor_to_ground[i]
end
```

Constraints per part: sum of outgoing view factors (to all other parts, to
sky, to ground) equals one. Reciprocity:
`area_i * view_factor_i_to_j == area_j * view_factor_j_to_i`.

Computation: precomputed at organism-construction. Analytical formulas for
simple pairs (disc-disc, parallel-cylinders); ray-cast / hemicube on the
composite geometry for the general case, reusing the rasterization
infrastructure that already backs `silhouette_rasterized`. Belongs in
`BiophysicalGeometry` — same layer as silhouette.

Default until this ships: `parts` matrix is zero; sky + ground sum to one per
part (today's model).

**Impact on today's design:**
- The Tier 2 radiation cache must have room for a per-part
  `incoming_infrared_from_parts` vector or an incoming-radiation pre-sum. Add
  as an optional field now (empty in the no-inter-part-radiation case) so
  consumers can already look for it.
- Introduces cross-part coupling in the outer loop: part `i`'s
  `environmental_heat` depends on part `j`'s `skin_temperature`. Either
  linearize (add off-diagonal terms to the compartment system) or iterate.
  Design the compartment solve to be extensible from block-diagonal
  (independent parts) to full linear system.

### 8.3 Local convection per part

**Motivation:** Wind is not uniform across the body. Leeward parts see less
convection; the dog's legs shelter the belly from horizontal wind. Effect
scales with the wind field's directional coherence and posture.

**Sketch:**

```julia
struct FlowOcclusion{NumParts, ElementType}
    per_part::SVector{NumParts, ElementType}   # multiplier on incident wind speed, 0..1
end
```

Precomputed per (posture, wind direction) — same cadence as posture-dependent
radiation cache. Each part's Nusselt-number calculation uses
`local_wind = ambient_wind * occlusion_factor[part]`. Behavioral effectors
(`TurnToWind`, `FaceHeadIntoWind`) invalidate this cache the same way `Roll`
invalidates the radiation cache.

**Impact on today's design:**
- Tier 2 cache design already accommodates per-part scalars keyed by part
  name. Add `flow_factor` as a Tier 2 entry (default 1.0 per part until
  computed). No structural change; just an extra field.

### 8.4 Selective brain cooling

**Motivation:** Ungulates, canids, felids, some birds. Nasal evaporation
cools venous return blood; carotid rete transfers cool from vein to artery
bound for the brain; brain temperature stays below core during heat stress.

**Sketch:**

Represent brain as either a distinct part (`parts.brain = Body(Sphere(...),
Naked())` inside the head) or as a layer inside a `LayeredShape` head. Couple
to the body's arterial system with a new `HeatCoupling`:

```julia
struct BloodPerfusionCoupling{Effectiveness, EvaporationSource} <: HeatCoupling
    rete_effectiveness::Effectiveness            # 0..1, may be effector-controlled
    evaporation_source::EvaporationSource        # part hosting the cooling evaporation,
                                                 # e.g. :head
end
```

At compartment-solve time, `BloodPerfusionCoupling` contributes a modified
row/column: the brain's `core_temperature` equals arterial core temperature
minus `rete_effectiveness * evaporation_source_heat`. The evaporation-source
heat is the *nasal* evaporative flux from a `nasal_part` (typically the head)
— **not** the pulmonary flux from `lung_part` (section 3.9). Nasal
countercurrent evaporation cools venous return specifically; alveolar
evaporation cools mixed venous return via general circulation and does not
contribute to selective brain cooling.

Adds a companion `OrganismTraits.nasal_part::Symbol` at the same time.

**Impact on today's design:**
- The `HeatCoupling` abstraction in `HeatExchange` (section 3.3) must be open
  to new types — `SharedCore`, `ConductiveCoupling`, `InsulatedJoin` today;
  `BloodPerfusionCoupling` and other perfusion models later. `HeatBalanceSolver`
  dispatches on coupling type when building `conductance_matrix` and
  `heat_load_vector`. Keep the interface
  `contribution_to_conductance_matrix` / `contribution_to_heat_load` for each
  coupling type, so adding new physics is a new method, not a rewrite.
- The `lung_part` tag from section 3.9 is a prerequisite for correct
  respiration attribution. Build it now. `nasal_part` is deferred until this
  extension.

---

**Cross-cutting design principle for section 8:** the state schema,
compartment graph, radiation cache, and coupling interface are all designed to
be *additive* — extensions add new node types, new coupling types, new cache
fields, without rewriting the outer loop. Every phase-8 exit criterion should
be checked against "would this survive adding one of section 8's extensions?"

---

## 9. Design notes worth preserving

- **Cache lives on the organism** (not threaded as a separate argument). Kept
  simple; consistent with the existing `@set organism.body = ...` pattern.
- **Physics annotations (`HeatCoupling`) live in `HeatExchange`**, not on
  `Join` in `BiophysicalGeometry` and not in `BiophysicalBehaviour`. Coupling
  types are heat-physics concepts; `BiophysicalGeometry` stays pure geometry
  (usable by aerodynamic / mechanical consumers); `BiophysicalBehaviour` only
  *stores* couplings in `OrganismTraits.couplings` and forwards them into
  `HeatBalanceProblem`.
- **`HeatExchange` implements `CommonSolve`.** `solve` for one-shot users,
  `init` / `reinit!` / `solve!` for `BiophysicalBehaviour`'s thermoregulation
  control loop, which warm-starts across effector iterations. Same interface
  the SciML ecosystem uses; no bespoke drive API.
- **Iterative and NLP paths share physics.** The `IterativeStrategy`
  (Picard + compartment `SMatrix` solve) and the `IPOPTControl` driver
  (Optimization.jl + OptimizationIpopt.jl, in a `BiophysicalBehaviour`
  extension) both consume the same `solve_part_heat_balance`,
  `contribution_to_conductance_matrix`, and `contribution_to_heat_load`
  primitives in `HeatExchange`. `HeatExchange` never depends on
  `Optimization.jl`.
- **Optimization.jl already implements `CommonSolve`.** No JuMP model to
  persist, no bespoke driver — `IPOPTControl` composes with the SciML
  stack directly.
- **Preserve numerical invariants from the current IPOPT branch:** the
  Float64 / Unitful boundary (now consolidated at the single
  `Flatten.reconstruct` call inside `nlp_residuals!`), the
  `log(metabolic_heat_flow)` transform (a `log_metabolic_heat_flow` leaf
  on the template), the per-effector quadratic-penalty objective
  normalised by `_range` values, and the single whole-organism Q10
  inequality.
- **NLP variable layout uses the DEBTOOL_J template pattern.** A nested
  `NamedTuple` template with `Unitful` values expresses the entire
  decision-variable structure; `Flatten.flatten` / `Flatten.reconstruct`
  (as in `DEBTOOL_J`'s `StateReconstructor`,
  `DEBTOOL_J/src/animals/ode/utils.jl:39`) do all pack / unpack work.
  Adding a part, compartment, or per-part decision variable is a one-line
  change to the template — no generated function, no compile-time
  `VariableLayout{Names, Sizes, Ranges}` type, no per-formulation `_pack`
  / `_unpack` helpers with hardcoded `effectors[i]` indices. Bounds and
  objective weights follow the same template pattern, so a new variable
  adds one field in each parallel table rather than N places across an
  indexed layout.
- **One home per fact — the user authors nothing per-part twice.**
  Per-part information (insulation, insulation piloerection range,
  skin wetness range, panting capacity, flesh-conductivity range) all
  live inside each part's own physiology entry — the same struct that
  owns the fur owns the piloerection range for that fur. Whole-organism
  scalars (setpoint, Q10, penalty weights, the four magic-number
  replacements) stay on `ThermoregulationLimits`, mostly identical to
  today. Every "parallel template" the optimizer needs (bounds,
  weights, references, ranges) is machine-derived by walking parts and
  reading their own physiology plus the whole-organism scalars — the
  user writes none of them. Per-part penalty variation is opt-in by
  swapping a scalar field for a `NamedTuple{PartNames}` on that field.
  For the common case (uniform tuning) the user surface is the same as
  today; adding a part adds no new user surface. See 3.7 "Coefficient
  plumbing" for the full mechanism.
- **Family-abstract dispatch beats mapping tables.** Adding a shape is one
  line in `BiophysicalGeometry`; physics works via `<:` automatically.
- **`ventral_fraction` is a modelling assumption smuggled into a data
  structure.** Replacement uses actual pose + silhouette. Small physics
  upgrade, one real change to the radiation code.
- **DEBTOOL_J's `foldl` over life stages** is the direct analogue for
  `foldl_parts` over body parts. Same immutable-state discipline.
- **Zero-allocation compartment solve is possible** because topology is
  type-parameter-encoded and small —
  `SMatrix{NumCompartments,NumCompartments} \ SVector{NumCompartments}` is
  unrolled and stack-allocated.
- **Spelled-out identifiers throughout.** No `T_`, `Q_`, `k_`, `Ra`, `Pr`,
  `Nu`, `G`, `b`, `IR`, `SBC`, `evap`, `resp`, `geom`, `env`, `pars`. New
  code introduced by this refactor uses fully-spelled names; existing
  abbreviated names are renamed as part of the phase that touches them.
