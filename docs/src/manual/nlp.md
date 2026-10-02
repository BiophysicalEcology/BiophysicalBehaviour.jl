# How the optimisation is built

[Thermoregulation by optimisation](optimisation.md) describes the problem [`IPOPTControl`](@ref) solves. This
page describes how it is put together: how the heat budget of HeatExchange.jl becomes the constraints of a
nonlinear program, where its derivatives come from, and what each design choice is for. It is the counterpart
of [Differentiability and the NLP interface](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/autodiff)
in the documentation of HeatExchange.jl, which describes the same boundary from the other side.

```@setup nlp
using Main.FigureHelpers
using CairoMakie
```

## What an interior-point solver needs

IPOPT minimises an objective ``f(x)`` subject to constraints ``g_L \le g(x) \le g_U`` and bounds
``x_L \le x \le x_U``. At each iteration it asks for:

| Callback | Returns | Here |
|:--|:--|:--|
| objective | ``f(x)`` | the weighted sum of squared departures from the targets |
| constraints | ``g(x)`` | the residuals of the heat budget |
| objective gradient | ``\nabla f`` | Enzyme, reverse mode |
| constraint Jacobian | ``\partial g_i / \partial x_j`` | Enzyme, reverse mode, one row for each constraint |
| Hessian of the Lagrangian | ``\nabla^2 (\sigma f + \lambda \cdot g)`` | Enzyme, forward mode over reverse mode, one column for each variable |

The Lagrangian combines objective and constraints, each constraint weighted by a multiplier ``\lambda_i``. The
multipliers are the dual variables. Each measures how much the objective would improve if its constraint were
relaxed.

## The physics as residuals

`solve_metabolic_rate` finds skin and surface temperatures by iteration, inside the function. An optimiser
cannot use it as a constraint: the temperatures would be found for the organism as the optimiser's present guess
describes it, with no derivative through the iteration worth having.

So the temperatures are made variables, and the physics is evaluated *without iteration* at whatever values the
optimiser proposes. `part_surface_residuals` of HeatExchange.jl takes the core, skin and surface temperatures
and the metabolic rate of one part and returns how far its heat budget is from balance. It is the same function
the iterative solver drives to zero, so the two cannot disagree about the physics.

```@example nlp
solver_paths_diagram() # hide
```

The residual function of the problem, `_heat_balance_residuals!`, walks the parts, calls
`part_surface_residuals` for each, closes the respiration balance once for the whole organism, and adds the
``Q_{10}`` inequality: two constraints for each part, then the whole-organism balance, then the inequality.

## The extension point

HeatExchange.jl knows nothing of IPOPT. It exports an abstract type and one function:

```julia
abstract type NLPStrategy end
nlp_pack(strategy, organism, environment, skin_temperature, insulation_temperature; smoothing)
```

This package defines the one concrete strategy, [`MultipartNLP`](@ref), and its method of `nlp_pack`. Packing
does once, before the solve, everything that does not depend on the variables: the geometry of each part, its
view factors, the properties of the air, the setup of each part's surface. The callbacks close over the
`MultipartNLPPacked` it returns.

| Owned by HeatExchange.jl | Owned by BiophysicalBehaviour.jl |
|:--|:--|
| the residuals of one part | which quantities are variables |
| smoothing of the kinks in the physics | their bounds and first guesses |
| the abstract `NLPStrategy` and `nlp_pack` | the objective and its weights |
| | scaling, the solver and its options |
| | assembling the output |

## One structure, flattened many ways

IPOPT works on a flat vector of `Float64`. The problem is naturally nested: some variables for the whole
organism, and a set for each named part. Indexing a flat vector by hand, `x[3 + 4i]`, is how a bound comes to be
applied to the wrong variable.

Instead the layout is defined once, as a nested NamedTuple:

```julia
(; core_temperature, log_metabolic_heat_flow, panting_rate,
   parts = (; torso = (; skin_temperature, insulation_temperature, flesh_conductivity, skin_wetness),
              head  = (; ...), ...))
```

and [Flatten.jl](https://github.com/rafaqz/Flatten.jl) converts between it and the flat vector. Variables, lower
bounds, upper bounds, first guesses and scale factors are all structures of this one shape, flattened in the same
order, so they are aligned by construction. Adding a part, or a variable to each part, changes the structure,
and every flat vector follows. The numbers of variables and constraints, ``3 + 4N`` and ``2N + 2`` for ``N``
parts, are counted from the structures the same way.

## Units at the boundary

Every quantity in the physics carries units. The optimiser's vector cannot. The boundary is the same structure:
its leaves are `Unitful` quantities taken from the organism's own traits, so the units live in the *types* of
the leaves. Flattening extracts the numbers. Reconstructing puts each back into its typed slot as a quantity.
No unit is stripped or attached by hand, and none is a literal that could drift from the units the model uses.

The residuals come back the other way: each is divided by one unit of itself at the last line of the residual
function. Between those two points everything is dimensioned, and a mistake of dimension in the physics is
still an error, see
[Units, dimensions and functional traits](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/units_traits#Where-units-are-stripped-in-this-package).

Metabolic heat production is carried as its logarithm, which keeps it positive and brings a variable that
ranges over an order of magnitude to the scale of the others.

## Smoothing

The physics has kinks: an `abs` in free convection, a `max` where a flow cannot be negative, a step where fur
is or is not present. Each is exact for a forward solve and gives an undefined or misleading derivative at the
kink. HeatExchange.jl passes a `SmoothingStrategy` through every such call. The rule-based controller uses
`HardBound()`, the exact functions. [`IPOPTControl`](@ref) carries a `SmoothBound(1e-5)`, which replaces each
kink with a smooth function that differs from it only within a narrow band:

```@example nlp
smoothing_figure() # hide
```

The band is narrow enough that the two controllers solve the same physics within their tolerances.

## Derivatives through three packages

[Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) differentiates compiled Julia code. It follows the
calculation from the residual function here, through the heat budget in HeatExchange.jl, into the properties of
air in FluidProperties.jl and the geometry of BiophysicalGeometry.jl, with no derivative written by hand. What
those packages owe it is code that is type-stable and does not allocate in the path differentiated, which is why
the residual function uses plain loops over tuples and writes into a buffer passed in.

The Hessian is exact. For each variable, one forward-mode pass over the reverse-mode gradient of the Lagrangian
gives one column. With exact second derivatives IPOPT takes full Newton steps and converges in tens of
iterations. The Lagrangian is a plain loop, not a call to `dot`, because nested differentiation of the BLAS
routine fails.

The working arrays of the callbacks are allocated once, in `IpoptCallbackBuffers`, and the callbacks are functor
structs, not closures, so their types are concrete.

## Scaling

Temperatures are near 300, skin wetness near 0.01. IPOPT is told the scale of each variable, through the same
structure flattened once more, instead of estimating scales from the gradient at the first point: 1/300 for
temperatures, 1 for the logarithm of metabolic rate, flesh conductivity and panting, 50 for skin wetness.

## Talking to IPOPT directly

The package calls [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) itself, not through a modelling layer, for the
warm start: through the direct interface the multipliers of the constraints and bounds at the solution can be
read and, for the next solve, written back.

An [`IPOPTSolverCache`](@ref) holds the callbacks, the buffers and the previous solution, primal and dual:

```julia
cache = IPOPTSolverCache(control, organism, environment, init)
out = thermoregulate(Endotherm(), control, organism, environment, init; cache)
```

Each call rebuilds the small C-side problem, since IPOPT fixes the bounds when a problem is created, and
restores the previous variables and multipliers into it. The callbacks read their parameters from a mutable
container in the cache, so a new environment is a matter of writing new values into it.

A cache is specific to the types it was built with: a new one is needed if the shape of the body, the number of
parts or the smoothing changes. `BiophysicalBehaviour.reset_warm_start!(cache)` discards the stored solution
after a large jump in conditions. A cache must not be shared between threads.

## Options

| Option | Value | Why |
|:--|:--|:--|
| `tol` | 1e-4 | |
| `acceptable_tol`, `acceptable_iter` | 1e-3, 5 | stop when progress stalls close to a solution |
| `max_iter` | 300 | |
| `mu_strategy` | `adaptive` | fewer iterations on a smooth problem than the monotone default |
| `nlp_scaling_method` | `user-scaling` | the scales above |
| `warm_start_init_point` | `yes` when the cache holds a solution | |

Passing `verbose = true` to `thermoregulate` prints IPOPT's iteration log.

## What is not yet done

- Fur depth and posture are not variables. They rebuild the body, and with it the setup of each part, inside the
  function to be differentiated.
- The status with which IPOPT stopped is not returned.
- The Jacobian and Hessian are passed as dense. With many parts they are sparse, since a part's residuals do not
  depend on another part's surface.
- The scale factors are constants, not derived from the bounds.
