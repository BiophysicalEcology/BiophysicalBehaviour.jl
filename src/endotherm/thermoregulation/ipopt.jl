# =============================================================================
# IPOPT-based endotherm thermoregulation control — the shared solver machinery.
#
# This file holds everything that is *independent of the NLP formulation*: the
# direct-Ipopt + Enzyme callback plumbing (`IPOPTSolverCache`, `IPOPTParameters`,
# the five `Evaluate*` functors, `LagrangianGradient`, `_lagrangian`, the warm-
# start solve loop `_ipopt_solve!`), and the `thermoregulate(::IPOPTControl)`
# dispatch. It is generic over the packed problem type: every strategy-specific
# quantity — `_objective_value`, `_heat_balance_residuals!`, `_inputs!`, `_scaling`,
# `_problem_size`, `_assemble` — is a method dispatched on that type.
#
# The one concrete strategy, `MultipartNLP` (multipart_nlp.jl), supplies those
# methods; the machinery here never mentions it. A single `Body` is the one-part
# case, so there is no separate single-body formulation.
#
# Execution talks to Ipopt.jl directly — no Optimization.jl / OptimizationIpopt
# wrappers. Going direct exposes `mult_g`, `mult_x_L`, `mult_x_U` on `IpoptProblem`,
# which is what makes dual warm-starting possible.
#
# Two execution paths:
#   1. `thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, init)`
#      — fresh per-call: allocates a single-use cache, runs one cold solve,
#      throws everything away. Use when each call has a structurally
#      different problem (different shape, different environment_vars layout).
#   2. `IPOPTSolverCache(control, organism, environment, init)`
#      then `thermoregulate(...; cache)` — warm-started: callback closures and
#      scratch buffers are reused, and primal+dual from the previous solve
#      are restored into the next `IpoptProblem` before `IpoptSolve` runs
#      with `warm_start_init_point = "yes"`. Use across a sweep where only
#      environment values change.
#
# Boundary between the flat optimiser vector and Unitful physics: the Ipopt
# callbacks operate purely on `Vector{Float64}`; each strategy method attaches
# units at the top (via `Flatten.reconstruct` against a Unitful template) and
# strips them at the bottom (residuals in W/K), so units never leak into the
# solver. See multipart_nlp.jl for that round-trip.
# =============================================================================

# =============================================================================
# Lagrangian — `sigma * f(x) + dot(lambda, g(x))`. IPOPT's `evaluate_hessian` callback
# wants the Hessian of *this*, not of the objective alone. We define it as a
# single scalar function so a forward-over-reverse Enzyme pass (one per
# decision-variable direction) yields one Hessian column at a time.
#
# The residual buffer `residual_buffer` is passed in explicitly rather than allocated
# inside. Under the inner reverse pass we mark it as `Duplicated(residual_buffer, residual_adjoint)`
# so Enzyme writes adjoints into the supplied `residual_adjoint` rather than allocating
# a fresh shadow each call. Under the outer forward pass, both `residual_buffer` and
# `residual_adjoint` are themselves marked `Duplicated` so their forward tangents
# propagate. All four arrays live on the cache; the callbacks just thread
# them through.
#
# `p` is a NamedTuple with fields `sigma::Float64`, `lambda::Vector{Float64}`,
# `ipopt_parameters::IPOPTParameters`. The mutable reference lets the Lagrangian
# see the *current* `objective_parameters` / `nlp_parameters` set by `_inputs!`
# — no per-call bookkeeping inside the Hessian callback.
# =============================================================================

function _lagrangian(nlp_packed, x, residual_buffer, p)
    _heat_balance_residuals!(nlp_packed, residual_buffer, x, p.ipopt_parameters.nlp_parameters)
    # Hand-rolled λ·g rather than `LinearAlgebra.dot`: under the nested
    # forward-over-reverse Enzyme pass, the BLAS `dot` dispatches to Enzyme's
    # fallback BLAS replacement, which segfaults. A plain loop keeps the whole
    # Lagrangian on the differentiable Julia path.
    constraint_term = zero(eltype(residual_buffer))
    @inbounds for i in eachindex(p.lambda, residual_buffer)
        constraint_term += p.lambda[i] * residual_buffer[i]
    end
    return p.sigma * _objective_value(nlp_packed, x, p.ipopt_parameters.objective_parameters) +
           constraint_term
end

# Per-call inputs (bounds, initial values, objective + NLP parameters) are supplied
# by the strategy-specific `_inputs!(::MultipartNLPPacked, …)` method (multipart_nlp.jl):
# it authors the bound/initial templates and the objective's `RegulationTarget` tuple,
# and populates the mutable `IPOPTParameters` container the callbacks read.

# =============================================================================
# Warm-start cache
#
# `IPOPTParameters` is a mutable struct that carries the two NamedTuples the
# Enzyme-AD callbacks consume (objective weights + NLP physics parameters).
# The Ipopt callbacks hold a reference to this container and read its fields on
# every invocation. To set up a new solve we just write fresh NamedTuples into
# the two fields — same concrete types, no reboxing — and the callbacks pick
# them up automatically without rebuilding.
#
# Both fields are parametrically typed so accesses stay type-stable: as long
# as the per-call `objective_parameters` and `nlp_parameters` NamedTuples have the same key/type
# signature (guaranteed by the `_inputs!` helpers above), assignment is a
# bit-for-bit field overwrite. If a structurally different organism or
# environment is fed in (e.g. a different body shape), construct a new cache.
# =============================================================================

mutable struct IPOPTParameters{ObjectiveParametersType,NlpParametersType}
    objective_parameters::ObjectiveParametersType
    nlp_parameters::NlpParametersType
end

# Pre-allocated workspaces consumed by the Enzyme-AD calls inside the Ipopt
# callbacks. One instance per cache, sized at construction from (n_variables, n_constraints).
#
# Jacobian section (used by `evaluate_constraint_jacobian`):
#   jacobian_residual           length n_constraints, residuals primal
#   jacobian_residual_seed          length n_constraints, reverse-mode seed (unit vector for row j)
#   jacobian_x_derivative          length n_variables, gradient output → one Jacobian row
#
# Hessian section (used by `evaluate_hessian`, forward-over-reverse):
#   hessian_gradient         length n_variables, inner-Reverse gradient
#   hessian_gradient_tangent         length n_variables, outer-Forward dual → one Hessian column
#   hessian_seed       length n_variables, outer-Forward seed (unit vector for column j)
#   hessian_multipliers     length n_constraints, current IPOPT λ_g (copied per `evaluate_hessian` call)
#   hessian_residual          length n_constraints, Lagrangian residual buffer
#   hessian_residual_adjoint         length n_constraints, inner-Reverse adjoint of `hessian_residual`
#   hessian_residual_tangent      length n_constraints, outer-Forward tangent of `hessian_residual`
#   hessian_residual_adjoint_tangent     length n_constraints, outer-Forward tangent of `hessian_residual_adjoint`
struct IpoptCallbackBuffers
    jacobian_residual::Vector{Float64}
    jacobian_residual_seed::Vector{Float64}
    jacobian_x_derivative::Vector{Float64}
    hessian_gradient::Vector{Float64}
    hessian_gradient_tangent::Vector{Float64}
    hessian_seed::Vector{Float64}
    hessian_multipliers::Vector{Float64}
    hessian_residual::Vector{Float64}
    hessian_residual_adjoint::Vector{Float64}
    hessian_residual_tangent::Vector{Float64}
    hessian_residual_adjoint_tangent::Vector{Float64}
end

function IpoptCallbackBuffers(n_variables::Int, n_constraints::Int)
    IpoptCallbackBuffers(
        zeros(n_constraints), zeros(n_constraints), zeros(n_variables),               # jacobian buffers
        zeros(n_variables), zeros(n_variables), zeros(n_variables), zeros(n_constraints),     # hessian gradient / tangent / seed / multipliers
        zeros(n_constraints), zeros(n_constraints), zeros(n_constraints), zeros(n_constraints),     # hessian residual + adjoint + their forward tangents
    )
end

"""
    IPOPTSolverCache

    IPOPTSolverCache(control, organism, environment, init)

Reusable solver cache for the IPOPT thermoregulation path. Holds the
Ipopt callback functors, the per-call bound and `initial_values` buffers, the
reverse-mode Jacobian scratch vectors, the Julia-side primal+dual warm-start
vectors, and a mutable `IPOPTParameters` whose fields the functors read on
every solve.

Construct once per organism + smoothing + NLP-strategy combination:

    cache = IPOPTSolverCache(control, organism, environment, init)

where `init` is a NamedTuple of initial guesses (`metabolic_heat_flow`,
`skin_temperature`, `insulation_temperature`).

Then pass the cache as the trailing argument of `thermoregulate` to skip
the closure / scratch-buffer rebuilding and to enable primal+dual warm-
starting from the previous solve. The cache must not be shared between
threads.

Invalidate the cache (build a new one) if the organism's body shape, NLP
strategy, or smoothing strategy changes — the closures bake in those types.
Routine per-iteration changes (environment temperatures, q10, etc.) reuse
the cache without rebuilding.

## Warm-start scope

The cache warm-starts **both** primal and dual:

- *Primal*: the previous solve's effector vector (clamped into the new
  per-call bounds) is fed in as the next solve's initial `x`. The fresh
  path seeds only from the `init` fields; the cache seeds every decision
  variable from the previous solution, a strictly better starting point
  in smooth regimes.
- *Dual*: the Lagrange multipliers from the previous solve (`mult_g` on
  constraints, `mult_x_L` / `mult_x_U` on bounds) are restored into the
  next `IpoptProblem` before `IpoptSolve`, and Ipopt is told to use them
  via `warm_start_init_point=yes`. Bound multipliers for inactive bounds
  carry over harmlessly; for active bounds they accelerate convergence.

We go directly through `Ipopt.jl` rather than `OptimizationIpopt` because
the latter creates a fresh `IpoptProblem` on every `solve` and discards
the final multipliers — no Julia-side hook for them. Bounds still need to
be set at `CreateIpoptProblem` time, so we rebuild the C-side problem each
call with the tight per-call bounds; that allocation is small and Julia-
side state (`previous_primal`, `prev_mult_*`) carries the warm-start info across
rebuilds.
"""
struct IPOPTSolverCache{ParametersType<:IPOPTParameters, ObjectiveCallback, ConstraintsCallback, ObjectiveGradientCallback, ConstraintJacobianCallback, HessianCallback}
    ipopt_parameters::ParametersType
    # Cached callback functors. These read `ipopt_parameters.objective_parameters`
    # and `ipopt_parameters.nlp_parameters` at solve time, so updating those
    # NamedTuples retargets every callback without rebuilding the functors.
    # Stored on the cache so the same callable objects are passed to every
    # CreateIpoptProblem call — that way Julia compiles a single specialised
    # IpoptProblem type once and reuses it.
    evaluate_objective::ObjectiveCallback
    evaluate_constraints::ConstraintsCallback
    evaluate_objective_gradient::ObjectiveGradientCallback
    evaluate_constraint_jacobian::ConstraintJacobianCallback
    evaluate_hessian::HessianCallback
    # Per-call tight bounds — written by `_inputs!` each solve to match the
    # fresh `thermoregulate(..., organism, env, init)` path. Ipopt bakes
    # bounds into the C-side problem at CreateIpoptProblem time, so we rebuild
    # the C handle every solve (cheap — a single malloc + struct init) rather
    # than loosening the bounds to keep one handle alive.
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    initial_values::Vector{Float64}
    constraint_lower_bounds::Vector{Float64}
    constraint_upper_bounds::Vector{Float64}
    n_variables::Int                              # decision variables (from `_problem_size`)
    n_constraints::Int                              # constraints: per-part equalities + whole-organism balance + 1 Q10
    # All AD scratch buffers consumed by the Ipopt callbacks. See the
    # `IpoptCallbackBuffers` definition above.
    buffers::IpoptCallbackBuffers
    # Warm-start state, persisted Julia-side across IpoptProblem rebuilds.
    # On the first solve `has_previous_solve = false` so we feed `initial_values` as the primal and
    # let Ipopt initialise its own multipliers. After IpoptSolve we copy
    # `prob.x`, `prob.mult_g`, `prob.mult_x_L`, `prob.mult_x_U` into these
    # vectors. On subsequent solves we restore them before calling IpoptSolve
    # with `warm_start_init_point = "yes"`.
    previous_primal::Vector{Float64}
    previous_constraint_multipliers::Vector{Float64}
    previous_lower_bound_multipliers::Vector{Float64}
    previous_upper_bound_multipliers::Vector{Float64}
    has_previous_solve::Base.RefValue{Bool}
    # User-provided variable scaling (one entry per decision variable). See
    # the `_scaling` methods for the values. Passed to
    # Ipopt via `SetIpoptProblemScaling` together with `obj_scaling = 1.0`
    # and an all-ones `g_scaling`, after `CreateIpoptProblem` but before
    # `IpoptSolve`. Allocated once at cache construction and reused for
    # every solve.
    x_scaling::Vector{Float64}
    g_scaling::Vector{Float64}
end
function IPOPTSolverCache(control::IPOPTControl, organism::Organism, environment::NamedTuple, init::NamedTuple)
    nlp_packed = HeatExchange.nlp_pack(control.nlp_strategy, organism, environment,
                                       init.skin_temperature, init.insulation_temperature;
                                       smoothing = control.smoothing)
    return IPOPTSolverCache(nlp_packed, organism, environment, init)
end
function IPOPTSolverCache(nlp_packed, organism::Organism, environment::NamedTuple, init::NamedTuple)
    n_variables, n_constraints = _problem_size(nlp_packed)
    lower_bounds   = zeros(n_variables)
    upper_bounds   = zeros(n_variables)
    initial_values = zeros(n_variables)
    objective_parameters, nlp_parameters = _inputs!(lower_bounds, upper_bounds, initial_values,
        nlp_packed, organism, environment, init)
    ipopt_parameters = IPOPTParameters(objective_parameters, nlp_parameters)
    buffers = IpoptCallbackBuffers(n_variables, n_constraints)
    evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_constraint_jacobian, evaluate_hessian =
        _build_ipopt_callbacks(ipopt_parameters, n_variables, n_constraints, buffers)
    constraint_lower_bounds, constraint_upper_bounds = _constraint_bounds(n_constraints)
    return IPOPTSolverCache(ipopt_parameters,
        evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_constraint_jacobian, evaluate_hessian,
        lower_bounds, upper_bounds, initial_values,
        constraint_lower_bounds, constraint_upper_bounds,
        n_variables, n_constraints, buffers,
        zeros(n_variables), zeros(n_constraints), zeros(n_variables), zeros(n_variables), Ref(false),
        _scaling(nlp_packed), ones(n_constraints))
end

# Direct-Ipopt callbacks. Strategy-specific behaviour comes from method
# dispatch on `ipopt_parameters.nlp_parameters.nlp_packed` — passed as `Const`
# first arg to every Enzyme.autodiff call so Julia picks the correct
# `_objective_value`, `_heat_balance_residuals!`, `_lagrangian` method for the
# packed type (`MultipartNLPPacked`). The rest of the wiring (Enzyme plumbing,
# sparsity, lower-triangle indexing) is shared.
#
# Functors read `ipopt_parameters.objective_parameters` /
# `ipopt_parameters.nlp_parameters` at solve time, so refreshing the parameters
# container retargets every callback without rebuilding them.
#
# Jacobian is treated as fully dense — for n_variables=9/11 sparse bookkeeping costs
# more than the saved entries. Sparsity convention: row-major dense, 1-based,
# entry k = (i-1)*n_variables + j ↔ (rows[k], cols[k]) = (i, j).
#
# Hessian uses forward-over-reverse Enzyme: the inner Reverse functor
# `LagrangianGradient` computes ∇ₓL into `buffers.hessian_gradient`, the outer
# Forward seeds e_j and pulls back column j of H into
# `buffers.hessian_gradient_tangent`. Sparsity is dense lower triangle
# (i ≥ j), entry k = i(i-1)/2 + j ↔ (i, j).
# --- Functor types for the five Ipopt callbacks ------------------------------
#
# Each callback is a `struct` with a call method. The struct fields stand in
# for what closures would capture (`ipopt_parameters`, `buffers`, problem dimensions).
# Two practical benefits over anonymous closures:
#   - the callable's type is named — easier to read in stack traces / `@code_*`.
#   - the captured fields are first-class struct fields with explicit types,
#     so type-stability is obvious by inspection rather than implied by capture.
# `LagrangianGradient` is the inner reverse-mode pass that `EvaluateHessian`
# differentiates forward over; it's its own functor for the same reasons.

struct EvaluateObjective{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
(f::EvaluateObjective)(x) =
    _objective_value(f.ipopt_parameters.nlp_parameters.nlp_packed, x, f.ipopt_parameters.objective_parameters)

struct EvaluateConstraints{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
(f::EvaluateConstraints)(x, g) =
    _heat_balance_residuals!(f.ipopt_parameters.nlp_parameters.nlp_packed, g, x, f.ipopt_parameters.nlp_parameters)

struct EvaluateObjectiveGradient{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
function (f::EvaluateObjectiveGradient)(x, grad_f)
    fill!(grad_f, 0)
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse), _objective_value,
        Enzyme.Active,
        Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
        Enzyme.Duplicated(x, grad_f),
        Enzyme.Const(f.ipopt_parameters.objective_parameters)
    )
    return nothing
end

struct EvaluateConstraintJacobian{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
    buffers::IpoptCallbackBuffers
    n_variables::Int
    n_constraints::Int
end
function (f::EvaluateConstraintJacobian)(x, rows, cols, values)
    if values === nothing
        k = 1
        @inbounds for i in 1:f.n_constraints, j in 1:f.n_variables
            rows[k] = i; cols[k] = j; k += 1
        end
    else
        b = f.buffers
        @inbounds for i in 1:f.n_constraints
            fill!(b.jacobian_residual, 0); fill!(b.jacobian_residual_seed, 0); fill!(b.jacobian_x_derivative, 0); b.jacobian_residual_seed[i] = 1.0
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse), _heat_balance_residuals!,
                Enzyme.Const,
                Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
                Enzyme.Duplicated(b.jacobian_residual, b.jacobian_residual_seed),
                Enzyme.Duplicated(x, b.jacobian_x_derivative),
                Enzyme.Const(f.ipopt_parameters.nlp_parameters)
            )
            for j in 1:f.n_variables
                values[(i-1)*f.n_variables + j] = b.jacobian_x_derivative[j]
            end
        end
    end
    return nothing
end

# Inner reverse-mode pass: Enzyme.Reverse on `_lagrangian`. `x` and the
# residual buffer are both Duplicated so Enzyme writes adjoints into the
# supplied buffers rather than allocating shadows. Dispatched to the right
# `_lagrangian` method via the concrete packed type held on `ipopt_parameters`.
struct LagrangianGradient{P<:IPOPTParameters} <: Function
    ipopt_parameters::P
end
function (f::LagrangianGradient)(g, x, residual_buffer, residual_adjoint, p)
    fill!(g, 0)
    fill!(residual_adjoint, 0)
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse), _lagrangian,
        Enzyme.Active,
        Enzyme.Const(f.ipopt_parameters.nlp_parameters.nlp_packed),
        Enzyme.Duplicated(x, g),
        Enzyme.Duplicated(residual_buffer, residual_adjoint),
        Enzyme.Const(p)
    )
    return nothing
end

struct EvaluateHessian{P<:IPOPTParameters, LG<:LagrangianGradient} <: Function
    ipopt_parameters::P
    buffers::IpoptCallbackBuffers
    n_variables::Int
    lagrangian_gradient::LG
end
function (f::EvaluateHessian)(x, rows, cols, sigma, lambda, values)
    if values === nothing
        k = 1
        @inbounds for i in 1:f.n_variables, j in 1:i
            rows[k] = i; cols[k] = j; k += 1
        end
    else
        b = f.buffers
        # IPOPT passes the current objective factor (sigma) and constraint
        # multipliers; copy lambda into the cache so the Lagrangian functor
        # (via `p.sigma` / `p.lambda`) sees them under nested AD.
        copyto!(b.hessian_multipliers, lambda)
        p = (; sigma, lambda = b.hessian_multipliers, ipopt_parameters = f.ipopt_parameters)
        @inbounds for j in 1:f.n_variables
            fill!(b.hessian_seed, 0); b.hessian_seed[j] = 1.0
            fill!(b.hessian_gradient,    0); fill!(b.hessian_gradient_tangent,     0)
            fill!(b.hessian_residual,    0); fill!(b.hessian_residual_tangent,  0)
            fill!(b.hessian_residual_adjoint,   0); fill!(b.hessian_residual_adjoint_tangent, 0)
            # Outer Forward propagates tangents through everything
            # `lagrangian_gradient` reads/writes — including `residual_buffer`
            # and `residual_adjoint` — so every mutable buffer needs its
            # forward dual here.
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Forward), Enzyme.Const(f.lagrangian_gradient),
                Enzyme.Const,
                Enzyme.Duplicated(b.hessian_gradient, b.hessian_gradient_tangent),
                Enzyme.Duplicated(x, b.hessian_seed),
                Enzyme.Duplicated(b.hessian_residual,  b.hessian_residual_tangent),
                Enzyme.Duplicated(b.hessian_residual_adjoint, b.hessian_residual_adjoint_tangent),
                Enzyme.Const(p)
            )
            for i in j:f.n_variables
                k = (i * (i - 1)) ÷ 2 + j
                values[k] = b.hessian_gradient_tangent[i]
            end
        end
    end
    return nothing
end

function _build_ipopt_callbacks(ipopt_parameters::IPOPTParameters, n_variables::Int, n_constraints::Int,
                                 buffers::IpoptCallbackBuffers)
    return (
        EvaluateObjective(ipopt_parameters),
        EvaluateConstraints(ipopt_parameters),
        EvaluateObjectiveGradient(ipopt_parameters),
        EvaluateConstraintJacobian(ipopt_parameters, buffers, n_variables, n_constraints),
        EvaluateHessian(ipopt_parameters, buffers, n_variables, LagrangianGradient(ipopt_parameters)),
    )
end

# Apply the project-wide Ipopt options to a fresh problem. Called once per
# `CreateIpoptProblem` — i.e. per cached solve — alongside the warm-start
# toggle. `print_level` is set per call because callers can flip `verbose`
# between solves; everything else is constant for the cache's lifetime.
function _apply_ipopt_options!(prob::Ipopt.IpoptProblem, verbose::Bool, warm_start::Bool)
    # `hessian_approximation` defaults to "exact" — we register a real
    # `evaluate_hessian` (Lagrangian Hessian via Enzyme forward-over-reverse), so no
    # override needed.
    Ipopt.AddIpoptNumOption(prob, "acceptable_tol",        1e-3)
    Ipopt.AddIpoptIntOption(prob, "acceptable_iter",       5)
    Ipopt.AddIpoptNumOption(prob, "tol",                   1e-4)
    Ipopt.AddIpoptIntOption(prob, "max_iter",              300)
    Ipopt.AddIpoptIntOption(prob, "print_level",           verbose ? 5 : 0)
    Ipopt.AddIpoptStrOption(prob, "warm_start_init_point", warm_start ? "yes" : "no")
    # Adaptive barrier-parameter strategy — Ipopt picks `mu` per-iteration
    # based on the optimality-error / complementarity ratio instead of the
    # default `monotone` schedule. On smooth NLPs (ours: ellipsoid heat
    # balance with continuous physical relationships) this typically
    # converges in fewer iterations because the barrier collapses faster
    # once the optimiser is in the right basin.
    Ipopt.AddIpoptStrOption(prob, "mu_strategy",           "adaptive")
    # User-provided variable scaling — set the actual factors via
    # `SetIpoptProblemScaling` after CreateIpoptProblem (see the `_scaling`
    # methods for values). Telling Ipopt
    # we'll supply our own scales skips its `gradient-based` heuristic pass
    # (an extra Jacobian eval at the initial point) and gives the solver a
    # well-conditioned internal problem from iteration zero.
    Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method",    "user-scaling")
    return prob
end

# Per-variable Ipopt scaling (`_scaling`) and the variable/constraint counts
# (`_problem_size`) are strategy-specific — both live in multipart_nlp.jl, authored
# and counted off the same variable/residual templates as the bounds, so they cannot
# drift out of alignment with the decision-variable layout.

# Constraint bounds: residuals[1:end-1] == 0 (equalities), residuals[end] >= 0
# (the single Q10 inequality, written last by `_heat_balance_residuals!`).
function _constraint_bounds(n_constraints::Int)
    constraint_lower_bounds = zeros(n_constraints)
    constraint_upper_bounds = zeros(n_constraints)
    constraint_upper_bounds[end] = Inf
    return constraint_lower_bounds, constraint_upper_bounds
end

"""
    reset_warm_start!(cache::IPOPTSolverCache)

Discard the previous solution stored in the cache. The next solve will use
the cold initial guess derived from the `init` NamedTuple. Use this when the
trajectory takes a large discontinuous step (e.g. switching organism,
changing q10 regime) where seeding from the last solution would land the
optimiser in a bad basin.
"""
reset_warm_start!(cache::IPOPTSolverCache) = (cache.has_previous_solve[] = false; cache)

# =============================================================================
# Per-call input refresh.
# Overwrite the cache's tight-bound buffers and the two NamedTuples on
# `ipopt_parameters` in place. Dispatch on the packed type happens inside
# `_inputs!`, so the unused branch is dead-code-eliminated and no Union return
# appears on the hot path.
# =============================================================================

function _refresh_cache!(cache::IPOPTSolverCache,
                   nlp_packed, organism, environment, init::NamedTuple)
    objective_parameters, nlp_parameters = _inputs!(cache.lower_bounds, cache.upper_bounds, cache.initial_values, nlp_packed,
        organism, environment, init)
    cache.ipopt_parameters.objective_parameters = objective_parameters
    cache.ipopt_parameters.nlp_parameters = nlp_parameters
    return cache
end
_refresh_cache!(cache::Nothing, args...) = IPOPTSolverCache(args...)

# Direct-Ipopt solve with primal+dual warm-start.
#
# Each call rebuilds the C-side IpoptProblem so the new tight bounds in
# `cache.lower_bounds` / `cache.upper_bounds` take effect — same bounds the fresh
# `thermoregulate(..., organism, env, init)` path uses, so solutions
# stay bit-identical to that path's first iterate. CreateIpoptProblem is a
# malloc + a few struct inits; well under a percent of an IPOPT iteration's
# cost.
#
# Returns the IpoptProblem so the caller can pull the final `x` (and
# optionally inspect the multipliers via `cache.prev_*`).
function _ipopt_solve!(cache::IPOPTSolverCache, verbose::Bool)
    prob = Ipopt.CreateIpoptProblem(
        cache.n_variables, cache.lower_bounds, cache.upper_bounds,
        cache.n_constraints, cache.constraint_lower_bounds, cache.constraint_upper_bounds,
        cache.n_constraints * cache.n_variables,                       # nele_jac: dense
        (cache.n_variables * (cache.n_variables + 1)) ÷ 2,           # nele_hess: dense lower triangle
        cache.evaluate_objective, cache.evaluate_constraints, cache.evaluate_objective_gradient, cache.evaluate_constraint_jacobian,
        cache.evaluate_hessian,
    )
    _apply_ipopt_options!(prob, verbose, cache.has_previous_solve[])
    # User-scaling pairs with `nlp_scaling_method = "user-scaling"` in the
    # options above. `obj_scaling = 1.0` is fine — the objective is already
    # a sum of dimensionless squared penalties so it's intrinsically O(1).
    Ipopt.SetIpoptProblemScaling(prob, 1.0, cache.x_scaling, cache.g_scaling)

    if cache.has_previous_solve[]
        # Restore the previous primal+dual, with the primal clamped into the
        # new tight bounds. `skin_T_min = air_T − 5` shifts between solves;
        # an out-of-range initial x makes Ipopt spend its first iteration
        # snapping back to the feasible box. Dual multipliers can carry over
        # as-is: inactive bounds have ≈0 multipliers (no information lost),
        # active bounds get the previous shadow value as a good guess.
        @inbounds for i in 1:cache.n_variables
            prob.x[i] = clamp(cache.previous_primal[i], cache.lower_bounds[i], cache.upper_bounds[i])
        end
        copyto!(prob.mult_g,   cache.previous_constraint_multipliers)
        copyto!(prob.mult_x_L, cache.previous_lower_bound_multipliers)
        copyto!(prob.mult_x_U, cache.previous_upper_bound_multipliers)
    else
        # Cold start — seed primal from the inputs-helper defaults, let Ipopt
        # initialise its own multipliers (warm_start_init_point = "no" above).
        @inbounds for i in 1:cache.n_variables
            prob.x[i] = cache.initial_values[i]
        end
    end

    Ipopt.IpoptSolve(prob)

    # Persist primal+dual for the next solve.
    copyto!(cache.previous_primal,       prob.x)
    copyto!(cache.previous_constraint_multipliers,  prob.mult_g)
    copyto!(cache.previous_lower_bound_multipliers, prob.mult_x_L)
    copyto!(cache.previous_upper_bound_multipliers, prob.mult_x_U)
    cache.has_previous_solve[] = true
    return prob
end

# Output assembly (`_assemble`) reconstructs the converged `x` into named per-part
# fields and evaluates the physics once at the solution; it is strategy-specific and
# lives in multipart_nlp.jl.

# =============================================================================
# Main thermoregulate dispatch
# =============================================================================

"""
    thermoregulate(::Endotherm, ::IPOPTControl, organism, environment, init;
                   cache=nothing, verbose=false)

Solve endotherm heat balance as a nonlinear program via IPOPT.

`init` is a NamedTuple of initial guesses with fields `metabolic_heat_flow`,
`skin_temperature`, and `insulation_temperature`.

The NLP formulation is `IPOPTControl.nlp_strategy`, `MultipartNLP()` by default: a
genuine multi-part organism with one regulated core shared across all parts, per-part
surface temperatures and effectors, and a single whole-organism respiration balance
closed at the lung. A single `Body` is solved as the one-part case; a `CompositeBody`
adds two heat-balance equality constraints per part. Geometry effectors (piloerection,
uncurling) are deferred, matching the multi-part rule-based ladder.

The objective penalises deviation from the setpoint core temperature and enforces
`metabolic_heat_flow >= heat_flow_min * Q10^((core_temperature − setpoint)/10)`.
Objective weights in `ThermoregulationLimits`:
  - `core_temperature_weight`  — lower → core_temperature rises sooner before effectors are exhausted
  - `metabolic_heat_weight`   — small regularisation (default 0.1) preventing high-panting/high-metabolic_heat_flow
                                  degeneracy at cold temperatures; overridden at hot temperatures by Q10 constraint
  - `gradient_weight`         — non-zero → penalise deviation from `target_core_skin_gradient`
  - `panting_weight`          — lower → panting activates sooner
  - `skin_wetness_weight`     — higher than `panting_weight` → panting activates before sweating

# Arguments
- `organism`: must have `OrganismTraits` with `ThermoregulationLimits`
- `environment`: NamedTuple with `environment_pars` and `environment_vars`
- `init`: NamedTuple of initial guesses — `metabolic_heat_flow` (W),
  `skin_temperature` (K), `insulation_temperature` (K). Use
  `initial_physiological_state(organism, environment.environment_vars)` for defaults.
- `cache` (keyword, optional): an `IPOPTSolverCache` produced by
  `IPOPTSolverCache(control, …)`. Reuse across calls in a sweep to skip per-solve
  allocation and enable primal+dual warm-starting. Must match the organism's body-
  shape type, NLP strategy and smoothing strategy — invalidate the cache and rebuild
  if any of those change.
"""
function thermoregulate(
    ::Endotherm,
    control::IPOPTControl,
    organism::Organism,
    environment::NamedTuple,
    init::NamedTuple;
    cache::Union{Nothing,IPOPTSolverCache} = nothing,
    verbose::Bool = false,
)
    nlp_packed = HeatExchange.nlp_pack(
        control.nlp_strategy, organism, environment, init.skin_temperature, init.insulation_temperature;
        smoothing = control.smoothing
    )
    # Either branch leaves the cache fully populated for `_ipopt_solve!`, 
    # with `_inputs!` called exactly once.
    cache = _refresh_cache!(cache, nlp_packed, organism, environment, init)
    prob = _ipopt_solve!(cache, verbose)
    return _assemble(nlp_packed, organism, environment, prob.x)
end
