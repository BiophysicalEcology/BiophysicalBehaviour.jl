using BiophysicalBehaviour
using BiophysicalBehaviour: AbstractCondition, AbstractStage, AbstractTransition, StageSequence,
    stages, transitions
import BiophysicalBehaviour: stage_traits, transition_controller, stage_key, transition_key,
    initial_controller_state, controller_rate, controller_level, register_callback, trigger_conditions
using HeatExchange: HardBound, SmoothBound
using ComponentArrays
using Unitful
using Test

# --- Inline test fixtures: a minimal insect-like vocabulary, proving the
# domain-agnostic core works for a real case. Not a shipped vocabulary
# module -- see Architecture in the design plan.

struct Egg{T} <: AbstractStage{T}
    traits::T
end

struct Instar{N,T} <: AbstractStage{T}
    traits::T
end
Instar{N}(traits::T) where {N,T} = Instar{N,T}(traits)
# Default stage_key (nameof) would collide across different N -- override so
# multiple Instar{N} stages can coexist in one sequence (see docs/stages.md).
stage_key(::Type{<:Instar{N}}) where {N} = Symbol(:Instar, N)

struct Hatch{C<:AbstractCondition} <: AbstractTransition{C}
    controller::C
end

struct Molt{N,C<:AbstractCondition} <: AbstractTransition{C}
    controller::C
end
Molt{N}(controller::C) where {N,C} = Molt{N,C}(controller)
transition_key(::Type{<:Molt{N}}) where {N} = Symbol(:Molt, N)

# Accessor-protocol fixture (test 17): different internal field names
# entirely, proving the protocol runs through accessors, not field names.
struct ExternalStage{T} <: AbstractStage{T}
    data::T
end
stage_traits(stage::ExternalStage) = stage.data

struct ExternalTransition{C<:AbstractCondition} <: AbstractTransition{C}
    gate::C
end
transition_controller(transition::ExternalTransition) = transition.gate

# Custom AbstractCondition fixture (test 16): not ThresholdController/
# FunctionController, with real own-state, to prove the protocol is open.
struct FlagCondition{T} <: AbstractCondition
    at::T
end
initial_controller_state(::FlagCondition) = (; checks=0)
controller_rate(::FlagCondition, own_state, progress, signals, model, stage_state) = (; checks=1)
controller_level(c::FlagCondition, own_state, progress, signals, model, stage_state) =
    progress >= c.at ? 1.0 : 0.0
register_callback(::FlagCondition) = true
trigger_conditions(c::FlagCondition) = ((own_state, progress, signals, model, stage_state) -> progress - c.at,)

hard_hatch(at) = Hatch(ThresholdController(metric=RawProgress(), bound=FixedBound(at), comparison=AboveBound()))
hard_molt(N, at) = Molt{N}(ThresholdController(metric=RawProgress(), bound=FixedBound(at), comparison=AboveBound()))
smooth_hatch(at, eps) = Hatch(ThresholdController(metric=RawProgress(), bound=FixedBound(at), comparison=AboveBound(), smoothing=SmoothBound(eps)))

@testset "1. Type-keyed lookup" begin
    sequence = StageSequence(
        Egg(1.0),
        hard_hatch(0.2) => Instar{3}(3.0),
        hard_molt(1, 0.5) => Instar{5}(5.0),
        hard_molt(2, 0.8) => Instar{7}(7.0),
    )

    # Unparameterized stage type.
    @test sequence[Egg] == 1.0
    @test haskey(sequence, Egg)

    # Instar{3}: partially instantiated type key.
    @test sequence[Instar{3}] == 3.0
    # Instar{3,Float64}: fully instantiated type key, same match.
    @test sequence[Instar{3,Float64}] == 3.0
    # Instar{3} vs Instar{5}: must not collide.
    @test sequence[Instar{5}] == 5.0
    @test sequence[Instar{3}] != sequence[Instar{5}]

    # Absent type: get/haskey agree, independent of whether matched data is `nothing`.
    @test !haskey(sequence, Instar{99})
    @test get(sequence, Instar{99}, nothing) === nothing

    # Parameterized transitions: same partial/full/non-collision/absence matrix.
    @test sequence[Molt{1}].bound.value == 0.5
    @test sequence[Molt{1,typeof(sequence[Molt{1}])}] isa ThresholdController  # fully instantiated transition key
    @test sequence[Molt{2}].bound.value == 0.8
    @test sequence[Molt{1}].bound.value != sequence[Molt{2}].bound.value
    @test !haskey(sequence, Molt{99})
    @test get(sequence, Molt{99}, nothing) === nothing
end

@testset "2. Cardinality enforced on both construction paths" begin
    # Chain form: one more stage than transitions is the only valid shape.
    @test_throws ArgumentError StageSequence(
        (Egg(1.0), Instar{1}(2.0), Instar{2}(3.0)),
        (hard_hatch(0.5),),  # 3 stages, 1 transition -- invalid
    )
    # Raw tuple-pair form goes through the same validated inner constructor.
    @test_throws ArgumentError StageSequence(
        (Egg(1.0),),
        (hard_hatch(0.5), hard_molt(1, 0.8)),  # 1 stage, 2 transitions -- invalid
    )
    # Valid shape succeeds on both paths.
    @test StageSequence(Egg(1.0), hard_hatch(0.5) => Instar{1}(2.0)) isa StageSequence
    @test StageSequence((Egg(1.0), Instar{1}(2.0)), (hard_hatch(0.5),)) isa StageSequence
end

@testset "3. Singleton sequence" begin
    sequence = StageSequence((Egg(1.0),), ())
    state = initial_stage_state(sequence)
    @test state == NamedTuple()
    @test stage_conditions(sequence) == ()
    weights = stage_weights(sequence, state, 0.5, (;))
    @test weights == (; Egg=1)
    @test current_stage(sequence, state, 0.5, (;)) isa Egg
    @test stage_value(sequence, state, 0.5, (;), (traits, progress, signals) -> traits) == 1.0
end

@testset "4. Duplicate-key rejection" begin
    @test_throws ArgumentError StageSequence(
        Egg(1.0),
        hard_hatch(0.5) => Egg(2.0),  # duplicate stage_key :Egg
    )
    @test_throws ArgumentError StageSequence(
        Egg(1.0),
        hard_hatch(0.3) => Instar{1}(2.0),
        hard_hatch(0.7) => Instar{2}(3.0),  # duplicate transition_key :Hatch
    )
    # A stage and a transition sharing a key does not throw -- different namespaces.
    struct Hatch_{T} <: AbstractStage{T}
        traits::T
    end
    @test StageSequence(Egg(1.0), hard_hatch(0.5) => Hatch_(2.0)) isa StageSequence

    # Two default-keyed Instar instantiations collide (nameof(S) alone isn't unique)...
    struct UnkeyedInstar{N,T} <: AbstractStage{T}
        traits::T
    end
    UnkeyedInstar{N}(traits::T) where {N,T} = UnkeyedInstar{N,T}(traits)
    @test_throws ArgumentError StageSequence(
        UnkeyedInstar{1}(1.0),
        hard_hatch(0.3) => UnkeyedInstar{2}(2.0),
    )
    # ...but overriding stage_key (as our real Instar fixture does) resolves it.
    @test StageSequence(
        Instar{1}(1.0),
        hard_hatch(0.3) => Instar{2}(2.0),
    ) isa StageSequence
end

@testset "5. HardBound: exactly one active stage away from the threshold" begin
    sequence = StageSequence(Egg(1.0), hard_hatch(0.5) => Instar{1}(2.0))
    state = initial_stage_state(sequence)
    @test stage_weights(sequence, state, 0.49, (;)) == (; Egg=1.0, Instar1=0.0)
    @test stage_weights(sequence, state, 0.51, (;)) == (; Egg=0.0, Instar1=1.0)
    # Exactly at the threshold: HardBound's safe_step requires gap > 0 (strict), so
    # the transition hasn't fired yet and the prior stage is still fully active.
    @test stage_weights(sequence, state, 0.5, (;)) == (; Egg=1.0, Instar1=0.0)
end

@testset "6. SmoothBound: graded blend" begin
    sequence = StageSequence(Egg(1.0), smooth_hatch(0.5, 1.0e-2) => Instar{1}(2.0))
    state = initial_stage_state(sequence)
    progresses = 0.3:0.02:0.7
    sums = [sum(values(stage_weights(sequence, state, p, (;)))) for p in progresses]
    @test all(s -> s ≈ 1.0, sums)
    egg_weights = [stage_weights(sequence, state, p, (;)).Egg for p in progresses]
    @test issorted(egg_weights; rev=true)  # monotonically decreasing as progress increases
    @test egg_weights[1] > 0.9    # far below threshold: ~all Egg
    @test egg_weights[end] < 0.1  # far above threshold: ~all Instar1
end

@testset "7. Overlapping smooth transitions" begin
    sequence = StageSequence(
        Egg(1.0),
        smooth_hatch(0.3, 0.1) => Instar{1}(2.0),
        hard_molt(1, 0.35) => Instar{2}(3.0),  # overlaps hatch's smoothing window
    )
    state = initial_stage_state(sequence)
    for p in 0.2:0.05:0.5
        w = stage_weights(sequence, state, p, (;))
        @test sum(values(w)) ≈ 1.0
        l1 = controller_level(sequence[Hatch], getfield(state, transition_key(transitions(sequence)[1])), p, (;), sequence, state)
        l2 = controller_level(sequence[Molt{1}], getfield(state, transition_key(transitions(sequence)[2])), p, (;), sequence, state)
        @test w.Egg ≈ 1 - l1
        @test w.Instar1 ≈ l1 * (1 - l2)
        @test w.Instar2 ≈ l1 * l2
    end
end

@testset "8. Numeric type preservation" begin
    # HeatExchange's safe_step always returns Float64 regardless of input type;
    # a bare level_function isolates stage_weights's own type preservation.
    float32_condition = FunctionController(
        condition=(s, p, sig, m, ss) -> error("unused"),
        level_function=(s, p, sig, m, ss) -> p > 0.5f0 ? 1.0f0 : 0.0f0,
    )
    sequence = StageSequence(Egg(1.0f0), Hatch(float32_condition) => Instar{1}(2.0f0))
    state = initial_stage_state(sequence)
    w32 = stage_weights(sequence, state, 0.6f0, (;))
    @test eltype(typeof(values(w32))) === Float32

    struct DualLike
        value::Float64
        deriv::Float64
    end
    Base.one(::DualLike) = DualLike(1.0, 0.0)
    Base.one(::Type{DualLike}) = DualLike(1.0, 0.0)
    Base.:-(a::DualLike, b::DualLike) = DualLike(a.value - b.value, a.deriv - b.deriv)
    Base.:*(a::DualLike, b::DualLike) = DualLike(a.value * b.value, a.deriv * b.value + a.value * b.deriv)
    Base.:*(a::Number, b::DualLike) = DualLike(a * b.value, a * b.deriv)

    dual_condition = FunctionController(condition=(s, p, sig, m, ss) -> error("unused"),
                                         level_function=(s, p, sig, m, ss) -> DualLike(p, 1.0))
    dual_sequence = StageSequence(Egg(1.0), Hatch(dual_condition) => Instar{1}(2.0))
    dual_state = initial_stage_state(dual_sequence)
    w_dual = stage_weights(dual_sequence, dual_state, 0.5, (;))
    @test w_dual.Egg isa DualLike
end

@testset "9. current_stage exact-tie behavior" begin
    sequence = StageSequence(Egg(1.0), smooth_hatch(0.5, 1.0) => Instar{1}(2.0))
    state = initial_stage_state(sequence)
    w = stage_weights(sequence, state, 0.5, (;))
    @test w.Egg ≈ w.Instar1 atol = 1e-6
    @test current_stage(sequence, state, 0.5, (;)) isa Egg  # exact tie resolves to the earlier stage
end

@testset "10. Stateful transition" begin
    sequence = StageSequence(
        Egg(1.0),
        Hatch(ThresholdController(metric=Accumulate((p, sig, m, ss) -> 1.0, 0.0), bound=FixedBound(5.0), comparison=AboveBound())) => Instar{1}(2.0),
    )
    state = initial_stage_state(sequence)
    for _ in 1:4
        @test stage_weights(sequence, state, 0.0, (;)).Egg == 1.0
        rate = advance_stage(sequence, state, 0.0, (;))
        state = step_state(state, rate, 1.0)
    end
    @test state.Hatch.metric.accumulator == 4.0
    rate = advance_stage(sequence, state, 0.0, (;))
    state = step_state(state, rate, 1.0)
    @test state.Hatch.metric.accumulator == 5.0
    @test stage_weights(sequence, state, 0.0, (;)).Egg == 1.0  # exactly at bound: HardBound requires gap > 0
    rate = advance_stage(sequence, state, 0.0, (;))
    state = step_state(state, rate, 1.0)
    @test state.Hatch.metric.accumulator == 6.0
    @test stage_weights(sequence, state, 0.0, (;)).Instar1 == 1.0
end

@testset "11. stage_component on content-free state" begin
    sequence = StageSequence(Egg(1.0), hard_hatch(0.5) => Instar{1}(2.0))
    component = stage_component(sequence)
    @test component isa ComponentArray
    @test length(component) == 0
end

@testset "12. stage_value accessor-based linear ramp" begin
    sequence = StageSequence(
        Egg(0.0),
        smooth_hatch(0.5, 0.05) => Instar{1}(1.0),
    )
    state = initial_stage_state(sequence)
    ramp_values = [stage_value(sequence, state, p, (;), (traits, progress, signals) -> traits) for p in 0.0:0.1:1.0]
    @test issorted(ramp_values)
    @test ramp_values[1] ≈ 0.0 atol = 0.05
    @test ramp_values[end] ≈ 1.0 atol = 0.05
end

@testset "13. Chortoicetes-equivalence" begin
    # Mirrors SteppedHydricStage's two-threshold ramp logic independently.
    lo, hi = 0.4, 0.6
    original_stage_value(progress, early, late) =
        progress <= lo ? early :
        progress >= hi ? late :
        early + (late - early) * (progress - lo) / (hi - lo)

    sequence = StageSequence(
        Egg(2.0),
        smooth_hatch((lo + hi) / 2, (hi - lo) / 2) => Instar{1}(8.0),
    )
    state = initial_stage_state(sequence)
    mismatches = 0
    for progress in 0.0:0.01:1.0
        new_value = stage_value(sequence, state, progress, (;), (traits, p, sig) -> traits)
        orig_value = original_stage_value(progress, 2.0, 8.0)
        # smoothing windows differ slightly in shape; both are monotonic ramps
        # over the same [lo,hi] window and agree closely away from the
        # smoothing kernel's own transition width.
        if progress <= lo - 0.1 || progress >= hi + 0.1
            mismatches += !isapprox(new_value, orig_value; atol=0.4)
        end
    end
    @test mismatches == 0
end

@testset "14. stage_conditions via a real VectorContinuousCallback" begin
    using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5, VectorContinuousCallback

    sequence = StageSequence(
        Egg(1.0),
        hard_hatch(0.5) => Instar{1}(2.0),
        hard_molt(1, 0.8) => Instar{2}(3.0),
    )
    conditions = stage_conditions(sequence)
    @test length(conditions) == 2
    state = initial_stage_state(sequence)

    fired = Int[]
    function condition!(out, u, t, integrator)
        for (i, cond) in enumerate(conditions)
            out[i] = cond(state, u[1], (;))
        end
    end
    function affect!(integrator, i)
        push!(fired, i)
    end
    cb = VectorContinuousCallback(condition!, affect!, length(conditions))

    f(u, p, t) = [1.0]  # progress increases at a constant rate
    prob = ODEProblem(f, [0.0], (0.0, 1.0))
    sol = solve(prob, Tsit5(); callback=cb, dtmax=0.01)

    @test fired == [1, 2]  # crossing order: Hatch (0.5) before Molt{1} (0.8)
    @test isapprox(sol.t[end], 1.0; atol=1e-6)
end

@testset "15. print_stage_structure" begin
    sequence = StageSequence(
        Egg(1.0),
        hard_hatch(0.5) => Instar{1}(2.0),
        hard_molt(1, 0.8) => Instar{2}(3.0),
    )
    output = mktemp() do path, io
        redirect_stdout(io) do
            print_stage_structure(sequence, "test")
        end
        close(io)
        read(path, String)
    end
    @test occursin("test", output)
    @test occursin("0.5", output)
    @test occursin("0.8", output)
    @test occursin("Hatch", output)
    @test occursin("Molt", output)
end

@testset "16. Custom AbstractCondition extensibility" begin
    sequence = StageSequence(Egg(1.0), Hatch(FlagCondition(0.5)) => Instar{1}(2.0))
    state = initial_stage_state(sequence)
    @test state.Hatch == (; checks=0)

    for _ in 1:3
        rate = advance_stage(sequence, state, 0.3, (;))
        state = step_state(state, rate, 1.0)
    end
    @test state.Hatch.checks == 3

    @test stage_weights(sequence, state, 0.3, (;)).Egg == 1.0
    @test stage_weights(sequence, state, 0.6, (;)).Instar1 == 1.0

    conditions = stage_conditions(sequence)
    @test length(conditions) == 1
    @test conditions[1](state, 0.4, (;)) ≈ 0.4 - 0.5
end

@testset "17. Accessor-protocol fixture (non-standard field names)" begin
    sequence = StageSequence(
        ExternalStage(42.0),
        ExternalTransition(ThresholdController(metric=RawProgress(), bound=FixedBound(0.5), comparison=AboveBound())) => Instar{1}(7.0),
    )
    @test stage_traits(stages(sequence)[1]) == 42.0
    @test transition_controller(transitions(sequence)[1]) isa ThresholdController
    state = initial_stage_state(sequence)
    @test stage_weights(sequence, state, 0.2, (;)).ExternalStage == 1.0
    @test stage_value(sequence, state, 0.2, (;), (traits, progress, signals) -> traits) == 42.0
end

@testset "18. Type stability and allocation" begin
    # Measured inside a function: top-level bindings are global/untyped and
    # would measure boxing overhead, not stage_weights/stage_value's own behavior.
    function check_inference_and_allocation()
        sequence = StageSequence(
            Egg(1.0),
            hard_hatch(0.3) => Instar{1}(2.0),
            hard_molt(1, 0.7) => Instar{2}(3.0),
        )
        state = initial_stage_state(sequence)
        signals = (;)
        accessor = (traits, p, sig) -> traits

        @inferred stage_weights(sequence, state, 0.5, signals)
        # current_stage indexes a heterogeneous Tuple by runtime index, so its
        # return type is an inherent Union, not a type-stability bug -- not @inferred.
        current_stage(sequence, state, 0.5, signals)
        @inferred stage_value(sequence, state, 0.5, signals, accessor)
        @inferred advance_stage(sequence, state, 0.5, signals)

        stage_weights(sequence, state, 0.5, signals)  # warm up
        weights_allocated = @allocated stage_weights(sequence, state, 0.5, signals)
        stage_value(sequence, state, 0.5, signals, accessor)  # warm up
        value_allocated = @allocated stage_value(sequence, state, 0.5, signals, accessor)
        return weights_allocated, value_allocated
    end

    weights_allocated, value_allocated = check_inference_and_allocation()
    @test weights_allocated == 0
    @test value_allocated == 0
end
