using BiophysicalBehaviour
using HeatExchange: HardBound, SmoothBound
using ComponentArrays
using Unitful
using Test

@testset "ThresholdController: HardBound step vs SmoothBound graded transition" begin
    hard = ThresholdController(metric=RawSignal(:x), bound=FixedBound(5.0), comparison=AboveBound())
    hs = initial_controller_state(hard)
    @test controller_level(hard, hs, 0.0, (; x=3.0), hard, hs) == 0.0
    @test controller_level(hard, hs, 0.0, (; x=7.0), hard, hs) == 1.0
    @test register_callback(hard)

    soft = ThresholdController(metric=RawSignal(:x), bound=FixedBound(5.0), comparison=AboveBound(), smoothing=SmoothBound(1.0e-3))
    ss = initial_controller_state(soft)
    levels = [controller_level(soft, ss, 0.0, (; x=v), soft, ss) for v in 3.0:0.5:7.0]
    @test levels[1] < 0.1 && levels[end] > 0.9 && issorted(levels)
    @test levels[5] ≈ 0.5  # exactly at the bound (x=5.0)
    @test !register_callback(soft)
end

@testset "SmoothBound with a Unitful metric/bound and matching scale" begin
    c = ThresholdController(metric=Accumulate((p, sig, m, as) -> 1.0, 0.0u"hr"), bound=FixedBound(10.0u"hr"),
                             comparison=BelowBound(), smoothing=SmoothBound(1.0e-1), scale=1.0u"hr")
    s = initial_controller_state(c)
    levels = Float64[]
    for _ in 1:20
        push!(levels, controller_level(c, s, 0.0, (;), c, s))
        r = controller_rate(c, s, 0.0, (;), c, s)
        s = step_state(s, r, 1.0u"hr")
    end
    @test levels[1] > 0.9 && levels[end] < 0.1 && issorted(levels; rev=true)
end

@testset "Direction gating (RisingDirection/FallingDirection)" begin
    c = ThresholdController(metric=RawSignal(:x), bound=FixedBound(5.0), comparison=BelowBound(), direction=FallingDirection())
    s = initial_controller_state(c)
    @test controller_level(c, s, 0.0, (; x=(value=3.0, rate=-1.0)), c, s) == 1.0
    @test controller_level(c, s, 0.0, (; x=(value=3.0, rate=1.0)), c, s) == 0.0
    @test controller_level(c, s, 0.0, (; x=(value=7.0, rate=-1.0)), c, s) == 0.0
end

@testset "RawProgress has no rate: RisingDirection/FallingDirection errors clearly" begin
    c = ThresholdController(metric=RawProgress(), bound=FixedBound(0.5), comparison=AboveBound(), direction=RisingDirection())
    s = initial_controller_state(c)
    @test_throws "RawProgress has no notion of rate" controller_level(c, s, 0.6, (;), c, s)
end

@testset "Window via AllController (two ThresholdControllers on RawProgress)" begin
    window = AllController(
        lo=ThresholdController(metric=RawProgress(), bound=FixedBound(0.25), comparison=AboveBound()),
        hi=ThresholdController(metric=RawProgress(), bound=FixedBound(0.30), comparison=BelowBound()),
    )
    s = initial_controller_state(window)
    @test controller_level(window, s, 0.20, (;), window, s) == 0.0
    @test controller_level(window, s, 0.27, (;), window, s) == 1.0
    @test controller_level(window, s, 0.35, (;), window, s) == 0.0
end

@testset "Accumulate + FixedBound (cumulation), BelowBound comparison" begin
    c = ThresholdController(metric=Accumulate((p, sig, m, as) -> 1.0, 0.0u"hr"), bound=FixedBound(10.0u"hr"), comparison=BelowBound())
    s = initial_controller_state(c)
    for _ in 1:5
        r = controller_rate(c, s, 0.0, (;), c, s)
        s = step_state(s, r, 1.0u"hr")
    end
    @test s.metric.accumulator == 5.0u"hr"
    @test controller_level(c, s, 0.0, (;), c, s) == 1.0
    for _ in 1:10
        r = controller_rate(c, s, 0.0, (;), c, s)
        s = step_state(s, r, 1.0u"hr")
    end
    @test controller_level(c, s, 0.0, (;), c, s) == 0.0
end

@testset "FunctionController: fires exactly at closure zero-crossing" begin
    c = FunctionController(condition=(s, progress, signals, m, as) -> signals.x - 5.0)
    s = initial_controller_state(c)
    @test controller_level(c, s, 0.0, (; x=3.0), c, s) == 0.0
    @test controller_level(c, s, 0.0, (; x=5.0), c, s) == 1.0
    @test controller_level(c, s, 0.0, (; x=7.0), c, s) == 1.0
end

@testset "Dual-accumulator: second stays at zero rate until first crosses its threshold" begin
    first_ctrl = ThresholdController(metric=Accumulate((p, sig, m, as) -> 1.0, 0.0), bound=FixedBound(5.0), comparison=AboveBound())
    second_rate = (p, sig, model, arrest_state) -> controller_level(model.controllers.first, arrest_state.first, p, sig, model, arrest_state) >= 1.0 ? 1.0 : 0.0
    second_ctrl = ThresholdController(metric=Accumulate(second_rate, 0.0), bound=FixedBound(3.0), comparison=AboveBound())
    dual = AllController(first=first_ctrl, second=second_ctrl)
    s = initial_controller_state(dual)
    hist = Float64[]
    for _ in 1:10
        push!(hist, s.second.metric.accumulator)
        r = controller_rate(dual, s, 0.0, (;), dual, s)
        s = step_state(s, r, 1.0)
    end
    @test hist[1:5] == zeros(5)
    @test hist[end] > 0.0
end

@testset "Gradual breaking: SmoothBound Accumulate gives a time-continuous transition" begin
    temp_rate = (progress, signals, m, as) -> ustrip(u"K", signals.temperature) > 285.0 ? 0.05 : 0.0
    c = ThresholdController(metric=Accumulate(temp_rate, 0.0), bound=FixedBound(1.0), comparison=AboveBound(), smoothing=SmoothBound(0.05))
    s = initial_controller_state(c)
    levels = Float64[]
    for t in 1:60
        temp = 280.0u"K" + t * 0.3u"K"
        push!(levels, controller_level(c, s, 0.0, (; temperature=temp), c, s))
        r = controller_rate(c, s, 0.0, (; temperature=temp), c, s)
        s = step_state(s, r, 1.0)
    end
    @test count(l -> 0.05 < l < 0.95, levels) >= 2
    @test issorted(levels)
end

@testset "Chortoicetes equivalence (discrete-loop, no callbacks)" begin
    # egg_model/params/chortoicetes.jl values, hours throughout for simplicity.
    cold_temp_K, diapause_lo, diapause_hi = 285.65, 0.45, 0.50
    quiescence_windows = ((0.25, 0.30), (0.45, 0.50))
    cold_hour_threshold, diapause_hour_threshold, desiccation_tolerance = 720.0, 1240.0, 0.6

    original_in_diapause(progress, chill, dur) =
        diapause_lo < progress < diapause_hi && chill <= cold_hour_threshold && dur <= diapause_hour_threshold
    original_quiescence(progress, hydration) =
        any(w -> w[1] < progress < w[2], quiescence_windows) && hydration < desiccation_tolerance
    original_arrested(progress, chill, dur, hydration) = original_in_diapause(progress, chill, dur) || original_quiescence(progress, hydration)

    duration_rate(progress, signals, model, arrest_state) =
        (diapause_lo < progress < diapause_hi &&
         arrest_state.induction.chill.metric.accumulator <= cold_hour_threshold &&
         arrest_state.induction.duration.metric.accumulator <= diapause_hour_threshold) ? 1.0 : 0.0
    diapause_induction = AllController(
        lo=ThresholdController(metric=RawProgress(), bound=FixedBound(diapause_lo), comparison=AboveBound()),
        hi=ThresholdController(metric=RawProgress(), bound=FixedBound(diapause_hi), comparison=BelowBound()),
        chill=ThresholdController(metric=Accumulate((p, sig, m, as) -> ustrip(u"K", sig.temperature) < cold_temp_K ? 1.0 : 0.0, 0.0),
                                   bound=FixedBound(cold_hour_threshold), comparison=BelowBound()),
        duration=ThresholdController(metric=Accumulate(duration_rate, 0.0), bound=FixedBound(diapause_hour_threshold), comparison=BelowBound()),
    )
    diapause_model = ComposedArrest(induction=diapause_induction)

    window_controller(lo, hi) = AllController(
        lo=ThresholdController(metric=RawProgress(), bound=FixedBound(lo), comparison=AboveBound()),
        hi=ThresholdController(metric=RawProgress(), bound=FixedBound(hi), comparison=BelowBound()),
        wet=ThresholdController(metric=RawSignal(:hydration), bound=FixedBound(desiccation_tolerance), comparison=BelowBound()),
    )
    quiescence_model = ComposedArrest(induction=AnyController(w1=window_controller(quiescence_windows[1]...), w2=window_controller(quiescence_windows[2]...)))

    chortoicetes_model = AnyArrestModel(diapause=diapause_model, quiescence=quiescence_model)

    n_hours = 2000
    progress_series = range(0.0, 1.0; length=n_hours)
    temp_series = [284.0 + 4.0 * sin(2pi * t / 400) for t in 1:n_hours]
    hydration_series = [t < 1000 ? 0.8 : 0.3 for t in 1:n_hours]

    orig_chill, orig_dur = 0.0, 0.0
    state = initial_arrest_state(chortoicetes_model)
    mismatches = 0
    for t in 1:n_hours
        progress = progress_series[t]
        signals = (; temperature=temp_series[t] * u"K", hydration=hydration_series[t])

        orig = original_arrested(progress, orig_chill, orig_dur, hydration_series[t])
        new = arrest_level(chortoicetes_model, state, progress, signals) > 0.5
        mismatches += (orig != new)

        orig_chill += temp_series[t] < cold_temp_K ? 1.0 : 0.0
        orig_dur += original_in_diapause(progress, orig_chill, orig_dur) ? 1.0 : 0.0

        rate = advance_arrest(chortoicetes_model, state, progress, signals)
        state = step_state(state, rate, 1.0)
    end
    @test mismatches == 0
end

@testset "print_arrest_structure runs without error" begin
    c = ThresholdController(metric=RawSignal(:x), bound=FixedBound(5.0), comparison=AboveBound())
    m = ComposedArrest(induction=c)
    output = mktemp() do path, io
        redirect_stdout(io) do
            print_arrest_structure(m, "test")
        end
        close(io)
        read(path, String)
    end
    @test occursin("test", output)
end
