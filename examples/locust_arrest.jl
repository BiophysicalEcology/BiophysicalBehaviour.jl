# Chortoicetes terminifera (Australian plague locust) egg dormancy: a demonstration of
# src/arrest/'s diapause and quiescence pathways (see docs/arrest.md). These pathways
# may occur in parallel in the locust.
# Diapause is induced within a single progress window, only while accumulated chilling and 
# diapause duration stay under their own limits. 
# Quiescence is a desiccation-triggered pause during two separate progress windows, the second
# overlapping with the diapause window. 
# The pathways are combined via AnyArrestModel so that either arrests the egg).
# Chill/duration accumulate in `u"hr"`, so SmoothBound's `scale` carries matching units;
# progress/hydration are dimensionless, so their controllers use the default scale=1.0.
using BiophysicalBehaviour
using HeatExchange
using Unitful

# --- Chortoicetes-specific parameters ---
const COLD_TEMPERATURE = u"K"(12.5u"°C")
const DIAPAUSE_WINDOW = (0.45, 0.50)                # development-progress fraction
const QUIESCENCE_WINDOWS = ((0.25, 0.30), (0.45, 0.50))
const COLD_HOUR_THRESHOLD = 720.0u"hr"              # 30 days
const DIAPAUSE_HOUR_THRESHOLD = 1240.0u"hr"
const DESICCATION_TOLERANCE = 0.6                   # hydration index

# individual elements, shown separately before they're assembled below
lo_example = ThresholdController(metric=RawProgress(), bound=FixedBound(DIAPAUSE_WINDOW[1]), comparison=Above())
hi_example = ThresholdController(metric=RawProgress(), bound=FixedBound(DIAPAUSE_WINDOW[2]), comparison=Below())
chill_example = ThresholdController(metric=Accumulate((p, sig, m, as) -> sig.temperature < COLD_TEMPERATURE ? 1.0 : 0.0, 0.0u"hr"),
                                     bound=FixedBound(COLD_HOUR_THRESHOLD), comparison=Below())
induction_example = AllController(lo=lo_example, hi=hi_example, chill=chill_example)  # AND: all three must hold
demo_model = ComposedArrest(induction=induction_example)  # duration omitted -- its rate reads a sibling accumulator, only meaningful once nested
demo_state = initial_arrest_state(demo_model)
demo_signals = (; temperature=u"K"(5.0u"°C"))  # below COLD_TEMPERATURE, so chill accumulates
demo_level = arrest_level(demo_model, demo_state, 0.47, demo_signals)
demo_rate = advance_arrest(demo_model, demo_state, 0.47, demo_signals)
demo_state = step_state(demo_state, demo_rate, 1.0u"hr")

# build_locust_arrest_model assembles the same kind of elements -- for both the diapause and
# quiescence pathways -- into the full AnyArrestModel
function build_locust_arrest_model(; smoothing::SmoothingStrategy=HardBound())
    
    # diapause submodel
    
    lo, hi = DIAPAUSE_WINDOW
    duration_rate(progress, signals, model, arrest_state) =
        (lo < progress < hi &&
         arrest_state.induction.chill.metric.accumulator <= COLD_HOUR_THRESHOLD &&
         arrest_state.induction.duration.metric.accumulator <= DIAPAUSE_HOUR_THRESHOLD) ? 1.0 : 0.0
    
    # scale=2400.0u"hr" with SmoothBound(0.01) gives a 24-hour smoothing window at the limits.
    diapause_induction = AllController(
        lo=ThresholdController(metric=RawProgress(), bound=FixedBound(lo), comparison=Above(); smoothing),
        hi=ThresholdController(metric=RawProgress(), bound=FixedBound(hi), comparison=Below(); smoothing),
        chill=ThresholdController(metric=Accumulate((p, sig, m, as) -> sig.temperature < COLD_TEMPERATURE ? 1.0 : 0.0, 0.0u"hr"),
                                   bound=FixedBound(COLD_HOUR_THRESHOLD), comparison=Below(); smoothing, scale=2400.0u"hr"),
        duration=ThresholdController(metric=Accumulate(duration_rate, 0.0u"hr"), bound=FixedBound(DIAPAUSE_HOUR_THRESHOLD),
                                      comparison=Below(); smoothing, scale=2400.0u"hr"),
    )
    diapause_model = ComposedArrest(induction=diapause_induction)

    # quiescence submodel
    window_controller(w) = AllController(
        lo=ThresholdController(metric=RawProgress(), bound=FixedBound(w[1]), comparison=Above(); smoothing),
        hi=ThresholdController(metric=RawProgress(), bound=FixedBound(w[2]), comparison=Below(); smoothing),
        wet=ThresholdController(metric=RawSignal(:hydration), bound=FixedBound(DESICCATION_TOLERANCE), comparison=Below(); smoothing),
    )
    quiescence_model = ComposedArrest(induction=AnyController(w1=window_controller(QUIESCENCE_WINDOWS[1]),
                                                               w2=window_controller(QUIESCENCE_WINDOWS[2])))
    # combined arrest model: either pathway arrests the egg
    AnyArrestModel(diapause=diapause_model, quiescence=quiescence_model)
end

print_arrest_structure(build_locust_arrest_model(), "Chortoicetes_terminifera")

# --- example forcing, contrived to exercise both pathways: progress ramps 0->1 over 4000 hours, so
# the diapause window and second quiescence window (both 0.45-0.50) fall around hours 1800-2000,
# and the first quiescence window (0.25-0.30) falls around hours 1000-1200. A warm baseline with
# a cold snap centred on hour 1900 puts the only chilling inside the diapause window; a short dry
# spell inside the first quiescence window triggers that pathway alone.
n_hours = 4000
development_rate = 1/n_hours

# seasonal cycle (smooth half-sine dip over the run) plus a daily cycle, so temperature
# repeatedly crosses the diapause-breakage threshold each night near the seasonal minimum.
seasonal_temperature(t) = u"K"(30.0u"°C") - 15.0u"K" * sin(pi * t / n_hours)
daily_temperature_swing = 4.0u"K"
temperature_series = [seasonal_temperature(t) + daily_temperature_swing * sin(2pi * t / 24) for t in 1:n_hours]

hydration_series = [0.8 - 0.5 * exp(-((t - 1100) / 100)^2) for t in 1:n_hours]

# one manual model step at hour 1 -- levels are read from the state *before* step_state
# advances it, matching the evaluate-then-advance order simulate_locust_arrest loops below.
development_state_example = 0.46 # after threshold for diapause/quiescence, so both pathways are active
signals_example = (; temperature=u"K"(30.0u"°C"), hydration=0.5) # hydration below DESICCATION_TOLERANCE, so quiescence triggers
hard_model_example = build_locust_arrest_model(; smoothing=HardBound())
initial_state_example = initial_arrest_state(hard_model_example)
levels_example = arrest_level(hard_model_example, initial_state_example, development_state_example, signals_example)
diapause_level_example = arrest_level(hard_model_example.models.diapause, initial_state_example.diapause, development_state_example, signals_example)
quiescence_level_example = arrest_level(hard_model_example.models.quiescence, initial_state_example.quiescence, development_state_example, signals_example)
rate_example = advance_arrest(hard_model_example, initial_state_example, development_state_example, signals_example)
next_state_example = step_state(initial_state_example, rate_example, 1.0u"hr")
next_development_state_example = development_state_example + development_rate * (1 - levels_example)
println("hour 1: level=", levels_example, " (diapause=", diapause_level_example, ", quiescence=", quiescence_level_example, ")")
println("hour 1 -> 2: chill accumulator ", initial_state_example.diapause.induction.chill.metric.accumulator,
        " -> ", next_state_example.diapause.induction.chill.metric.accumulator)

# put them together in a loop to track each pathway's own level
function simulate_locust_arrest(model)
    state = initial_arrest_state(model)
    levels = Vector{Float64}(undef, n_hours)
    diapause_levels = Vector{Float64}(undef, n_hours)
    quiescence_levels = Vector{Float64}(undef, n_hours)
    developmental_progress_series = zeros(n_hours)
    for t in 1:n_hours
        progress = developmental_progress_series[t]
        signals = (; temperature=temperature_series[t], hydration=hydration_series[t])
        levels[t] = arrest_level(model, state, progress, signals)
        diapause_levels[t] = arrest_level(model.models.diapause, state.diapause, progress, signals)
        quiescence_levels[t] = arrest_level(model.models.quiescence, state.quiescence, progress, signals)
        t == n_hours && break
        rate = advance_arrest(model, state, progress, signals)
        state = step_state(state, rate, 1.0u"hr")
        # increment development progress for the next hour, scaled by the current arrest level (1 - level)
        developmental_progress_series[t + 1] = progress + development_rate * (1 - levels[t])
    end
    (; levels, diapause_levels, quiescence_levels, developmental_progress_series)
end

hard = simulate_locust_arrest(build_locust_arrest_model(; smoothing=HardBound()))
smooth_levels = simulate_locust_arrest(build_locust_arrest_model(; smoothing=SmoothBound(0.01))).levels

println("hours arrested overall (HardBound, either pathway): ", count(>(0.5), hard.levels), " / ", n_hours)
println("via diapause (cold snap inside the diapause window): ", count(>(0.5), hard.diapause_levels))
println("via quiescence (dry spell inside the first window): ", count(>(0.5), hard.quiescence_levels))

# Uncomment to plot temperature/hydration/progress forcing and arrest_level under both
# HardBound and SmoothBound side by side.
# using Plots
# p1 = plot(1:n_hours, u"°C".(temperature_series); label=nothing, title="Temperature")
# hline!(p1, [u"°C"(COLD_TEMPERATURE)]; label="cold_temperature", linestyle=:dash)

# p2 = plot(1:n_hours, hydration_series; label=nothing, title="Hydration index", ylabel="index")
# hline!(p2, [DESICCATION_TOLERANCE]; label="desiccation_tolerance", linestyle=:dash)

# p3 = plot(1:n_hours, hard.developmental_progress_series; label=nothing, title="Development progress", ylabel="fraction")
# hline!(p3, [DIAPAUSE_WINDOW...]; label="diapause_window", linestyle=:dash, color=:red)
# hline!(p3, [QUIESCENCE_WINDOWS[1]...]; label="quiescence_window_1", linestyle=:dot, color=:blue)

# p4 = plot(1:n_hours, hard.levels; label="HardBound (binary)", title="arrest_level", ylabel="level", fillrange=0, fillalpha=0.3)
# plot!(p4, 1:n_hours, smooth_levels; label="SmoothBound (gradual)", fillrange=0, fillalpha=0.3)

# panel = plot(p1, p2, p3, p4; layout=(4, 1), size=(1000, 1200), xlabel="hour",
#              plot_title="Locust (Chortoicetes) egg arrest model")
# display(panel)
