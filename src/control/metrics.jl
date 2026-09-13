# Composable metrics for ThresholdController: what's being compared to the bound.

abstract type AbstractMetric end

metric_state(::AbstractMetric) = NamedTuple()
metric_rate(::AbstractMetric, own_state, progress, signals, model, arrest_state) = NamedTuple()
metric_value(m::AbstractMetric, own_state, progress, signals) = error("no metric_value method for $(typeof(m))")
metric_rate_value(m::AbstractMetric, own_state, progress, signals) = error("no metric_rate_value method for $(typeof(m))")

# reads signals[signal] directly, no state of its own.
struct RawSignal{S<:Symbol} <: AbstractMetric
    signal::S
end

metric_value(m::RawSignal, own_state, progress, signals) = signal_value(getfield(signals, m.signal))
metric_rate_value(m::RawSignal, own_state, progress, signals) = signal_rate(getfield(signals, m.signal))

# reads `progress` directly (not a signal) -- e.g. a diapause/quiescence
# window is a threshold on the process's own progress, not an external cue.
struct RawProgress <: AbstractMetric end
metric_value(::RawProgress, own_state, progress, signals) = progress
# progress is driven externally by the host, not tracked as own-state -- no rate to report.
metric_rate_value(::RawProgress, own_state, progress, signals) =
    error("RawProgress has no notion of rate; RisingDirection/FallingDirection is not supported for RawProgress")

# integrates a rate functor into its own accumulator, from `init`. Generalizes
# chill_accumulation/diapause_duration and DualAccumulatorArrest -- two
# Accumulate metrics under AllController, the second's rate reading the
# first's accumulator via the full `arrest_state` arg.
struct Accumulate{R,I} <: AbstractMetric
    rate::R
    init::I
end
Accumulate(rate) = Accumulate(rate, 0.0)

metric_state(m::Accumulate) = (; accumulator=m.init)
metric_rate(m::Accumulate, own_state, progress, signals, model, arrest_state) =
    (; accumulator=m.rate(progress, signals, model, arrest_state))

metric_value(m::Accumulate, own_state, progress, signals) = own_state.accumulator
# an accumulating metric's "trend" is just the sign of its own rate.
metric_rate_value(m::Accumulate, own_state, progress, signals) = own_state.accumulator
