# Metric-vs-bound comparison: covers threshold and cumulation via composition.
Base.@kwdef struct ThresholdController{M<:AbstractMetric,B<:AbstractBound,D<:AbstractDirection,C<:AbstractComparison,S<:SmoothingStrategy,SC} <: AbstractArrestController
    metric::M
    bound::B
    direction::D = AnyDirection()
    comparison::C = BelowBound()
    smoothing::S = HardBound()
    scale::SC = 1.0
end

function initial_controller_state(c::ThresholdController)
    (; metric=metric_state(c.metric), bound=bound_state(c.bound))
end

function controller_rate(c::ThresholdController, own_state, progress, signals, model, arrest_state)
    (; metric=metric_rate(c.metric, own_state.metric, progress, signals, model, arrest_state), bound=NamedTuple())
end

function _gap(c::ThresholdController, own_state, progress, signals)
    value = metric_value(c.metric, own_state.metric, progress, signals)
    bound = bound_value(c.bound, own_state.bound)
    signed_gap(c.comparison, value, bound)
end

function controller_level(c::ThresholdController, own_state, progress, signals, model, arrest_state)
    gap = _gap(c, own_state, progress, signals)
    base = safe_step(c.smoothing, gap; scale=c.scale)
    if c.direction isa AnyDirection
        base
    else
        rate = metric_rate_value(c.metric, own_state.metric, progress, signals)
        direction_gate(c.direction, rate) ? base : 0.0
    end
end

# only HardBound has a real discontinuity worth root-finding.
register_callback(c::ThresholdController) = c.smoothing isa HardBound

function trigger_conditions(c::ThresholdController)
    register_callback(c) || return ()
    ((own_state, progress, signals, model, arrest_state) -> _gap(c, own_state, progress, signals),)
end
