# Proportional control: a metric's raw fraction toward a target, read directly as
# controller_level -- unlike ThresholdController (even under SmoothBound, which
# only smooths a narrow window right at the crossing), this ramps across the
# metric's whole range from zero to target.
Base.@kwdef struct ProportionalController{M<:AbstractMetric,T} <: AbstractCondition
    metric::M
    target::T
end

function initial_controller_state(c::ProportionalController)
    (; metric=metric_state(c.metric))
end

function controller_rate(c::ProportionalController, own_state, progress, signals, model, arrest_state)
    (; metric=metric_rate(c.metric, own_state.metric, progress, signals, model, arrest_state))
end

function controller_level(c::ProportionalController, own_state, progress, signals, model, arrest_state)
    value = metric_value(c.metric, own_state.metric, progress, signals)
    clamp(value / c.target, 0.0, 1.0)
end

# smooth, monotone ramp -- no crossing event to root-find.
register_callback(::ProportionalController) = false
