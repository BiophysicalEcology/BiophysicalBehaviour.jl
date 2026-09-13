# Composable bounds for ThresholdController: what the metric is compared against.

abstract type AbstractBound end

bound_state(::AbstractBound) = NamedTuple()
bound_value(b::AbstractBound, own_state) = error("no bound_value method for $(typeof(b))")

struct FixedBound{T} <: AbstractBound
    value::T
end
bound_value(b::FixedBound, own_state) = b.value
