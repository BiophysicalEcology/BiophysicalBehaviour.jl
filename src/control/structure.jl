# Human-readable strings for a metric, bound, comparison, direction, or leaf
# AbstractCondition; shared by print_arrest_structure/print_stage_structure.

describe_metric(m::AbstractMetric) = error("no describe_metric method for $(typeof(m))")
describe_metric(m::RawSignal) = "signal :$(m.signal)"
describe_metric(::RawProgress) = "progress"
describe_metric(::Accumulate) = "accumulator"

describe_bound(b::AbstractBound) = error("no describe_bound method for $(typeof(b))")
describe_bound(b::FixedBound) = string(b.value)

describe_comparison(c::AbstractComparison) = error("no describe_comparison method for $(typeof(c))")
describe_comparison(::BelowBound) = "<"
describe_comparison(::AboveBound) = ">"

describe_direction(d::AbstractDirection) = error("no describe_direction method for $(typeof(d))")
describe_direction(::AnyDirection) = ""
describe_direction(::RisingDirection) = " (rising)"
describe_direction(::FallingDirection) = " (falling)"

node_label(c::AbstractCondition) = error("no node_label method for $(typeof(c))")
function node_label(c::ThresholdController)
    label = "$(describe_metric(c.metric)) $(describe_comparison(c.comparison)) $(describe_bound(c.bound))$(describe_direction(c.direction))"
    c.smoothing isa HardBound ? label : label * " [smooth]"
end
node_label(c::ProportionalController) = "$(describe_metric(c.metric)) / $(c.target) (proportional)"
node_label(::FunctionController) = "<custom function>"
node_label(::NeverController) = "never"
