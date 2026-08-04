# Structural diagram of an arrest model -- the pathways/states it encodes
# (like a life-cycle flowchart), not a simulation result. Walks the
# composition tree and prints an indented text diagram: AND/OR gates at each
# composite, a human-readable condition at each leaf controller.

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

node_label(c::AbstractArrestController) = error("no node_label method for $(typeof(c))")
function node_label(c::ThresholdController)
    label = "$(describe_metric(c.metric)) $(describe_comparison(c.comparison)) $(describe_bound(c.bound))$(describe_direction(c.direction))"
    c.smoothing isa HardBound ? label : label * " [smooth]"
end
node_label(::FunctionController) = "<custom function>"
node_label(::NeverController) = "never"

"""
    print_arrest_structure(model, label="model")

Prints the composition tree of an arrest model or controller as an indented
text diagram: AND/OR gates at each composite, a human-readable condition
(metric, comparison, bound) at each leaf. Shows the pathways/states a model
encodes, not a simulation result.
"""
function print_arrest_structure(model, label="model"; prefix="", is_last::Bool=true)
    connector = is_last ? "└─ " : "├─ "
    child_prefix = prefix * (is_last ? "   " : "│  ")

    if model isa Union{AnyArrestModel,AllArrestModel}
        gate = model isa AnyArrestModel ? "OR" : "AND"
        println(prefix, connector, label, " (", gate, ")")
        ks = keys(model.models)
        for (i, k) in enumerate(ks)
            print_arrest_structure(getfield(model.models, k), string(k); prefix=child_prefix, is_last=i == length(ks))
        end
    elseif model isa ComposedArrest
        println(prefix, connector, label, " (induction AND NOT breakage)")
        breaks = !(model.breakage isa NeverController)
        print_arrest_structure(model.induction, "induction"; prefix=child_prefix, is_last=!breaks)
        breaks && print_arrest_structure(model.breakage, "breakage"; prefix=child_prefix, is_last=true)
    elseif model isa Union{AnyController,AllController}
        gate = model isa AnyController ? "OR" : "AND"
        println(prefix, connector, label, " (", gate, ")")
        ks = keys(model.controllers)
        for (i, k) in enumerate(ks)
            print_arrest_structure(getfield(model.controllers, k), string(k); prefix=child_prefix, is_last=i == length(ks))
        end
    else
        println(prefix, connector, label, ": ", node_label(model))
    end
end
