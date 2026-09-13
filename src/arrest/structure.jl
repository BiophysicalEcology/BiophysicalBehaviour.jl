# Structural diagram of an arrest model -- the pathways/states it encodes
# (like a life-cycle flowchart), not a simulation result. Walks the
# composition tree and prints an indented text diagram: AND/OR gates at each
# composite, a human-readable condition at each leaf controller.

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
