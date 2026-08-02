# Escape hatch for anything not expressible as metric-vs-bound -- genuinely
# custom/order-dependent logic stays plain code here (not a composable
# "Rule" type), matching RuleBasedSequentialControl's hardcoded if/elseif idiom.
#
# `condition`: (own_state, progress, signals, model, arrest_state) -> Real,
# same signed-distance convention as ThresholdController's `_gap`.
Base.@kwdef struct FunctionController{F,L} <: AbstractArrestController
    condition::F
    level_function::L = nothing
    needs_callback::Bool = true
end

function level(c::FunctionController, own_state, progress, signals, model, arrest_state)
    if c.level_function === nothing
        c.condition(own_state, progress, signals, model, arrest_state) >= 0 ? 1.0 : 0.0
    else
        c.level_function(own_state, progress, signals, model, arrest_state)
    end
end

register_callback(c::FunctionController) = c.needs_callback
function trigger_conditions(c::FunctionController)
    c.needs_callback || return ()
    (c.condition,)
end
