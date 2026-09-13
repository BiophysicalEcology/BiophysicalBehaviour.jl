# Bare condition closures with no attached transition identity -- mapping a
# fired callback index back to a transition is the host's responsibility.

function stage_conditions(sequence::StageSequence)
    _flatten_tuples(map(transitions(sequence)) do transition
        key = transition_key(transition)
        map(trigger_conditions(transition_controller(transition))) do trigger
            (stage_state, progress, signals) -> trigger(getfield(stage_state, key), progress, signals, sequence, stage_state)
        end
    end)
end
