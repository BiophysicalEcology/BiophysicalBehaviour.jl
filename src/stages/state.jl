# State-threading: a stage_state is a NamedTuple keyed by transition_key,
# mirroring arrest_state. Always fully populated -- a stateless RawProgress
# condition still gets a content-free, present entry, never a missing key.

function initial_stage_state(sequence::StageSequence{<:Any,TransitionKeys}) where {TransitionKeys}
    NamedTuple{TransitionKeys}(map(transitions(sequence)) do transition
        initial_controller_state(transition_controller(transition))
    end)
end

"""
    advance_stage(sequence, stage_state, progress, signals)

RATES, not advanced state -- apply with `step_state`.
"""
function advance_stage(sequence::StageSequence{<:Any,TransitionKeys}, stage_state, progress, signals) where {TransitionKeys}
    NamedTuple{TransitionKeys}(map(transitions(sequence), TransitionKeys) do transition, key
        controller_rate(transition_controller(transition), getfield(stage_state, key), progress, signals, sequence, stage_state)
    end)
end
