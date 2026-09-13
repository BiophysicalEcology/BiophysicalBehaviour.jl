# Linear diagram of a StageSequence's topology -- not a simulation result.
# Shows both the transition type and its threshold, so same-type transitions
# with different thresholds are distinguishable.

"""
    print_stage_structure(sequence, label="sequence")

Prints `stage --[TransitionType: condition]--> stage --> ...` for the given
`StageSequence`.
"""
function print_stage_structure(sequence::StageSequence, label="sequence")
    stage_list = stages(sequence)
    transition_list = transitions(sequence)
    print(label, ": ", stage_key(stage_list[1]))
    for i in eachindex(transition_list)
        transition = transition_list[i]
        print(" --[", nameof(typeof(transition)), ": ", node_label(transition_controller(transition)), "]--> ")
        print(stage_key(stage_list[i+1]))
    end
    println()
end
