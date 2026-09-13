# Wraps initial_stage_state in a ComponentArray for named ODE-state access
# (u.stages.<transition_key>... instead of fixed positional indices).

stage_component(sequence::StageSequence) = ComponentArray(initial_stage_state(sequence))
