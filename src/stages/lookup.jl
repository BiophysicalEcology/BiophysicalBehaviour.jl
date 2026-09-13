# Type-keyed lookup via `isa` (matches partial parametric types, e.g.
# Instar{3} against any concrete Instar{3,...}). haskey tests the matched
# wrapper's presence, not its unwrapped value, so a present-but-`nothing`
# value is never confused with absence.

function _get_wrapper(members::Tuple, ::Type{K}) where {K}
    for member in members
        member isa K && return member
    end
    return nothing
end

Base.haskey(sequence::StageSequence, ::Type{K}) where {K<:AbstractStage} =
    _get_wrapper(stages(sequence), K) !== nothing
Base.haskey(sequence::StageSequence, ::Type{K}) where {K<:AbstractTransition} =
    _get_wrapper(transitions(sequence), K) !== nothing

function Base.getindex(sequence::StageSequence, ::Type{K}) where {K<:AbstractStage}
    wrapper = _get_wrapper(stages(sequence), K)
    wrapper === nothing && throw(ArgumentError("no stage matching $K in this StageSequence"))
    stage_traits(wrapper)
end
function Base.getindex(sequence::StageSequence, ::Type{K}) where {K<:AbstractTransition}
    wrapper = _get_wrapper(transitions(sequence), K)
    wrapper === nothing && throw(ArgumentError("no transition matching $K in this StageSequence"))
    transition_controller(wrapper)
end

function Base.get(sequence::StageSequence, ::Type{K}, default) where {K<:AbstractStage}
    wrapper = _get_wrapper(stages(sequence), K)
    wrapper === nothing ? default : stage_traits(wrapper)
end
function Base.get(sequence::StageSequence, ::Type{K}, default) where {K<:AbstractTransition}
    wrapper = _get_wrapper(transitions(sequence), K)
    wrapper === nothing ? default : transition_controller(wrapper)
end
