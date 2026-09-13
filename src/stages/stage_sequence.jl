# Ordered life-cycle topology: N stages, N-1 transitions. See docs/stages.md.

abstract type AbstractStage{T} end
abstract type AbstractTransition{T<:AbstractCondition} end

stage_traits(stage::AbstractStage) = stage.traits
transition_controller(transition::AbstractTransition) = transition.controller

# Override the Type{S} method, not the instance method, when nameof(S) alone
# isn't a unique key (e.g. distinct parametric instantiations of one type).
stage_key(::Type{S}) where {S<:AbstractStage} = nameof(S)
stage_key(stage::AbstractStage) = stage_key(typeof(stage))

transition_key(::Type{S}) where {S<:AbstractTransition} = nameof(S)
transition_key(transition::AbstractTransition) = transition_key(typeof(transition))

function _check_unique_keys(stage_keys::Tuple{Vararg{Symbol}}, transition_keys::Tuple{Vararg{Symbol}})
    allunique(stage_keys) || throw(ArgumentError("duplicate stage_key among stages: $stage_keys"))
    allunique(transition_keys) || throw(ArgumentError("duplicate transition_key among transitions: $transition_keys"))
    nothing
end

# StageKeys/TransitionKeys: resolved once at construction, carried as type
# parameters for type-stable NamedTuple{...} construction in hot paths.
struct StageSequence{StageKeys,TransitionKeys,S<:Tuple{Vararg{AbstractStage}},T<:Tuple{Vararg{AbstractTransition}}}
    stages::S
    transitions::T

    function StageSequence(stages::S, transitions::T) where {S<:Tuple{Vararg{AbstractStage}},T<:Tuple{Vararg{AbstractTransition}}}
        length(stages) == length(transitions) + 1 ||
            throw(ArgumentError("a StageSequence requires exactly one more stage than transitions"))
        stage_keys = map(stage_key, stages)
        transition_keys = map(transition_key, transitions)
        _check_unique_keys(stage_keys, transition_keys)
        new{stage_keys,transition_keys,S,T}(stages, transitions)
    end
end

function StageSequence(initial::AbstractStage, rest::Vararg{Pair{<:AbstractTransition,<:AbstractStage}})
    StageSequence((initial, map(Base.last, rest)...), map(Base.first, rest))
end

stages(sequence::StageSequence) = sequence.stages
transitions(sequence::StageSequence) = sequence.transitions
