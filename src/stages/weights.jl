# Sequential gating: for transition levels l_1, ..., l_{n-1},
#   w_1 = 1 - l_1
#   w_i = l_1 * ... * l_{i-1} * (1 - l_i)   for 1 < i < n
#   w_n = l_1 * ... * l_{n-1}
# Nonnegative for levels in [0,1], sums to 1 algebraically. Numeric-type-
# preserving: no literal Float64 anywhere, so Float32/Unitful/dual-number
# controller levels pass through unchanged.

function stage_weights(sequence::StageSequence{StageKeys,TransitionKeys}, stage_state, progress, signals) where {StageKeys,TransitionKeys}
    n = length(StageKeys)
    n == 1 && return NamedTuple{StageKeys}((1,))  # singleton: documented Int exception, no level to derive a type from

    levels = map(transitions(sequence), TransitionKeys) do transition, key
        controller_level(transition_controller(transition), getfield(stage_state, key), progress, signals, sequence, stage_state)
    end
    identity_level = one(first(levels))
    cumulative = accumulate(*, levels)                  # (l1, l1*l2, ..., prod(levels)) -- inclusive scan
    remaining = (identity_level, cumulative...)         # remaining[i] == prod(levels[1:i-1])
    weights = ntuple(n) do i
        i == n ? remaining[n] : remaining[i] * (one(levels[i]) - levels[i])
    end
    NamedTuple{StageKeys}(weights)
end

"""
    current_stage(sequence, stage_state, progress, signals)

The single stage with the highest weight (ties resolve to the earlier
stage). For a blend, use `stage_value` instead.
"""
current_stage(sequence, stage_state, progress, signals) =
    stages(sequence)[argmax(values(stage_weights(sequence, stage_state, progress, signals)))]

"""
    stage_value(sequence, stage_state, progress, signals, accessor)

Weighted sum of `accessor(stage_traits(stage), progress, signals)` over
`stage_weights`. `accessor`'s return type must support `*` with the weight
type and `+` with itself, and must be callable for every stage in the
sequence (not only the nominal "current" one -- a SmoothBound transition's
tail is never exactly zero far from threshold). For non-blendable stage data
(a whole model object), use `stage_traits(current_stage(...))` instead.
"""
function stage_value(sequence::StageSequence, stage_state, progress, signals, accessor)
    weights = values(stage_weights(sequence, stage_state, progress, signals))
    terms = map(weights, stages(sequence)) do weight, stage
        weight * accessor(stage_traits(stage), progress, signals)
    end
    reduce(+, terms)
end
