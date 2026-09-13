# Recursively advances a nested NamedTuple state by rate*dt; an empty rate
# NamedTuple leaves the matching state branch unchanged.
step_state(state::NamedTuple, rate::NamedTuple{(),Tuple{}}, dt) = state
step_state(state::NamedTuple, rate::NamedTuple, dt) =
    NamedTuple{keys(state)}(map((s, r) -> step_state(s, r, dt), values(state), values(rate)))
step_state(state::Number, rate::Number, dt) = state + rate * dt
