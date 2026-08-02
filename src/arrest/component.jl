# Wraps initial_arrest_state in a ComponentArray for named ODE-state access
# (u.arrest.induction... instead of fixed positional indices).

arrest_component(model::AbstractArrestModel) =
    ComponentArray(initial_arrest_state(model))
