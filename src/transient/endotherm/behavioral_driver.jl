# Event-driven activity/rest cycling for the transient endotherm model: an animal alternates
# between an active phase (elevated metabolic rate, e.g. flight or running - and typically
# elevated wind exposure) and resting (basal metabolic rate) purely as a function of core
# temperature, using ContinuousCallback the same way simulate_transient_behavior does for
# ectotherms. No day/night cycle - see docs/transient_body_temperature.md for scope notes.
# Insulated bodies only: metabolic_heat_flow must be an explicit per-phase input, which
# HeatExchange.onelump only supports for Insulated (not Naked) organisms.

abstract type EndothermActivityPhase end
struct EndothermActivePhase <: EndothermActivityPhase end
struct EndothermRestingPhase <: EndothermActivityPhase end

endotherm_phase_forcing(::EndothermActivePhase, active, resting) = active
endotherm_phase_forcing(::EndothermRestingPhase, active, resting) = resting

endotherm_phase_metabolic_heat_flow(::EndothermActivePhase, active_mhf, resting_mhf) = active_mhf
endotherm_phase_metabolic_heat_flow(::EndothermRestingPhase, active_mhf, resting_mhf) = resting_mhf

endotherm_phase_state(::EndothermActivePhase) = Active()
endotherm_phase_state(::EndothermRestingPhase) = Resting()

endotherm_next_phase(::EndothermActivePhase) = EndothermRestingPhase()
endotherm_next_phase(::EndothermRestingPhase) = EndothermActivePhase()

# Conditions take plain Float64 (u in K), negative while in-phase, crossing zero when the
# phase should end (ContinuousCallback root) - same convention as the ectotherm driver.
endotherm_phase_condition(::EndothermActivePhase, active_temperature_max, resume_temperature) =
    (u, t) -> u - ustrip(u"K", active_temperature_max)
endotherm_phase_condition(::EndothermRestingPhase, active_temperature_max, resume_temperature) =
    (u, t) -> ustrip(u"K", resume_temperature) - u

"""
    simulate_endotherm_activity_cycle(times, core_temperature_init, organism::Organism, environment_pars,
                                       active_forcing::EnvironmentForcing, resting_forcing::EnvironmentForcing;
                                       active_metabolic_heat_flow, active_temperature_max, resume_temperature,
                                       resting_metabolic_heat_flow=metabolism_pars(organism).metabolic_heat_flow,
                                       solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
                                       smoothing=HardBound(), max_bouts=50*length(times), bout_chunk=3600.0)

Event-driven active/resting cycling for an insulated endotherm: active (elevated metabolic
rate, e.g. flight or running, typically under `active_forcing`'s higher wind speed) until core
temperature reaches `active_temperature_max`, then resting (`resting_metabolic_heat_flow`,
`resting_forcing`) until it cools to `resume_temperature`, repeating. Effectors (insulation,
panting, skin wetness, flesh conductivity) are held fixed throughout, as in
`HeatExchange.onelump`. Starts in the active phase.

# Returns
NamedTuple with `t` (s), `core_temperature` (K), `core_temperature_rate` (K/s), `state`
(`Vector{OrganismState}`, `Active`/`Resting`), one entry per accepted solver step across all bouts.
"""
function simulate_endotherm_activity_cycle(
    times, core_temperature_init, organism::Organism, environment_pars,
    active_forcing::EnvironmentForcing, resting_forcing::EnvironmentForcing;
    active_metabolic_heat_flow, active_temperature_max, resume_temperature,
    resting_metabolic_heat_flow=metabolism_pars(organism).metabolic_heat_flow,
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), max_bouts=50 * length(times),
    bout_chunk=3600.0,
)
    t0, tend = ustrip(u"s", first(times)), ustrip(u"s", last(times))
    phase = EndothermActivePhase()
    u_current = ustrip(u"K", core_temperature_init)
    t_current = t0

    all_t = Float64[]
    all_u = Float64[]
    all_state = OrganismState[]

    for _ in 1:max_bouts
        t_current >= tend && break

        condition = endotherm_phase_condition(phase, active_temperature_max, resume_temperature)
        condition(u_current, t_current) >= 0 && (phase = endotherm_next_phase(phase))

        forcing = endotherm_phase_forcing(phase, active_forcing, resting_forcing)
        mhf = endotherm_phase_metabolic_heat_flow(phase, active_metabolic_heat_flow, resting_metabolic_heat_flow)
        f = (u, _, t) -> ustrip(u"K/s", HeatExchange.onelump(
            u * u"K", t * u"s", organism, (; environment_pars, environment_vars=forcing(t * u"s"; smoothing));
            metabolic_heat_flow=mhf, smoothing,
        ).core_temperature_rate)

        condition = endotherm_phase_condition(phase, active_temperature_max, resume_temperature)
        # affect_neg!=nothing: only the upward crossing of the (phase-specific) signal ends a phase.
        callback = OrdinaryDiffEqTsit5.ContinuousCallback(
            (u, t, integrator) -> condition(u, t), OrdinaryDiffEqTsit5.terminate!;
            affect_neg! = nothing,
        )
        t_chunk_end = min(t_current + bout_chunk, tend)
        problem = OrdinaryDiffEqTsit5.ODEProblem(f, u_current, (t_current, t_chunk_end))
        sol = OrdinaryDiffEqTsit5.solve(problem, solver; callback, solver_kwargs...)

        append!(all_t, sol.t)
        append!(all_u, sol.u)
        append!(all_state, fill(endotherm_phase_state(phase), length(sol.t)))

        t_next, u_next = sol.t[end], sol.u[end]
        t_next < t_chunk_end && (phase = endotherm_next_phase(phase))
        t_current, u_current = t_next, u_next
    end

    return (; t=all_t * u"s", core_temperature=all_u * u"K", state=all_state)
end
