# Event-driven diurnal behavior for the transient body-temperature models. Port of NicheMapR's
# trans_behav.R (sleep -> bask -> forage <-> cool), using ContinuousCallback in place of
# R's rootfun events. Simplifications relative to R noted in docs/transient_body_temperature.md.

abstract type TransientBehavioralPhase end
struct SleepPhase <: TransientBehavioralPhase end
struct BaskPhase <: TransientBehavioralPhase end
struct ForagePhase <: TransientBehavioralPhase end
struct CoolPhase <: TransientBehavioralPhase end

_phase_forcing(::SleepPhase, sun, shade) = shade
_phase_forcing(::BaskPhase, sun, shade) = sun
_phase_forcing(::ForagePhase, sun, shade) = sun
_phase_forcing(::CoolPhase, sun, shade) = shade

_phase_posture(::SleepPhase) = Intermediate()
_phase_posture(::BaskPhase) = NormalToSun()
_phase_posture(::ForagePhase) = Intermediate()
_phase_posture(::CoolPhase) = Intermediate()

_phase_state(::SleepPhase) = Resting()
_phase_state(::BaskPhase) = Basking()
_phase_state(::ForagePhase) = Active()
_phase_state(::CoolPhase) = Resting()

# State is Float64 (one-lump, K) or Vector{Float64} (two-lump, [core, shell] K) - dispatch
# picks out the core temperature that drives all behavioral thresholds.
_core_temperature_of(u::Real) = u
_core_temperature_of(u::AbstractVector) = u[1]

# Conditions take plain Float64 (u in K, t in s), negative while in-phase, crossing zero
# when the phase should end (ContinuousCallback root). Sub-signals keep their own unit
# scale — min/max only need matching sign, giving AND=`min`, OR=`max`.
_basking_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.basking_temperature_min)
_active_min_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.active_temperature_min)
_active_max_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.active_temperature_max)

function _phase_condition(::SleepPhase, zenith_signal, limits)
    (u, t) -> min(-zenith_signal(t), _basking_signal(u, limits))
end
function _phase_condition(::BaskPhase, zenith_signal, limits)
    (u, t) -> max(_active_min_signal(u, limits), zenith_signal(t))
end
function _phase_condition(::ForagePhase, zenith_signal, limits)
    (u, t) -> max(_active_max_signal(u, limits), -_active_min_signal(u, limits), zenith_signal(t))
end
# trans_behav.R's `forage` event: resume at T_F_min, or shade air temp + 1K if that's
# higher (otherwise unreachable in warm shade), falling back to T_F_min if shade is
# already near T_F_max.
function _cool_resume_threshold(t, limits, shade_air_temperature)
    active_min = ustrip(u"K", limits.active_temperature_min)
    active_max = ustrip(u"K", limits.active_temperature_max)
    shade_air = shade_air_temperature(t)
    shade_target = shade_air > active_max - 2.0 ? 0.0 : shade_air + 1.0
    max(active_min, shade_target)
end

function _phase_condition(::CoolPhase, zenith_signal, limits, shade_air_temperature)
    (u, t) -> max(_cool_resume_threshold(t, limits, shade_air_temperature) - _core_temperature_of(u), zenith_signal(t))
end

_next_phase(::SleepPhase, u, zenith_signal, t, limits, shade_air_temperature) = BaskPhase()
function _next_phase(::BaskPhase, u, zenith_signal, t, limits, shade_air_temperature)
    zenith_signal(t) >= 0.0 ? SleepPhase() : ForagePhase()
end
function _next_phase(::ForagePhase, u, zenith_signal, t, limits, shade_air_temperature)
    if zenith_signal(t) >= 0.0
        SleepPhase()
    # Compare the two signals to each other, not each to zero independently, so root-finding
    # noise at either boundary can't flip this to the wrong phase.
    elseif _active_max_signal(u, limits) >= -_active_min_signal(u, limits)
        CoolPhase()
    else
        BaskPhase()
    end
end
function _next_phase(::CoolPhase, u, zenith_signal, t, limits, shade_air_temperature)
    zenith_signal(t) >= 0.0 ? SleepPhase() : ForagePhase()
end

# RHS dispatch: one-lump (scalar state, HeatExchange.ectotherm_onelump) vs two-lump (2-vector
# [core, shell], HeatExchange.ectotherm_twolump). Both methods take the same keywords so callers
# don't need to know which model is in play - `shell_thickness` is simply unused for ectotherm_onelump.
function _body_temperature_rate(
    u::Real, t, body, environment_pars, environment_vars;
    internal_conduction, shell_thickness=nothing, posture, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
)
    ustrip(u"K/s", HeatExchange.ectotherm_onelump(
        u * u"K", t * u"s", body, environment_pars, environment_vars;
        internal_conduction, posture, body_absorptivity, emissivity,
        sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
    ))
end
function _body_temperature_rate(
    u::AbstractVector, t, body, environment_pars, environment_vars;
    internal_conduction, shell_thickness, posture, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
)
    out = HeatExchange.ectotherm_twolump(
        (core_temperature=u[1] * u"K", shell_temperature=u[2] * u"K"), t * u"s", body, environment_pars, environment_vars;
        internal_conduction, shell_thickness, posture, body_absorptivity, emissivity,
        sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
    )
    Float64[ustrip(u"K/s", out.core_temperature_rate), ustrip(u"K/s", out.shell_temperature_rate)]
end

_package_trajectory(all_t, all_u::Vector{Float64}, all_state) =
    (; t=all_t * u"s", core_temperature=all_u * u"K", state=all_state)
function _package_trajectory(all_t, all_u::Vector{Vector{Float64}}, all_state)
    (;
        t=all_t * u"s",
        core_temperature=[uu[1] for uu in all_u] * u"K",
        shell_temperature=[uu[2] for uu in all_u] * u"K",
        state=all_state,
    )
end

"""
    simulate_diurnal_behavior(times, core_temperature_init, body, environment_pars,
                               sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
                               limits::EctothermBehavioralLimits;
                               internal_conduction, body_absorptivity, emissivity,
                               sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
                               solver=OrdinaryDiffEqTsit5.Tsit5(),
                               solver_kwargs=(;), smoothing=HardBound(), max_bouts=100*length(times),
                               bout_chunk=3600.0)
    simulate_diurnal_behavior(times, (; core_temperature, shell_temperature), body, environment_pars,
                               sun_forcing, shade_forcing, limits; internal_conduction, shell_thickness,
                               body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
                               metabolic_heat_volumetric, solver=OrdinaryDiffEqTsit5.Tsit5(),
                               solver_kwargs=(;), smoothing=HardBound(), max_bouts=100*length(times),
                               bout_chunk=3600.0)

Event-driven diurnal behavioral thermoregulation (sleep → bask → forage ⇄ cool → sleep).
Dispatches on the initial state: a plain temperature runs the one-lump model
(`HeatExchange.ectotherm_onelump`); a `(; core_temperature, shell_temperature)` NamedTuple runs the
two-lump model (`HeatExchange.ectotherm_twolump`, needs a `shell_thickness` keyword). Both reuse
`limits`' `active_temperature_min/max`/`basking_temperature_min` thresholds. Port of
NicheMapR's `trans_behav.R` (see file banner for what's simplified relative to R).

# Returns
NamedTuple with `t` (s), `core_temperature` (K) [, `shell_temperature` (K) for the two-lump
method], `state` (`Vector{OrganismState}`), one entry per accepted solver step across all bouts.
"""
function simulate_diurnal_behavior(
    times, core_temperature_init::Unitful.Quantity, body, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    internal_conduction, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), max_bouts=100 * length(times),
    bout_chunk=3600.0,
)
    _simulate_diurnal_behavior(
        times, ustrip(u"K", core_temperature_init), body, environment_pars, sun_forcing, shade_forcing, limits;
        internal_conduction, shell_thickness=nothing, body_absorptivity, emissivity,
        sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
        solver, solver_kwargs, smoothing, max_bouts, bout_chunk,
    )
end
function simulate_diurnal_behavior(
    times, state::NamedTuple{(:core_temperature, :shell_temperature)}, body, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    internal_conduction, shell_thickness, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), max_bouts=100 * length(times),
    bout_chunk=3600.0,
)
    u0 = Float64[ustrip(u"K", state.core_temperature), ustrip(u"K", state.shell_temperature)]
    _simulate_diurnal_behavior(
        times, u0, body, environment_pars, sun_forcing, shade_forcing, limits;
        internal_conduction, shell_thickness, body_absorptivity, emissivity,
        sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
        solver, solver_kwargs, smoothing, max_bouts, bout_chunk,
    )
end

function _simulate_diurnal_behavior(
    times, u0, body, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    internal_conduction, shell_thickness, body_absorptivity, emissivity,
    sky_view_factor, ground_view_factor, metabolic_heat_volumetric,
    solver, solver_kwargs, smoothing::SmoothingStrategy, max_bouts, bout_chunk,
)
    t0, tend = ustrip(u"s", first(times)), ustrip(u"s", last(times))
    zenith_signal(t) = ustrip(u"°", sun_forcing(t * u"s"; smoothing).zenith_angle - 90.0u"°")
    shade_air_temperature(t) = ustrip(u"K", shade_forcing(t * u"s"; smoothing).air_temperature)

    phase = zenith_signal(t0) < 0.0 ? ForagePhase() : SleepPhase()
    u_current = u0
    t_current = t0

    all_t = Float64[]
    all_u = typeof(u0)[]
    all_state = OrganismState[]

    for _ in 1:max_bouts
        t_current >= tend && break

        for _ in 1:8
            condition = phase isa CoolPhase ?
                _phase_condition(phase, zenith_signal, limits, shade_air_temperature) :
                _phase_condition(phase, zenith_signal, limits)
            condition(u_current, t_current) < 0 && break
            phase = _next_phase(phase, u_current, zenith_signal, t_current, limits, shade_air_temperature)
        end

        forcing = _phase_forcing(phase, sun_forcing, shade_forcing)
        posture = _phase_posture(phase)
        f = (u, _, t) -> _body_temperature_rate(
            u, t, body, environment_pars, forcing(t * u"s"; smoothing);
            internal_conduction, shell_thickness, posture, body_absorptivity, emissivity,
            sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
        )
        condition = phase isa CoolPhase ?
            _phase_condition(phase, zenith_signal, limits, shade_air_temperature) :
            _phase_condition(phase, zenith_signal, limits)
        # affect_neg!=nothing: only the upward crossing ends a phase.
        callback = OrdinaryDiffEqTsit5.ContinuousCallback(
            (u, t, integrator) -> condition(u, t), OrdinaryDiffEqTsit5.terminate!;
            affect_neg! = nothing,
        )
        t_chunk_end = min(t_current + bout_chunk, tend)
        problem = OrdinaryDiffEqTsit5.ODEProblem(f, u_current, (t_current, t_chunk_end))
        sol = OrdinaryDiffEqTsit5.solve(problem, solver; callback, solver_kwargs...)

        append!(all_t, sol.t)
        append!(all_u, sol.u)
        append!(all_state, fill(_phase_state(phase), length(sol.t)))

        t_next = sol.t[end]
        u_next = sol.u[end]
        if t_next < t_chunk_end
            phase = _next_phase(phase, u_next, zenith_signal, t_next, limits, shade_air_temperature)
        end
        t_current, u_current = t_next, u_next
    end

    return _package_trajectory(all_t, all_u, all_state)
end
