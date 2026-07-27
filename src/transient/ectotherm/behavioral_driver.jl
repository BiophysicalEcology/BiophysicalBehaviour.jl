# Event-driven diurnal behavior for the transient body-temperature models. Port of NicheMapR's
# trans_behav.R (sleep -> bask -> forage <-> cool), using ContinuousCallback in place of
# R's rootfun events. Simplifications relative to R noted in docs/transient_body_temperature.md.

"""
    TransientBehavioralPhase

Abstract supertype for a phase in `simulate_transient_behavior`'s event-driven transient model.

To add a custom phase, define a singleton subtype and implement methods for: `phase_forcing`,
`phase_posture`, `phase_state`, `phase_key`, `phase_condition`, `next_phase` (all ordinary
dispatch, no registration needed). `SleepPhase`/`BaskPhase`/`ForagePhase`/`CoolPhase`/
`ClimbPhase`/`BurrowPhase`/`RefugePhase` below are worked examples.

`phase_condition`/`next_phase` take an internal `phase_context::NamedTuple` with fields:
`activity_signal` (function of `t`, negative during the activity window), `limits`
(`EctothermBehavioralLimits`), `active_min_hysteresis`, `activity_hysteresis`, `forcings`
(`(; sun, shade, sleep, climb, underground)`, `climb`/`underground` may be `nothing`),
`smoothing`, `night_phase`, `escape_phase`, `refuge_phase`, `cool_resume_margin`,
`cool_resume_offset`.
"""
abstract type TransientBehavioralPhase end
struct SleepPhase <: TransientBehavioralPhase end
struct BaskPhase <: TransientBehavioralPhase end
struct ForagePhase <: TransientBehavioralPhase end
struct CoolPhase <: TransientBehavioralPhase end
struct ClimbPhase <: TransientBehavioralPhase end
struct BurrowPhase <: TransientBehavioralPhase end
struct RefugePhase <: TransientBehavioralPhase end

"""
    phase_forcing(phase::TransientBehavioralPhase, forcings) -> EnvironmentForcing

`EnvironmentForcing` driving the ODE while `phase` is active. Part of the
`TransientBehavioralPhase` interface. `forcings = (; sun, shade, sleep, climb, underground)`.
"""
phase_forcing(::SleepPhase, forcings) = forcings.sleep
phase_forcing(::BaskPhase, forcings) = forcings.sun
phase_forcing(::ForagePhase, forcings) = forcings.sun
phase_forcing(::CoolPhase, forcings) = forcings.shade
phase_forcing(::ClimbPhase, forcings) = forcings.climb
phase_forcing(::BurrowPhase, forcings) = forcings.underground
phase_forcing(::RefugePhase, forcings) = forcings.underground

"Solar posture used for the heat-balance solve while `phase` is active. Part of the `TransientBehavioralPhase` interface."
phase_posture(::SleepPhase) = Intermediate()
phase_posture(::BaskPhase) = NormalToSun()
phase_posture(::ForagePhase) = Intermediate()
phase_posture(::CoolPhase) = Intermediate()
phase_posture(::ClimbPhase) = Intermediate()
phase_posture(::BurrowPhase) = Intermediate()
phase_posture(::RefugePhase) = Intermediate()

"OrganismState reported in the output while `phase` is active. Reporting only, no effect on the physics. Part of the `TransientBehavioralPhase` interface."
phase_state(::SleepPhase) = Resting()
phase_state(::BaskPhase) = Basking()
phase_state(::ForagePhase) = Active()
phase_state(::CoolPhase) = Resting()
phase_state(::ClimbPhase) = Resting()
phase_state(::BurrowPhase) = Resting()
phase_state(::RefugePhase) = Resting()

"Symbol key for `metabolic_multipliers` lookup (e.g. `(; forage=3.0)`); unlisted phases default to 1.0. Part of the `TransientBehavioralPhase` interface."
phase_key(::SleepPhase) = :sleep
phase_key(::BaskPhase) = :bask
phase_key(::ForagePhase) = :forage
phase_key(::CoolPhase) = :cool
phase_key(::ClimbPhase) = :climb
phase_key(::BurrowPhase) = :burrow
phase_key(::RefugePhase) = :refuge
_metabolic_multiplier(multipliers::NamedTuple, phase) = Float64(get(multipliers, phase_key(phase), 1.0))

# Wraps the organism's metabolic rate model so it scales by `multiplier`, via HeatExchange's
# existing "any callable works as MetabolismParameters.model" fallback - no change needed there.
function _scale_metabolism(organism::Organism, multiplier)
    base_model = HeatExchange.metabolism_pars(organism).model
    scaled_model(mass, core_temperature) = HeatExchange.metabolic_rate(base_model, mass, core_temperature) * multiplier
    @set organism.traits.heat_exchange.metabolism_pars.model = scaled_model
end
_phase_organism(organism, multiplier) = multiplier == 1.0 ? organism : _scale_metabolism(organism, multiplier)

# State is Float64 (one-lump, K) or Vector{Float64} (two-lump, [core, shell] K) - dispatch
# picks out the core temperature that drives all behavioral thresholds.
_core_temperature_of(u::Real) = u
_core_temperature_of(u::AbstractVector) = u[1]

# Conditions take plain Float64 (u in K, t in s), negative while in-phase, crossing zero
# when the phase should end (ContinuousCallback root). Sub-signals keep their own unit
# scale — min/max only need matching sign, giving AND=`min` (all signals must be >0), 
# and OR=`max` (any signal is > 0).
_basking_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.basking_temperature_min)
_active_min_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.active_temperature_min)
_active_max_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.active_temperature_max)
_emerge_signal(u, limits) = _core_temperature_of(u) - ustrip(u"K", limits.emerge_temperature_min)

_air_temperature(forcing, t, smoothing) = ustrip(u"K", forcing(t * u"s"; smoothing).air_temperature)

# Continuous analogue of is_active (thermoregulation.jl): negative in the active window.
# zenith is plain Float64 degrees; zenith=90 is sunset/sunrise.
_activity_signal(::Diurnal, zenith) = zenith - 90.0
_activity_signal(::Nocturnal, zenith) = 90.0 - zenith
_activity_signal(::Crepuscular, zenith) = abs(zenith - 90.0) - 5.0
_activity_signal(a::CombinedActivity, zenith) = minimum(p -> _activity_signal(p, zenith), a.periods)

"""
    phase_condition(phase::TransientBehavioralPhase, phase_context) -> (u, t) -> Float64

`ContinuousCallback` root function for `phase`'s exit: negative while `phase` continues,
crossing zero when it should end. Combine competing exit reasons with `max` (any one reason
ends the phase) or `min` (all must agree). Part of the `TransientBehavioralPhase` interface.
"""
# Sleep/Burrow only exit `activity_hysteresis` past the day/night boundary, not exactly at it -
# Bask/Forage/Cool/Climb's own event check (below) sits at the boundary. Any ambiguity is
# resolved by comparing competing reasons within one phase - see next_phase).
function phase_condition(::SleepPhase, phase_context)
    (u, t) -> min(-phase_context.activity_signal(t) - phase_context.activity_hysteresis, _basking_signal(u, phase_context.limits))
end
# Bask/Forage hysteresis: exit thresholds sit `active_min_hysteresis` above/below
# active_temperature_min so they never coincide (see docs).
function phase_condition(::BaskPhase, phase_context)
    (u, t) -> max(_active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis, phase_context.activity_signal(t))
end
function phase_condition(::ForagePhase, phase_context)
    (u, t) -> max(_active_max_signal(u, phase_context.limits), -_active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis, phase_context.activity_signal(t))
end
# Resume at active_temperature_min, or the candidate site's air temp + cool_resume_offset 
# if that's higher (otherwise unreachable in a warm site), falling  back to active_temperature_min
# if the site is within cool_resume_margin of active_temperature_max.
# Shared by CoolPhase (shade) and ClimbPhase (climb site).
function _cool_resume_threshold(limits, candidate_air_temperature, cool_resume_margin, cool_resume_offset)
    active_min = ustrip(u"K", limits.active_temperature_min)
    active_max = ustrip(u"K", limits.active_temperature_max)
    target = candidate_air_temperature > active_max - cool_resume_margin ? 0.0 : candidate_air_temperature + cool_resume_offset
    max(active_min, target)
end
_cool_resume_signal(phase_context, forcing, u, t) =
    _cool_resume_threshold(
        phase_context.limits, _air_temperature(forcing, t, phase_context.smoothing),
        phase_context.cool_resume_margin, phase_context.cool_resume_offset,
    ) - _core_temperature_of(u)

# -Inf when refuge_phase is nothing (unavailable), so this branch can never win a max/relative
# comparison - Cool/Climb behave exactly as before when can_retreat_underground/underground_forcing
# aren't both given.
_escape_max_signal(u, phase_context) =
    phase_context.refuge_phase === nothing ? -Inf : _core_temperature_of(u) - ustrip(u"K", phase_context.limits.escape_temperature_max)

function phase_condition(::CoolPhase, phase_context)
    (u, t) -> max(_cool_resume_signal(phase_context, phase_context.forcings.shade, u, t), phase_context.activity_signal(t), _escape_max_signal(u, phase_context))
end
function phase_condition(::ClimbPhase, phase_context)
    (u, t) -> max(_cool_resume_signal(phase_context, phase_context.forcings.climb, u, t), phase_context.activity_signal(t), _escape_max_signal(u, phase_context))
end
function phase_condition(::BurrowPhase, phase_context)
    (u, t) -> min(-phase_context.activity_signal(t) - phase_context.activity_hysteresis, _emerge_signal(u, phase_context.limits))
end
# Resumes once core temperature (not the underground site's temperature - same reasoning as
# BurrowPhase's own emergence check) drops back below active_temperature_max, offset by
# active_min_hysteresis so it doesn't land exactly on ForagePhase's own hot-exit threshold.
function phase_condition(::RefugePhase, phase_context)
    (u, t) -> (ustrip(u"K", phase_context.limits.active_temperature_max) - phase_context.active_min_hysteresis) - _core_temperature_of(u)
end

"""
    next_phase(phase::TransientBehavioralPhase, u, t, phase_context) -> TransientBehavioralPhase

Phase to transition to once `phase_condition(phase, ...)` signals exit. Route by comparing the
same signals used in `phase_condition` relatively against each other, not independently against
zero. Part of the `TransientBehavioralPhase` interface.
"""
# Sleep and Burrow play the identical downstream role (the night phase always hands off to
# Bask unconditionally at its own exit) - one method covers both.
next_phase(::Union{SleepPhase,BurrowPhase}, u, t, phase_context) = BaskPhase()

# Compare competing exit signals to each other rather than testing each against zero
# independently (same principle as the ForagePhase fix below). Root finding may stop a
# few floating-point representable values (ULPs) from the exact threshold, so relative 
# comparisons are more robust.
function next_phase(::BaskPhase, u, t, phase_context)
    day = phase_context.activity_signal(t)
    warm = _active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis
    day >= warm ? phase_context.night_phase : ForagePhase()
end
function next_phase(::ForagePhase, u, t, phase_context)
    day = phase_context.activity_signal(t)
    hot = _active_max_signal(u, phase_context.limits)
    cold = -_active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis
    if day >= max(hot, cold)
        phase_context.night_phase
    elseif hot >= cold
        phase_context.escape_phase
    else
        BaskPhase()
    end
end
# Cool and Climb additionally escalate to RefugePhase if shade/climbing alone isn't enough
# (core temperature keeps climbing toward escape_temperature_max) - a third exit reason
# compared relatively against the other two, same principle as ForagePhase's hot/cold branches.
function next_phase(::CoolPhase, u, t, phase_context)
    day = phase_context.activity_signal(t)
    resume = _cool_resume_signal(phase_context, phase_context.forcings.shade, u, t)
    escape = _escape_max_signal(u, phase_context)
    if escape >= max(day, resume)
        phase_context.refuge_phase
    elseif day >= resume
        phase_context.night_phase
    else
        ForagePhase()
    end
end
function next_phase(::ClimbPhase, u, t, phase_context)
    day = phase_context.activity_signal(t)
    resume = _cool_resume_signal(phase_context, phase_context.forcings.climb, u, t)
    escape = _escape_max_signal(u, phase_context)
    if escape >= max(day, resume)
        phase_context.refuge_phase
    elseif day >= resume
        phase_context.night_phase
    else
        ForagePhase()
    end
end
# RefugePhase only has one exit reason (cooled down enough), so it always resumes straight to
# ForagePhase - the animal was already active before fleeing the heat, no warm-up phase needed.
next_phase(::RefugePhase, u, t, phase_context) = ForagePhase()

# RHS dispatch: one-lump (scalar state, HeatExchange.onelump) vs two-lump (2-vector
# [core, shell], HeatExchange.twolump). Both methods take the same keywords so callers
# don't need to know which model is in play - `shell_thickness` is simply unused for onelump.
function _body_temperature_rate(
    u::Real, t, organism::Organism, environment_pars, environment_vars;
    shell_thickness=nothing, posture, surface_solve, smoothing,
)
    ustrip(u"K/s", HeatExchange.onelump(
        u * u"K", t * u"s", organism, (; environment_pars, environment_vars); posture, smoothing,
    ).core_temperature_rate)
end
function _body_temperature_rate(
    u::AbstractVector, t, organism::Organism, environment_pars, environment_vars;
    shell_thickness, posture, surface_solve, smoothing,
)
    out = HeatExchange.twolump(
        (core_temperature=u[1] * u"K", shell_temperature=u[2] * u"K"), t * u"s", organism, (; environment_pars, environment_vars);
        shell_thickness, posture, surface_solve, smoothing,
    )
    Float64[ustrip(u"K/s", out.core_temperature_rate), ustrip(u"K/s", out.shell_temperature_rate)]
end

_package_trajectory(all_t, all_u::Vector{Float64}, all_state, all_phase) =
    (; t=all_t * u"s", core_temperature=all_u * u"K", state=all_state, phase=all_phase)
function _package_trajectory(all_t, all_u::Vector{Vector{Float64}}, all_state, all_phase)
    (;
        t=all_t * u"s",
        core_temperature=[uu[1] for uu in all_u] * u"K",
        shell_temperature=[uu[2] for uu in all_u] * u"K",
        state=all_state,
        phase=all_phase,
    )
end

# Catches the class of bug fixed in BaskPhase/ForagePhase (two phases sharing one exact
# threshold, causing a zero-duration stall) before it can resurface via a misconfigured
# EctothermBehavioralLimits. The base ordering is shared with the steady-state loop
# (_validate_ectotherm_thresholds); this adds the extra active_min_hysteresis margin the
# transient driver's root-finding needs on top of it.
function _validate_thresholds(limits, active_min_hysteresis)
    _validate_ectotherm_thresholds(limits)
    limits.basking_temperature_min + active_min_hysteresis * u"K" < limits.active_temperature_min ||
        throw(ArgumentError("basking_temperature_min must be below active_temperature_min by more than active_min_hysteresis"))
    return nothing
end

"""
    simulate_transient_behavior(times, 
                               core_temperature_init, 
                               organism::Organism, 
                               environment_pars,
                               sun_forcing::EnvironmentForcing, 
                               shade_forcing::EnvironmentForcing,
                               limits::EctothermBehavioralLimits;
                               solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;), smoothing=HardBound(), 
                               max_bouts=100*length(times),
                               bout_chunk=3600.0, 
                               active_min_hysteresis=0.15u"K",
                               climb_forcing=nothing, 
                               underground_forcing=nothing,
                               sleep_forcing=shade_forcing, 
                               initial_phase=nothing,
                               activity_period=Diurnal(), 
                               activity_hysteresis=0.1u"°",
                               metabolic_multipliers=NamedTuple(),
                               cool_resume_margin=2.0u"K", 
                               cool_resume_offset=1.0u"K"
                               )
    simulate_transient_behavior(times, 
                               (; core_temperature, shell_temperature), 
                               organism, 
                               environment_pars,
                               sun_forcing, 
                               shade_forcing, 
                               limits; 
                               shell_thickness,
                               surface_solve=LinearizedSurface(), solver=OrdinaryDiffEqTsit5.Tsit5(),
                               solver_kwargs=(;), smoothing=HardBound(), 
                               max_bouts=100*length(times),
                               bout_chunk=3600.0, 
                               active_min_hysteresis=0.15u"K",
                               climb_forcing=nothing, 
                               underground_forcing=nothing,
                               sleep_forcing=shade_forcing, 
                               initial_phase=nothing,
                               activity_period=Diurnal(), 
                               activity_hysteresis=0.1u"°",
                               metabolic_multipliers=NamedTuple(),
                               cool_resume_margin=2.0u"K", 
                               cool_resume_offset=1.0u"K"
                               )

Event-driven diurnal behavioral thermoregulation for an organism with a transient
heat budget (thermal inertia).

Behavior is represented as a sequence of phases:

- `SleepPhase` (or `BurrowPhase`)
- `BaskPhase`
- `ForagePhase`
- `CoolPhase` (or `ClimbPhase`)
- `RefugePhase`

The model switches between phases using continuous event detection rather than fixed
time steps. This is a Julia port of NicheMapR's `trans_behav.R`.

# Heat budget

Dispatches on the initial thermal state:

- `initial_temperature::Quantity` → one-lump model (`HeatExchange.onelump`),
  suitable for small animals (≈ <5 kg).
- `(; core_temperature, shell_temperature)` → two-lump model
  (`HeatExchange.twolump`; requires `shell_thickness`), suitable for larger
  animals or shelled organisms.

# Habitat options

Providing additional environmental forcings enables alternative behavioural phases:

- `underground_forcing` + `limits.can_retreat_underground`
  enables `BurrowPhase`.
- `climb_forcing` + `limits.can_climb`
  enables `ClimbPhase`.

`sleep_forcing` specifies the inactive location (defaults to `shade_forcing`).

# Activity

`activity_period` (`Diurnal()` by default) defines the active window
(`BaskPhase`/`ForagePhase`) and inactive window (`SleepPhase`/`BurrowPhase`).
Supported values are `Diurnal`, `Nocturnal`, `Crepuscular`, and
`CombinedActivity(...)`.

`initial_phase` overrides the default day/night starting phase.

# Temperature limits

Behavior is controlled by the thresholds in `limits`.

Thresholds are validated on entry to ensure

    escape_min < emerge_min < bask_min < active_min < active_max < escape_max

with the required hysteresis around `active_min`.

`active_min_hysteresis` prevents oscillation between `BaskPhase` and
`ForagePhase` when body temperature sits exactly at `active_temperature_min`.

`activity_hysteresis` similarly separates day/night transitions.

# Metabolism

`metabolic_multipliers` scales metabolic heat production independently for each
phase (e.g. `(; forage=3, bask=1.5)`).

# Cool-phase resume

`cool_resume_margin` and `cool_resume_offset` determine when
`CoolPhase`/`ClimbPhase` return to `ForagePhase`.

# Returns

A `NamedTuple` containing

- `t`
- `core_temperature`
- `shell_temperature` (two-lump model only)
- `state`
- `phase`

with one entry per accepted solver step.
"""
function simulate_transient_behavior(
    times, core_temperature_init::Unitful.Quantity, organism::Organism, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), max_bouts=100 * length(times),
    bout_chunk=3600.0, active_min_hysteresis=0.15u"K",
    climb_forcing::Union{Nothing,EnvironmentForcing}=nothing,
    underground_forcing::Union{Nothing,EnvironmentForcing}=nothing,
    sleep_forcing::EnvironmentForcing=shade_forcing,
    initial_phase::Union{Nothing,TransientBehavioralPhase}=nothing,
    activity_period::ActivityPeriod=Diurnal(),
    activity_hysteresis=0.1u"°",
    metabolic_multipliers::NamedTuple=NamedTuple(),
    cool_resume_margin=2.0u"K",
    cool_resume_offset=1.0u"K",
)
    _simulate_transient_behavior(
        times, ustrip(u"K", core_temperature_init), organism, environment_pars, sun_forcing, shade_forcing, limits;
        shell_thickness=nothing, surface_solve=HeatExchange.LinearizedSurface(),
        solver, solver_kwargs, smoothing, max_bouts, bout_chunk,
        active_min_hysteresis=ustrip(u"K", active_min_hysteresis),
        climb_forcing, underground_forcing, sleep_forcing, initial_phase, activity_period,
        activity_hysteresis=ustrip(u"°", activity_hysteresis),
        metabolic_multipliers,
        cool_resume_margin=ustrip(u"K", cool_resume_margin),
        cool_resume_offset=ustrip(u"K", cool_resume_offset),
    )
end
function simulate_transient_behavior(
    times, state::NamedTuple{(:core_temperature, :shell_temperature)}, organism::Organism, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    shell_thickness, surface_solve::SurfaceSolveStrategy=HeatExchange.LinearizedSurface(),
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), max_bouts=100 * length(times),
    bout_chunk=3600.0, active_min_hysteresis=0.15u"K",
    climb_forcing::Union{Nothing,EnvironmentForcing}=nothing,
    underground_forcing::Union{Nothing,EnvironmentForcing}=nothing,
    sleep_forcing::EnvironmentForcing=shade_forcing,
    initial_phase::Union{Nothing,TransientBehavioralPhase}=nothing,
    activity_period::ActivityPeriod=Diurnal(),
    activity_hysteresis=0.1u"°",
    metabolic_multipliers::NamedTuple=NamedTuple(),
    cool_resume_margin=2.0u"K",
    cool_resume_offset=1.0u"K",
)
    u0 = Float64[ustrip(u"K", state.core_temperature), ustrip(u"K", state.shell_temperature)]
    _simulate_transient_behavior(
        times, u0, organism, environment_pars, sun_forcing, shade_forcing, limits;
        shell_thickness, surface_solve,
        solver, solver_kwargs, smoothing, max_bouts, bout_chunk,
        active_min_hysteresis=ustrip(u"K", active_min_hysteresis),
        climb_forcing, underground_forcing, sleep_forcing, initial_phase, activity_period,
        activity_hysteresis=ustrip(u"°", activity_hysteresis),
        metabolic_multipliers,
        cool_resume_margin=ustrip(u"K", cool_resume_margin),
        cool_resume_offset=ustrip(u"K", cool_resume_offset),
    )
end

function _simulate_transient_behavior(
    times, u0, organism::Organism, environment_pars,
    sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
    limits::EctothermBehavioralLimits;
    shell_thickness, surface_solve,
    solver, solver_kwargs, smoothing::SmoothingStrategy, max_bouts, bout_chunk, active_min_hysteresis,
    climb_forcing, underground_forcing, sleep_forcing, initial_phase, activity_period, activity_hysteresis,
    metabolic_multipliers, cool_resume_margin, cool_resume_offset,
)
    _validate_thresholds(limits, active_min_hysteresis)

    t0, tend = ustrip(u"s", first(times)), ustrip(u"s", last(times))
    activity_signal(t) = _activity_signal(activity_period, ustrip(u"°", sun_forcing(t * u"s"; smoothing).zenith_angle))

    has_climb = limits.can_climb && climb_forcing !== nothing
    has_burrow = limits.can_retreat_underground && underground_forcing !== nothing
    night_phase = has_burrow ? BurrowPhase() : SleepPhase()
    escape_phase = has_climb ? ClimbPhase() : CoolPhase()
    refuge_phase = has_burrow ? RefugePhase() : nothing

    if initial_phase isa ClimbPhase && !has_climb
        throw(ArgumentError("initial_phase=ClimbPhase() requires limits.can_climb=true and a non-nothing climb_forcing"))
    end
    if initial_phase isa BurrowPhase && !has_burrow
        throw(ArgumentError("initial_phase=BurrowPhase() requires limits.can_retreat_underground=true and a non-nothing underground_forcing"))
    end
    if initial_phase isa RefugePhase && !has_burrow
        throw(ArgumentError("initial_phase=RefugePhase() requires limits.can_retreat_underground=true and a non-nothing underground_forcing"))
    end

    forcings = (; sun=sun_forcing, shade=shade_forcing, sleep=sleep_forcing, climb=climb_forcing, underground=underground_forcing)
    phase_context = (;
        activity_signal, limits, active_min_hysteresis, activity_hysteresis, forcings, smoothing,
        night_phase, escape_phase, refuge_phase, cool_resume_margin, cool_resume_offset,
    )

    phase = initial_phase !== nothing ? initial_phase : (activity_signal(t0) < 0.0 ? ForagePhase() : night_phase)
    u_current = u0
    t_current = t0

    all_t = Float64[]
    all_u = typeof(u0)[]
    all_state = OrganismState[]
    all_phase = TransientBehavioralPhase[]

    for _ in 1:max_bouts
        t_current >= tend && break

        for _ in 1:8
            condition = phase_condition(phase, phase_context)
            condition(u_current, t_current) < 0 && break
            phase = next_phase(phase, u_current, t_current, phase_context)
        end

        forcing = phase_forcing(phase, forcings)
        posture = phase_posture(phase)
        bout_organism = _phase_organism(organism, _metabolic_multiplier(metabolic_multipliers, phase))
        f = (u, _, t) -> _body_temperature_rate(
            u, t, bout_organism, environment_pars, forcing(t * u"s"; smoothing);
            shell_thickness, posture, surface_solve, smoothing,
        )
        condition = phase_condition(phase, phase_context)
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
        append!(all_state, fill(phase_state(phase), length(sol.t)))
        append!(all_phase, fill(phase, length(sol.t)))

        t_next = sol.t[end]
        u_next = sol.u[end]
        if t_next < t_chunk_end
            phase = next_phase(phase, u_next, t_next, phase_context)
        end
        t_current, u_current = t_next, u_next
    end

    return _package_trajectory(all_t, all_u, all_state, all_phase)
end
