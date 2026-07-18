# Event-driven diurnal behavior for the transient body-temperature models. Port of NicheMapR's
# trans_behav.R (sleep -> bask -> forage <-> cool), using ContinuousCallback in place of
# R's rootfun events. Simplifications relative to R noted in docs/transient_body_temperature.md.

abstract type TransientBehavioralPhase end
struct SleepPhase <: TransientBehavioralPhase end
struct BaskPhase <: TransientBehavioralPhase end
struct ForagePhase <: TransientBehavioralPhase end
struct CoolPhase <: TransientBehavioralPhase end
struct ClimbPhase <: TransientBehavioralPhase end
struct BurrowPhase <: TransientBehavioralPhase end
struct RefugePhase <: TransientBehavioralPhase end

# forcings = (; sun, shade, sleep, climb, underground) - climb/underground may be `nothing`,
# but are only ever looked up when ClimbPhase/BurrowPhase is reachable (see has_climb/has_burrow).
_phase_forcing(::SleepPhase, forcings) = forcings.sleep
_phase_forcing(::BaskPhase, forcings) = forcings.sun
_phase_forcing(::ForagePhase, forcings) = forcings.sun
_phase_forcing(::CoolPhase, forcings) = forcings.shade
_phase_forcing(::ClimbPhase, forcings) = forcings.climb
_phase_forcing(::BurrowPhase, forcings) = forcings.underground
_phase_forcing(::RefugePhase, forcings) = forcings.underground

_phase_posture(::SleepPhase) = Intermediate()
_phase_posture(::BaskPhase) = NormalToSun()
_phase_posture(::ForagePhase) = Intermediate()
_phase_posture(::CoolPhase) = Intermediate()
_phase_posture(::ClimbPhase) = Intermediate()
_phase_posture(::BurrowPhase) = Intermediate()
_phase_posture(::RefugePhase) = Intermediate()

_phase_state(::SleepPhase) = Resting()
_phase_state(::BaskPhase) = Basking()
_phase_state(::ForagePhase) = Active()
_phase_state(::CoolPhase) = Resting()
_phase_state(::ClimbPhase) = Resting()
_phase_state(::BurrowPhase) = Resting()
_phase_state(::RefugePhase) = Resting()

# Keys into the user-supplied `metabolic_multipliers` NamedTuple (e.g. `(; forage=3.0, climb=2.5)`).
# Phases not present in that NamedTuple default to a multiplier of 1.0 (no scaling).
_phase_key(::SleepPhase) = :sleep
_phase_key(::BaskPhase) = :bask
_phase_key(::ForagePhase) = :forage
_phase_key(::CoolPhase) = :cool
_phase_key(::ClimbPhase) = :climb
_phase_key(::BurrowPhase) = :burrow
_phase_key(::RefugePhase) = :refuge
_metabolic_multiplier(multipliers::NamedTuple, phase) = Float64(get(multipliers, _phase_key(phase), 1.0))

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
# scale — min/max only need matching sign, giving AND=`min`, OR=`max`.
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

# phase_context = (; activity_signal, limits, active_min_hysteresis, activity_hysteresis,
#                    forcings, smoothing, night_phase, escape_phase, refuge_phase,
#                    cool_resume_margin, cool_resume_offset)
# Sleep/Burrow only exit `activity_hysteresis` past the day/night boundary, not exactly at it -
# Bask/Forage/Cool/Climb's own day-arrived check (below) sits right at the boundary, so this
# is genuine separation between two different phases' thresholds, not a self-comparison (that
# ambiguity is instead resolved by comparing competing reasons within one phase - see _next_phase).
function _phase_condition(::SleepPhase, phase_context)
    (u, t) -> min(-phase_context.activity_signal(t) - phase_context.activity_hysteresis, _basking_signal(u, phase_context.limits))
end
# Bask/Forage hysteresis: exit thresholds sit `active_min_hysteresis` above/below
# active_temperature_min so they never coincide (see docs).
function _phase_condition(::BaskPhase, phase_context)
    (u, t) -> max(_active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis, phase_context.activity_signal(t))
end
function _phase_condition(::ForagePhase, phase_context)
    (u, t) -> max(_active_max_signal(u, phase_context.limits), -_active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis, phase_context.activity_signal(t))
end
# trans_behav.R's `forage` event: resume at active_temperature_min, or the candidate site's air
# temp + cool_resume_offset if that's higher (otherwise unreachable in a warm site), falling
# back to active_temperature_min if the site is within cool_resume_margin of
# active_temperature_max. Shared by CoolPhase (shade) and ClimbPhase (climb site).
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

function _phase_condition(::CoolPhase, phase_context)
    (u, t) -> max(_cool_resume_signal(phase_context, phase_context.forcings.shade, u, t), phase_context.activity_signal(t), _escape_max_signal(u, phase_context))
end
function _phase_condition(::ClimbPhase, phase_context)
    (u, t) -> max(_cool_resume_signal(phase_context, phase_context.forcings.climb, u, t), phase_context.activity_signal(t), _escape_max_signal(u, phase_context))
end
function _phase_condition(::BurrowPhase, phase_context)
    (u, t) -> min(-phase_context.activity_signal(t) - phase_context.activity_hysteresis, _emerge_signal(u, phase_context.limits))
end
# Resumes once core temperature (not the underground site's temperature - same reasoning as
# BurrowPhase's own emergence check) drops back below active_temperature_max, offset by
# active_min_hysteresis so it doesn't land exactly on ForagePhase's own hot-exit threshold.
function _phase_condition(::RefugePhase, phase_context)
    (u, t) -> (ustrip(u"K", phase_context.limits.active_temperature_max) - phase_context.active_min_hysteresis) - _core_temperature_of(u)
end

# Sleep and Burrow play the identical downstream role (the night phase always hands off to
# Bask unconditionally at its own exit) - one method covers both.
_next_phase(::Union{SleepPhase,BurrowPhase}, u, t, phase_context) = BaskPhase()

# Routing compares each phase's own competing exit reasons to each other, not any one of them
# to zero independently (same principle as the ForagePhase fix below) - root-finding can land
# a root a few ULPs off an exact threshold, which an independent >=0 check can misread, but a
# relative comparison against a clearly-negative-or-positive competitor can't.
function _next_phase(::BaskPhase, u, t, phase_context)
    day = phase_context.activity_signal(t)
    warm = _active_min_signal(u, phase_context.limits) - phase_context.active_min_hysteresis
    day >= warm ? phase_context.night_phase : ForagePhase()
end
function _next_phase(::ForagePhase, u, t, phase_context)
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
function _next_phase(::CoolPhase, u, t, phase_context)
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
function _next_phase(::ClimbPhase, u, t, phase_context)
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
_next_phase(::RefugePhase, u, t, phase_context) = ForagePhase()

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
    simulate_diurnal_behavior(times, core_temperature_init, organism::Organism, environment_pars,
                               sun_forcing::EnvironmentForcing, shade_forcing::EnvironmentForcing,
                               limits::EctothermBehavioralLimits;
                               solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
                               smoothing=HardBound(), max_bouts=100*length(times),
                               bout_chunk=3600.0, active_min_hysteresis=0.15u"K",
                               climb_forcing=nothing, underground_forcing=nothing,
                               sleep_forcing=shade_forcing, initial_phase=nothing,
                               activity_period=Diurnal(), activity_hysteresis=0.1u"°",
                               metabolic_multipliers=NamedTuple(),
                               cool_resume_margin=2.0u"K", cool_resume_offset=1.0u"K")
    simulate_diurnal_behavior(times, (; core_temperature, shell_temperature), organism, environment_pars,
                               sun_forcing, shade_forcing, limits; shell_thickness,
                               surface_solve=LinearizedSurface(), solver=OrdinaryDiffEqTsit5.Tsit5(),
                               solver_kwargs=(;), smoothing=HardBound(), max_bouts=100*length(times),
                               bout_chunk=3600.0, active_min_hysteresis=0.15u"K",
                               climb_forcing=nothing, underground_forcing=nothing,
                               sleep_forcing=shade_forcing, initial_phase=nothing,
                               activity_period=Diurnal(), activity_hysteresis=0.1u"°",
                               metabolic_multipliers=NamedTuple(),
                               cool_resume_margin=2.0u"K", cool_resume_offset=1.0u"K")

Event-driven diurnal behavioral thermoregulation. Seven phases: `SleepPhase`/`BurrowPhase`
(inactive site), `BaskPhase`, `ForagePhase`, `CoolPhase`/`ClimbPhase` (too-hot escape site),
`RefugePhase` (deeper too-hot escalation from `CoolPhase`/`ClimbPhase`).
Dispatches on the initial state: a plain temperature runs the one-lump model
(`HeatExchange.onelump`); a `(; core_temperature, shell_temperature)` NamedTuple runs the
two-lump model (`HeatExchange.twolump`, needs a `shell_thickness` keyword). Port of NicheMapR's
`trans_behav.R` (see file banner for what's simplified relative to R).

`BurrowPhase` substitutes for `SleepPhase` whenever `limits.can_retreat_underground` and
`underground_forcing` is given; `ClimbPhase` substitutes for `CoolPhase` whenever
`limits.can_climb` and `climb_forcing` is given — unconditional substitution, not a competing
choice. `sleep_forcing` sets the inactive-period site (default `shade_forcing`; pass
`sun_forcing` for an unshaded sleeper). `initial_phase` overrides the default day/night
starting logic; passing `ClimbPhase()`/`BurrowPhase()` without the matching capability +
forcing throws `ArgumentError`. `activity_period` (`Diurnal()` default, or `Nocturnal()`/
`Crepuscular()`/`CombinedActivity(...)`) sets which part of the day counts as the active
window (`Bask`/`Forage`) vs. inactive (`night_phase`) — the continuous-time analogue of
`is_active` (`thermoregulation.jl`); `ResponsiveActivity` is not supported.

`active_min_hysteresis` separates Bask's and Forage's exit thresholds around
`active_temperature_min` (± 0.15 K by default), preventing an indefinite zero-duration
oscillation if core temperature settles exactly on that boundary.
`limits`' thresholds are validated at call time (`escape_temperature_min <
emerge_temperature_min < basking_temperature_min < active_temperature_min <
active_temperature_max < escape_temperature_max`, with the `active_min_hysteresis` margin
below `active_temperature_min`) to catch the same class of bug up front.

Two complementary fixes protect the day/night boundary: `activity_hysteresis` (default
`0.1u"°"`) separates `SleepPhase`/`BurrowPhase`'s exit from `BaskPhase`/`ForagePhase`/
`CoolPhase`/`ClimbPhase`'s day-arrived check — genuinely different phases sharing one exact
`activity_signal=0` threshold, the same failure class as `active_min_hysteresis` above. Within
each of those phases, day-arrived is *also* compared against that phase's own temperature exit
reason (not against zero independently), the same relative-comparison principle already used
for `ForagePhase`'s hot/cold branches — this handles one phase's own routing ambiguity, which a
margin alone can't fix.

`limits.escape_temperature_max` also gates `CoolPhase`/`ClimbPhase`'s further escalation to
`RefugePhase` (underground) if shade/climbing alone can't keep core temperature down, compared
relatively against those phases' other exit reasons for the same reason as above. `RefugePhase`
resumes to `ForagePhase` once core temperature drops back to `active_temperature_max -
active_min_hysteresis`.

`metabolic_multipliers` scales metabolic heat production per phase, e.g.
`(; forage=3.0, climb=2.5, bask=1.5)` — unlisted phases (keys: `sleep`, `bask`, `forage`, `cool`,
`climb`, `burrow`, `refuge`) default to `1.0` (no scaling). Applied by wrapping the organism's
metabolic rate model in a closure for the duration of each bout (constant while a phase is
active, since phase never changes mid-bout), so it composes with any `MetabolicRateEquation`.

`cool_resume_margin` (default `2.0u"K"`) and `cool_resume_offset` (default `1.0u"K"`) tune
`CoolPhase`/`ClimbPhase`'s resume threshold: resume at `active_temperature_min`, or the
candidate site's air temperature + `cool_resume_offset` if that's higher — unless the site's
air temperature is within `cool_resume_margin` of `active_temperature_max`, in which case the
site is too warm for that offset to be meaningful and the threshold falls back to
`active_temperature_min`.

# Returns
NamedTuple with `t` (s), `core_temperature` (K) [, `shell_temperature` (K) for the two-lump
method], `state` (`Vector{OrganismState}`), `phase` (`Vector{TransientBehavioralPhase}`), one
entry per accepted solver step across all bouts.
"""
function simulate_diurnal_behavior(
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
    _simulate_diurnal_behavior(
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
function simulate_diurnal_behavior(
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
    _simulate_diurnal_behavior(
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

function _simulate_diurnal_behavior(
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
            condition = _phase_condition(phase, phase_context)
            condition(u_current, t_current) < 0 && break
            phase = _next_phase(phase, u_current, t_current, phase_context)
        end

        forcing = _phase_forcing(phase, forcings)
        posture = _phase_posture(phase)
        bout_organism = _phase_organism(organism, _metabolic_multiplier(metabolic_multipliers, phase))
        f = (u, _, t) -> _body_temperature_rate(
            u, t, bout_organism, environment_pars, forcing(t * u"s"; smoothing);
            shell_thickness, posture, surface_solve, smoothing,
        )
        condition = _phase_condition(phase, phase_context)
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
        append!(all_phase, fill(phase, length(sol.t)))

        t_next = sol.t[end]
        u_next = sol.u[end]
        if t_next < t_chunk_end
            phase = _next_phase(phase, u_next, t_next, phase_context)
        end
        t_current, u_current = t_next, u_next
    end

    return _package_trajectory(all_t, all_u, all_state, all_phase)
end
