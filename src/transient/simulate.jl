"""
    simulate_onelump(times, core_temperature_init, organism::Organism, environment_pars, forcing::EnvironmentForcing;
                      solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;), smoothing=HardBound(), kw...)

Integrate `HeatExchange.onelump`'s derivative forward through time-varying environmental
forcing, holding `organism`'s thermoregulatory effectors fixed. `kw...` is forwarded to
`onelump` — pass `metabolic_heat_flow` there for insulated bodies (a `Quantity` or
`Function(core_temperature)`); `Naked` bodies need no extra keywords.

# Returns
NamedTuple with `t` (s), `core_temperature` (K), `core_temperature_rate` (K/s),
`skin_temperature` (K), `diagnostics` (one `onelump` output per saved timestep,
reconstructed post-hoc), and the raw `solution` (`ODESolution`, Unitful-stripped).
"""
function simulate_onelump(
    times, core_temperature_init, organism::Organism, environment_pars, forcing::EnvironmentForcing;
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(), kw...,
)
    u0 = ustrip(u"K", core_temperature_init)
    tspan = (ustrip(u"s", first(times)), ustrip(u"s", last(times)))
    f = (u, _, t) -> ustrip(u"K/s", HeatExchange.onelump(
        u * u"K", t * u"s", organism, (; environment_pars, environment_vars=forcing(t * u"s"; smoothing)); smoothing, kw...,
    ).core_temperature_rate)
    problem = OrdinaryDiffEqTsit5.ODEProblem(f, u0, tspan)
    solution = OrdinaryDiffEqTsit5.solve(problem, solver; saveat=ustrip.(u"s", times), solver_kwargs...)
    core_temperature = solution.u * u"K"

    diagnostics = [
        HeatExchange.onelump(
            Tc, ti * u"s", organism, (; environment_pars, environment_vars=forcing(ti * u"s"; smoothing)); smoothing, kw...,
        )
        for (Tc, ti) in zip(core_temperature, solution.t)
    ]
    core_temperature_rate = [ustrip(u"K/s", d.core_temperature_rate) for d in diagnostics] * u"K/s"
    skin_temperature = [d.skin_temperature for d in diagnostics]

    return (; t=solution.t * u"s", core_temperature, core_temperature_rate, skin_temperature, diagnostics, solution)
end

"""
    simulate_twolump(times, core_temperature_init, shell_temperature_init, organism::Organism, environment_pars,
                      forcing::EnvironmentForcing; shell_thickness, posture=radiation_pars(organism).solar_orientation,
                      surface_solve=LinearizedSurface(), solver=OrdinaryDiffEqTsit5.Tsit5(),
                      solver_kwargs=(;), smoothing=HardBound())

Integrate `HeatExchange.twolump`'s derivatives forward through time-varying environmental
forcing.

# Returns
NamedTuple with `t` (s), `core_temperature`, `shell_temperature` (K), and the raw
`solution` (`ODESolution`, Unitful-stripped, state `[core_temperature, shell_temperature]`)
for callers who need dense interpolation or `retcode`.
"""
function simulate_twolump(
    times, core_temperature_init, shell_temperature_init, organism::Organism, environment_pars, forcing::EnvironmentForcing;
    shell_thickness, posture=radiation_pars(organism).solar_orientation,
    surface_solve::SurfaceSolveStrategy=HeatExchange.LinearizedSurface(),
    solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(),
)
    u0 = [ustrip(u"K", core_temperature_init), ustrip(u"K", shell_temperature_init)]
    tspan = (ustrip(u"s", first(times)), ustrip(u"s", last(times)))
    f = (u, _, t) -> begin
        out = HeatExchange.twolump(
            (core_temperature=u[1] * u"K", shell_temperature=u[2] * u"K"), t * u"s", organism,
            (; environment_pars, environment_vars=forcing(t * u"s"; smoothing));
            shell_thickness, posture, surface_solve, smoothing,
        )
        [ustrip(u"K/s", out.core_temperature_rate), ustrip(u"K/s", out.shell_temperature_rate)]
    end
    problem = OrdinaryDiffEqTsit5.ODEProblem(f, u0, tspan)
    solution = OrdinaryDiffEqTsit5.solve(problem, solver; saveat=ustrip.(u"s", times), solver_kwargs...)
    return (;
        t=solution.t * u"s",
        core_temperature=[u[1] for u in solution.u] * u"K",
        shell_temperature=[u[2] for u in solution.u] * u"K",
        solution,
    )
end
