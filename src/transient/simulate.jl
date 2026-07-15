"""
    simulate_onelump(times, core_temperature_init, body, environment_pars, forcing::EnvironmentForcing;
                      internal_conduction, posture=Intermediate(),
                      body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
                      metabolic_heat_volumetric, solver=OrdinaryDiffEqTsit5.Tsit5(),
                      solver_kwargs=(;), smoothing=HardBound())

Integrate `HeatExchange.onelump`'s derivative form forward through time-varying
environmental forcing.

# Returns
NamedTuple with `t` (s), `core_temperature` (K), and the raw `solution` (`ODESolution`,
Unitful-stripped) for callers who need dense interpolation or `retcode`.
"""
function simulate_onelump(
    times, core_temperature_init, body, environment_pars, forcing::EnvironmentForcing;
    internal_conduction, posture=Intermediate(),
    body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    metabolic_heat_volumetric, solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(),
)
    u0 = ustrip(u"K", core_temperature_init)
    tspan = (ustrip(u"s", first(times)), ustrip(u"s", last(times)))
    f = (u, _, t) -> ustrip(u"K/s", HeatExchange.onelump(
        u * u"K", t * u"s", body, environment_pars, forcing(t * u"s"; smoothing);
        internal_conduction, posture, body_absorptivity, emissivity,
        sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
    ))
    problem = OrdinaryDiffEqTsit5.ODEProblem(f, u0, tspan)
    solution = OrdinaryDiffEqTsit5.solve(problem, solver; saveat=ustrip.(u"s", times), solver_kwargs...)
    return (; t=solution.t * u"s", core_temperature=solution.u * u"K", solution)
end

"""
    simulate_twolump(times, core_temperature_init, shell_temperature_init, body, environment_pars,
                      forcing::EnvironmentForcing; internal_conduction, shell_thickness,
                      posture=Intermediate(), body_absorptivity, emissivity, sky_view_factor,
                      ground_view_factor, metabolic_heat_volumetric, solver=OrdinaryDiffEqTsit5.Tsit5(),
                      solver_kwargs=(;), smoothing=HardBound())

Integrate `HeatExchange.twolump`'s derivatives forward through time-varying
environmental forcing.

# Returns
NamedTuple with `t` (s), `core_temperature`, `shell_temperature` (K), and the raw
`solution` (`ODESolution`, Unitful-stripped, state `[core_temperature, shell_temperature]`)
for callers who need dense interpolation or `retcode`.
"""
function simulate_twolump(
    times, core_temperature_init, shell_temperature_init, body, environment_pars, forcing::EnvironmentForcing;
    internal_conduction, shell_thickness, posture=Intermediate(),
    body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    metabolic_heat_volumetric, solver=OrdinaryDiffEqTsit5.Tsit5(), solver_kwargs=(;),
    smoothing::SmoothingStrategy=HeatExchange.HardBound(),
)
    u0 = [ustrip(u"K", core_temperature_init), ustrip(u"K", shell_temperature_init)]
    tspan = (ustrip(u"s", first(times)), ustrip(u"s", last(times)))
    f = (u, _, t) -> begin
        out = HeatExchange.twolump(
            (core_temperature=u[1] * u"K", shell_temperature=u[2] * u"K"), t * u"s", body, environment_pars, forcing(t * u"s"; smoothing);
            internal_conduction, shell_thickness, posture, body_absorptivity, emissivity,
            sky_view_factor, ground_view_factor, metabolic_heat_volumetric, smoothing,
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
