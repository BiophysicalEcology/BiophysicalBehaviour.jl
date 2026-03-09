"""
    torpor(organism, environment, T_skin, T_insulation, limits) → (endotherm_out, organism)

Iteratively lower core temperature to find torpor equilibrium.

Ports NicheMapR SOLVENDO.f torpor loop (TORPOR=1 path, lines 204–519).

## Algorithm

1. Solve heat balance at the current (normothermic) T_core.
2. Compute Q_basal = Q_minimum_ref × q10^((T_core − T_core_ref) / 10).
3. If Q_gen ≤ Q_basal, return immediately — no torpor entry needed.
4. Otherwise, iteratively decrement T_core by `limits.torpid.T_core_step`, re-scaling
   Q_basal via Q10 and re-solving the heat balance at each step, until any stopping
   condition is met:
   - **Equilibrium** (`Q_gen ≤ Q_basal`): metabolic demand has dropped to the Q10-scaled
     minimum — torpor depth is sufficient.
   - **Floor** (`T_core ≤ T_core_min`): minimum torpor temperature reached.
   - **Environment too warm** (`Q_gen / Q_basal < −torptol`): Q_gen has gone negative,
     meaning the environment supplies more heat than the animal needs even at minimum MR;
     revert one step and exit (mirrors NicheMapR RESPFN/RESPGEN < −TORPTOL check).

## Returns

`(endotherm_out, organism_torpid)` — the heat-balance result and organism at the torpor
equilibrium T_core. If no torpor was entered, both reflect the normothermic state.
"""
function torpor(
    organism::Organism,
    environment::NamedTuple,
    T_skin,
    T_insulation,
    limits::ThermoregulationLimits,
)
    torpid_limits  = limits.torpid
    T_core_min     = torpid_limits.T_core_min
    T_core_step    = torpid_limits.T_core_step
    torptol        = torpid_limits.torptol
    Q_basal_ref    = limits.Q_minimum_ref
    T_core_ref     = limits.T_core.reference   # normothermic reference (same as hyperthermia())
    q10            = HeatExchange.metabolismpars(organism).q10

    # Initial solve at current (normothermic) T_core
    endotherm_out  = solve_metabolic_rate(organism, environment, T_skin, T_insulation)
    T_skin         = endotherm_out.thermoregulation.T_skin
    T_insulation   = endotherm_out.thermoregulation.T_insulation
    Q_gen          = endotherm_out.energy_fluxes.Q_gen

    T_core_current = HeatExchange.metabolismpars(organism).T_core
    q10mult        = q10 ^ (ustrip(u"K", T_core_current - T_core_ref) / 10)
    Q_basal        = Q_basal_ref * q10mult

    # Fast exit: already at or below basal rate — no torpor entry needed
    Q_gen <= Q_basal && return endotherm_out, organism

    # Torpor descent loop (SOLVENDO.f: do while QGEN > QBASAL .AND. TC > TC_MIN)
    organism_prev      = organism
    endotherm_out_prev = endotherm_out

    while Q_gen > Q_basal && T_core_current > T_core_min
        organism_prev      = organism
        endotherm_out_prev = endotherm_out

        # Decrement T_core (SOLVENDO.f: TC = TC - TC_INC)
        T_core_new  = max(T_core_current - T_core_step, T_core_min)

        # Q10-scale Q_basal at new T_core (SOLVENDO.f: Q10mult = Q10**((TC - TC_REF)/10))
        q10mult_new = q10 ^ (ustrip(u"K", T_core_new - T_core_ref) / 10)
        Q_basal_new = Q_basal_ref * q10mult_new

        # Apply new T_core and Q10-scaled minimum MR to organism
        organism = @set organism.traits.heat_exchange.metabolism_pars.T_core       = T_core_new
        organism = @set organism.traits.heat_exchange.metabolism_pars.Q_metabolism = Q_basal_new

        # Solve heat balance at new T_core
        endotherm_out = solve_metabolic_rate(organism, environment, T_skin, T_insulation)
        T_skin        = endotherm_out.thermoregulation.T_skin
        T_insulation  = endotherm_out.thermoregulation.T_insulation
        Q_gen         = endotherm_out.energy_fluxes.Q_gen

        # Environment too warm: Q_gen strongly negative → revert one step and exit
        # (SOLVENDO.f: IF (RESPFN/RESPGEN) < -TORPTOL → TC = TC_LAST, GO TO 222)
        if Q_basal_new > zero(Q_basal_new) &&
                ustrip(Q_gen / Q_basal_new) < -torptol
            organism      = organism_prev
            endotherm_out = endotherm_out_prev
            break
        end

        T_core_current = T_core_new
        Q_basal        = Q_basal_new
    end

    return endotherm_out, organism
end

"""
    thermoregulate(::Heterotherm, organism, environment, Q_gen, T_skin, T_insulation)

Run heterotherm thermoregulation, dispatching on the organism's control strategy.
"""
function thermoregulate(
    ::Heterotherm,
    organism::Organism,
    environment::NamedTuple,
    Q_gen,
    T_skin,
    T_insulation,
)
    thermoregulate(
        Heterotherm(),
        control_strategy(organism),
        organism,
        environment,
        Q_gen,
        T_skin,
        T_insulation,
    )
end

"""
    thermoregulate(::Heterotherm, ::RuleBasedSequentialControl, organism, environment,
                   Q_gen, T_skin, T_insulation)

Heterotherm thermoregulation using rule-based sequential control.

## Algorithm

1. Run the full endotherm physiological loop (piloerection → uncurl → vasodilate →
   hyperthermia → pant → sweat). If this achieves heat balance at or above Q_minimum,
   the animal is normothermic.
2. Attempt torpor via `torpor()`, which fast-exits if Q_gen is already ≤ Q_basal at
   normothermic T_core.

## Returns

The standard `solve_metabolic_rate` NamedTuple (`thermoregulation`, `morphology`,
`energy_fluxes`, `mass_fluxes`) merged with `torpor_T_core` — the core temperature
achieved. Equals normothermic T_core when no torpor was entered.
"""
function thermoregulate(
    ::Heterotherm,
    ::RuleBasedSequentialControl,
    organism::Organism,
    environment::NamedTuple,
    Q_gen,
    T_skin,
    T_insulation,
)
    limits = thermoregulation(organism)

    # 1. Run the normothermic endotherm physiological loop
    endotherm_out = thermoregulate(
        Endotherm(),
        RuleBasedSequentialControl(),
        organism,
        environment,
        Q_gen,
        T_skin,
        T_insulation,
    )

    # 2. Attempt torpor — torpor() fast-exits if Q_gen ≤ Q_basal at normothermic T_core
    T_skin_out  = endotherm_out.thermoregulation.T_skin
    T_insul_out = endotherm_out.thermoregulation.T_insulation
    torpor_out, organism_torpid = torpor(organism, environment, T_skin_out, T_insul_out, limits)

    torpor_T_core = HeatExchange.metabolismpars(organism_torpid).T_core
    return merge(torpor_out, (; torpor_T_core))
end
