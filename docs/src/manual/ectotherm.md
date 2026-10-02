# Ectotherm thermoregulation

The body temperature of an ectotherm is whatever balances its heat budget where it is. It regulates by changing
where it is, how it is oriented and what colour it is. This page describes the sequence of decisions, that of
the NicheMapR ectotherm model (Kearney and Porter 2020), going back to Porter et al. (1973).

```@setup ectotherm
using Main.FigureHelpers
using Main.ExampleEnvironments
using CairoMakie
using BiophysicalBehaviour, HeatExchange, BiophysicalGeometry, Unitful
```

## The aim

For each hour, find a position and posture at which the steady-state body temperature ``T_b`` lies between the
minimum temperature for activity and the target temperature,

```math
T_{F,\text{min}} \le T_b \le T_{\text{pref}}
```

or, failing that, the best that can be had. The state is ``T_b``, from [`solve_body_temperature`](@ref). The
reference is the pair of threshold traits. The error is which side of the range ``T_b`` falls on, and the
controller moves one actuator one step, in a fixed order, and solves again, see
[Behaviour as control](control.md#Kinds-of-controller).

One call does this for one hour:

```julia
thermoregulate(organism, environments, limits, site, step, previous_depth; activity_commenced)
```

| Argument | Meaning |
|:--|:--|
| `organism` | an `Organism` whose traits have the [`Ectotherm`](@ref) strategy |
| `environments` | the [`AvailableEnvironments`](@ref), see [Activity and available environments](environments.md) |
| `limits` | the [`EctothermBehavioralLimits`](@ref) of the organism, from `thermoregulation(organism)` |
| `site` | the `EnvironmentalPars` of HeatExchange.jl |
| `step` | the hour, as an index into the microclimate |
| `previous_depth` | the soil node occupied in the previous hour. 1 is the surface |
| `activity_commenced` | whether the animal has already been out today |

## The sequence

**1. Start from the reference position.** Each hour begins on the surface in the least shade, at the lowest
height, darkest colour and a neutral posture, with the target temperature at its starting value
([`reset_position`](@ref)). The animal carries no memory of the previous hour except its depth.

**2. Is this the time of day for activity?** The [`ActivityPeriod`](@ref) is tested against the zenith angle
and sunlight ([`is_active`](@ref)). If not, and the animal can retreat underground, [`select_depth`](@ref)
finds it a depth, and the hour is done.

**3. Can it emerge?** After an hour underground the animal comes out only if the soil where it is has reached
`emerge_temperature_min`. With a non-zero `emerge_signal` it also waits for the soil to be warming, or
cooling, at that rate, as an animal deep in a burrow has no other cue to the time of day. Otherwise it stays
down, at a depth chosen again for this hour.

**4. Regulate.** On the surface, body temperature is solved and compared with the thresholds. One response is
made, body temperature is solved again, and the loop repeats until it is in range or nothing is left to try.

```@example ectotherm
ladder_diagram(["lighten", "seek shade", "raise target", "climb", "pant", "retreat\nunderground"]; # hide
    title = "Too hot: body temperature above the target", colour = RGBf(0.95, 0.80, 0.78)) # hide
```

```@example ectotherm
ladder_diagram(["darken", "face the sun", "press to\nground", "leave shade", "retreat\nunderground"]; # hide
    title = "Too cold: body temperature below the basking minimum", colour = RGBf(0.78, 0.87, 0.95)) # hide
```

**Too hot**, ``T_b > T_{\text{pref}}``. The first of these that is allowed and not exhausted:

| | Response | Function | Until |
|:--|:--|:--|:--|
| 1 | become paler | [`lighten`](@ref) | absorptivity reaches its minimum |
| 2 | move into more shade | [`seek_shade`](@ref) | shade reaches the maximum available |
| 3 | tolerate a warmer body | [`increment_target_temperature`](@ref) | the target reaches `active_temperature_max` |
| 4 | climb | [`climb`](@ref) | the highest node |
| 5 | pant | [`pant`](@ref) | the maximum panting rate |
| 6 | go underground | [`select_depth`](@ref) | the loop ends |

**Too cold**, ``T_b < T_{B,\text{min}}``:

| | Response | Function | Until |
|:--|:--|:--|:--|
| 1 | become darker | [`darken`](@ref) | absorptivity reaches its maximum |
| 2 | turn broadside to the sun | [`orient_perpendicular`](@ref) | done in one step |
| 3 | press against the ground | [`press_to_ground`](@ref) | done in one step |
| 4 | by day, leave the shade | [`avoid_shade`](@ref) | shade reaches the minimum |
| 5 | by night, move under cover, away from the cold sky | [`seek_shade`](@ref) | shade reaches the maximum |
| 6 | if below the critical minimum, climb | [`climb`](@ref) | the highest node |
| 7 | go underground | [`select_depth`](@ref) | the loop ends |

**Basking**, ``T_{B,\text{min}} \le T_b < T_{F,\text{min}}``: turn broadside to the sun if not already, and
otherwise accept the state. A lizard that has basked broadside into the activity range returns to a neutral
posture, and its body temperature is solved once more.

Each response is switched on or off by a capability flag: `can_change_absorptivity`, `can_seek_shade`,
`can_climb`, `can_pant`, `can_solar_orient`, `can_press_to_ground`, `can_retreat_underground`. These are model
traits, see [States, thresholds and traits](states_traits.md#The-four-classes-of-functional-trait).

**5. Classify.** The final body temperature decides the [`OrganismState`](@ref): [`Active`](@ref),
[`Basking`](@ref) or [`Resting`](@ref). Underground, by default, body temperature is that of the soil, as in
NicheMapR. Set `solve_underground = true` to solve the heat budget there too.

## What the order means

The order is that of cost. Colour and posture are free and immediate. Shade costs little. Raising the target
accepts a body temperature above the preferred one before giving up the surface. Climbing and panting come
late, and the burrow last, because an animal underground is not feeding.

The third step deserves a note. The target temperature starts at the preferred temperature and rises in steps
to the maximum for activity. It is the reference of the controller that moves, not an actuator. The animal uses
all the shade it has to stay at its preferred temperature, and tolerates more only when the shade runs out.

## Choosing a depth

[`select_depth`](@ref) returns the shallowest node, from `depth_min_underground` down, at which the soil is
warmer than the critical minimum and cooler than a point half way between `active_temperature_max` and the
critical maximum. If none qualifies it returns the deepest. The animal is as near the surface as is safe, where
it will be first to detect that conditions above have changed.

## An hour at a time

The lizard of [Get started](../get_started.md), on the middle day of July at Madison:

```@example ectotherm
environments = madison_environments()
site = example_environment_pars(; elevation = 226.0u"m")
lizard(; kw...) = Organism(Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()),
    example_ectotherm_organism_traits(; kw...))

july = (7 - 1) * 24 .+ (1:24)
day = simulate(lizard(), environments, site, july)
state_name(s) = lowercase(string(nameof(typeof(s)))) # hide
markdown_table(["Hour", "Body temperature", "State", "Shade", "Depth", "Posture"], [ # hide
    (h - 1, celsius(s.core_temperature), state_name(s.state), s.depth_node == 1 ? s.shade : "–", # hide
     environments.depths[s.depth_node], s.depth_node > 1 ? "–" : s.sun_orientation == 90.0 ? "facing sun" : "neutral") # hide
    for (h, s) in enumerate(day) if isodd(h)]) # hide
```

The output of each hour is a NamedTuple:

| Field | Content |
|:--|:--|
| `core_temperature` | body temperature |
| `state` | [`Resting`](@ref), [`Basking`](@ref) or [`Active`](@ref) |
| `shade` | the shade fraction chosen |
| `depth_node` | the soil node. 1 is the surface |
| `height` | the height above the ground, or the depth as a negative height |
| `absorptivity` | the solar absorptivity chosen |
| `sun_orientation` | 90 facing the sun, 45 neutral, 0 parallel to it |
| `pressed_to_ground`, `pant_rate` | the remaining responses |
| `ectotherm_out` | the full heat budget of HeatExchange.jl in the final state |

## What each behaviour is worth

Because each response is a flag, its value can be measured by taking it away. Hours of activity on the middle
day of each month, with all behaviours, without shade, without a burrow, and without posture:

```@example ectotherm
year = 1:288
active_hours(results) = [count(s -> s.state isa Active, results[(m - 1) * 24 .+ (1:24)]) for m in 1:12]
variants = (
    "all behaviours" => lizard(),
    "no shade" => lizard(; can_seek_shade = false),
    "no burrow" => lizard(; can_retreat_underground = false),
    "no posture" => lizard(; can_solar_orient = false),
)
hours = [label => active_hours(simulate(animal, environments, site, year)) for (label, animal) in variants]
fig, ax = figure_axis("Month", "Hours active on the middle day"; size = (720, 380), xticks = 1:12) # hide
for (label, h) in hours # hide
    scatterlines!(ax, 1:12, h; label, linewidth = 2) # hide
end # hide
axislegend(ax; position = :lt, labelsize = 11) # hide
fig # hide
```

The same loop gives the highest body temperature reached in each case:

```@example ectotherm
peak(animal) = maximum(s.core_temperature for s in simulate(animal, environments, site, july))
markdown_table(["Behaviours", "Highest body temperature in July"], # hide
    [(label, celsius(peak(animal))) for (label, animal) in variants]) # hide
```

[A lizard's day](../tutorials/lizard.md) does the same for a desert iguana through a year.

## An endotherm that chooses where to be

An endotherm has the same choices of position, for a different reason: every degree it avoids by moving is
water it need not evaporate. Given available environments and an `EctothermBehavioralLimits`,
`thermoregulate` for an [`Endotherm`](@ref) first runs a loop like the one above on its *operative
temperature*, the body temperature it would have as a passive object, then thermoregulates physiologically at
the position chosen, as in [Endotherm thermoregulation by rules](endotherm_rules.md). In the heat it tries a
paler colour, a posture parallel to the sun ([`orient_parallel`](@ref)), shade, height and the burrow in turn,
before any panting or sweating. See [A desert mammal through the year](../tutorials/endotherm_year.md).
