# For NicheMapR users

In [NicheMapR](https://github.com/mrke/NicheMapR) behaviour is inside the models it belongs to: the ectotherm
model moves its animal between shade, height and depth within the Fortran program `ECTOTHERM`, and `endoR` has
a thermoregulatory loop around its heat budget. Here the heat budget is in HeatExchange.jl, see its
[page for NicheMapR users](https://biophysicalecology.github.io/HeatExchange.jl/dev/manual/nichemapr), and
everything the animal does about it is in this package.

## What moved where

| NicheMapR | BiophysicalBehaviour.jl |
|:--|:--|
| `ectotherm`: the hourly loop of `ECTOTHERM.f` with `THERMO.f`, `SHADEADJUST.f`, `SELDEP.f`, `ABOVEGROUND.f`, `BELOWGROUND.f` | [`thermoregulate`](@ref) for an [`Ectotherm`](@ref), called once for each hour, see [Ectotherm thermoregulation](ectotherm.md) |
| `metout`, `shadmet`, `soil`, `shadsoil` from `micro_*` | [`AvailableEnvironments`](@ref) over two results of Microclimate.jl, see [Activity and available environments](environments.md) |
| `endoR` and `endoR_devel` with `THERMOREG = 1` | [`thermoregulate`](@ref) for an [`Endotherm`](@ref) with [`RuleBasedSequentialControl`](@ref), see [Endotherm thermoregulation by rules](endotherm_rules.md) |
| `TREGMODE` 1, 2, 3 | [`CoreFirst`](@ref), [`CoreAndPantingFirst`](@ref), [`CorePantingSweatingFirst`](@ref) |
| `HomoTherm`: `endoR` for each part inside a thermoregulatory loop | an organism on a `CompositeBody`, see [Bodies of many parts](multipart.md) and [A human that thermoregulates](../tutorials/human.md) |
| the `ACT` column of `environ`, 0, 1, 2 | [`Resting`](@ref), [`Basking`](@ref), [`Active`](@ref) |
| nothing | [`IPOPTControl`](@ref), see [Thermoregulation by optimisation](optimisation.md) |
| `onelump`, `onelump_var`, `twolump` | transient heat budgets, to be added to HeatExchange.jl: the residual of its heat balance tracked through time and turned into body temperature by the heat capacity |
| `trans_behav`, and the behaviour of the transient option of `ectotherm` | behaviour during a transient: in development here, not in this documentation |
| the Dynamic Energy Budget model within `ectotherm` | not here: to come through [AnimalMapper.jl](https://github.com/BiophysicalEcology/AnimalMapper.jl), by way of [DEBtool_J.jl](https://github.com/add-my-pet/DEBtool_J.jl) |

## What is different

**The loop over hours is the user's.** `ectotherm` takes a microclimate and returns tables for the whole
period. Here `thermoregulate` does one hour, and the loop around it is a few lines, see
[Get started](../get_started.md). Only the depth of the previous hour, and whether the animal has yet been
active today, carry from hour to hour.

**Capabilities are flags.** Whether an animal seeks shade, climbs, burrows, changes colour, pants or orients to
the sun is each a field of [`EctothermBehavioralLimits`](@ref).

**Every graded response has the same form.** In NicheMapR each has its own trio of names, such as `AK1`,
`AK1_MAX` and `AK1_INC`. Here each is a [`SteppedParameter`](@ref) with `current`, `reference`, `max` and
`step`.

**The controller is a value.** The order of responses is one choice of
[`AbstractControlStrategy`](@ref) among others.

**An endotherm can choose where to be.** `endoR` is given one environment. Here an endotherm can be given the
[`AvailableEnvironments`](@ref) of an ectotherm, and selects among them before it thermoregulates
physiologically.

## Parameter names: `ectotherm`

| NicheMapR | Here | In |
|:--|:--|:--|
| `T_F_min`, `T_F_max` | `active_temperature_min`, `active_temperature_max` | [`EctothermBehavioralLimits`](@ref) |
| `T_B_min` | `basking_temperature_min` | |
| `T_RB_min` | `emerge_temperature_min` | |
| `T_pref` | `target_temperature` | |
| `CT_min`, `CT_max` | `critical_temperature_min`, `critical_temperature_max` | |
| `minshade`, `maxshade` (%) | `shade.reference`, `shade.max` (fractions) | |
| `delta_shade` (%) | `shade.step` | |
| `shade_seek` | `can_seek_shade` | |
| `burrow` | `can_retreat_underground` | |
| `climb` | `can_climb` | |
| `mindepth`, `maxdepth` (soil nodes) | `depth_min_underground`, `depth.max` | |
| `shdburrow` 0, 1, 2 | `burrow_shade_mode`: [`MinShadeOnly`](@ref), [`AdaptiveBurrowShade`](@ref), [`MaxShadeOnly`](@ref) | |
| `warmsig` | `emerge_signal` | |
| `alpha_min`, `alpha_max` | `absorptivity.reference`, `absorptivity.max`, with `can_change_absorptivity` | |
| `postur` | `can_solar_orient`, and `solar_orientation` of the radiation parameters | |
| `pct_cond` | `can_press_to_ground`, and `conduction_fraction` of the conduction parameters | |
| `pantmax` | `pant_rate.max`, with `can_pant` | |
| `diurn`, `nocturn`, `crepus` | `activity_period`: [`Diurnal`](@ref), [`Nocturnal`](@ref), [`Crepuscular`](@ref), [`CombinedActivity`](@ref) | [`BehavioralTraits`](@ref) |

`example_ectotherm_behavioral_limits` has the defaults of `ectotherm`.

## Parameter names: `endoR`

| NicheMapR | Here | In |
|:--|:--|:--|
| `THERMOREG` | the controller acts whenever a response has a range | [`ThermoregulationLimits`](@ref) |
| `TREGMODE` | `control.mode` | [`RuleBasedSequentialControl`](@ref) |
| `QBASAL` | `minimum_heat_flow` | [`ThermoregulationLimits`](@ref) |
| `ZFURD`, `ZFURD_MAX`, `ZFURV`, `ZFURV_MAX` | `insulation.dorsal`, `insulation.ventral` | [`InsulationLimits`](@ref) |
| `SHAPE_B`, `SHAPE_B_MAX`, `UNCURL` | `axis_ratio_factor`: `current`, `max`, `step` | |
| `AK1`, `AK1_MAX`, `AK1_INC` | `flesh_conductivity` | |
| `TC`, `TC_MAX`, `TC_INC` | `core_temperature` | |
| `PANT`, `PANT_MAX`, `PANT_INC` | `panting.pant` | [`PantingLimits`](@ref) |
| `PANT_MULT` | `panting.multiplier` | |
| `PCTWET`, `PCTWET_MAX`, `PCTWET_INC` (%) | `skin_wetness` (fractions) | |
| `Q10` | `q10` of the metabolism parameters of HeatExchange.jl | |
| `DIFTOL`, `BRENTOL` | the options of HeatExchange.jl | |

`example_thermoregulation_limits` has the defaults of `endoR`, except that panting can reach 10 times the
resting rate and all of the skin can be wet.

## How closely they agree

The tests compare both models with NicheMapR.

- **Ectotherm.** The microclimate tables of NicheMapR are read in as the available environments, and hourly body
  temperature, shade, depth and activity are compared for four combinations of behaviours.
- **Endotherm.** [A mammal across air temperatures](../tutorials/mammal.md#Comparison-with-NicheMapR) computes the
  comparison on the page: each response begins at the same air temperature and reaches the same value as in
  `endoR_devel` from 0 to 50 °C, and metabolic rate agrees to 1 % on average and 5 % at worst. The differences
  are those of the heat budget, see the documentation of HeatExchange.jl, and of where each loop stops within
  its tolerance.
- **Human.** [A human that thermoregulates](../tutorials/human.md) compares a six-part person with `HomoTherm`,
  closely to 34 °C, and says where and why the two loops part above it.

`example_respiration_pars` exhales air at lung temperature. NicheMapR's default, `DELTAR = 0`, exhales at air
temperature, and the comparison uses `DELTAR = 100` to match.
