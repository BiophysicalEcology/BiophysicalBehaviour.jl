```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "BiophysicalBehaviour.jl"
  text: "Behaviour as the control of heat and water budgets"
  tagline: "what an organism does about its body temperature, metabolic rate and water loss: where it goes, how it holds itself and what it does with its fur, blood and breath, by rules or by optimisation."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🎛️ Behaviour as control
    details: The heat budget is the plant, body temperature the state, a threshold trait the setpoint, and <a class="highlight-link">thermoregulation the controller</a>.
    link: /manual/control
  - title: 🦎 Ectotherms
    details: Emerge, bask, change colour and posture, <a class="highlight-link">seek shade, climb and retreat underground</a>, hour by hour through a microclimate.
    link: /manual/ectotherm
  - title: 🐇 Endotherms
    details: Raise and flatten fur, uncurl, vasodilate, <a class="highlight-link">let core temperature rise, pant and sweat</a>, in a fixed order until the heat budget balances.
    link: /manual/endotherm_rules
  - title: 📉 Optimisation
    details: Or state the costs and let <a class="highlight-link">IPOPT, with derivatives from Enzyme</a>, find the combination of responses.
    link: /manual/optimisation
  - title: 🐕 Bodies of many parts
    details: A trunk, head and limbs with <a class="highlight-link">their own physiology</a>, one core and one pair of lungs.
    link: /manual/multipart
  - title: 🌡️ Thresholds are the traits
    details: Body temperature is a state. <a class="highlight-link">The temperatures at which an animal acts</a> are its traits.
    link: /manual/states_traits
  - title: ↔️ Gradients and control
    details: Heat flows down gradients of potential. Behaviour follows <a class="highlight-link">gradients of information</a>, between a state and its target.
    link: /manual/gradients
  - title: 📦 For NicheMapR users
    details: How the behaviour of <a class="highlight-link">ectotherm, endoR and HomoTherm</a> maps onto this package.
    link: /manual/nichemapr
---
```

## How to install BiophysicalBehaviour.jl?

BiophysicalBehaviour.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl")
```

## Manual

BiophysicalBehaviour.jl computes what an organism does to regulate its body temperature. It takes the heat
budget of [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) and wraps a controller around
it: the organism changes its position, its posture or its physiology, the heat budget is solved again, and the
result is compared with a target. The package could as well have been called BiophysicalControl.jl, see
[Behaviour as control](manual/control.md).

It is part of the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche
modelling. The bodies come from
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl), the heat budgets from
HeatExchange.jl, and the environments to choose between from
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl). See the
[Introduction](manual/introduction.md) for the design of the package and
[For NicheMapR users](manual/nichemapr.md) for its origin in [NicheMapR](https://github.com/mrke/NicheMapR).
