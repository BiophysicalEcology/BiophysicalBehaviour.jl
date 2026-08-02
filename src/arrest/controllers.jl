# Controller interface: (own_state, progress, signals, model, arrest_state)
# -> level in [0,1]. `level` is the primary quantity -- read directly as a
# continuous multiplier on a process rate; hosts derive their own booleans
# from it if they need one. SmoothingStrategy picks exact step (HardBound)
# vs smooth, differentiable transition (SmoothBound) -- see threshold.jl.

abstract type AbstractArrestController end

# flattens a Tuple-of-Tuples via recursive tuple destructuring. Used to 
# trigger_conditions/arrest_conditions concretely typed for type-stable 
# root-finding.
_flatten_tuples(::Tuple{}) = ()
_flatten_tuples(t::Tuple) = (first(t)..., _flatten_tuples(Base.tail(t))...)

# comparison/direction as dispatchable types, not Symbols (zero-cost per concrete controller type).
abstract type AbstractComparison end
struct Below <: AbstractComparison end
struct Above <: AbstractComparison end

# positive when condition holds -- signed gap fed to safe_step.
signed_gap(c::AbstractComparison, value, bound) = error("no signed_gap method for $(typeof(c))")
signed_gap(::Below, value, bound) = bound - value
signed_gap(::Above, value, bound) = value - bound

abstract type AbstractDirection end
struct AnyDirection <: AbstractDirection end
struct Rising <: AbstractDirection end
struct Falling <: AbstractDirection end

direction_gate(d::AbstractDirection, rate) = error("no direction_gate method for $(typeof(d))")
direction_gate(::AnyDirection, rate) = true
direction_gate(::Rising, rate) = rate > zero(rate)
direction_gate(::Falling, rate) = rate < zero(rate)

# signals is a plain NamedTuple; bare value or (;value,rate) for directional controllers.
signal_value(x) = x
signal_value(x::NamedTuple{(:value, :rate)}) = x.value

signal_rate(x::NamedTuple{(:value, :rate)}) = x.rate
signal_rate(x) = error("directional controller needs a `(; value, rate)` signal, got $(typeof(x))")

initial_controller_state(::AbstractArrestController) = NamedTuple()
controller_rate(::AbstractArrestController, own_state, progress, signals, model, arrest_state) = NamedTuple()

function level end

# continuous-ODE hosts only; only HardBound controllers need one (SmoothBound has no discontinuity to root-find).
trigger_conditions(::AbstractArrestController) = ()
register_callback(::AbstractArrestController) = false

struct NeverController <: AbstractArrestController end
level(::NeverController, own_state, progress, signals, model, arrest_state) = 0.0
