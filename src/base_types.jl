# Base abstract types and aliases needed by all other type definition files.

abstract type Simulation end
const Parameters = Dict{String,Any}

# Math errors from constitutive models (e.g. J2 plasticity raising a negative
# number to a fractional power) that should be treated as evaluation failures
# rather than hard crashes.
const _MATH_ERRORS = Union{DomainError, InexactError, OverflowError}
