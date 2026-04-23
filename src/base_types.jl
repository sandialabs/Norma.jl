# Base abstract types and aliases needed by all other type definition files.

abstract type Simulation end
const Parameters = Dict{String,Any}

# Stable identifier for a sub-simulation slot within a MultiDomainSimulation.
# Used by Schwarz BCs so they don't hold direct object pointers to partner subsims.
struct DomainHandle
    id::Int64
end

# Math errors from constitutive models (e.g. J2 plasticity raising a negative
# number to a fractional power) that should be treated as evaluation failures
# rather than hard crashes.
const _MATH_ERRORS = Union{DomainError, InexactError, OverflowError}
