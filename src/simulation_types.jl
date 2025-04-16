# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type Simulation end
const Parameters = Dict{String,Any}

include("constitutive_types.jl")
include("ics_bcs_types.jl")
include("model_types.jl")
include("time_integrator_types.jl")
include("solver_types.jl")
include("schwarz_types.jl")

mutable struct SingleDomainSimulation <: Simulation
    name::String
    params::Parameters
    controller::SolidSingleController
    integrator::TimeIntegrator
    solver::Solver
    model::Model
    failed::Bool
end

mutable struct MultiDomainSimulation <: Simulation
    name::String
    params::Parameters
    controller::MultiDomainController
    subsims::Vector{SingleDomainSimulation}
    subsim_name_index_map::Dict{String,Int64}
    failed::Bool
end
