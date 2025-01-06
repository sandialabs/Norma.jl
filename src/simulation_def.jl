abstract type Simulation end

include("constitutive_def.jl")
include("ics_bcs_def.jl")
include("model_def.jl")
include("time_integrator_def.jl")
include("solver_def.jl")
include("schwarz_def.jl")

mutable struct SingleDomainSimulation <: Simulation
    name::String
    params::Dict{String,Any}
    integrator::TimeIntegrator
    solver::Solver
    model::Model
    failed::Bool
end

mutable struct MultiDomainSimulation <: Simulation
    name::String
    params::Dict{String,Any}
    schwarz_controller::SchwarzController
    subsims::Vector{SingleDomainSimulation}
    subsim_name_index_map::Dict{String,Int64}
    failed::Bool
end
