include("rbf.jl")

mutable struct RBFKernelROM <: OpInfModel
    interpolant::RBFKernelInterpolant
    forcing::Vector{Float64}
    basis::Array{Float64}
    reduced_state::Vector{Float64}
    reduced_velocity::Vector{Float64}
    reduced_boundary_forcing::Vector{Float64}
    internal_force::Vector{Float64}
    free_dofs::BitVector
    boundary_conditions::Vector{BoundaryCondition}
    time::Float64
    failed::Bool
    fom_model::SolidMechanics
    reference::Matrix{Float64}
    inclined_support::Bool
end
