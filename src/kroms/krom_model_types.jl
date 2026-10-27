# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

include("rbf.jl")

mutable struct RBFKernelROM <: OpInfModel
    interpolant::RBFKernelInterpolant
    opinf_rom::Dict{Any,Any}
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
