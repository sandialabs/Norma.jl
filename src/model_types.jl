# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type Model end
abstract type RomModel <: Model end
abstract type OpInfModel <: RomModel end
using SparseArrays

@enum Kinematics begin
    Undefined
    Infinitesimal
    Finite
end

mutable struct COOVector
    index::Vector{Int64}
    vals::Vector{Float64}
    len::Int64  # logical length
end

mutable struct COOMatrix
    rows::Vector{Int64}
    cols::Vector{Int64}
    vals::Vector{Float64}
    len::Int64 # logical length
end

mutable struct SolidMechanics <: Model
    mesh::ExodusDatabase
    materials::Vector{Solid}
    reference::Matrix{Float64}
    current::Matrix{Float64}
    velocity::Matrix{Float64}
    acceleration::Matrix{Float64}
    internal_force::Vector{Float64}
    boundary_force::Vector{Float64}
    boundary_conditions::Vector{BoundaryCondition}
    stress::Vector{Vector{Vector{Vector{Float64}}}}
    stored_energy::Vector{Vector{Float64}}
    strain_energy::Float64
    stiffness::SparseMatrixCSC{Float64,Int64}
    mass::SparseMatrixCSC{Float64,Int64}
    lumped_mass::Vector{Float64}
    body_force::Vector{Float64}
    free_dofs::BitVector
    time::Float64
    compute_stiffness::Bool
    compute_mass::Bool
    compute_lumped_mass::Bool
    failed::Bool
    mesh_smoothing::Bool
    smooth_reference::String
    inclined_support::Bool
    global_transform::SparseMatrixCSC{Float64,Int64}
    kinematics::Kinematics
end

# TODO: Add potential energy as in the above
mutable struct HeatConduction <: Model
    mesh::ExodusDatabase
    materials::Vector{Thermal}
    reference::Matrix{Float64}
    temperature::Vector{Float64}
    rate::Vector{Float64}
    internal_heat_flux::Vector{Float64}
    boundary_heat_flux::Vector{Float64}
    boundary_conditions::Vector{BoundaryCondition}
    flux::Vector{Vector{Vector{Vector{Float64}}}}
    stored_energy::Vector{Vector{Float64}}
    free_dofs::BitVector
    time::Float64
    failed::Bool
end

mutable struct QuadraticOpInfRom <: OpInfModel
    opinf_rom::Dict{Any,Any}
    basis::Array{Float64}
    reduced_state::Vector{Float64}
    reduced_boundary_forcing::Vector{Float64}
    #internal_force not used, but include to ease interfacing in Schwarz
    internal_force::Vector{Float64}
    free_dofs::BitVector
    boundary_conditions::Vector{BoundaryCondition}
    time::Float64
    failed::Bool
    fom_model::SolidMechanics
    reference::Matrix{Float64}
    inclined_support::Bool
end

mutable struct LinearOpInfRom <: OpInfModel
    opinf_rom::Dict{Any,Any}
    basis::Array{Float64}
    reduced_state::Vector{Float64}
    reduced_boundary_forcing::Vector{Float64}
    #internal_force not used, but include to ease interfacing in Schwarz
    internal_force::Vector{Float64}
    free_dofs::BitVector
    boundary_conditions::Vector{BoundaryCondition}
    time::Float64
    failed::Bool
    fom_model::SolidMechanics
    reference::Matrix{Float64}
    inclined_support::Bool
end
