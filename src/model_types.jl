# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type Model end
using Exodus
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

struct EvaluationFlags
    is_dynamic::Bool
    is_implicit::Bool
    is_hessian_opt::Bool
    is_matrix_free::Bool
    need_lumped_mass::Bool
    need_stiffness::Bool
    need_mass::Bool
    compute_lumped_mass::Bool
    compute_stiffness::Bool
    compute_mass::Bool
    mesh_smoothing::Bool
end

struct SMThreadLocalArrays{V,M}
    energy::Vector{Float64}
    internal_force::Vector{V}
    lumped_mass::Vector{V}
    stiffness::Vector{M}
    mass::Vector{M}
end

struct SMElementThreadLocalArrays{T,DOFV,IFV,LMV,SM,MM}
    energy::Vector{T}
    dofs::Vector{DOFV}
    internal_force::Vector{IFV}
    lumped_mass::Vector{LMV}
    stiffness::Vector{SM}
    mass::Vector{MM}
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


include("opinf/opinf_model_types.jl")
