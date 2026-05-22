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
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    acceleration::Matrix{Float64}
    internal_force::Vector{Float64}
    boundary_force::Vector{Float64}
    boundary_conditions::Vector{BoundaryCondition}
    state_old::Vector{Vector{Vector{Vector{Float64}}}}
    state::Vector{Vector{Vector{Vector{Float64}}}}
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
    kinematics::Kinematics
    recovery_data::AbstractRecoveryData
    recovered_stress::Matrix{Float64}
    recovered_internal_variables::Matrix{Float64}
    # Populated only when recovery_data isa BothRecovery; otherwise zeros(0,0).
    lumped_recovered_stress::Matrix{Float64}
    consistent_recovered_stress::Matrix{Float64}
    lumped_recovered_internal_variables::Matrix{Float64}
    consistent_recovered_internal_variables::Matrix{Float64}
    # Nodal von Mises stress derived from the recovered nodal stress tensor.
    # Populated after compute_nodal_von_mises! is called (which must follow recover_stress!).
    # For single recovery mode (lumped or consistent): recovered_von_mises is used.
    # For BothRecovery: lumped_recovered_von_mises and consistent_recovered_von_mises are used.
    # Empty (Float64[]) when the corresponding recovery mode is not active.
    recovered_von_mises::Vector{Float64}
    lumped_recovered_von_mises::Vector{Float64}
    consistent_recovered_von_mises::Vector{Float64}
    num_int_pts::Vector{Int}
end

include("opinf/opinf_model_types.jl")
include("kroms/krom_model_types.jl")