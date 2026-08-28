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
    prev_state_old::Vector{Vector{Vector{Vector{Float64}}}}  # saved copy for step-failure rollback
    stop_state_old::Vector{Vector{Vector{Vector{Float64}}}}  # saved copy for Schwarz-stop rollback
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
    size_field::Union{Function,Nothing}  # compiled s(t,x,y,z) target edge length, or nothing
    kinematics::Kinematics
    recovery_data::AbstractRecoveryData
    # Single-mode recovery (recovery_data isa LumpedRecovery or ConsistentRecovery).
    # Each is zeros(0, 0) when the quantity is not enabled OR when in BothRecovery mode.
    recovered_stress::Matrix{Float64}                   # 6 × n_nodes
    recovered_von_mises::Matrix{Float64}                # 1 × n_nodes
    recovered_F::Matrix{Float64}                        # 9 × n_nodes
    recovered_internal_variables::Matrix{Float64}       # n_iv × n_nodes
    # BothRecovery variants.  Each is zeros(0, 0) outside both-mode or when
    # the quantity is not enabled.
    lumped_recovered_stress::Matrix{Float64}
    consistent_recovered_stress::Matrix{Float64}
    lumped_recovered_von_mises::Matrix{Float64}
    consistent_recovered_von_mises::Matrix{Float64}
    lumped_recovered_F::Matrix{Float64}
    consistent_recovered_F::Matrix{Float64}
    lumped_recovered_internal_variables::Matrix{Float64}
    consistent_recovered_internal_variables::Matrix{Float64}
    num_int_pts::Vector{Int}
    # Tracks whether this model was constructed from restart
    # snapshot data 
    restarted::Bool
end

# Whether a Model subtype's restart snapshot (nodal displacement/velocity)
# fully captures the state needed to resume a run. Defaults to false; a
# model type opts in by overriding this trait next to its own struct
# definition, rather than being enumerated in a separately maintained list.
# process_restart!() (simulation.jl) calls
# supports_restart(model_type_for(model_type)) (model.jl) instead of
# checking membership in a hand-kept table, so a new Model subtype can't
# silently fall out of sync with create_model() the way the old
# RESTART_SUPPORTED_MODEL_TYPES string table could.
supports_restart(::Type{<:Model}) = false
supports_restart(::Type{SolidMechanics}) = true

include("opinf/opinf_model_types.jl")
include("kroms/krom_model_types.jl")

# Every ROM type is layered on the FOM restart machinery via its internal
# fom_model::SolidMechanics (see process_restart!() in simulation.jl and
# apply_ics(::Parameters, ::RomModel, ...) in opinf_ics_bcs.jl), so restart
# support is a property of RomModel itself, not of each individual ROM
# subtype. This one method covers LinearOpInfRom, QuadraticOpInfRom,
# CubicOpInfRom, NeuralNetworkOpInfRom, RBFKernelROM, and any future
# RomModel subtype automatically -- no separate entry needed per type.
supports_restart(::Type{<:RomModel}) = true
