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

# One value segment per element for the assembled vectors, per block. The
# threaded element loop writes disjoint segments, so no thread-local copies or
# merges are needed; the reduction into the global vector is one pass in
# element order, which also makes the result independent of the thread count.
struct ElementValueBuffers
    force::Vector{Vector{Float64}}        # per block, num_elements * dofs per element
    lumped_mass::Vector{Vector{Float64}}
end

# Fixed sparsity pattern of the global matrices with, per block and element,
# the position in the value array of every entry of the element matrix in the
# column-major order of the element matrix, and the per-element value segments
# of the stiffness and the mass. Built once on first use; the matrices are
# rebuilt each evaluation from the pattern without sorting.
struct MatrixPattern
    colptr::Vector{Int64}
    rowval::Vector{Int64}
    slots::Vector{Vector{Int64}}
    stiffness::Vector{Vector{Float64}}
    mass::Vector{Vector{Float64}}
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

struct SMElementThreadLocalArrays{T,DOFV,IFV,LMV,SM,MM}
    energy::Vector{T}
    dofs::Vector{DOFV}
    internal_force::Vector{IFV}
    lumped_mass::Vector{LMV}
    stiffness::Vector{SM}
    mass::Vector{MM}
end

# Per-block mesh data read once at construction: the connectivity and the
# shape function tables were read from the Exodus file and rebuilt on every
# evaluation before. The tables are static arrays whose type depends on the
# element type, so they are stored untyped and passed through a function
# barrier once per block.
struct ElementBlockData
    id::Int64
    element_type::ElementType
    num_points::Int64
    num_elements::Int64
    num_nodes_per_element::Int64
    connectivity::Matrix{Int64}
    N::Any
    dN::Any
    weights::Any
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
    blocks::Vector{ElementBlockData}
    value_buffers::ElementValueBuffers
    matrix_pattern::Union{Nothing,MatrixPattern}
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
