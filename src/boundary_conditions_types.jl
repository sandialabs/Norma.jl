# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type BoundaryCondition end
abstract type InitialCondition end
abstract type SolidMechanicsBoundaryCondition <: BoundaryCondition end
abstract type SolidMechanicsRegularBoundaryCondition <: SolidMechanicsBoundaryCondition end
abstract type SolidMechanicsNeumannRobinBoundaryCondition <: SolidMechanicsRegularBoundaryCondition end
abstract type SolidMechanicsSchwarzBoundaryCondition <: SolidMechanicsBoundaryCondition end
abstract type SolidMechanicsCouplingSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition end

using Symbolics

mutable struct SolidMechanicsDirichletBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_fun::Function
    velo_fun::Function
    acce_fun::Function
end

mutable struct SolidMechanicsNeumannBoundaryCondition <: SolidMechanicsNeumannRobinBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_fun::Function
end

mutable struct SolidMechanicsRobinBoundaryCondition <: SolidMechanicsNeumannRobinBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_fun::Function
    robin_parameter::Float64
end

mutable struct SolidMechanicsNeumannPressureBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    pressure_fun::Function
end

mutable struct SolidMechanicsContactSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    rotation_matrix::Matrix{Float64}
    active_contact::Bool
    friction_type::Int64
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

mutable struct SolidMechanicsOverlapSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    coupled_block_name::String
    search_tolerance::Float64
    dirichlet_projector::Matrix{Float64}
    use_weak::Bool
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# Impedance-matching overlap Schwarz: replaces DBC-DBC with absorbing
# conditions on the overlap boundaries. Same strong interpolation as
# regular overlap, but applies t + Z u̇ = g as a force (not a constraint).
mutable struct SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    square_projector::Matrix{Float64}
    impedance::Float64
    robin_parameter::Float64     # α for displacement penalty (0 = pure impedance)
    impedance_scale::Vector{Float64}  # multiplier on Z per step (default [1.0])
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

mutable struct SolidMechanicsNonOverlapSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

mutable struct SolidMechanicsRobinSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    robin_parameter::Float64
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# Impedance-matching Robin-Robin Schwarz: t + Z u̇ + α W u = g
# Z = ρ c_p (characteristic impedance) absorbs outgoing waves at the interface,
# preventing reflections that cause energy growth with mixed integrators.
# The displacement penalty α W u provides quasi-static stability.
mutable struct SolidMechanicsImpedanceSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    impedance::Float64           # Z = ρ c_p = √(ρ(λ + 2μ))
    robin_parameter::Float64     # α for displacement penalty (0 = pure impedance)
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

include("opinf/opinf_ics_bcs_types.jl")
