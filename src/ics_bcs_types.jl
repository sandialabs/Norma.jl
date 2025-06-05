# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type BoundaryCondition end
abstract type InitialCondition end
abstract type SolidMechanicsBoundaryCondition <: BoundaryCondition end
abstract type SolidMechanicsRegularBoundaryCondition <: SolidMechanicsBoundaryCondition end
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

mutable struct SolidMechanicsInclinedDirichletBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_fun::Function
    velo_fun::Function
    acce_fun::Function
    reference_funs::Vector{Function}
end

mutable struct SolidMechanicsNeumannBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_fun::Function
end

mutable struct SolidMechanicsRobinBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    #IKT 6/5/2025: new struct for regular Robin BCs.  TODO fill in.
    #Code below is cut/paste from SolidMechanicsNeumannBoundaryCondition
    #name::String
    #offset::Int64
    #side_set_id::Int64
    #num_nodes_per_side::Vector{Int64}
    #side_set_node_indices::Vector{Int64}
    #traction_fun::Function
end

mutable struct SolidMechanicsContactSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_subsim::Simulation
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    rotation_matrix::Matrix{Float64}
    active_contact::Bool
    friction_type::Int64
    variational::Bool
end

mutable struct SolidMechanicsOverlapSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_node_indices::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    coupled_subsim::Simulation
    subsim::Simulation
    variational::Bool
end

#IKT 6/5/2025: Q for Alejandro.  Do we want to rename this SolidMechanicsDNSchwarzBoundaryCondition or something like this,
#now that we are planning to have a Robin-Robin option too?
mutable struct SolidMechanicsNonOverlapSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_subsim::Simulation
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    variational::Bool
end
