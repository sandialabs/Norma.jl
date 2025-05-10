# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type BoundaryCondition end
abstract type SchwarzBoundaryCondition <: BoundaryCondition end
abstract type RegularBoundaryCondition <: BoundaryCondition end
abstract type ContactSchwarzBoundaryCondition <: SchwarzBoundaryCondition end
abstract type CouplingSchwarzBoundaryCondition <: SchwarzBoundaryCondition end
abstract type OverlapSchwarzBoundaryCondition <: CouplingSchwarzBoundaryCondition end
abstract type NonOverlapSchwarzBoundaryCondition <: CouplingSchwarzBoundaryCondition end
abstract type InitialCondition end

using Symbolics

mutable struct SMDirichletBC <: RegularBoundaryCondition
    name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_num::Num
    velo_num::Num
    acce_num::Num
end

mutable struct SMDirichletInclined <: RegularBoundaryCondition
    name::String
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_expression::Vector{Num}
    velo_expression::Vector{Num}
    acce_expression::Vector{Num}
    reference_normal::Vector{Num}
    off_axis_free::Bool
end

mutable struct SMNeumannBC <: RegularBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_num::Num
end

mutable struct SMContactSchwarzBC <: ContactSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64, Int64}
    global_from_local_map::Vector{Int64}
    coupled_subsim::Simulation
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichelt_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    rotation_matrix::Matrix{Float64}
    active_contact::Bool
    friction_type::Int64
end

mutable struct SMOverlapSchwarzBC <: OverlapSchwarzBoundaryCondition
    name::String
    side_set_node_indices::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    coupled_subsim::Simulation
    subsim::Simulation
end

mutable struct SMNonOverlapSchwarzBC <: NonOverlapSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64, Int64}
    global_from_local_map::Vector{Int64}
    coupled_subsim::Simulation
    coupled_bc_name::String
    coupled_bc_index::Simulation
    dirichelt_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
end
