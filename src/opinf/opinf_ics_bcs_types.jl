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

using Exodus
using Symbolics

mutable struct SMOpInfDirichletBC <: RegularBoundaryCondition 
    node_set_name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_num::Num
    velo_num::Num
    acce_num::Num
    fom_bc::SMDirichletBC
    nn_model::Any
    basis::Array{Float64}
end


mutable struct SMOpInfOverlapSchwarzBC <: OverlapSchwarzBoundaryCondition
    side_set_name::String
    side_set_node_indices::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    coupled_subsim::Simulation
    subsim::Simulation
    is_dirichlet::Bool
    swap_bcs::Bool
    fom_bc::SMOverlapSchwarzBC
    nn_model::Any
    basis::Array{Float64}
end


