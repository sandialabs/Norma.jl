# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

mutable struct SolidMechanicsOpInfOverlapSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition
    name::String
    side_set_node_indices::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    coupled_subsim::Simulation
    subsim::Simulation
    variational::Bool
    fom_bc::SolidMechanicsOverlapSchwarzBoundaryCondition
    nn_model::Any
    basis::Array{Float64}
end

mutable struct SolidMechanicsOpInfDirichletBC <: SolidMechanicsRegularBoundaryCondition
    name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_num::Num
    velo_num::Num
    acce_num::Num
    fom_bc::SolidMechanicsDirichletBoundaryCondition
    nn_model::Any
    basis::Array{Float64}
end
