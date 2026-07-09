# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
# An overlapping coupling, so it belongs under the overlap-coupling abstract
# alongside the FOM overlap BCs (it previously subtyped the grandparent
# SolidMechanicsSchwarzBoundaryCondition directly, so `isa CouplingSchwarz`
# checks silently skipped it).  apply_bc_detail for the neural-network ROM is
# dispatched on the concrete type in ext/NormaPyTorchExt.jl and stays strictly
# more specific than the OpInfModel/CouplingSchwarz fallback, so this reparenting
# does not change method resolution for it.
mutable struct SolidMechanicsOpInfOverlapSchwarzBoundaryCondition <: SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_node_indices::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    fom_bc::SolidMechanicsOverlapSchwarzBoundaryCondition
    nn_model::Any
    basis::Array{Float64}
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
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
