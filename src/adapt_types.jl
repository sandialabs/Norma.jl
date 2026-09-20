# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Options of the topological phase of the adaptivity loop
# (docs/notes/ems-adaptivity), from the `adaptivity` block of the input.
struct AdaptivityOptions
    desired_density::Float64     # elements above this energy density are candidates
    allowed_density::Float64     # no accepted operation may create an element above this
    minimum_scaled_jacobian::Float64  # geometric floor for the worst element an operation creates
    minimum_decrease::Float64    # relative decrease of the cavity energy an operation must achieve
    adjacency_layers::Int        # rings of adjacent elements added to the candidate set
    maximum_passes::Int          # passes of a topology phase
    outer_iterations::Int        # alternations of smoothing and topology
    swaps::Bool
    collapses::Bool
    splits::Bool
    size_by_length::Bool         # size operations accepted on the edge length alone, not the energy
end

function AdaptivityOptions(
    desired_density::Real,
    allowed_density::Real,
    minimum_scaled_jacobian::Real,
    minimum_decrease::Real,
    adjacency_layers::Integer,
    maximum_passes::Integer,
    outer_iterations::Integer,
    swaps::Bool,
    collapses::Bool,
    splits::Bool,
)
    return AdaptivityOptions(
        desired_density,
        allowed_density,
        minimum_scaled_jacobian,
        minimum_decrease,
        adjacency_layers,
        maximum_passes,
        outer_iterations,
        swaps,
        collapses,
        splits,
        false,
    )
end

# The node an edge split adds: its position after placement and local
# relaxation, the sets it inherits, the edge it splits, and the nodal metric
# data it carries (nothing unless the metric is carried by the nodes).
struct SplitNode
    position::SVector{3,Float64}
    node_sets::Vector{Int}
    side_sets::Vector{Int}
    edge::Tuple{Int,Int}
    metric::Union{Nothing,Vector{Float64}}
end

# Result of one topological operation: the cavity replaced, the energies
# before and after, for a collapse the node removed and the node it was moved
# onto (zero otherwise), and for a split the node added (nothing otherwise),
# which the new connectivity refers to by the index it will receive.
struct CavityProposal
    old_elements::Vector{Int}
    new_connectivity::Matrix{Int}
    block::Int
    energy_before::Float64
    energy_after::Float64
    removed_node::Int
    surviving_node::Int
    split::Union{Nothing,SplitNode}
end
