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
    shape_by_quality::Bool       # shape operations accepted on the cavity minimum of the scaled Jacobian
    face_swaps::Bool             # try the face swap (two elements into three) besides the edge swap
    boundary_swaps::Bool         # try the swap of boundary edges between faces of one flat patch
    boundary_swap_angle::Float64 # largest angle in radians between the two boundary faces of such an edge
    desired_quality::Float64     # under the scaled Jacobian criterion, elements below this are candidates
    size_first::Bool             # collapse and split before swapping in every pass, rather than after
    shape_every_pass::Bool       # try the shape-driven collapses and splits in every pass
end

# The options before the criteria and the extra operators, with those as
# keywords, for the tests.
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
    splits::Bool;
    size_by_length::Bool=false,
    shape_by_quality::Bool=false,
    face_swaps::Bool=false,
    boundary_swaps::Bool=false,
    boundary_swap_angle::Real=0.0,
    desired_quality::Real=0.9,
    size_first::Bool=false,
    shape_every_pass::Bool=false,
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
        size_by_length,
        shape_by_quality,
        face_swaps,
        boundary_swaps,
        boundary_swap_angle,
        desired_quality,
        size_first,
        shape_every_pass,
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

# The change of a side set by the swap of a boundary edge: the two faces
# that contained the edge are replaced by the two that contain the new edge.
struct SurfaceSwap
    side_set::Int
    removed::Vector{NTuple{3,Int}}
    added::Vector{NTuple{3,Int}}
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
    surface::Union{Nothing,SurfaceSwap}
end

function CavityProposal(
    old_elements, new_connectivity, block, energy_before, energy_after, removed_node, surviving_node, split
)
    return CavityProposal(
        old_elements, new_connectivity, block, energy_before, energy_after, removed_node, surviving_node, split, nothing
    )
end

# The edges and faces whose operation a topology phase has refused, kept
# until an element around them changes: nodes do not move within a phase,
# so a refused proposal stays refused until its cavity does.
struct PhaseMemory
    swaps::Set{Tuple{Int,Int}}
    boundary_swaps::Set{Tuple{Int,Int}}
    face_swaps::Set{NTuple{3,Int}}
    collapses::Set{Tuple{Int,Int}}
    splits::Set{Tuple{Int,Int}}
end
function PhaseMemory()
    return PhaseMemory(
        Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int}}(), Set{NTuple{3,Int}}(), Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int}}()
    )
end

# Edges and faces created by the operations of a pass, which the adjacency
# does not know until the pass ends.
struct CreatedEntities
    edges::Set{Tuple{Int,Int}}
    faces::Set{NTuple{3,Int}}
end
CreatedEntities() = CreatedEntities(Set{Tuple{Int,Int}}(), Set{NTuple{3,Int}}())
