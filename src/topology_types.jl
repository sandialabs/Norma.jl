# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# In-memory topology of a four-node tetrahedral mesh for the adaptivity loop
# (docs/notes/ems-adaptivity).  The Exodus database stays the source of truth
# for the simulation; this structure is built from it at the start of a
# topology phase, modified through tombstones (dead nodes and elements are
# marked, new ones appended), compacted and written back at the end.  The
# adjacency arrays describe the alive elements at the last call of
# build_adjacency!; within a pass an operation may make them stale around the
# cavity it changed, which the driver handles by deferring cavities that touch
# modified elements to the next pass.
mutable struct MeshTopology
    positions::Matrix{Float64}      # 3 × num_nodes, current coordinates
    connectivity::Matrix{Int}       # 4 × num_elements, positively oriented
    block::Vector{Int}              # block index of every element
    block_ids::Vector{Int}          # Exodus block id per block index
    block_names::Vector{String}
    node_alive::BitVector
    element_alive::BitVector
    # Node-to-element adjacency in compressed sparse row form
    node_to_elements_offsets::Vector{Int}
    node_to_elements::Vector{Int}
    # Edge (sorted node pair) and face (sorted node triple) to incident elements
    edges::Dict{Tuple{Int,Int},Vector{Int}}
    faces::Dict{NTuple{3,Int},Vector{Int}}
    boundary_edges::Set{Tuple{Int,Int}}
    # Sets: membership flags per node, boundary faces per side set, and the
    # nodes of each side set derived from its faces
    node_sets::Dict{Int,BitVector}
    node_set_names::Dict{Int,String}
    side_sets::Dict{Int,Set{NTuple{3,Int}}}
    side_set_names::Dict{Int,String}
    node_side_sets::Dict{Int,BitVector}
end

# A matrix of nodal columns extended by one column for a node that is not
# yet added, so that a split can be evaluated without copying the data: the
# positions of a topology, or the nodal values of a metric.
struct MatrixWithColumn{N} <: AbstractMatrix{Float64}
    base::Matrix{Float64}
    extra::SVector{N,Float64}
end
Base.size(m::MatrixWithColumn{N}) where {N} = (N, size(m.base, 2) + 1)
function Base.getindex(m::MatrixWithColumn, i::Int, j::Int)
    return j ≤ size(m.base, 2) ? m.base[i, j] : m.extra[i]
end
const PositionsWithNode = MatrixWithColumn{3}
