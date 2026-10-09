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
#
# The positions and the connectivity grow by one node or a few elements per
# accepted operation.  They are stored with spare columns, whose number is
# doubled when they run out, and read as `topology.positions` (3 × num_nodes)
# and `topology.connectivity` (4 × num_elements, positively oriented), views
# of the columns in use.  Appending a column by `hcat` copied the whole
# matrix, so the cost of a pass grew with the product of the mesh size and
# the number of accepted operations (issue #231).
mutable struct MeshTopology
    position_storage::Matrix{Float64}
    num_positions::Int
    connectivity_storage::Matrix{Int}
    num_connectivity::Int
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

function MeshTopology(positions::AbstractMatrix{Float64}, connectivity::AbstractMatrix{Int}, rest...)
    return MeshTopology(
        Matrix(positions), size(positions, 2), Matrix(connectivity), size(connectivity, 2), rest...
    )
end

function Base.getproperty(topology::MeshTopology, name::Symbol)
    if name === :positions
        return view(getfield(topology, :position_storage), :, 1:getfield(topology, :num_positions))
    elseif name === :connectivity
        return view(getfield(topology, :connectivity_storage), :, 1:getfield(topology, :num_connectivity))
    end
    return getfield(topology, name)
end

function Base.setproperty!(topology::MeshTopology, name::Symbol, value)
    if name === :positions
        setfield!(topology, :position_storage, Matrix{Float64}(value))
        return setfield!(topology, :num_positions, size(value, 2))
    elseif name === :connectivity
        setfield!(topology, :connectivity_storage, Matrix{Int}(value))
        return setfield!(topology, :num_connectivity, size(value, 2))
    end
    return setfield!(topology, name, convert(fieldtype(MeshTopology, name), value))
end

# A matrix of nodal columns extended by one column for a node that is not
# yet added, so that a split can be evaluated without copying the data: the
# positions of a topology, or the nodal values of a metric.
struct MatrixWithColumn{N,M<:AbstractMatrix{Float64}} <: AbstractMatrix{Float64}
    base::M
    extra::SVector{N,Float64}
end
function MatrixWithColumn{N}(base::AbstractMatrix{Float64}, extra::AbstractVector{Float64}) where {N}
    return MatrixWithColumn{N,typeof(base)}(base, SVector{N,Float64}(extra))
end
Base.size(m::MatrixWithColumn{N}) where {N} = (N, size(m.base, 2) + 1)
function Base.getindex(m::MatrixWithColumn, i::Int, j::Int)
    return j ≤ size(m.base, 2) ? m.base[i, j] : m.extra[i]
end
const PositionsWithNode = MatrixWithColumn{3}
