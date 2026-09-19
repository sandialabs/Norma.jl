# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Exodus side numbering of a four-node tetrahedron: local nodes of each side.
const TETRA4_SIDES = ((1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2))

sorted_edge(a::Integer, b::Integer) = a < b ? (Int(a), Int(b)) : (Int(b), Int(a))

function sorted_face(a::Integer, b::Integer, c::Integer)
    x, y, z = Int(a), Int(b), Int(c)
    if x > y
        x, y = y, x
    end
    if y > z
        y, z = z, y
    end
    if x > y
        x, y = y, x
    end
    return (x, y, z)
end

# Faces of element `e` as sorted triples, in the Exodus side order.
function element_faces(topology::MeshTopology, e::Int)
    c = view(topology.connectivity, :, e)
    return ntuple(s -> sorted_face(c[TETRA4_SIDES[s][1]], c[TETRA4_SIDES[s][2]], c[TETRA4_SIDES[s][3]]), 4)
end

# Local Exodus side number of the face `face` of element `e`, or 0.
function local_side(topology::MeshTopology, e::Int, face::NTuple{3,Int})
    faces = element_faces(topology, e)
    for s in 1:4
        faces[s] == face && return s
    end
    return 0
end

function tetrahedron_volume(x::AbstractMatrix{Float64})
    u = SVector{3,Float64}(x[1, 2] - x[1, 1], x[2, 2] - x[2, 1], x[3, 2] - x[3, 1])
    v = SVector{3,Float64}(x[1, 3] - x[1, 1], x[2, 3] - x[2, 1], x[3, 3] - x[3, 1])
    w = SVector{3,Float64}(x[1, 4] - x[1, 1], x[2, 4] - x[2, 1], x[3, 4] - x[3, 1])
    return dot(u, cross(v, w)) / 6.0
end

function element_volume(topology::MeshTopology, e::Int)
    return tetrahedron_volume(view(topology.positions, :, view(topology.connectivity, :, e)))
end

# Build the topology of a smoothing model from its mesh at the current
# positions.  Every block must be TETRA4.
function build_topology(model::SolidMechanics)
    # The model's handle may already be closed when the topology is built
    # after a smoothing run, so the sets are read from a fresh read-only
    # handle on the same file.
    mesh = Exodus.ExodusDatabase(model.mesh.file_name, "r")
    positions = model.reference + model.displacement
    num_nodes = size(positions, 2)
    connectivity = zeros(Int, 4, 0)
    block = Int[]
    block_ids = Int[]
    block_names = String[]
    for (block_index, block_data) in enumerate(model.blocks)
        block_data.element_type == TETRA4 || norma_abort("The adaptivity loop supports four-node tetrahedra only")
        connectivity = hcat(connectivity, Int.(block_data.connectivity))
        append!(block, fill(block_index, block_data.num_elements))
        push!(block_ids, Int(block_data.id))
        push!(block_names, Exodus.read_name(mesh, Block, block_data.id))
    end
    num_elements = size(connectivity, 2)
    node_sets = Dict{Int,BitVector}()
    node_set_names = Dict{Int,String}()
    for id in Exodus.read_ids(mesh, NodeSet)
        flags = falses(num_nodes)
        flags[Int.(Exodus.read_node_set_nodes(mesh, id))] .= true
        node_sets[Int(id)] = flags
        node_set_names[Int(id)] = Exodus.read_name(mesh, NodeSet, id)
    end
    side_sets = Dict{Int,Set{NTuple{3,Int}}}()
    side_set_names = Dict{Int,String}()
    for id in Exodus.read_ids(mesh, SideSet)
        counts, nodes = Exodus.read_side_set_node_list(mesh, id)
        faces = Set{NTuple{3,Int}}()
        offset = 0
        for count in counts
            count == 3 || norma_abort("Side set $id has a side with $count nodes; only triangular sides are supported")
            push!(faces, sorted_face(nodes[offset + 1], nodes[offset + 2], nodes[offset + 3]))
            offset += count
        end
        side_sets[Int(id)] = faces
        side_set_names[Int(id)] = Exodus.read_name(mesh, SideSet, id)
    end
    Exodus.close(mesh)
    topology = MeshTopology(
        positions,
        connectivity,
        block,
        block_ids,
        block_names,
        trues(num_nodes),
        trues(num_elements),
        Int[],
        Int[],
        Dict{Tuple{Int,Int},Vector{Int}}(),
        Dict{NTuple{3,Int},Vector{Int}}(),
        Set{Tuple{Int,Int}}(),
        node_sets,
        node_set_names,
        side_sets,
        side_set_names,
        Dict{Int,BitVector}(),
    )
    for e in 1:num_elements
        element_volume(topology, e) > 0.0 || norma_abort("Element $e of the input mesh is inverted or degenerate")
    end
    build_adjacency!(topology)
    return topology
end

# Topology from arrays, without sets: one block, for tests and prototypes.
function build_topology(
    positions::Matrix{Float64}, connectivity::Matrix{Int}; block_id::Int=1, block_name::String="block"
)
    num_nodes = size(positions, 2)
    num_elements = size(connectivity, 2)
    topology = MeshTopology(
        copy(positions),
        copy(connectivity),
        fill(1, num_elements),
        [block_id],
        [block_name],
        trues(num_nodes),
        trues(num_elements),
        Int[],
        Int[],
        Dict{Tuple{Int,Int},Vector{Int}}(),
        Dict{NTuple{3,Int},Vector{Int}}(),
        Set{Tuple{Int,Int}}(),
        Dict{Int,BitVector}(),
        Dict{Int,String}(),
        Dict{Int,Set{NTuple{3,Int}}}(),
        Dict{Int,String}(),
        Dict{Int,BitVector}(),
    )
    for e in 1:num_elements
        element_volume(topology, e) > 0.0 || norma_abort("Element $e is inverted or degenerate")
    end
    build_adjacency!(topology)
    return topology
end

# Rebuild every adjacency from the alive elements.
function build_adjacency!(topology::MeshTopology)
    num_nodes = size(topology.positions, 2)
    num_elements = size(topology.connectivity, 2)
    counts = zeros(Int, num_nodes)
    for e in 1:num_elements
        topology.element_alive[e] || continue
        for i in 1:4
            counts[topology.connectivity[i, e]] += 1
        end
    end
    offsets = Vector{Int}(undef, num_nodes + 1)
    offsets[1] = 1
    for n in 1:num_nodes
        offsets[n + 1] = offsets[n] + counts[n]
    end
    fill!(counts, 0)
    elements = Vector{Int}(undef, offsets[end] - 1)
    for e in 1:num_elements
        topology.element_alive[e] || continue
        for i in 1:4
            n = topology.connectivity[i, e]
            elements[offsets[n] + counts[n]] = e
            counts[n] += 1
        end
    end
    topology.node_to_elements_offsets = offsets
    topology.node_to_elements = elements
    edges = Dict{Tuple{Int,Int},Vector{Int}}()
    faces = Dict{NTuple{3,Int},Vector{Int}}()
    sizehint!(edges, 7 * num_nodes)
    sizehint!(faces, 2 * num_elements)
    for e in 1:num_elements
        topology.element_alive[e] || continue
        c = view(topology.connectivity, :, e)
        for i in 1:4, j in (i + 1):4
            push!(get!(() -> Int[], edges, sorted_edge(c[i], c[j])), e)
        end
        for face in element_faces(topology, e)
            push!(get!(() -> Int[], faces, face), e)
        end
    end
    topology.edges = edges
    topology.faces = faces
    boundary_edges = Set{Tuple{Int,Int}}()
    node_side_sets = Dict{Int,BitVector}(id => falses(num_nodes) for id in keys(topology.side_sets))
    for (face, incident) in faces
        length(incident) ≤ 2 || norma_abort("Face $face is shared by $(length(incident)) elements")
        length(incident) == 1 || continue
        push!(boundary_edges, sorted_edge(face[1], face[2]))
        push!(boundary_edges, sorted_edge(face[1], face[3]))
        push!(boundary_edges, sorted_edge(face[2], face[3]))
        for (id, side_faces) in topology.side_sets
            if face in side_faces
                flags = node_side_sets[id]
                flags[face[1]] = flags[face[2]] = flags[face[3]] = true
            end
        end
    end
    topology.boundary_edges = boundary_edges
    topology.node_side_sets = node_side_sets
    return topology
end

num_alive_nodes(topology::MeshTopology) = count(topology.node_alive)
num_alive_elements(topology::MeshTopology) = count(topology.element_alive)

# Elements around a node and around an edge, from the adjacency of the last build.
function node_elements(topology::MeshTopology, n::Int)
    offsets = topology.node_to_elements_offsets
    return view(topology.node_to_elements, offsets[n]:(offsets[n + 1] - 1))
end

function edge_elements(topology::MeshTopology, a::Int, b::Int)
    return get(topology.edges, sorted_edge(a, b), Int[])
end

is_boundary_edge(topology::MeshTopology, a::Int, b::Int) = sorted_edge(a, b) in topology.boundary_edges

function boundary_faces(topology::MeshTopology)
    return [face for (face, incident) in topology.faces if length(incident) == 1]
end

# Euler characteristic V - E + F - T of the alive mesh from the adjacency of
# the last build; one for a mesh of a ball.
function euler_characteristic(topology::MeshTopology)
    return num_alive_nodes(topology) - length(topology.edges) + length(topology.faces) - num_alive_elements(topology)
end

# Append a node; it inherits the set memberships given.  Returns its index.
function add_node!(
    topology::MeshTopology,
    position::AbstractVector{Float64};
    node_sets::Vector{Int}=Int[],
    side_sets::Vector{Int}=Int[],
)
    topology.positions = hcat(topology.positions, Vector{Float64}(position))
    push!(topology.node_alive, true)
    n = size(topology.positions, 2)
    for (id, flags) in topology.node_sets
        push!(flags, id in node_sets)
    end
    for (id, flags) in topology.node_side_sets
        push!(flags, id in side_sets)
    end
    return n
end

# Append elements (4 × m connectivity into the current nodes) to a block.
# Returns their indices.  Side-set faces of the new elements are added by the
# caller through add_side_set_face!, since only the operator knows which of
# its faces lie on the boundary.
function add_elements!(topology::MeshTopology, connectivity::AbstractMatrix{<:Integer}, block::Int)
    size(connectivity, 1) == 4 || norma_abort("Elements must be four-node tetrahedra")
    first = size(topology.connectivity, 2) + 1
    topology.connectivity = hcat(topology.connectivity, Int.(connectivity))
    m = size(connectivity, 2)
    append!(topology.block, fill(block, m))
    append!(topology.element_alive, trues(m))
    for e in first:(first + m - 1)
        element_volume(topology, e) > 0.0 || norma_abort("Element added to the topology is inverted or degenerate")
    end
    return collect(first:(first + m - 1))
end

function remove_elements!(topology::MeshTopology, elements::AbstractVector{<:Integer})
    for e in elements
        topology.element_alive[e] = false
    end
    return topology
end

function remove_node!(topology::MeshTopology, n::Int)
    topology.node_alive[n] = false
    return topology
end

function add_side_set_face!(topology::MeshTopology, id::Int, face::NTuple{3,Int})
    push!(topology.side_sets[id], face)
    return topology
end

function remove_side_set_face!(topology::MeshTopology, id::Int, face::NTuple{3,Int})
    delete!(topology.side_sets[id], face)
    return topology
end

# Renumber the alive nodes and elements contiguously, drop the dead ones, and
# rebuild the adjacency.  Returns the maps from old to new indices (zero for a
# dropped entity).
function compact!(topology::MeshTopology)
    node_map = zeros(Int, length(topology.node_alive))
    k = 0
    for n in 1:length(topology.node_alive)
        if topology.node_alive[n]
            k += 1
            node_map[n] = k
        end
    end
    element_map = zeros(Int, length(topology.element_alive))
    k = 0
    for e in 1:length(topology.element_alive)
        if topology.element_alive[e]
            k += 1
            element_map[e] = k
        end
    end
    alive_nodes = findall(topology.node_alive)
    alive_elements = findall(topology.element_alive)
    topology.positions = topology.positions[:, alive_nodes]
    connectivity = topology.connectivity[:, alive_elements]
    for i in eachindex(connectivity)
        connectivity[i] = node_map[connectivity[i]]
        connectivity[i] > 0 || norma_abort("An alive element refers to a dead node")
    end
    topology.connectivity = connectivity
    topology.block = topology.block[alive_elements]
    topology.node_alive = trues(length(alive_nodes))
    topology.element_alive = trues(length(alive_elements))
    for (id, flags) in topology.node_sets
        topology.node_sets[id] = flags[alive_nodes]
    end
    for (id, faces) in topology.side_sets
        renumbered = Set{NTuple{3,Int}}()
        for face in faces
            a, b, c = node_map[face[1]], node_map[face[2]], node_map[face[3]]
            (a > 0 && b > 0 && c > 0) && push!(renumbered, sorted_face(a, b, c))
        end
        topology.side_sets[id] = renumbered
    end
    build_adjacency!(topology)
    return node_map, element_map
end

# Write the compacted topology as a new Exodus mesh: blocks, node sets, and
# side sets, with their names.  Boundary faces of a side set that no longer
# exist are dropped.
function write_topology(topology::MeshTopology, file_name::String)
    all(topology.node_alive) && all(topology.element_alive) || norma_abort("Compact the topology before writing it")
    num_nodes = size(topology.positions, 2)
    num_elements = size(topology.connectivity, 2)
    num_blocks = length(topology.block_ids)
    init = Exodus.Initialization{Int32}(
        Int32(3),
        Int32(num_nodes),
        Int32(num_elements),
        Int32(num_blocks),
        Int32(length(topology.node_sets)),
        Int32(length(topology.side_sets)),
    )
    isfile(file_name) && rm(file_name; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(file_name, "w", init)
    Exodus.write_coordinates(exo, topology.positions)
    # Elements are written block by block; Exodus numbers them globally in
    # that order, which the side sets refer to.
    global_index = zeros(Int, num_elements)
    next = 0
    for (block_index, block_id) in enumerate(topology.block_ids)
        elements = findall(==(block_index), topology.block)
        for e in elements
            next += 1
            global_index[e] = next
        end
        Exodus.write_block(exo, block_id, "TETRA4", Matrix{Int32}(topology.connectivity[:, elements]))
        Exodus.write_name(exo, Block, block_id, topology.block_names[block_index])
    end
    for id in sort(collect(keys(topology.node_sets)))
        node_set = Exodus.NodeSet(Int32(id), Vector{Int32}(findall(topology.node_sets[id])))
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, topology.node_set_names[id])
    end
    for id in sort(collect(keys(topology.side_sets)))
        elements = Int32[]
        sides = Int32[]
        for face in sort(collect(topology.side_sets[id]))
            incident = get(topology.faces, face, Int[])
            length(incident) == 1 || continue
            e = incident[1]
            push!(elements, Int32(global_index[e]))
            push!(sides, Int32(local_side(topology, e, face)))
        end
        side_set = Exodus.SideSet{Int32,Vector{Int32}}(Int32(id), elements, sides, Int32[], Int32[])
        Exodus.write_set(exo, side_set)
        Exodus.write_name(exo, side_set, topology.side_set_names[id])
    end
    Exodus.close(exo)
    return file_name
end
