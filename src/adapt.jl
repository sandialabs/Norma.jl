# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Topological operations of the adaptivity loop (docs/notes/ems-adaptivity):
# cavity proposals accepted by one test on the smoothing energy, or, for the
# size operations under `size criterion: length`, on the edge length in the
# prescribed target with the energy as a validity check only, and, for the
# shape operations under `shape criterion: scaled Jacobian`, on the minimum
# scaled Jacobian of the cavity.

function AdaptivityOptions(params::Parameters)
    size_criterion = get(params, "size criterion", "energy")
    if !(size_criterion in ("energy", "length"))
        norma_abort("\"size criterion\" in \"adaptivity\" must be \"energy\" or \"length\"; got \"$size_criterion\"")
    end
    shape_criterion = get(params, "shape criterion", "energy")
    if !(shape_criterion in ("energy", "scaled Jacobian"))
        norma_abort(
            "\"shape criterion\" in \"adaptivity\" must be \"energy\" or \"scaled Jacobian\"; got \"$shape_criterion\"",
        )
    end
    return AdaptivityOptions(
        Float64(get(params, "desired energy density", 0.1)),
        Float64(get(params, "allowed energy density", Inf)),
        Float64(get(params, "minimum scaled Jacobian", 0.0)),
        Float64(get(params, "minimum decrease", 1.0e-08)),
        Int(get(params, "adjacency layers", 4)),
        Int(get(params, "maximum passes", 20)),
        Int(get(params, "outer iterations", 5)),
        Bool(get(params, "swaps", true)),
        Bool(get(params, "collapses", true)),
        Bool(get(params, "splits", true)),
        size_criterion == "length",
        shape_criterion == "scaled Jacobian",
        Bool(get(params, "face swaps", false)),
        Bool(get(params, "boundary swaps", false)),
        deg2rad(Float64(get(params, "boundary swap angle", 20.0))),
        Float64(get(params, "desired scaled Jacobian", 0.9)),
        Bool(get(params, "size operations first", false)),
        Bool(get(params, "shape operations every pass", false)),
    )
end

# Forget the refusals of every edge and face of the elements added since
# `first_new`, whose cavities have changed, and drop the dead ones.
function forget_changed!(memory::PhaseMemory, topology::MeshTopology, first_new::Int)
    for e in first_new:length(topology.element_alive)
        c = view(topology.connectivity, :, e)
        for i in 1:4, j in (i + 1):4
            edge = sorted_edge(c[i], c[j])
            delete!(memory.swaps, edge)
            delete!(memory.boundary_swaps, edge)
            delete!(memory.collapses, edge)
            delete!(memory.splits, edge)
        end
        for face in element_faces(topology, e)
            delete!(memory.face_swaps, face)
        end
    end
    return memory
end

# Renumber the memory after a compaction.
function remap!(memory::PhaseMemory, node_map::Vector{Int})
    all(i -> node_map[i] == i, eachindex(node_map)) && return memory
    for set in (memory.swaps, memory.boundary_swaps, memory.collapses, memory.splits)
        kept = Tuple{Int,Int}[]
        for (a, b) in set
            (node_map[a] > 0 && node_map[b] > 0) && push!(kept, sorted_edge(node_map[a], node_map[b]))
        end
        empty!(set)
        union!(set, kept)
    end
    kept_faces = NTuple{3,Int}[]
    for face in memory.face_swaps
        all(n -> node_map[n] > 0, face) || continue
        push!(kept_faces, sorted_face(node_map[face[1]], node_map[face[2]], node_map[face[3]]))
    end
    empty!(memory.face_swaps)
    union!(memory.face_swaps, kept_faces)
    return memory
end

# Relative increase of the cavity minimum of the scaled Jacobian that an
# operation must achieve under `shape criterion: scaled Jacobian`.
const MINIMUM_QUALITY_INCREASE = 1.0e-6

# Scaled Jacobians of the elements of a connectivity, in the metric space of
# the target when a metric field is prescribed: each element is mapped by
# the factor F_M of the metric sampled on it, so that an element that
# matches an anisotropic target is regular and measures one, and the shape
# criteria of the loop cannot fight the smoother over the elongation the
# target asks for.  Without a metric field the measure is the geometric one
# (a scalar size scales all edges alike and leaves it unchanged).
function quality_jacobians(
    model::SolidMechanics, positions::AbstractMatrix{Float64}, connectivity::AbstractMatrix{<:Integer}
)
    model.metric_field === nothing && return scaled_jacobians(positions, connectivity)
    quality = Vector{Float64}(undef, size(connectivity, 2))
    for e in 1:size(connectivity, 2)
        node_indices = view(connectivity, :, e)
        X = SMatrix{3,4,Float64,12}(positions[i, node_indices[j]] for i in 1:3, j in 1:4)
        _, F_M, _ = create_metric_reference(model.metric_field, X, model.time; node_indices)
        quality[e] = tetrahedron_scaled_jacobian(F_M * X)
    end
    return quality
end

# Whether the new elements of a cavity raise its minimum scaled Jacobian.
function raises_minimum_quality(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
)
    old_minimum = minimum(quality_jacobians(model, topology.positions, topology.connectivity[:, old_elements]))
    new_minimum = minimum(quality_jacobians(model, positions, new_connectivity))
    return new_minimum > old_minimum * (1.0 + MINIMUM_QUALITY_INCREASE)
end

# Volumes of the ideal elements of the given connectivity, with the target
# sampled at `sample_positions`.
function ideal_element_volumes(
    model::SolidMechanics,
    block_index::Int,
    connectivity::AbstractMatrix{<:Integer},
    sample_positions::AbstractMatrix{Float64};
    metric::Union{MetricField,Nothing}=model.metric_field,
)
    element_type = model.blocks[block_index].element_type
    volumes = Vector{Float64}(undef, size(connectivity, 2))
    for e in 1:size(connectivity, 2)
        node_indices = view(connectivity, :, e)
        sample = SMatrix{3,4,Float64,12}(sample_positions[i, node_indices[j]] for i in 1:3, j in 1:4)
        X, _, _ = smoothing_reference(model, element_type, sample; node_indices, metric)
        volumes[e] = abs(tetrahedron_volume(X))
    end
    return volumes
end

# Energy densities (energy per unit ideal volume) of the alive elements of a
# topology at its current positions; dead elements get NaN.
function energy_densities(model::SolidMechanics, topology::MeshTopology)
    alive = findall(topology.element_alive)
    densities = fill(NaN, length(topology.element_alive))
    for block_index in unique(topology.block[alive])
        elements = [e for e in alive if topology.block[e] == block_index]
        connectivity = topology.connectivity[:, elements]
        energies = element_energies(model, block_index, connectivity, topology.positions)
        volumes = ideal_element_volumes(model, block_index, connectivity, topology.positions)
        densities[elements] = energies ./ volumes
    end
    return densities
end

# The geometric floor: an operation may not create an element whose scaled
# Jacobian is below the floor, unless the worst element of the old cavity was
# already below it and the new worst is no worse.  The energy sums the cavity,
# so this is what keeps a single element from being sacrificed for the sum.
function passes_scaled_jacobian_floor(old_minimum::Float64, new_minimum::Float64, floor::Float64)
    new_minimum ≥ floor && return true
    return new_minimum ≥ old_minimum
end

function passes_scaled_jacobian_floor(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
    floor::Float64,
)
    floor ≤ 0.0 && return true
    old_minimum = minimum(quality_jacobians(model, topology.positions, topology.connectivity[:, old_elements]))
    new_minimum = minimum(quality_jacobians(model, positions, new_connectivity))
    return passes_scaled_jacobian_floor(old_minimum, new_minimum, floor)
end

# The acceptance test: a proposal is accepted when the energy of the new
# elements is below that of the old ones by the relative margin, no new
# element exceeds the allowed density, and the geometric floor holds.
# Returns the proposal or nothing.
function accept_proposal(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    options::AdaptivityOptions;
    require_decrease::Bool=true,
)
    positions = topology.positions
    old_energy = sum(element_energies(model, block, topology.connectivity[:, old_elements], positions))
    new_energies = element_energies(model, block, new_connectivity, positions)
    all(isfinite, new_energies) || return nothing
    new_energy = sum(new_energies)
    if require_decrease
        if options.shape_by_quality
            raises_minimum_quality(model, topology, old_elements, new_connectivity, positions) || return nothing
        else
            new_energy ≤ old_energy - options.minimum_decrease * old_energy || return nothing
        end
    end
    if isfinite(options.allowed_density)
        volumes = ideal_element_volumes(model, block, new_connectivity, positions)
        maximum(new_energies ./ volumes) ≤ options.allowed_density || return nothing
    end
    floor = options.minimum_scaled_jacobian
    passes_scaled_jacobian_floor(model, topology, old_elements, new_connectivity, positions, floor) || return nothing
    return CavityProposal(old_elements, new_connectivity, block, old_energy, new_energy, 0, 0, nothing)
end

function apply!(topology::MeshTopology, proposal::CavityProposal; metric::Union{MetricField,Nothing}=nothing)
    if proposal.split !== nothing
        split = proposal.split
        m = add_node!(topology, split.position; node_sets=split.node_sets, side_sets=split.side_sets)
        m == size(topology.positions, 2) || norma_abort("The split node did not receive the expected index")
        split_sets!(topology, split.edge[1], split.edge[2], m)
        add_metric_node!(metric, split.metric)
    end
    remove_elements!(topology, proposal.old_elements)
    add_elements!(topology, proposal.new_connectivity, proposal.block)
    if proposal.surface !== nothing
        for face in proposal.surface.removed
            remove_side_set_face!(topology, proposal.surface.side_set, face)
        end
        for face in proposal.surface.added
            add_side_set_face!(topology, proposal.surface.side_set, face)
        end
    end
    if proposal.removed_node > 0
        collapse_sets!(topology, proposal.removed_node, proposal.surviving_node)
        remove_node!(topology, proposal.removed_node)
    end
    return topology
end

# Elements around an interior edge in cyclic order, and the ring nodes p_1,
# ..., p_n such that element i has nodes {a, b, p_i, p_(i+1)}.  Returns
# nothing when the edge is on the boundary, when any incident element is dead,
# or when the elements do not form a single closed ring.
function edge_ring(topology::MeshTopology, a::Int, b::Int)
    is_boundary_edge(topology, a, b) && return nothing
    elements = edge_elements(topology, a, b)
    n = length(elements)
    n ≥ 3 || return nothing
    all(topology.element_alive[e] for e in elements) || return nothing
    others = Vector{Tuple{Int,Int}}(undef, n)
    for (i, e) in enumerate(elements)
        c = view(topology.connectivity, :, e)
        pq = Int[]
        for node in c
            (node == a || node == b) || push!(pq, node)
        end
        length(pq) == 2 || return nothing
        others[i] = (pq[1], pq[2])
    end
    used = falses(n)
    ordered = Int[elements[1]]
    ring = Int[others[1][1], others[1][2]]
    used[1] = true
    for _ in 2:n
        last = ring[end]
        found = 0
        for i in 1:n
            used[i] && continue
            if others[i][1] == last || others[i][2] == last
                found = i
                break
            end
        end
        found == 0 && return nothing
        used[found] = true
        push!(ordered, elements[found])
        push!(ring, others[found][1] == last ? others[found][2] : others[found][1])
    end
    ring[end] == ring[1] || return nothing
    pop!(ring)
    length(unique(ring)) == n || return nothing
    return ordered, ring
end

# All triangulations of the polygon with vertices 1..n, as lists of triangles
# of vertex indices (Catalan number of the polygon).
function polygon_triangulations(n::Int)
    return chain_triangulations(1, n)
end

function chain_triangulations(i::Int, j::Int)
    j - i < 2 && return [Vector{NTuple{3,Int}}()]
    result = Vector{Vector{NTuple{3,Int}}}()
    for k in (i + 1):(j - 1)
        for left in chain_triangulations(i, k), right in chain_triangulations(k, j)
            push!(result, vcat(left, [(i, k, j)], right))
        end
    end
    return result
end

# Connectivity of the elements that replace the ring of edge (a, b) for one
# triangulation of the ring polygon: two tetrahedra per triangle, one with
# each edge node as apex.  A triangulation tiles the cavity only if every
# triangle, taken in the cyclic order of the ring, has the orientation of the
# polygon in the projection along the edge (an ear at a reflex vertex lies
# outside the cavity even though its tetrahedra can be positively oriented),
# and if the edge nodes lie on opposite sides of every triangle.  Returns
# nothing when either fails.
function swapped_connectivity(
    topology::MeshTopology, a::Int, b::Int, ring::Vector{Int}, triangles::Vector{NTuple{3,Int}}
)
    positions = topology.positions
    xa = SVector{3,Float64}(view(positions, :, a))
    xb = SVector{3,Float64}(view(positions, :, b))
    d = xb - xa
    x(n) = SVector{3,Float64}(view(positions, :, n))
    # Orientation of the ring polygon in the projection along the edge, from
    # its first sector, which the old element guarantees to be well formed.
    s = sign(dot(cross(x(ring[1]) - xa, x(ring[2]) - xa), d))
    s == 0.0 && return nothing
    connectivity = zeros(Int, 4, 2 * length(triangles))
    column = 0
    for (i, j, k) in triangles
        p, q, r = ring[i], ring[j], ring[k]
        xp, xq, xr = x(p), x(q), x(r)
        normal = cross(xq - xp, xr - xp)
        sign(dot(normal, d)) == s || return nothing
        # With the triangle oriented toward b, (p, q, r, b) is positive and
        # (p, r, q, a) is positive when a lies on the other side.
        if s < 0.0
            q, r = r, q
            xq, xr = xr, xq
            normal = -normal
        end
        dot(normal, xb - xp) > 0.0 || return nothing
        dot(normal, xa - xp) < 0.0 || return nothing
        column += 1
        connectivity[:, column] .= (p, q, r, b)
        column += 1
        connectivity[:, column] .= (p, r, q, a)
    end
    return connectivity
end

# Try to swap the interior edge (a, b): evaluate every triangulation of its
# ring, keep the one of least energy, and submit it to the acceptance test.  Returns
# the accepted proposal or nothing.
function try_edge_swap(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    options::AdaptivityOptions;
    created::CreatedEntities=CreatedEntities(),
)
    ring = edge_ring(topology, a, b)
    ring === nothing && return nothing
    elements, ring_nodes = ring
    n = length(ring_nodes)
    block = topology.block[elements[1]]
    all(topology.block[e] == block for e in elements) || return nothing
    best = nothing
    best_score = Inf
    for triangles in polygon_triangulations(n)
        # A chord of the ring polygon becomes an edge; it may not exist
        # already, in the mesh or from an operation of this pass, since the
        # link of an edge must be a single ring.
        chords_are_new = true
        for (i, j, k) in triangles, (u, v) in ((i, j), (j, k), (k, i))
            mod(u - v, n) in (1, n - 1) && continue
            edge = sorted_edge(ring_nodes[u], ring_nodes[v])
            if haskey(topology.edges, edge) || edge in created.edges
                chords_are_new = false
                break
            end
        end
        chords_are_new || continue
        connectivity = swapped_connectivity(topology, a, b, ring_nodes, triangles)
        connectivity === nothing && continue
        faces_are_new(topology, elements, connectivity, created) || continue
        score = configuration_score(model, topology, block, connectivity, options)
        if score < best_score
            best_score = score
            best = connectivity
        end
    end
    best === nothing && return nothing
    return accept_proposal(model, topology, elements, best, block, options)
end

# Score of a candidate configuration of a swap, lower is better: the cavity
# energy, or under `shape criterion: scaled Jacobian` the negative of its
# minimum scaled Jacobian; infinite for an invalid configuration.
function configuration_score(
    model::SolidMechanics, topology::MeshTopology, block::Int, connectivity::Matrix{Int}, options::AdaptivityOptions
)
    # Under the scaled Jacobian criterion the energies of the chosen
    # configuration are computed once by the acceptance test, as its
    # validity check.
    options.shape_by_quality && return -minimum(quality_jacobians(model, topology.positions, connectivity))
    energies = element_energies(model, block, connectivity, topology.positions)
    all(isfinite, energies) || return Inf
    return sum(energies)
end


# Elements around a boundary edge (a, b) as an open chain from one boundary
# face to the other, with the chain nodes p_1, ..., p_n such that element i
# has nodes {a, b, p_i, p_(i+1)} and the faces (a, b, p_1) and (a, b, p_n) are
# on the boundary.  Returns nothing when the edge is interior, when an
# element is dead, or when the chain is not simple.
function boundary_edge_chain(topology::MeshTopology, a::Int, b::Int)
    is_boundary_edge(topology, a, b) || return nothing
    elements = edge_elements(topology, a, b)
    n = length(elements)
    n ≥ 2 || return nothing
    all(topology.element_alive[e] for e in elements) || return nothing
    others = Vector{Tuple{Int,Int}}(undef, n)
    for (i, e) in enumerate(elements)
        pq = Int[]
        for node in view(topology.connectivity, :, e)
            (node == a || node == b) || push!(pq, node)
        end
        length(pq) == 2 || return nothing
        others[i] = (pq[1], pq[2])
    end
    # The chain starts at an element with a boundary face through the edge.
    start = 0
    first_node = 0
    for (i, (p, q)) in enumerate(others)
        for node in (p, q)
            if length(get(topology.faces, sorted_face(a, b, node), Int[])) == 1
                start, first_node = i, node
                break
            end
        end
        start > 0 && break
    end
    start > 0 || return nothing
    used = falses(n)
    used[start] = true
    ordered = Int[elements[start]]
    chain = Int[first_node, others[start][1] == first_node ? others[start][2] : others[start][1]]
    for _ in 2:n
        last = chain[end]
        found = 0
        for i in 1:n
            used[i] && continue
            if others[i][1] == last || others[i][2] == last
                found = i
                break
            end
        end
        found == 0 && return nothing
        used[found] = true
        push!(ordered, elements[found])
        push!(chain, others[found][1] == last ? others[found][2] : others[found][1])
    end
    length(unique(chain)) == n + 1 || return nothing
    length(get(topology.faces, sorted_face(a, b, chain[end]), Int[])) == 1 || return nothing
    return ordered, chain
end

# Outward unit normal of the boundary face `face` of element `e`.
function outward_normal(topology::MeshTopology, e::Int, face::NTuple{3,Int})
    x(n) = SVector{3,Float64}(view(topology.positions, :, n))
    normal = cross(x(face[2]) - x(face[1]), x(face[3]) - x(face[1]))
    apex = first(n for n in view(topology.connectivity, :, e) if !(n in face))
    dot(normal, x(apex) - x(face[1])) > 0.0 && (normal = -normal)
    return normal / norm(normal)
end

# Try to swap the boundary edge (a, b) whose two boundary faces belong to one
# side set and make an angle below the tolerance, so that they form a flat
# patch: the edge is replaced by the edge (p_1, p_n) between the far nodes
# of the two faces, the faces (a, b, p_1) and (a, b, p_n) by (a, p_1, p_n)
# and (b, p_1, p_n), and the chain of elements behind them by the elements of
# a triangulation of the closed polygon p_1, ..., p_n.  Returns the accepted
# proposal or nothing.
function try_boundary_edge_swap(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    options::AdaptivityOptions;
    created::CreatedEntities=CreatedEntities(),
)
    chain = boundary_edge_chain(topology, a, b)
    chain === nothing && return nothing
    elements, nodes = chain
    n = length(nodes)
    n ≥ 3 || return nothing
    block = topology.block[elements[1]]
    all(topology.block[e] == block for e in elements) || return nothing
    old_faces = (sorted_face(a, b, nodes[1]), sorted_face(a, b, nodes[n]))
    # Both faces in one side set, or both in none.
    side_set = 0
    for (id, faces) in topology.side_sets
        in_set = (old_faces[1] in faces, old_faces[2] in faces)
        in_set[1] == in_set[2] || return nothing
        in_set[1] && (side_set = id)
    end
    normal_1 = outward_normal(topology, elements[1], old_faces[1])
    normal_2 = outward_normal(topology, elements[end], old_faces[2])
    cos_angle = clamp(dot(normal_1, normal_2), -1.0, 1.0)
    acos(cos_angle) ≤ options.boundary_swap_angle || return nothing
    new_edge = sorted_edge(nodes[1], nodes[n])
    (haskey(topology.edges, new_edge) || new_edge in created.edges) && return nothing
    best = nothing
    best_score = Inf
    for triangles in polygon_triangulations(n)
        chords_are_new = true
        for (i, j, k) in triangles, (u, v) in ((i, j), (j, k), (k, i))
            abs(u - v) == 1 && continue
            edge = sorted_edge(nodes[u], nodes[v])
            if haskey(topology.edges, edge) || edge in created.edges
                chords_are_new = false
                break
            end
        end
        chords_are_new || continue
        connectivity = swapped_connectivity(topology, a, b, nodes, triangles)
        connectivity === nothing && continue
        faces_are_new(topology, elements, connectivity, created) || continue
        score = configuration_score(model, topology, block, connectivity, options)
        if score < best_score
            best_score = score
            best = connectivity
        end
    end
    best === nothing && return nothing
    proposal = accept_proposal(model, topology, elements, best, block, options)
    proposal === nothing && return nothing
    new_faces = [sorted_face(a, nodes[1], nodes[n]), sorted_face(b, nodes[1], nodes[n])]
    surface = side_set > 0 ? SurfaceSwap(side_set, collect(old_faces), new_faces) : nothing
    return CavityProposal(
        proposal.old_elements,
        proposal.new_connectivity,
        proposal.block,
        proposal.energy_before,
        proposal.energy_after,
        0,
        0,
        nothing,
        surface,
    )
end

function record!(created::CreatedEntities, connectivity::AbstractMatrix{<:Integer})
    for column in eachcol(connectivity)
        for i in 1:4, j in (i + 1):4
            push!(created.edges, sorted_edge(column[i], column[j]))
        end
        for (i, j, k) in TETRA4_SIDES
            push!(created.faces, sorted_face(column[i], column[j], column[k]))
        end
    end
    return created
end

# Whether every face of the new elements that is not a face of the old ones
# is new to the mesh, so that no face ends up shared by more than two
# elements: the link of an edge must be a single ring.
function faces_are_new(
    topology::MeshTopology,
    old_elements::AbstractVector{<:Integer},
    connectivity::AbstractMatrix{<:Integer},
    created::CreatedEntities,
)
    old_faces = Set{NTuple{3,Int}}()
    for e in old_elements, face in element_faces(topology, e)
        push!(old_faces, face)
    end
    for column in eachcol(connectivity), (i, j, k) in TETRA4_SIDES
        face = sorted_face(column[i], column[j], column[k])
        face in old_faces && continue
        (haskey(topology.faces, face) || face in created.faces) && return false
    end
    return true
end

# Try to swap the interior face (p, q, r) shared by two elements with apexes
# d and e: the two elements are replaced by the three around the new edge
# (d, e), one for each edge of the face.  Valid when the segment (d, e)
# crosses the face, so that the three new elements are positive, and when
# the edge (d, e) does not exist already.  Returns the accepted proposal or
# nothing.
function try_face_swap(
    model::SolidMechanics,
    topology::MeshTopology,
    face::NTuple{3,Int},
    options::AdaptivityOptions;
    created::CreatedEntities=CreatedEntities(),
)
    incident = get(topology.faces, face, Int[])
    length(incident) == 2 || return nothing
    all(topology.element_alive[e] for e in incident) || return nothing
    block = topology.block[incident[1]]
    topology.block[incident[2]] == block || return nothing
    apex(e) = first(n for n in view(topology.connectivity, :, e) if !(n in face))
    d, e = apex(incident[1]), apex(incident[2])
    d == e && return nothing
    new_edge = sorted_edge(d, e)
    (haskey(topology.edges, new_edge) || new_edge in created.edges) && return nothing
    positions = topology.positions
    connectivity = zeros(Int, 4, 3)
    p, q, r = face
    for (k, (u, v)) in enumerate(((p, q), (q, r), (r, p)))
        tetra = [d, e, u, v]
        volume = tetrahedron_volume(positions[:, tetra])
        if volume < 0.0
            tetra[3], tetra[4] = tetra[4], tetra[3]
            volume = -volume
        end
        volume > 0.0 || return nothing
        connectivity[:, k] = tetra
    end
    faces_are_new(topology, incident, connectivity, created) || return nothing
    isfinite(configuration_score(model, topology, block, connectivity, options)) || return nothing
    return accept_proposal(model, topology, collect(incident), connectivity, block, options)
end

# Whether node b may be removed by collapsing it onto node a.  The survivor
# keeps its position, so every set membership of b must be one of a: a node
# of a node set cannot vanish unless a carries the same node sets, and a node
# on the boundary can only slide along the boundary onto a node of the same
# surfaces, along a boundary edge, so that the surface and its feature lines
# are preserved.  Nodes on three surfaces (corners) are never removed.
function may_collapse(topology::MeshTopology, b::Int, a::Int)
    for (id, flags) in topology.node_sets
        flags[b] && !flags[a] && return false
    end
    on_boundary = false
    for (id, flags) in topology.node_side_sets
        flags[b] || continue
        on_boundary = true
        flags[a] || return false
    end
    if on_boundary
        is_boundary_edge(topology, a, b) || return false
    end
    return true
end

# Try to remove node b by collapsing the edge (a, b) onto a: the elements
# that contain both nodes are removed and the others around b receive a in
# its place.  Returns the accepted proposal or nothing.
function try_edge_collapse(
    model::SolidMechanics,
    topology::MeshTopology,
    b::Int,
    a::Int,
    options::AdaptivityOptions;
    by_length::Bool=false,
    created::CreatedEntities=CreatedEntities(),
)
    (topology.node_alive[a] && topology.node_alive[b]) || return nothing
    may_collapse(topology, b, a) || return nothing
    star = collect(node_elements(topology, b))
    isempty(star) && return nothing
    all(topology.element_alive[e] for e in star) || return nothing
    block = topology.block[star[1]]
    all(topology.block[e] == block for e in star) || return nothing
    kept = Int[]
    for e in star
        a in view(topology.connectivity, :, e) || push!(kept, e)
    end
    isempty(kept) && return nothing
    new_connectivity = topology.connectivity[:, kept]
    for i in eachindex(new_connectivity)
        new_connectivity[i] == b && (new_connectivity[i] = a)
    end
    faces_are_new(topology, star, new_connectivity, created) || return nothing
    if options.size_by_length
        # No edge from the surviving node may become longer than the band.
        for q in unique(vec(new_connectivity))
            q == a && continue
            metric_edge_length(model, topology, a, q) ≤ LENGTH_BAND[2] || return nothing
        end
    end
    return accept_proposal(model, topology, star, new_connectivity, block, options, b, a; require_decrease=!by_length)
end

function accept_proposal(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    options::AdaptivityOptions,
    removed_node::Int,
    surviving_node::Int;
    require_decrease::Bool=true,
)
    proposal = accept_proposal(model, topology, old_elements, new_connectivity, block, options; require_decrease)
    proposal === nothing && return nothing
    return CavityProposal(
        proposal.old_elements,
        proposal.new_connectivity,
        proposal.block,
        proposal.energy_before,
        proposal.energy_after,
        removed_node,
        surviving_node,
        nothing,
    )
end

# The level-set closures of the Surface boundary conditions on the given side
# sets, for placing a new boundary node on its surfaces.
function surface_constraints(model::SolidMechanics, side_set_ids::Vector{Int})
    constraints = Tuple{Function,Function}[]
    for bc in model.boundary_conditions
        bc isa SolidMechanicsSurfaceBoundaryCondition || continue
        Int(bc.side_set_id) in side_set_ids || continue
        push!(constraints, (bc.level_set_fun, bc.level_set_grad))
    end
    return constraints
end

# Closest-point projection of a point onto its surfaces, as return_to_surface!.
function project_to_surfaces(x::SVector{3,Float64}, constraints::Vector{Tuple{Function,Function}}, time::Float64)
    isempty(constraints) && return x
    y = Vector{Float64}(x)
    for _ in 1:SURFACE_RETURN_MAX_ITERS
        _, g, A = surface_normal_frame(constraints, y, time)
        maximum(abs, g) ≤ SURFACE_RETURN_TOL && break
        y .-= permutedims(A) * ((A * permutedims(A)) \ g)
    end
    return SVector{3,Float64}(y)
end

# Energy of the elements around a new node as a function of its position.
function split_cavity_energy(
    model::SolidMechanics,
    topology::MeshTopology,
    block::Int,
    connectivity::Matrix{Int},
    position::SVector{3,Float64},
    metric::Union{MetricField,Nothing},
)
    positions = PositionsWithNode(topology.positions, position)
    return sum(element_energies(model, block, connectivity, positions; metric))
end

# Local relaxation of a new node with the rest of the cavity fixed: a few
# steps of descent on the energy of the new elements, the gradient by central
# differences of the position (three degrees of freedom), each step returned
# to the node's surfaces.  The operation is then judged at its locally
# relaxed position rather than at the midpoint.
const SPLIT_RELAXATION_STEPS = 5

function relax_split_node(
    model::SolidMechanics,
    topology::MeshTopology,
    block::Int,
    connectivity::Matrix{Int},
    position::SVector{3,Float64},
    constraints::Vector{Tuple{Function,Function}},
    scale::Float64,
    metric::Union{MetricField,Nothing},
)
    energy(x) = split_cavity_energy(model, topology, block, connectivity, x, metric)
    x = position
    value = energy(x)
    isfinite(value) || return x, value
    δ = 1.0e-6 * scale
    for _ in 1:SPLIT_RELAXATION_STEPS
        gradient = SVector{3,Float64}(
            (energy(x + SVector(δ, 0.0, 0.0)) - energy(x - SVector(δ, 0.0, 0.0))) / (2δ),
            (energy(x + SVector(0.0, δ, 0.0)) - energy(x - SVector(0.0, δ, 0.0))) / (2δ),
            (energy(x + SVector(0.0, 0.0, δ)) - energy(x - SVector(0.0, 0.0, δ))) / (2δ),
        )
        all(isfinite, gradient) || break
        step_length = 0.1 * scale / max(norm(gradient), floatmin(Float64))
        improved = false
        for _ in 1:8
            trial = project_to_surfaces(x - step_length * gradient, constraints, model.time)
            trial_value = energy(trial)
            if isfinite(trial_value) && trial_value < value
                x, value = trial, trial_value
                improved = true
                break
            end
            step_length *= 0.5
        end
        improved || break
    end
    return x, value
end

# Try to split the edge (a, b) at a new node: every element around the edge
# is bisected.  The node starts at the midpoint, is returned to the surfaces
# of the boundary faces that contain the edge, inherits the side sets of
# those faces and the node sets common to both ends of a boundary edge, is
# relaxed locally, and the result is submitted to the acceptance test.
function try_edge_split(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    options::AdaptivityOptions;
    by_length::Bool=false,
)
    (topology.node_alive[a] && topology.node_alive[b]) || return nothing
    ring = edge_elements(topology, a, b)
    isempty(ring) && return nothing
    all(topology.element_alive[e] for e in ring) || return nothing
    block = topology.block[ring[1]]
    all(topology.block[e] == block for e in ring) || return nothing
    side_sets = Int[]
    on_boundary = false
    for e in ring, face in element_faces(topology, e)
        (a in face && b in face) || continue
        length(get(topology.faces, face, Int[])) == 1 || continue
        on_boundary = true
        for (id, faces) in topology.side_sets
            face in faces && !(id in side_sets) && push!(side_sets, id)
        end
    end
    node_sets = Int[]
    if on_boundary
        for (id, flags) in topology.node_sets
            flags[a] && flags[b] && push!(node_sets, id)
        end
    end
    m = size(topology.positions, 2) + 1
    connectivity = zeros(Int, 4, 2 * length(ring))
    for (k, e) in enumerate(ring)
        c = topology.connectivity[:, e]
        connectivity[:, 2k - 1] = replace(c, b => m)
        connectivity[:, 2k] = replace(c, a => m)
    end
    xa = SVector{3,Float64}(view(topology.positions, :, a))
    xb = SVector{3,Float64}(view(topology.positions, :, b))
    # A metric carried by the nodes gives the new node the mean of the ends.
    node_metric = metric_node_data(model.metric_field, a, b)
    metric = metric_with_node(model.metric_field, node_metric)
    constraints = surface_constraints(model, side_sets)
    position = project_to_surfaces(0.5 * (xa + xb), constraints, model.time)
    if on_boundary && isempty(constraints)
        # A boundary edge whose surfaces have no analytic description: the
        # node stays at the midpoint, which lies on the boundary facets,
        # since a relaxation could move it off the boundary.
        energy_after = split_cavity_energy(model, topology, block, connectivity, position, metric)
    else
        position, energy_after =
            relax_split_node(model, topology, block, connectivity, position, constraints, norm(xb - xa), metric)
    end
    isfinite(energy_after) || return nothing
    energy_before = sum(element_energies(model, block, topology.connectivity[:, ring], topology.positions))
    if options.size_by_length
        # No new edge from the split node may be shorter than the band.
        for p in unique(vec(topology.connectivity[:, ring]))
            (p == a || p == b) && continue
            xp = SVector{3,Float64}(view(topology.positions, :, p))
            metric_length(model, position, xp, (a, b, p, p)) ≥ LENGTH_BAND[1] || return nothing
        end
    end
    positions = PositionsWithNode(topology.positions, position)
    if !by_length
        if options.shape_by_quality
            raises_minimum_quality(model, topology, collect(ring), connectivity, positions) || return nothing
        else
            energy_after ≤ energy_before - options.minimum_decrease * energy_before || return nothing
        end
    end
    if isfinite(options.allowed_density)
        energies = element_energies(model, block, connectivity, positions; metric)
        volumes = ideal_element_volumes(model, block, connectivity, positions; metric)
        maximum(energies ./ volumes) ≤ options.allowed_density || return nothing
    end
    floor = options.minimum_scaled_jacobian
    passes_scaled_jacobian_floor(model, topology, collect(ring), connectivity, positions, floor) || return nothing
    split = SplitNode(position, node_sets, side_sets, (a, b), node_metric)
    return CavityProposal(collect(ring), connectivity, block, energy_before, energy_after, 0, 0, split)
end

# Edges longer than sqrt(2) in the prescribed target, the split candidates of
# the size phase; empty without a target.
function long_edges(model::SolidMechanics, topology::MeshTopology)
    edges = Tuple{Int,Int}[]
    (model.metric_field === nothing && model.size_field === nothing) && return edges
    for edge in keys(topology.edges)
        metric_edge_length(model, topology, edge[1], edge[2]) > LENGTH_BAND[2] && push!(edges, edge)
    end
    return edges
end

# One pass of edge splits over the candidate edges, in decreasing order of
# the energy around them: the long edges of the target always, and the edges
# of the elements above the desired density when the shape-driven operations
# are on.
function split_pass!(
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions,
    shape::Bool;
    memory::PhaseMemory=PhaseMemory(),
)
    densities = energy_densities(model, topology)
    edges = Set{Tuple{Int,Int}}()
    if shape
        for e in candidate_elements(model, topology, densities, options)
            c = view(topology.connectivity, :, e)
            for i in 1:4, j in (i + 1):4
                push!(edges, sorted_edge(c[i], c[j]))
            end
        end
    end
    size_edges = Set(long_edges(model, topology))
    union!(edges, size_edges)
    setdiff!(edges, memory.splits)
    ranked = collect(edges)
    # The edges beyond the band are split longest first, in the target, as
    # in longest-edge bisection, which is what bounds the shape of the
    # children; the shape-driven candidates follow by their cavity rank.
    function split_rank(edge)
        (a, b) = edge
        if options.size_by_length && edge in size_edges
            return (1, metric_edge_length(model, topology, a, b))
        end
        return (0, cavity_rank(model, topology, densities, edge_elements(topology, a, b), options))
    end
    sort!(ranked; by=split_rank, rev=true)
    accepted = 0
    decrease = 0.0
    for (a, b) in ranked
        by_length = options.size_by_length && (a, b) in size_edges
        ring_intact = all(topology.element_alive[e] for e in edge_elements(topology, a, b))
        proposal = try_edge_split(model, topology, a, b, options; by_length)
        if proposal === nothing
            ring_intact && push!(memory.splits, (a, b))
            continue
        end
        apply!(topology, proposal; metric=model.metric_field)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# Length of the segment between two points measured in the prescribed
# target, or NaN without one: the scalar 1/h at the midpoint for a size
# field, F_M for a metric field, sampled at the midpoint by the function
# sources and averaged over `nodes` by the nodal sources (the nodes whose
# mean is the midpoint: the two ends of an edge, or (a, b, p, p) for the
# segment from the midpoint of edge (a, b) to a node p).
function metric_length(model::SolidMechanics, xa::SVector{3,Float64}, xb::SVector{3,Float64}, nodes)
    midpoint = 0.5 * (xa + xb)
    if model.metric_field !== nothing
        h, R = principal_metric(model.metric_field.source, nodes, midpoint, model.time)
        F_M = SMatrix{3,3,Float64,9}(Diagonal(SVector{3,Float64}(1.0 / h[1], 1.0 / h[2], 1.0 / h[3]))) * R'
        return norm(F_M * (xb - xa))
    elseif model.size_field !== nothing
        return norm(xb - xa) / model.size_field((model.time, midpoint[1], midpoint[2], midpoint[3]))
    end
    return NaN
end

# Length of an edge measured in the prescribed target, or NaN without one.
function metric_edge_length(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int)
    xa = SVector{3,Float64}(view(topology.positions, :, a))
    xb = SVector{3,Float64}(view(topology.positions, :, b))
    return metric_length(model, xa, xb, (a, b))
end

# The length band of the prescribed target: edges shorter than the lower
# bound are collapsed, edges longer than the upper bound are split, and under
# `size criterion: length` no collapse or split may create an edge outside
# the band, so that the operations cannot undo each other and the energy
# decides the shape within the band only.
const LENGTH_BAND = (1.0 / sqrt(2.0), sqrt(2.0))

# Edges shorter than 1/sqrt(2) in the prescribed target, the collapse
# candidates of the size phase; empty without a target.
function short_edges(model::SolidMechanics, topology::MeshTopology)
    edges = Tuple{Int,Int}[]
    (model.metric_field === nothing && model.size_field === nothing) && return edges
    for edge in keys(topology.edges)
        metric_edge_length(model, topology, edge[1], edge[2]) < LENGTH_BAND[1] && push!(edges, edge)
    end
    return edges
end

# One pass of edge collapses over the candidate edges, in decreasing order
# of the energy around them: the short edges of the target always, and the
# edges of the elements above the desired density when the shape-driven
# operations are on.  Both directions of each edge are tried.  Returns the
# number of accepted collapses and the decrease.
function collapse_pass!(
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions,
    shape::Bool;
    memory::PhaseMemory=PhaseMemory(),
)
    densities = energy_densities(model, topology)
    edges = Set{Tuple{Int,Int}}()
    if shape
        for e in candidate_elements(model, topology, densities, options)
            c = view(topology.connectivity, :, e)
            for i in 1:4, j in (i + 1):4
                push!(edges, sorted_edge(c[i], c[j]))
            end
        end
    end
    size_edges = Set(short_edges(model, topology))
    union!(edges, size_edges)
    setdiff!(edges, memory.collapses)
    ranked = collect(edges)
    edge_energy(edge) = cavity_rank(model, topology, densities, edge_elements(topology, edge[1], edge[2]), options)
    sort!(ranked; by=edge_energy, rev=true)
    accepted = 0
    decrease = 0.0
    created = CreatedEntities()
    for (a, b) in ranked
        (topology.node_alive[a] && topology.node_alive[b]) || continue
        by_length = options.size_by_length && (a, b) in size_edges
        stars_intact = all(topology.element_alive[e] for n in (a, b) for e in node_elements(topology, n))
        proposal = try_edge_collapse(model, topology, b, a, options; by_length, created)
        proposal === nothing && (proposal = try_edge_collapse(model, topology, a, b, options; by_length, created))
        if proposal === nothing
            stars_intact && push!(memory.collapses, (a, b))
            continue
        end
        apply!(topology, proposal; metric=model.metric_field)
        record!(created, proposal.new_connectivity)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# Candidate elements: those above the desired density, dilated by the given
# number of layers of node adjacency.
function candidate_elements(
    model::SolidMechanics, topology::MeshTopology, densities::Vector{Float64}, options::AdaptivityOptions
)
    candidates = falses(length(densities))
    if options.shape_by_quality
        # An operation accepted on the cavity minimum can only help the
        # elements below the desired quality, and every edge of such an
        # element is incident to it, so no dilation is needed.
        alive = findall(topology.element_alive)
        sj = quality_jacobians(model, topology.positions, topology.connectivity[:, alive])
        for (k, e) in enumerate(alive)
            sj[k] < options.desired_quality && (candidates[e] = true)
        end
        return findall(candidates)
    end
    for e in eachindex(densities)
        topology.element_alive[e] && densities[e] > options.desired_density && (candidates[e] = true)
    end
    for _ in 1:options.adjacency_layers
        dilated = copy(candidates)
        for e in findall(candidates)
            for node in view(topology.connectivity, :, e), f in node_elements(topology, node)
                dilated[f] = true
            end
        end
        candidates = dilated
    end
    return findall(candidates)
end

# Score by which the operations of a pass are ranked: the energy of the
# elements around an edge or a face, or, under the scaled Jacobian
# criterion, the negative of their minimum scaled Jacobian, so that the
# worst cavities are tried first in either case.
function cavity_rank(
    model::SolidMechanics, topology::MeshTopology, densities::Vector{Float64}, elements, options::AdaptivityOptions
)
    if options.shape_by_quality
        isempty(elements) && return -1.0
        return -minimum(quality_jacobians(model, topology.positions, topology.connectivity[:, collect(elements)]))
    end
    return sum(densities[e] for e in elements; init=0.0)
end

# One pass of swaps over the candidate elements in decreasing order of
# cavity energy: the boundary edges (when enabled), then the interior edges,
# then the faces (when enabled).  An edge or face whose cavity was changed
# earlier in the pass is deferred to the next pass.  Returns the number of
# accepted swaps and the total decrease of the energy.
function swap_pass!(
    model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions; memory::PhaseMemory=PhaseMemory()
)
    densities = energy_densities(model, topology)
    candidates = candidate_elements(model, topology, densities, options)
    edges = Set{Tuple{Int,Int}}()
    boundary_edges = Set{Tuple{Int,Int}}()
    faces = Set{NTuple{3,Int}}()
    for e in candidates
        c = view(topology.connectivity, :, e)
        for i in 1:4, j in (i + 1):4
            edge = sorted_edge(c[i], c[j])
            if is_boundary_edge(topology, c[i], c[j])
                options.boundary_swaps && !(edge in memory.boundary_swaps) && push!(boundary_edges, edge)
            elseif !(edge in memory.swaps)
                push!(edges, edge)
            end
        end
        if options.face_swaps
            for face in element_faces(topology, e)
                length(get(topology.faces, face, Int[])) == 2 && !(face in memory.face_swaps) && push!(faces, face)
            end
        end
    end
    # Rank the edges and faces by the energy of the elements around them, in
    # decreasing order.
    ring_energy(edge) = cavity_rank(model, topology, densities, edge_elements(topology, edge[1], edge[2]), options)
    face_energy(face) = cavity_rank(model, topology, densities, get(topology.faces, face, Int[]), options)
    ranked_boundary = sort!(collect(boundary_edges); by=ring_energy, rev=true)
    ranked = sort!(collect(edges); by=ring_energy, rev=true)
    ranked_faces = sort!(collect(faces); by=face_energy, rev=true)
    accepted = 0
    decrease = 0.0
    created = CreatedEntities()
    # A proposal is refused for the phase when its cavity is intact and the
    # test refuses it; when the cavity was changed earlier in the pass the
    # edge or face is deferred to the next pass instead.
    function accept!(proposal, key, refused, cavity_intact)
        if proposal === nothing
            cavity_intact && push!(refused, key)
            return nothing
        end
        apply!(topology, proposal; metric=model.metric_field)
        record!(created, proposal.new_connectivity)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
        return nothing
    end
    intact(a, b) = all(topology.element_alive[e] for e in edge_elements(topology, a, b))
    for (a, b) in ranked_boundary
        proposal = try_boundary_edge_swap(model, topology, a, b, options; created)
        accept!(proposal, (a, b), memory.boundary_swaps, intact(a, b))
    end
    for (a, b) in ranked
        proposal = try_edge_swap(model, topology, a, b, options; created)
        accept!(proposal, (a, b), memory.swaps, intact(a, b))
    end
    for face in ranked_faces
        proposal = try_face_swap(model, topology, face, options; created)
        accept!(proposal, face, memory.face_swaps, all(topology.element_alive[e] for e in topology.faces[face]))
    end
    return accepted, decrease
end

# The topology phase: passes of operations until none is accepted or the
# cap is reached.  Compacts the topology after every operator, so the
# adjacency is current at the start of the next.  Returns the number of
# accepted operations and the total decrease of the energy.
function topology_phase!(model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions)
    total_accepted = 0
    total_decrease = 0.0
    memory = PhaseMemory()
    previous_swaps = 1
    for pass in 1:options.maximum_passes
        accepted_swaps = accepted_collapses = accepted_splits = 0
        decrease = 0.0
        # The topology is compacted after every operator, so that each one
        # starts from a current adjacency: the operations of one operator
        # are kept independent by the dead elements of their cavities, but
        # the next operator would otherwise not see the elements they made.
        function compact_all!(first_new)
            forget_changed!(memory, topology, first_new)
            node_map, _ = compact!(topology)
            compact_metric!(model.metric_field, findall(>(0), node_map))
            remap!(memory, node_map)
            return nothing
        end
        function swaps!()
            options.swaps || return nothing
            first_new = length(topology.element_alive) + 1
            accepted_swaps, decrease_swaps = swap_pass!(model, topology, options; memory)
            decrease += decrease_swaps
            compact_all!(first_new)
            return nothing
        end
        # The edges outside the length band of a prescribed target are
        # collapsed or split in every pass.  Collapses and splits driven by
        # the shape alone remove or add resolution, so they are tried only in
        # a pass where no swap was accepted (in the pass before, when the size
        # operations come first).
        function size_operations!(shape)
            if options.collapses
                first_new = length(topology.element_alive) + 1
                accepted_collapses, decrease_collapses = collapse_pass!(model, topology, options, shape; memory)
                decrease += decrease_collapses
                compact_all!(first_new)
            end
            if options.splits
                first_new = length(topology.element_alive) + 1
                accepted_splits, decrease_splits = split_pass!(model, topology, options, shape; memory)
                decrease += decrease_splits
                compact_all!(first_new)
            end
            return nothing
        end
        if options.size_first
            size_operations!(options.shape_every_pass || previous_swaps == 0)
            swaps!()
        else
            swaps!()
            size_operations!(options.shape_every_pass || accepted_swaps == 0)
        end
        previous_swaps = accepted_swaps
        accepted = accepted_swaps + accepted_collapses + accepted_splits
        norma_logf(
            0,
            :info,
            "Topology pass %d: %d operations accepted (%d swaps, %d collapses, %d splits), energy decrease %.6e",
            pass,
            accepted,
            accepted_swaps,
            accepted_collapses,
            accepted_splits,
            decrease,
        )
        total_accepted += accepted
        total_decrease += decrease
        accepted == 0 && break
    end
    return total_accepted, total_decrease
end

# Optional observer of the topology phases, called as
# observer(model, topology, accepted, decrease) after each phase, before the
# adapted mesh is written; nothing by default.
const TOPOLOGY_PHASE_OBSERVER = Ref{Any}(nothing)

# The coupled loop: smoothing on the current mesh, a topology phase, a new
# mesh written to disk, and smoothing again on it, until a topology phase
# accepts nothing or the outer iterations are exhausted.  Each smoothed and
# adapted mesh is a separate Exodus file, numbered after the input name.
function run_adaptive(params::Parameters)
    options = AdaptivityOptions(get(params, "adaptivity", Parameters()))
    name = stripped_name(params["output mesh file"])
    sim = create_simulation(params)
    run(sim)
    model = sim.model
    model isa SolidMechanics && model.mesh_smoothing || norma_abort("Adaptivity requires a mesh smoothing model")
    for iteration in 1:options.outer_iterations
        topology = build_topology(model)
        energy_before = sum(energy_densities(model, topology) .* ideal_element_volumes(
            model, 1, topology.connectivity, topology.positions
        ))
        accepted, decrease = topology_phase!(model, topology, options)
        observer = TOPOLOGY_PHASE_OBSERVER[]
        observer === nothing || observer(model, topology, accepted, decrease)
        norma_logf(
            0, :info, "Adaptivity iteration %d: %d operations accepted, energy %.6e -> %.6e",
            iteration, accepted, energy_before, energy_before - decrease,
        )
        accepted == 0 && break
        mesh_file = "$name-adapted-$iteration.g"
        write_topology(topology, mesh_file; nodal_variables=metric_nodal_variables(model.metric_field))
        next_params = deepcopy(params)
        next_params["input mesh file"] = mesh_file
        next_params["output mesh file"] = "$name-adapted-$iteration.e"
        sim = create_simulation(next_params)
        run(sim)
        model = sim.model
    end
    return sim
end
