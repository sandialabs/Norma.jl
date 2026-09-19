# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Topological operations of the adaptivity loop (docs/notes/ems-adaptivity):
# cavity proposals accepted by one test on the smoothing energy.

function AdaptivityOptions(params::Parameters)
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
    )
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
    block = model.blocks[block_index]
    volumes = Vector{Float64}(undef, size(connectivity, 2))
    for e in 1:size(connectivity, 2)
        node_indices = view(connectivity, :, e)
        X, _, _ = smoothing_reference(
            model, block.element_type, Matrix(sample_positions[:, node_indices]); node_indices, metric
        )
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
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
    floor::Float64,
)
    floor ≤ 0.0 && return true
    old_minimum = minimum(scaled_jacobians(topology.positions, topology.connectivity[:, old_elements]))
    new_minimum = minimum(scaled_jacobians(positions, new_connectivity))
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
    options::AdaptivityOptions,
)
    positions = topology.positions
    old_energy = sum(element_energies(model, block, topology.connectivity[:, old_elements], positions))
    new_energies = element_energies(model, block, new_connectivity, positions)
    all(isfinite, new_energies) || return nothing
    new_energy = sum(new_energies)
    new_energy ≤ old_energy - options.minimum_decrease * old_energy || return nothing
    if isfinite(options.allowed_density)
        volumes = ideal_element_volumes(model, block, new_connectivity, positions)
        maximum(new_energies ./ volumes) ≤ options.allowed_density || return nothing
    end
    floor = options.minimum_scaled_jacobian
    passes_scaled_jacobian_floor(topology, old_elements, new_connectivity, positions, floor) || return nothing
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
function try_edge_swap(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int, options::AdaptivityOptions)
    ring = edge_ring(topology, a, b)
    ring === nothing && return nothing
    elements, ring_nodes = ring
    n = length(ring_nodes)
    block = topology.block[elements[1]]
    all(topology.block[e] == block for e in elements) || return nothing
    best = nothing
    best_energy = Inf
    for triangles in polygon_triangulations(n)
        connectivity = swapped_connectivity(topology, a, b, ring_nodes, triangles)
        connectivity === nothing && continue
        energies = element_energies(model, block, connectivity, topology.positions)
        all(isfinite, energies) || continue
        energy = sum(energies)
        if energy < best_energy
            best_energy = energy
            best = connectivity
        end
    end
    best === nothing && return nothing
    return accept_proposal(model, topology, elements, best, block, options)
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
function try_edge_collapse(model::SolidMechanics, topology::MeshTopology, b::Int, a::Int, options::AdaptivityOptions)
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
    return accept_proposal(model, topology, star, new_connectivity, block, options, b, a)
end

function accept_proposal(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    options::AdaptivityOptions,
    removed_node::Int,
    surviving_node::Int,
)
    proposal = accept_proposal(model, topology, old_elements, new_connectivity, block, options)
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
function try_edge_split(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int, options::AdaptivityOptions)
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
    energy_after ≤ energy_before - options.minimum_decrease * energy_before || return nothing
    positions = PositionsWithNode(topology.positions, position)
    if isfinite(options.allowed_density)
        energies = element_energies(model, block, connectivity, positions; metric)
        volumes = ideal_element_volumes(model, block, connectivity, positions; metric)
        maximum(energies ./ volumes) ≤ options.allowed_density || return nothing
    end
    floor = options.minimum_scaled_jacobian
    passes_scaled_jacobian_floor(topology, collect(ring), connectivity, positions, floor) || return nothing
    split = SplitNode(position, node_sets, side_sets, (a, b), node_metric)
    return CavityProposal(collect(ring), connectivity, block, energy_before, energy_after, 0, 0, split)
end

# Edges longer than sqrt(2) in the prescribed target, the split candidates of
# the size phase; empty without a target.
function long_edges(model::SolidMechanics, topology::MeshTopology)
    edges = Tuple{Int,Int}[]
    (model.metric_field === nothing && model.size_field === nothing) && return edges
    for edge in keys(topology.edges)
        metric_edge_length(model, topology, edge[1], edge[2]) > sqrt(2.0) && push!(edges, edge)
    end
    return edges
end

# One pass of edge splits over the candidate edges, in decreasing order of
# the energy around them: the long edges of the target always, and the edges
# of the elements above the desired density when the shape-driven operations
# are on.
function split_pass!(model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions, shape::Bool)
    densities = energy_densities(model, topology)
    edges = Set{Tuple{Int,Int}}()
    if shape
        for e in candidate_elements(topology, densities, options)
            c = view(topology.connectivity, :, e)
            for i in 1:4, j in (i + 1):4
                push!(edges, sorted_edge(c[i], c[j]))
            end
        end
    end
    union!(edges, long_edges(model, topology))
    ranked = collect(edges)
    edge_energy(edge) = sum(densities[e] for e in edge_elements(topology, edge[1], edge[2]); init=0.0)
    sort!(ranked; by=edge_energy, rev=true)
    accepted = 0
    decrease = 0.0
    for (a, b) in ranked
        proposal = try_edge_split(model, topology, a, b, options)
        proposal === nothing && continue
        apply!(topology, proposal; metric=model.metric_field)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# Metric factor of the prescribed target on an edge: the scalar 1/h at the
# midpoint for a size field, F_M for a metric field (sampled at the midpoint
# by the function sources, averaged over the two nodes by the nodal sources),
# nothing for the legacy rules, which prescribe no length.
function metric_factor_on_edge(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int)
    xa = SVector{3,Float64}(view(topology.positions, :, a))
    xb = SVector{3,Float64}(view(topology.positions, :, b))
    midpoint = 0.5 * (xa + xb)
    if model.metric_field !== nothing
        h, R = principal_metric(model.metric_field.source, (a, b), midpoint, model.time)
        return SMatrix{3,3,Float64,9}(Diagonal(SVector{3,Float64}(1.0 / h[1], 1.0 / h[2], 1.0 / h[3]))) * R'
    elseif model.size_field !== nothing
        return SMatrix{3,3,Float64,9}(I) / model.size_field((model.time, midpoint[1], midpoint[2], midpoint[3]))
    end
    return nothing
end

# Length of an edge measured in the prescribed target, or NaN without one.
function metric_edge_length(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int)
    F_M = metric_factor_on_edge(model, topology, a, b)
    F_M === nothing && return NaN
    xa = SVector{3,Float64}(view(topology.positions, :, a))
    xb = SVector{3,Float64}(view(topology.positions, :, b))
    return norm(F_M * (xb - xa))
end

# Edges shorter than 1/sqrt(2) in the prescribed target, the collapse
# candidates of the size phase; empty without a target.
function short_edges(model::SolidMechanics, topology::MeshTopology)
    edges = Tuple{Int,Int}[]
    (model.metric_field === nothing && model.size_field === nothing) && return edges
    for edge in keys(topology.edges)
        metric_edge_length(model, topology, edge[1], edge[2]) < 1.0 / sqrt(2.0) && push!(edges, edge)
    end
    return edges
end

# One pass of edge collapses over the candidate edges, in decreasing order
# of the energy around them: the short edges of the target always, and the
# edges of the elements above the desired density when the shape-driven
# operations are on.  Both directions of each edge are tried.  Returns the
# number of accepted collapses and the decrease.
function collapse_pass!(model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions, shape::Bool)
    densities = energy_densities(model, topology)
    edges = Set{Tuple{Int,Int}}()
    if shape
        for e in candidate_elements(topology, densities, options)
            c = view(topology.connectivity, :, e)
            for i in 1:4, j in (i + 1):4
                push!(edges, sorted_edge(c[i], c[j]))
            end
        end
    end
    union!(edges, short_edges(model, topology))
    ranked = collect(edges)
    edge_energy(edge) = sum(densities[e] for e in edge_elements(topology, edge[1], edge[2]); init=0.0)
    sort!(ranked; by=edge_energy, rev=true)
    accepted = 0
    decrease = 0.0
    for (a, b) in ranked
        (topology.node_alive[a] && topology.node_alive[b]) || continue
        proposal = try_edge_collapse(model, topology, b, a, options)
        proposal === nothing && (proposal = try_edge_collapse(model, topology, a, b, options))
        proposal === nothing && continue
        apply!(topology, proposal; metric=model.metric_field)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# Candidate elements: those above the desired density, dilated by the given
# number of layers of node adjacency.
function candidate_elements(topology::MeshTopology, densities::Vector{Float64}, options::AdaptivityOptions)
    candidates = falses(length(densities))
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

# One pass of edge swaps over the candidate edges in decreasing order of
# cavity energy.  An edge whose ring was changed earlier in the pass is
# deferred to the next pass.  Returns the number of accepted swaps and the
# total decrease of the energy.
function swap_pass!(model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions)
    densities = energy_densities(model, topology)
    candidates = candidate_elements(topology, densities, options)
    edges = Set{Tuple{Int,Int}}()
    for e in candidates
        c = view(topology.connectivity, :, e)
        for i in 1:4, j in (i + 1):4
            edge = sorted_edge(c[i], c[j])
            is_boundary_edge(topology, c[i], c[j]) || push!(edges, edge)
        end
    end
    # Rank the edges by the energy of their rings, in decreasing order.
    ranked = collect(edges)
    ring_energy(edge) = sum(densities[e] for e in edge_elements(topology, edge[1], edge[2]); init=0.0)
    sort!(ranked; by=ring_energy, rev=true)
    accepted = 0
    decrease = 0.0
    for (a, b) in ranked
        proposal = try_edge_swap(model, topology, a, b, options)
        proposal === nothing && continue
        apply!(topology, proposal; metric=model.metric_field)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# The topology phase: passes of operations until none is accepted or the
# cap is reached.  Compacts the topology after every pass, so the adjacency
# is current at the start of the next.  Returns the number of accepted
# operations and the total decrease of the energy.
function topology_phase!(model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions)
    total_accepted = 0
    total_decrease = 0.0
    for pass in 1:options.maximum_passes
        accepted = 0
        decrease = 0.0
        if options.swaps
            accepted_swaps, decrease_swaps = swap_pass!(model, topology, options)
            accepted += accepted_swaps
            decrease += decrease_swaps
        end
        # The edges outside the length band of a prescribed target are
        # collapsed or split in every pass.  Collapses and splits driven by
        # the shape alone remove or add resolution, so they are tried only in
        # a pass where no swap was accepted.
        shape = accepted == 0
        if options.collapses
            accepted_collapses, decrease_collapses = collapse_pass!(model, topology, options, shape)
            accepted += accepted_collapses
            decrease += decrease_collapses
        end
        if options.splits
            accepted_splits, decrease_splits = split_pass!(model, topology, options, shape)
            accepted += accepted_splits
            decrease += decrease_splits
        end
        node_map, _ = compact!(topology)
        compact_metric!(model.metric_field, findall(>(0), node_map))
        norma_logf(0, :info, "Topology pass %d: %d operations accepted, energy decrease %.6e", pass, accepted, decrease)
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
