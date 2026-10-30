# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Topological operations of the adaptivity loop (docs/notes/ems-adaptivity):
# cavity proposals accepted by one test on the smoothing energy, or, for the
# size operations under `size criterion: length`, on the edge length in the
# prescribed target with the energy as a validity check only, and, for the
# shape operations under `shape criterion: scaled Jacobian`, on the minimum
# scaled Jacobian of the cavity.  Under `cavity relaxation` a swap or
# collapse refused at fixed positions is judged again after the interior
# nodes of its cavity are relaxed against the new connectivity
# (relaxed_proposal).

function AdaptivityOptions(params::Parameters)
    size_criterion = get(params, "size criterion", "energy")
    if !(size_criterion in ("energy", "length"))
        norma_abort("\"size criterion\" in \"adaptivity\" must be \"energy\" or \"length\"; got \"$size_criterion\"")
    end
    shape_criterion = get(params, "shape criterion", "energy")
    if !(shape_criterion in ("energy", "scaled Jacobian"))
        norma_abort(
            "\"shape criterion\" in \"adaptivity\" must be \"energy\" or \"scaled Jacobian\"; got \"$shape_criterion\""
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
        Bool(get(params, "cavity relaxation", false)),
    )
end

# Forget the refusals of every edge and face of the elements added since
# `first_new`, whose cavities have changed, and drop the dead ones.
function forget_changed!(memory::PhaseMemory, topology::MeshTopology, first_new::Int)
    return forget_elements!(memory, topology, first_new:length(topology.element_alive))
end

# Forget the refusals of every edge and face of the given elements.
function forget_elements!(memory::PhaseMemory, topology::MeshTopology, elements)
    for e in elements
        for edge in element_edges(topology, e)
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
        X = tetrahedron_coordinates(positions, node_indices)
        _, F_M, _ = create_metric_reference(model.metric_field, X, model.time; node_indices)
        quality[e] = tetrahedron_scaled_jacobian(F_M * X)
    end
    return quality
end

# The minimum scaled Jacobian of the old elements of a cavity and that of its
# new elements.
function cavity_quality_minima(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
)
    old_minimum = minimum(quality_jacobians(model, topology.positions, topology.connectivity[:, old_elements]))
    new_minimum = minimum(quality_jacobians(model, positions, new_connectivity))
    return old_minimum, new_minimum
end

# Whether the new elements of a cavity raise its minimum scaled Jacobian.
function raises_minimum_quality(old_minimum::Float64, new_minimum::Float64)
    return new_minimum > old_minimum * (1.0 + MINIMUM_QUALITY_INCREASE)
end

function raises_minimum_quality(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
)
    return raises_minimum_quality(cavity_quality_minima(model, topology, old_elements, new_connectivity, positions)...)
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
        sample = tetrahedron_coordinates(sample_positions, node_indices)
        X, _, _ = smoothing_reference(model, element_type, sample; node_indices, metric)
        volumes[e] = abs(tetrahedron_volume(X))
    end
    return volumes
end

# Energy densities (energy per unit ideal volume) of the alive elements of a
# topology at its current positions; dead elements get NaN.
function energy_densities(model::SolidMechanics, topology::MeshTopology)
    densities = fill(NaN, length(topology.element_alive))
    energy_densities!(densities, model, topology, findall(topology.element_alive))
    return densities
end

# Energy densities of the given elements, written into `densities`.
function energy_densities!(
    densities::Vector{Float64}, model::SolidMechanics, topology::MeshTopology, selected::AbstractVector{Int}
)
    for block_index in unique(topology.block[selected])
        elements = [e for e in selected if topology.block[e] == block_index]
        connectivity = topology.connectivity[:, elements]
        energies = element_energies(model, block_index, connectivity, topology.positions)
        volumes = ideal_element_volumes(model, block_index, connectivity, topology.positions)
        densities[elements] = energies ./ volumes
    end
    return densities
end

# Scaled Jacobians of the given elements, written into `quality`.
function element_qualities!(
    quality::Vector{Float64}, model::SolidMechanics, topology::MeshTopology, selected::AbstractVector{Int}
)
    quality[selected] = quality_jacobians(model, topology.positions, topology.connectivity[:, selected])
    return quality
end

# The energy densities and, under the scaled Jacobian criterion, the scaled
# Jacobians of the elements of the topology, from the memory of the phase
# when it matches the topology, and evaluated and stored in it otherwise.
function phase_element_values!(
    memory::PhaseMemory, model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions
)
    n = length(topology.element_alive)
    if length(memory.density) != n || (options.shape_by_quality && length(memory.quality) != n)
        alive = findall(topology.element_alive)
        resize!(memory.density, n)
        fill!(memory.density, NaN)
        energy_densities!(memory.density, model, topology, alive)
        if options.shape_by_quality
            resize!(memory.quality, n)
            fill!(memory.quality, NaN)
            element_qualities!(memory.quality, model, topology, alive)
        end
    end
    return memory.density, memory.quality
end

# Renumber the element values of the memory after a compaction and evaluate
# those of the elements added since `first_new`; empty them when they did not
# match the topology before it.
function remap_element_values!(
    memory::PhaseMemory,
    model::SolidMechanics,
    topology::MeshTopology,
    element_map::Vector{Int},
    first_new::Int,
    options::AdaptivityOptions,
)
    for values in (memory.density, memory.quality)
        if length(values) != first_new - 1
            empty!(values)
            continue
        end
        remapped = fill(NaN, length(topology.element_alive))
        for e in 1:(first_new - 1)
            element_map[e] > 0 && (remapped[element_map[e]] = values[e])
        end
        resize!(values, length(remapped))
        copyto!(values, remapped)
    end
    added = [element_map[e] for e in first_new:length(element_map) if element_map[e] > 0]
    isempty(memory.density) || energy_densities!(memory.density, model, topology, added)
    options.shape_by_quality && !isempty(memory.quality) && element_qualities!(memory.quality, model, topology, added)
    return memory
end

# The geometric floor: an operation may not create an element whose scaled
# Jacobian is below the floor, unless the worst element of the old cavity was
# already below it and the new worst is no worse.  The energy sums the cavity,
# so this is what keeps a single element from being sacrificed for the sum.
function passes_scaled_jacobian_floor(old_minimum::Float64, new_minimum::Float64, floor::Float64)
    new_minimum ≥ floor && return true
    return new_minimum ≥ old_minimum
end

# The tests of a proposal with finite energies: the decrease of the energy by
# the relative margin, or under the scaled Jacobian criterion the rise of the
# cavity minimum, when a decrease is required; the allowed density of the new
# elements; and the geometric floor.  The energies of the new elements are
# evaluated here when the density test needs them and none are given.  The
# minima of the scaled Jacobian are evaluated once, when a test first needs
# them.
function passes_acceptance_tests(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    positions::AbstractMatrix{Float64},
    metric::Union{MetricField,Nothing},
    new_energies::Union{Vector{Float64},Nothing},
    energy_before::Float64,
    energy_after::Float64,
    options::AdaptivityOptions;
    require_decrease::Bool=true,
)
    minima = nothing
    if require_decrease
        if options.shape_by_quality
            minima = cavity_quality_minima(model, topology, old_elements, new_connectivity, positions)
            raises_minimum_quality(minima...) || return false
        else
            energy_after ≤ energy_before - options.minimum_decrease * energy_before || return false
        end
    end
    if isfinite(options.allowed_density)
        if new_energies === nothing
            new_energies = element_energies(model, block, new_connectivity, positions; metric)
        end
        volumes = ideal_element_volumes(model, block, new_connectivity, positions; metric)
        maximum(new_energies ./ volumes) ≤ options.allowed_density || return false
    end
    floor = options.minimum_scaled_jacobian
    floor ≤ 0.0 && return true
    if minima === nothing
        minima = cavity_quality_minima(model, topology, old_elements, new_connectivity, positions)
    end
    return passes_scaled_jacobian_floor(minima..., floor)
end

# The acceptance test: a proposal is accepted when the energy of the new
# elements is below that of the old ones by the relative margin, no new
# element exceeds the allowed density, and the geometric floor holds.
# Returns the proposal, with the node a collapse removes and the node it
# moves onto, or the change of a side set by a boundary swap, when given; or
# nothing.  Under `cavity relaxation` a proposal refused on the energy at
# fixed positions is judged again after the relaxation of the interior nodes
# of its cavity (relaxed_proposal); the test on the length and the test on
# the scaled Jacobian are not relaxed, since the relaxation descends on the
# energy.
function accept_proposal(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    options::AdaptivityOptions;
    require_decrease::Bool=true,
    removed_node::Int=0,
    surviving_node::Int=0,
    surface::Union{SurfaceSwap,Nothing}=nothing,
)
    positions = topology.positions
    old_energy = sum(element_energies(model, block, topology.connectivity[:, old_elements], positions))
    new_energies = element_energies(model, block, new_connectivity, positions)
    all(isfinite, new_energies) || return nothing
    new_energy = sum(new_energies)
    if passes_acceptance_tests(
        model,
        topology,
        old_elements,
        new_connectivity,
        block,
        positions,
        model.metric_field,
        new_energies,
        old_energy,
        new_energy,
        options;
        require_decrease,
    )
        return CavityProposal(
            old_elements,
            new_connectivity,
            block,
            old_energy,
            new_energy,
            removed_node,
            surviving_node,
            nothing,
            surface,
            Int[],
            zeros(3, 0),
        )
    end
    (options.cavity_relaxation && require_decrease && !options.shape_by_quality) || return nothing
    return relaxed_proposal(
        model,
        topology,
        old_elements,
        new_connectivity,
        block,
        old_energy,
        new_energy,
        options,
        removed_node,
        surviving_node,
        surface,
    )
end

# Descent steps of the relaxation of each node of a cavity, and the ratio
# of the energy of the new and surrounding elements after the first step to
# that of the old and surrounding elements above which the relaxation stops:
# on the awful cube loop 98.2% of the operations the full relaxation accepts
# are below 1.5 after the first step, and 71% of the relaxations tried are
# above it (study of October 30, 2026, relaxation.md).
const CAVITY_RELAXATION_STEPS = 5
const CAVITY_RELAXATION_CONTINUE_RATIO = 1.5

# Diagnostic of the studies: when `RECORD_RELAXATION_RATIOS[]` is set, every
# proposal the relaxation judges appends to `RELAXATION_RATIOS` the ratio of
# the new to the old cavity energy at the fixed positions, the ratio of the
# energy of the new and surrounding elements after the first step over the
# nodes to that of the old and surrounding elements, and whether it was
# accepted.
const RECORD_RELAXATION_RATIOS = Ref(false)
const RELAXATION_RATIOS = Tuple{Float64,Float64,Bool}[]
const RELAXATION_RATIOS_LOCK = ReentrantLock()

# The nodes of a cavity that its relaxation may move: the nodes of its old
# and new elements, other than the node a collapse removes, that are in no
# node set and on no side set (so that the boundary and its features are
# kept), that the adjacency knows, and whose every element outside the
# cavity is alive and in the block of the cavity.  An element around the
# node that died in the pass means that an operation of the pass replaced
# it by elements the adjacency does not know, which the relaxation could
# not evaluate.
function movable_cavity_nodes(
    topology::MeshTopology, old_elements::Vector{Int}, new_connectivity::Matrix{Int}, block::Int, removed_node::Int
)
    nodes = Int[]
    in_set(sets, n) = any(n ≤ length(flags) && flags[n] for flags in values(sets))
    for n in unique(vcat(vec(topology.connectivity[:, old_elements]), vec(new_connectivity)))
        n == removed_node && continue
        n < length(topology.node_to_elements_offsets) || continue
        topology.node_alive[n] || continue
        (in_set(topology.node_sets, n) || in_set(topology.node_side_sets, n)) && continue
        movable = true
        for e in node_elements(topology, n)
            e in old_elements && continue
            if !topology.element_alive[e] || topology.block[e] != block
                movable = false
                break
            end
        end
        movable && push!(nodes, n)
    end
    return nodes
end

# The elements outside the cavity around the given nodes, without repetition.
function outer_elements(topology::MeshTopology, old_elements::Vector{Int}, nodes::Vector{Int})
    outer = Int[]
    for n in nodes, e in node_elements(topology, n)
        (e in old_elements || e in outer) || push!(outer, e)
    end
    return outer
end

# Relax the given nodes against the new elements of a cavity and the
# elements around the nodes outside it, one node at a time with the others
# fixed, by descent steps on the energy of the elements around the node, as
# relax_split_node: the gradient from element_energy_gradients, the step
# scaled by the mean length of the edges from the node, halved until the
# energy decreases, and the target sampled at the positions before the
# relaxation (see SPLIT_RELAXATION_STEPS).  After the first step over the nodes the energy of the
# elements is compared with `energy_before`, that of the old and surrounding
# elements, and the relaxation stops when the ratio is above
# CAVITY_RELAXATION_CONTINUE_RATIO.  Returns the positions with the nodes
# moved, or nothing when no node moved or the relaxation stopped, and the
# energy after the first step (NaN when no step was taken).
function relax_cavity_nodes(
    model::SolidMechanics,
    topology::MeshTopology,
    block::Int,
    new_connectivity::Matrix{Int},
    outer_connectivity::Matrix{Int},
    nodes::Vector{Int},
    energy_before::Float64,
)
    positions = topology.positions
    trial = MovedPositions(positions, nodes, Matrix{Float64}(positions[:, nodes]))
    around = hcat(new_connectivity, outer_connectivity)
    stars = [around[:, [n in column for column in eachcol(around)]] for n in nodes]
    scales = Float64[]
    for (k, n) in enumerate(nodes)
        x = SVector{3,Float64}(view(positions, :, n))
        neighbors = unique(m for m in vec(stars[k]) if m != n)
        push!(scales, sum(norm(SVector{3,Float64}(view(positions, :, m)) - x) for m in neighbors) / length(neighbors))
    end
    moved = false
    first_step_energy = NaN
    for step_index in 1:CAVITY_RELAXATION_STEPS
        improved = false
        for (k, n) in enumerate(nodes)
            star = stars[k]
            function energy(x)
                trial.moved[:, k] = x
                return sum(element_energies(model, block, star, trial; sample_positions=positions))
            end
            x = SVector{3,Float64}(view(trial.moved, :, k))
            energies, gradients = element_energy_gradients(model, block, star, trial; sample_positions=positions)
            value = sum(energies)
            isfinite(value) || return nothing
            step = node_descent_step(energy, identity, x, value, node_gradient(gradients, star, n), scales[k])
            if step === nothing
                trial.moved[:, k] = x
            else
                trial.moved[:, k] = step[1]
                improved = moved = true
            end
        end
        improved || break
        if step_index == 1
            first_step_energy = sum(element_energies(model, block, around, trial; sample_positions=positions))
            first_step_energy > CAVITY_RELAXATION_CONTINUE_RATIO * energy_before && return nothing, first_step_energy
        end
    end
    return moved ? trial : nothing, first_step_energy
end

# The acceptance test of a proposal after the relaxation of its cavity: the
# energy and the floor compare the old elements with the new ones at the
# relaxed positions, both together with the elements around the moved nodes
# that the operation keeps, since those change with the positions.  Returns
# the proposal with the moved nodes and their positions, or nothing.
function relaxed_proposal(
    model::SolidMechanics,
    topology::MeshTopology,
    old_elements::Vector{Int},
    new_connectivity::Matrix{Int},
    block::Int,
    old_energy::Float64,
    new_energy::Float64,
    options::AdaptivityOptions,
    removed_node::Int,
    surviving_node::Int,
    surface::Union{SurfaceSwap,Nothing},
)
    nodes = movable_cavity_nodes(topology, old_elements, new_connectivity, block, removed_node)
    isempty(nodes) && return nothing
    outer = outer_elements(topology, old_elements, nodes)
    outer_connectivity = topology.connectivity[:, outer]
    energy_before = old_energy + sum(element_energies(model, block, outer_connectivity, topology.positions))
    positions, first_step_energy = relax_cavity_nodes(
        model, topology, block, new_connectivity, outer_connectivity, nodes, energy_before
    )
    accepted = false
    energy_after = Inf
    extended = hcat(new_connectivity, outer_connectivity)
    energies = Float64[]
    if positions !== nothing
        energies = element_energies(model, block, extended, positions)
        if all(isfinite, energies)
            energy_after = sum(energies)
            accepted = passes_acceptance_tests(
                model,
                topology,
                vcat(old_elements, outer),
                extended,
                block,
                positions,
                model.metric_field,
                energies,
                energy_before,
                energy_after,
                options,
            )
        end
    end
    if RECORD_RELAXATION_RATIOS[]
        lock(RELAXATION_RATIOS_LOCK) do
            return push!(RELAXATION_RATIOS, (new_energy / old_energy, first_step_energy / energy_before, accepted))
        end
    end
    accepted || return nothing
    return CavityProposal(
        old_elements,
        new_connectivity,
        block,
        energy_before,
        energy_after,
        removed_node,
        surviving_node,
        nothing,
        surface,
        nodes,
        positions.moved,
    )
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
    for (k, n) in enumerate(proposal.moved_nodes)
        topology.positions[:, n] = proposal.moved_positions[:, k]
    end
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

# The block of the elements of a cavity; nothing when one of them is dead or
# when they do not all belong to one block.
function cavity_block(topology::MeshTopology, elements::AbstractVector{Int})
    all(topology.element_alive[e] for e in elements) || return nothing
    block = topology.block[elements[1]]
    all(topology.block[e] == block for e in elements) || return nothing
    return block
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
    others = opposite_node_pairs(topology, elements, a, b)
    others === nothing && return nothing
    walk = walk_edge_elements(elements, others, 1, others[1][1])
    walk === nothing && return nothing
    ordered, ring = walk
    ring[end] == ring[1] || return nothing
    pop!(ring)
    length(unique(ring)) == n || return nothing
    return ordered, ring
end

# The two nodes other than a and b of each element around the edge (a, b);
# nothing when an element does not have exactly two such nodes.
function opposite_node_pairs(topology::MeshTopology, elements::AbstractVector{Int}, a::Int, b::Int)
    pairs = Vector{Tuple{Int,Int}}(undef, length(elements))
    for (i, e) in enumerate(elements)
        pq = Int[]
        for node in view(topology.connectivity, :, e)
            (node == a || node == b) || push!(pq, node)
        end
        length(pq) == 2 || return nothing
        pairs[i] = (pq[1], pq[2])
    end
    return pairs
end

# The elements around an edge in the order of a walk that starts at element
# `start`, entered through its node `first_node`, and passes each time to an
# unvisited element that contains the last node reached; and the nodes
# reached, one more than the elements.  Returns nothing when the walk ends
# before it has visited every element.
function walk_edge_elements(elements::AbstractVector{Int}, pairs::Vector{Tuple{Int,Int}}, start::Int, first_node::Int)
    n = length(elements)
    used = falses(n)
    used[start] = true
    ordered = Int[elements[start]]
    nodes = Int[first_node, pairs[start][1] == first_node ? pairs[start][2] : pairs[start][1]]
    for _ in 2:n
        last = nodes[end]
        found = 0
        for i in 1:n
            used[i] && continue
            if pairs[i][1] == last || pairs[i][2] == last
                found = i
                break
            end
        end
        found == 0 && return nothing
        used[found] = true
        push!(ordered, elements[found])
        push!(nodes, pairs[found][1] == last ? pairs[found][2] : pairs[found][1])
    end
    return ordered, nodes
end

# Orientation of the ring polygon of the edge (a, b) in the projection along
# the edge, from its first sector, which the old element guarantees to be
# well formed: +1 or -1, or 0 for a degenerate sector.
function ring_orientation(topology::MeshTopology, a::Int, b::Int, ring::Vector{Int})
    x(n) = SVector{3,Float64}(view(topology.positions, :, n))
    xa = x(a)
    return sign(dot(cross(x(ring[1]) - xa, x(ring[2]) - xa), x(b) - xa))
end

# Connectivity of the two tetrahedra that replace the ring of edge (a, b)
# under the triangle (i, j, k) of the ring polygon, one with each edge node
# as apex, for the ring orientation s.  The triangle is admissible only if,
# taken in the cyclic order of the ring, it has the orientation of the
# polygon in the projection along the edge (an ear at a reflex vertex lies
# outside the cavity even though its tetrahedra can be positively oriented),
# and if the edge nodes lie on opposite sides of it.  Returns nothing when
# either fails.
function triangle_connectivity(
    topology::MeshTopology, a::Int, b::Int, ring::Vector{Int}, s::Float64, triangle::NTuple{3,Int}
)
    x(n) = SVector{3,Float64}(view(topology.positions, :, n))
    xa, xb = x(a), x(b)
    i, j, k = triangle
    p, q, r = ring[i], ring[j], ring[k]
    xp, xq, xr = x(p), x(q), x(r)
    normal = cross(xq - xp, xr - xp)
    sign(dot(normal, xb - xa)) == s || return nothing
    # With the triangle oriented toward b, (p, q, r, b) is positive and
    # (p, r, q, a) is positive when a lies on the other side.
    if s < 0.0
        q, r = r, q
        normal = -normal
    end
    dot(normal, xb - xp) > 0.0 || return nothing
    dot(normal, xa - xp) < 0.0 || return nothing
    return [p p; q r; r q; b a]
end

# Connectivity of the elements that replace the ring of edge (a, b) for one
# triangulation of the ring polygon: the two tetrahedra of each triangle, in
# the order of the triangles.  Returns nothing when a triangle is not
# admissible (triangle_connectivity).
function swapped_connectivity(
    topology::MeshTopology, a::Int, b::Int, ring::Vector{Int}, triangles::Vector{NTuple{3,Int}}
)
    s = ring_orientation(topology, a, b, ring)
    s == 0.0 && return nothing
    connectivity = zeros(Int, 4, 2 * length(triangles))
    for (t, triangle) in enumerate(triangles)
        columns = triangle_connectivity(topology, a, b, ring, s, triangle)
        columns === nothing && return nothing
        connectivity[:, (2t - 1):(2t)] .= columns
    end
    return connectivity
end

# Whether the edge (u, v) exists, in the mesh or from an operation of this
# pass.
function edge_exists(topology::MeshTopology, u::Int, v::Int, created::CreatedEntities)
    edge = sorted_edge(u, v)
    return haskey(topology.edges, edge) || edge in created.edges
end

# Whether a chord of the ring polygon of a triangulation is an edge already,
# in the mesh or from an operation of this pass.
function chords_exist(
    topology::MeshTopology, ring_nodes::Vector{Int}, triangles::Vector{NTuple{3,Int}}, created::CreatedEntities
)
    n = length(ring_nodes)
    for (i, j, k) in triangles, (u, v) in ((i, j), (j, k), (k, i))
        mod(u - v, n) in (1, n - 1) && continue
        edge_exists(topology, ring_nodes[u], ring_nodes[v], created) && return true
    end
    return false
end

# The ring of the interior edge (a, b) and the triangulation of least energy
# among those whose chords and faces are new: a named tuple with the ring
# elements and nodes, their block, and the best triangles and connectivity
# (nothing when no triangulation is admissible); nothing when the edge has no
# ring of elements of one block.
function best_edge_swap(
    model::SolidMechanics, topology::MeshTopology, a::Int, b::Int, options::AdaptivityOptions, created::CreatedEntities
)
    ring = edge_ring(topology, a, b)
    ring === nothing && return nothing
    elements, ring_nodes = ring
    block = cavity_block(topology, elements)
    block === nothing && return nothing
    best, best_triangles = best_triangulation(model, topology, a, b, elements, ring_nodes, block, options, created)
    return (; elements, ring_nodes, block, triangles=best_triangles, connectivity=best)
end

# The triangulation of least score of the polygon of the nodes p_1, ..., p_n
# around the edge (a, b), a closed ring or the open chain of a boundary edge,
# among those whose chords and faces are new, with its connectivity; nothing
# for both when no triangulation is admissible.
#
# The polygon has C(n - 2) triangulations (C the Catalan numbers, 3.5e7 for
# n = 18), so they are not enumerated (issue #247).  The score of a
# triangulation is the sum of the scores of its triangles under the energy
# criterion and their maximum under the scaled Jacobian criterion
# (configuration_score), and every admissibility condition holds or fails
# for one triangle or one chord at a time.  The least score is therefore
# found by dynamic programming over the sub-polygons p_i, ..., p_j (Klincsek
# 1980): the least score of a sub-polygon is the least, over its apexes k,
# of the combination of the least scores of p_i, ..., p_k and p_k, ..., p_j
# with that of the triangle (i, k, j), which takes C(n, 3) evaluations of a
# triangle and storage of order n³.
#
# Under the energy criterion each sub-polygon keeps its first apex of least
# score.  Under the scaled Jacobian criterion many triangulations share the
# least score, which is that of their worst triangle; the one returned is
# the first of them in the order in which the triangulations were
# enumerated before issue #247 (apex of p_1, ..., p_n ascending, then the
# triangulation of the left part, then that of the right part).  A
# sub-polygon has a triangulation of score at most v exactly when its least
# score is at most v, so that triangulation is found by descending from the
# whole polygon with the bound v the least score, taking at each
# sub-polygon the first apex whose triangle and two parts meet the bound.
# Under either rule, when the admissible triangulations are reduced and the
# returned one remains admissible, it is returned again
# (edge_swap_from_evaluation).
function best_triangulation(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    elements::Vector{Int},
    nodes::Vector{Int},
    block::Int,
    options::AdaptivityOptions,
    created::CreatedEntities,
)
    n = length(nodes)
    s = ring_orientation(topology, a, b, nodes)
    s == 0.0 && return nothing, nothing
    combine = options.shape_by_quality ? max : (+)
    # Least score of the sub-polygon p_i, ..., p_j and its apex; +Inf when
    # no triangulation of it is admissible.  A side of the polygon has the
    # neutral score of the combination.
    score = fill(Inf, n, n)
    apex = zeros(Int, n, n)
    # Scores of the triangles, kept under the scaled Jacobian criterion for
    # the descent; +Inf for a triangle not evaluated.
    triangle_scores = options.shape_by_quality ? fill(Inf, n, n, n) : fill(Inf, 0, 0, 0)
    for i in 1:(n - 1)
        score[i, i + 1] = options.shape_by_quality ? -Inf : 0.0
    end
    for width in 2:(n - 1), i in 1:(n - width)
        j = i + width
        # A chord of the polygon becomes an edge; it may not exist already,
        # since the link of an edge must be a single ring.  The side (p_1,
        # p_n) of the polygon is not a chord.
        (i, j) != (1, n) && edge_exists(topology, nodes[i], nodes[j], created) && continue
        for k in (i + 1):(j - 1)
            (score[i, k] == Inf || score[k, j] == Inf) && continue
            triangle = triangle_score(model, topology, a, b, elements, nodes, s, (i, k, j), block, options, created)
            options.shape_by_quality && (triangle_scores[i, k, j] = triangle)
            value = combine(combine(score[i, k], triangle), score[k, j])
            if value < score[i, j]
                score[i, j] = value
                apex[i, j] = k
            end
        end
    end
    score[1, n] == Inf && return nothing, nothing
    triangles = NTuple{3,Int}[]
    if options.shape_by_quality
        first_triangles!(triangles, score, triangle_scores, 1, n, score[1, n])
    else
        collect_triangles!(triangles, apex, 1, n)
    end
    return swapped_connectivity(topology, a, b, nodes, triangles), triangles
end

# Score of the triangle (i, j, k) of the polygon of best_triangulation: the
# score of its two tetrahedra (configuration_score), or +Inf when the
# triangle is not admissible or a face of its tetrahedra exists already.
function triangle_score(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    elements::Vector{Int},
    nodes::Vector{Int},
    s::Float64,
    triangle::NTuple{3,Int},
    block::Int,
    options::AdaptivityOptions,
    created::CreatedEntities,
)
    connectivity = triangle_connectivity(topology, a, b, nodes, s, triangle)
    connectivity === nothing && return Inf
    faces_are_new(topology, elements, connectivity, created) || return Inf
    value = configuration_score(model, topology, block, connectivity, options)
    return isnan(value) ? Inf : value
end

# The triangles of the sub-polygon p_i, ..., p_j from the apexes of
# best_triangulation, in the order left part, apex triangle, right part.
function collect_triangles!(triangles::Vector{NTuple{3,Int}}, apex::Matrix{Int}, i::Int, j::Int)
    j - i < 2 && return triangles
    k = apex[i, j]
    collect_triangles!(triangles, apex, i, k)
    push!(triangles, (i, k, j))
    collect_triangles!(triangles, apex, k, j)
    return triangles
end

# The triangles of the first triangulation of the sub-polygon p_i, ..., p_j
# whose triangles all have scores at most `bound`, in the order of
# best_triangulation, from its least scores of the sub-polygons and scores
# of the triangles; the sub-polygon must have such a triangulation.
function first_triangles!(
    triangles::Vector{NTuple{3,Int}},
    score::Matrix{Float64},
    triangle_scores::Array{Float64,3},
    i::Int,
    j::Int,
    bound::Float64,
)
    j - i < 2 && return triangles
    for k in (i + 1):(j - 1)
        (score[i, k] ≤ bound && triangle_scores[i, k, j] ≤ bound && score[k, j] ≤ bound) || continue
        first_triangles!(triangles, score, triangle_scores, i, k, bound)
        push!(triangles, (i, k, j))
        first_triangles!(triangles, score, triangle_scores, k, j, bound)
        return triangles
    end
    return norma_abort("No triangulation of a swap polygon meets the bound of its least score.")
end

# Try to swap the interior edge (a, b): find the triangulation of its ring
# of least score (best_triangulation) and submit it to the acceptance test.
# Returns the accepted proposal or nothing.
function try_edge_swap(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    options::AdaptivityOptions;
    created::CreatedEntities=CreatedEntities(),
)
    best = best_edge_swap(model, topology, a, b, options, created)
    (best === nothing || best.connectivity === nothing) && return nothing
    return accept_proposal(model, topology, best.elements, best.connectivity, best.block, options)
end

# The swap of an interior edge evaluated ahead of the sequential loop of a
# pass, on the topology at the start of its chunk and with no entity of the
# pass excluded: the best triangulation and the outcome of its acceptance
# test (issue #231).
struct EdgeSwapEvaluation
    best::Any
    proposal::Union{CavityProposal,Nothing}
end

function evaluate_edge_swap(model::SolidMechanics, topology::MeshTopology, a::Int, b::Int, options::AdaptivityOptions)
    best = best_edge_swap(model, topology, a, b, options, CreatedEntities())
    (best === nothing || best.connectivity === nothing) && return EdgeSwapEvaluation(best, nothing)
    return EdgeSwapEvaluation(
        best, accept_proposal(model, topology, best.elements, best.connectivity, best.block, options)
    )
end

# The result of try_edge_swap for the current state of the pass, from an
# evaluation made before the operations accepted since.  A dead ring element
# gives nothing, as edge_ring does.  The evaluation chose among a superset of
# the triangulations admissible now.  When its best is still admissible,
# best_triangulation returns it again.  Under the scaled Jacobian criterion
# it is the first of least score in a fixed order, which it remains in a
# subset.  Under the energy criterion the least scores of the sub-polygons
# of that triangulation are unchanged, those of the other sub-polygons can
# only have increased, and each sub-polygon keeps its first apex of least
# score.  The outcome of the acceptance test, which depends on the cavity
# alone, is then the same.  Otherwise the swap is evaluated again.  Under
# cavity relaxation the test also depends on the positions of the cavity
# nodes and on the elements around them, so the swap is evaluated again
# when a node of the cavity was moved in the pass or an element around one
# died since (cavity_surroundings_intact).
function edge_swap_from_evaluation(
    model::SolidMechanics,
    topology::MeshTopology,
    a::Int,
    b::Int,
    options::AdaptivityOptions,
    evaluation::EdgeSwapEvaluation,
    created::CreatedEntities,
)
    best = evaluation.best
    best === nothing && return nothing
    all(topology.element_alive[e] for e in best.elements) || return nothing
    # The admissible triangulations depend on the positions as well, so
    # this test precedes the one on the connectivity found.
    if options.cavity_relaxation &&
        !cavity_surroundings_intact(topology, created, best.elements, (a, b, best.ring_nodes...))
        return try_edge_swap(model, topology, a, b, options; created)
    end
    best.connectivity === nothing && return nothing
    if chords_exist(topology, best.ring_nodes, best.triangles, created) ||
        !faces_are_new(topology, best.elements, best.connectivity, created)
        return try_edge_swap(model, topology, a, b, options; created)
    end
    return evaluation.proposal
end

# Whether every element around the nodes of a cavity outside the cavity is
# alive and no node of those elements, or of the cavity, was moved in the
# pass, so that an evaluation of the cavity made earlier in the pass under
# cavity relaxation, which reads exactly those elements and positions, is
# still current.
function cavity_surroundings_intact(topology::MeshTopology, created::CreatedEntities, elements::Vector{Int}, nodes)
    for n in nodes
        n in created.moved && return false
        n < length(topology.node_to_elements_offsets) || continue
        for e in node_elements(topology, n)
            e in elements && continue
            topology.element_alive[e] || return false
            any(m in created.moved for m in view(topology.connectivity, :, e)) && return false
        end
    end
    return true
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
    others = opposite_node_pairs(topology, elements, a, b)
    others === nothing && return nothing
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
    walk = walk_edge_elements(elements, others, start, first_node)
    walk === nothing && return nothing
    ordered, chain = walk
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
    block = cavity_block(topology, elements)
    block === nothing && return nothing
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
    # The side (p_1, p_n) of the polygon is the new edge, absent as checked
    # above, so the chords that best_triangulation checks are all the others.
    best, _ = best_triangulation(model, topology, a, b, elements, nodes, block, options, created)
    best === nothing && return nothing
    new_faces = [sorted_face(a, nodes[1], nodes[n]), sorted_face(b, nodes[1], nodes[n])]
    surface = side_set > 0 ? SurfaceSwap(side_set, collect(old_faces), new_faces) : nothing
    return accept_proposal(model, topology, elements, best, block, options; surface)
end

function record!(created::CreatedEntities, connectivity::AbstractMatrix{<:Integer})
    for column in eachcol(connectivity)
        for edge in tetrahedron_edges(column)
            push!(created.edges, edge)
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
    # A cavity has a few tens of faces, which a vector searches faster than
    # a set hashes them.
    old_faces = NTuple{3,Int}[]
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
    block = cavity_block(topology, incident)
    block === nothing && return nothing
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
    block = cavity_block(topology, star)
    block === nothing && return nothing
    kept = Int[]
    for e in star
        a in view(topology.connectivity, :, e) || push!(kept, e)
    end
    isempty(kept) && return nothing
    new_connectivity = replace!(topology.connectivity[:, kept], b => a)
    faces_are_new(topology, star, new_connectivity, created) || return nothing
    if options.size_by_length
        # No edge from the surviving node may become longer than the band.
        for q in unique(vec(new_connectivity))
            q == a && continue
            metric_edge_length(model, topology, a, q) ≤ LENGTH_BAND[2] || return nothing
        end
    end
    return accept_proposal(
        model,
        topology,
        star,
        new_connectivity,
        block,
        options;
        require_decrease=(!by_length),
        removed_node=b,
        surviving_node=a,
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

# Energy of the elements around a new node as a function of its position,
# with the target sampled at `sample`, the positions with the node at its
# initial position during the relaxation, and at the positions themselves
# for the acceptance test.
function split_cavity_energy(
    model::SolidMechanics,
    topology::MeshTopology,
    block::Int,
    connectivity::Matrix{Int},
    position::SVector{3,Float64},
    metric::Union{MetricField,Nothing};
    sample::Union{PositionsWithNode,Nothing}=nothing,
)
    positions = PositionsWithNode(topology.positions, position)
    sample_positions = sample === nothing ? positions : sample
    return sum(element_energies(model, block, connectivity, positions; sample_positions, metric))
end

# Local relaxation of a new node with the rest of the cavity fixed: a few
# steps of descent on the energy of the new elements, with the gradient of
# the energy with respect to the node from element_energy_gradients (one
# evaluation of the cavity instead of the six of central differences), each
# step returned to the node's surfaces.  The operation is then judged at its
# locally relaxed position rather than at the midpoint.
#
# Within a relaxation the target is sampled at the positions before it (the
# new node at its initial position), so that the gradient is the exact
# derivative of the energy descended on; the acceptance test and the stored
# energies sample the target at the current positions, as everywhere in a
# topology phase, and the energy of the relaxed configuration is evaluated
# again in that convention after the relaxation.  The two differ only where
# the target depends on the positions (the equal-volume rule, a size field
# varying over an element, a metric field); the derivative of the energy
# through the target, the configurational force of the fictitious problem,
# is what the gradient within the relaxation leaves out.
const SPLIT_RELAXATION_STEPS = 5

# Gradient of the energy of the elements of `connectivity` with respect to
# node `n`, from their gradients.
function node_gradient(gradients::Vector{SMatrix{3,4,Float64,12}}, connectivity::AbstractMatrix{<:Integer}, n::Int)
    g = zero(SVector{3,Float64})
    for e in eachindex(gradients), a in 1:4
        connectivity[a, e] == n && (g += gradients[e][:, a])
    end
    return g
end

# One descent step of a node on the energy `energy(x)` of the elements
# around it from the gradient `g` at `x` with value `value`: the step length
# starts at a tenth of the node's scale and is halved up to eight times
# until the energy decreases, each trial returned by `project`.  Returns the
# new position and value, or nothing when no step decreased the energy.  (A
# second trial at the minimizer of the quadratic fitted through the value,
# the slope, and the first trial was measured on the awful cube loop: more
# operations accepted and a lower energy, 6.75 against 7.08, but a lower
# minimum scaled Jacobian, 0.363 against 0.382, and 47% more time per
# relaxation; not kept.)
function node_descent_step(
    energy::Function, project::Function, x::SVector{3,Float64}, value::Float64, g::SVector{3,Float64}, scale::Float64
)
    all(isfinite, g) || return nothing
    step_length = 0.1 * scale / max(norm(g), floatmin(Float64))
    for _ in 1:8
        trial = project(x - step_length * g)
        trial_value = energy(trial)
        isfinite(trial_value) && trial_value < value && return trial, trial_value
        step_length *= 0.5
    end
    return nothing
end

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
    sample = PositionsWithNode(topology.positions, position)
    energy(x) = split_cavity_energy(model, topology, block, connectivity, x, metric; sample)
    project(x) = project_to_surfaces(x, constraints, model.time)
    m = size(topology.positions, 2) + 1
    x = position
    value = energy(x)
    isfinite(value) || return x, value
    for _ in 1:SPLIT_RELAXATION_STEPS
        _, gradients = element_energy_gradients(
            model, block, connectivity, PositionsWithNode(topology.positions, x); sample_positions=sample, metric
        )
        step = node_descent_step(energy, project, x, value, node_gradient(gradients, connectivity, m), scale)
        step === nothing && break
        x, value = step
    end
    x == position && return x, value
    return x, split_cavity_energy(model, topology, block, connectivity, x, metric)
end

# Try to split the edge (a, b) at a new node: every element around the edge
# is bisected.  The node starts at the midpoint, is returned to the surfaces
# of the boundary faces that contain the edge, inherits the side sets of
# those faces and the node sets common to both ends of a boundary edge, is
# relaxed locally, and the result is submitted to the acceptance test.
function try_edge_split(
    model::SolidMechanics, topology::MeshTopology, a::Int, b::Int, options::AdaptivityOptions; by_length::Bool=false
)
    (topology.node_alive[a] && topology.node_alive[b]) || return nothing
    ring = edge_elements(topology, a, b)
    isempty(ring) && return nothing
    block = cavity_block(topology, ring)
    block === nothing && return nothing
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
        position, energy_after = relax_split_node(
            model, topology, block, connectivity, position, constraints, norm(xb - xa), metric
        )
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
    passes_acceptance_tests(
        model,
        topology,
        collect(ring),
        connectivity,
        block,
        positions,
        metric,
        nothing,
        energy_before,
        energy_after,
        options;
        require_decrease=(!by_length),
    ) || return nothing
    split = SplitNode(position, node_sets, side_sets, (a, b), node_metric)
    return CavityProposal(
        collect(ring), connectivity, block, energy_before, energy_after, 0, 0, split, nothing, Int[], zeros(3, 0)
    )
end

# Edges whose length in the prescribed target satisfies `outside`, a test
# against one bound of the length band; empty without a target.
function band_edges(outside::Function, model::SolidMechanics, topology::MeshTopology)
    edges = Tuple{Int,Int}[]
    (model.metric_field === nothing && model.size_field === nothing) && return edges
    for edge in keys(topology.edges)
        outside(metric_edge_length(model, topology, edge[1], edge[2])) && push!(edges, edge)
    end
    return edges
end

# Edges longer than sqrt(2) in the prescribed target, the split candidates of
# the size phase, and edges shorter than 1/sqrt(2), the collapse candidates.
long_edges(model::SolidMechanics, topology::MeshTopology) = band_edges(>(LENGTH_BAND[2]), model, topology)
short_edges(model::SolidMechanics, topology::MeshTopology) = band_edges(<(LENGTH_BAND[1]), model, topology)

# The candidate edges of a pass of collapses or splits: the edges of the
# candidate elements when the shape-driven operations are on, and the edges
# outside the length band always, less those refused earlier in the phase.
function candidate_edges(
    model::SolidMechanics,
    topology::MeshTopology,
    densities::Vector{Float64},
    quality::Vector{Float64},
    options::AdaptivityOptions,
    shape::Bool,
    size_edges::Set{Tuple{Int,Int}},
    refused::Set{Tuple{Int,Int}},
)
    edges = Set{Tuple{Int,Int}}()
    if shape
        for e in candidate_elements(model, topology, densities, options, quality)
            for edge in element_edges(topology, e)
                push!(edges, edge)
            end
        end
    end
    union!(edges, size_edges)
    setdiff!(edges, refused)
    return collect(edges)
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
    densities, quality = phase_element_values!(memory, model, topology, options)
    size_edges = Set(long_edges(model, topology))
    ranked = candidate_edges(model, topology, densities, quality, options, shape, size_edges, memory.splits)
    # The edges beyond the band are split longest first, in the target, as
    # in longest-edge bisection, which is what bounds the shape of the
    # children; the shape-driven candidates follow by their cavity rank.
    function split_rank(edge)
        (a, b) = edge
        if options.size_by_length && edge in size_edges
            return (1, metric_edge_length(model, topology, a, b))
        end
        return (0, cavity_rank(model, topology, densities, edge_elements(topology, a, b), options, quality))
    end
    ranked = ranked_by(split_rank, ranked)
    accepted = 0
    decrease = 0.0
    # The splits are evaluated in parallel, a chunk of the ranked edges at a
    # time, and accepted in rank order; the result is that of the sequential
    # loop (see edge_split_from_evaluation).
    for chunk in Iterators.partition(ranked, SPLIT_EVALUATION_CHUNK)
        evaluations = Vector{Union{CavityProposal,Nothing}}(undef, length(chunk))
        evaluated_node = size(topology.positions, 2) + 1
        Threads.@threads for k in eachindex(chunk)
            (a, b) = chunk[k]
            by_length = options.size_by_length && (a, b) in size_edges
            evaluations[k] = try_edge_split(model, topology, a, b, options; by_length)
        end
        for (k, (a, b)) in enumerate(chunk)
            ring_intact = all(topology.element_alive[e] for e in edge_elements(topology, a, b))
            proposal = edge_split_from_evaluation(topology, evaluations[k], evaluated_node, ring_intact)
            if proposal === nothing
                ring_intact && push!(memory.splits, (a, b))
                continue
            end
            apply!(topology, proposal; metric=model.metric_field)
            accepted += 1
            decrease += proposal.energy_before - proposal.energy_after
        end
    end
    return accepted, decrease
end

# Number of ranked edges whose splits are evaluated in parallel at a time.
const SPLIT_EVALUATION_CHUNK = 1024

# The result of try_edge_split for the current state of the pass, from an
# evaluation made before the splits accepted since (issue #231).  A split
# reads its ring, the positions of existing nodes, which do not move, and
# the set memberships of the faces of its ring, which an earlier split
# changes only for the faces of its own cavity, whose elements are dead.  So
# with its ring intact the evaluation is the sequential result, except for
# the index of the new node: the evaluation numbered it `evaluated_node`,
# and the splits accepted since have taken that index and those after it.
function edge_split_from_evaluation(
    topology::MeshTopology, evaluation::Union{CavityProposal,Nothing}, evaluated_node::Int, ring_intact::Bool
)
    (evaluation === nothing || !ring_intact) && return nothing
    node = size(topology.positions, 2) + 1
    if node != evaluated_node
        replace!(evaluation.new_connectivity, evaluated_node => node)
    end
    return evaluation
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
        return norm(metric_factor(h, R) * (xb - xa))
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
    densities, quality = phase_element_values!(memory, model, topology, options)
    size_edges = Set(short_edges(model, topology))
    ranked = candidate_edges(model, topology, densities, quality, options, shape, size_edges, memory.collapses)
    edge_energy(edge) =
        cavity_rank(model, topology, densities, edge_elements(topology, edge[1], edge[2]), options, quality)
    ranked = ranked_by(edge_energy, ranked)
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
        apply_accepted!(memory, model, topology, options, created, proposal)
        accepted += 1
        decrease += proposal.energy_before - proposal.energy_after
    end
    return accepted, decrease
end

# Candidate elements: those above the desired density, dilated by the given
# number of layers of node adjacency.
function candidate_elements(
    model::SolidMechanics,
    topology::MeshTopology,
    densities::Vector{Float64},
    options::AdaptivityOptions,
    quality::Vector{Float64}=Float64[],
)
    candidates = falses(length(densities))
    if options.shape_by_quality
        # An operation accepted on the cavity minimum can only help the
        # elements below the desired quality, and every edge of such an
        # element is incident to it, so no dilation is needed.
        alive = findall(topology.element_alive)
        if length(quality) == length(densities)
            for e in alive
                quality[e] < options.desired_quality && (candidates[e] = true)
            end
            return findall(candidates)
        end
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

# The items in decreasing order of their rank, each rank computed once.
# `sort!` with `by` evaluates the key at every comparison, which for the
# cavity ranks recomputed the energies or the scaled Jacobians of each
# cavity about log2(n) times (issue #231). Both sorts are stable, so the
# order, and with it the result of a pass, is unchanged.
function ranked_by(rank, items::Vector)
    ranks = map(rank, items)
    return items[sortperm(ranks; rev=true)]
end

# Score by which the operations of a pass are ranked: the energy of the
# elements around an edge or a face, or, under the scaled Jacobian
# criterion, the negative of their minimum scaled Jacobian, so that the
# worst cavities are tried first in either case.
function cavity_rank(
    model::SolidMechanics,
    topology::MeshTopology,
    densities::Vector{Float64},
    elements,
    options::AdaptivityOptions,
    quality::Vector{Float64}=Float64[],
)
    if options.shape_by_quality
        isempty(elements) && return -1.0
        length(quality) == length(densities) && return -minimum(quality[e] for e in elements)
        return -minimum(quality_jacobians(model, topology.positions, topology.connectivity[:, collect(elements)]))
    end
    return sum(densities[e] for e in elements; init=0.0)
end

# Number of ranked edges whose swaps are evaluated in parallel at a time.
const SWAP_EVALUATION_CHUNK = 4096

# One pass of swaps over the candidate elements in decreasing order of
# cavity energy: the boundary edges (when enabled), then the interior edges,
# then the faces (when enabled).  An edge or face whose cavity was changed
# earlier in the pass is deferred to the next pass.  Returns the number of
# accepted swaps and the total decrease of the energy.
function swap_pass!(
    model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions; memory::PhaseMemory=PhaseMemory()
)
    densities, quality = phase_element_values!(memory, model, topology, options)
    candidates = candidate_elements(model, topology, densities, options, quality)
    edges = Set{Tuple{Int,Int}}()
    boundary_edges = Set{Tuple{Int,Int}}()
    faces = Set{NTuple{3,Int}}()
    for e in candidates
        for edge in element_edges(topology, e)
            if edge in topology.boundary_edges
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
    ring_energy(edge) =
        cavity_rank(model, topology, densities, edge_elements(topology, edge[1], edge[2]), options, quality)
    face_energy(face) = cavity_rank(model, topology, densities, get(topology.faces, face, Int[]), options, quality)
    ranked_boundary = ranked_by(ring_energy, collect(boundary_edges))
    ranked = ranked_by(ring_energy, collect(edges))
    ranked_faces = ranked_by(face_energy, collect(faces))
    accepted = 0
    decrease = 0.0
    created = CreatedEntities()
    intact(a, b) = all(topology.element_alive[e] for e in edge_elements(topology, a, b))
    for (a, b) in ranked_boundary
        proposal = try_boundary_edge_swap(model, topology, a, b, options; created)
        if apply_swap!(memory, model, topology, options, created, proposal, (a, b), memory.boundary_swaps, intact(a, b))
            accepted += 1
            decrease += proposal.energy_before - proposal.energy_after
        end
    end
    # The swaps are evaluated in parallel, a chunk of the ranked edges at a
    # time, and accepted in rank order; the result is that of the sequential
    # loop (see edge_swap_from_evaluation).
    for chunk in Iterators.partition(ranked, SWAP_EVALUATION_CHUNK)
        evaluations = Vector{EdgeSwapEvaluation}(undef, length(chunk))
        Threads.@threads for k in eachindex(chunk)
            evaluations[k] = evaluate_edge_swap(model, topology, chunk[k][1], chunk[k][2], options)
        end
        for (k, (a, b)) in enumerate(chunk)
            proposal = edge_swap_from_evaluation(model, topology, a, b, options, evaluations[k], created)
            if apply_swap!(memory, model, topology, options, created, proposal, (a, b), memory.swaps, intact(a, b))
                accepted += 1
                decrease += proposal.energy_before - proposal.energy_after
            end
        end
    end
    for face in ranked_faces
        proposal = try_face_swap(model, topology, face, options; created)
        cavity_intact = all(topology.element_alive[e] for e in topology.faces[face])
        if apply_swap!(memory, model, topology, options, created, proposal, face, memory.face_swaps, cavity_intact)
            accepted += 1
            decrease += proposal.energy_before - proposal.energy_after
        end
    end
    return accepted, decrease
end

# Apply a swap proposal and record the entities it creates; returns whether
# it was applied.  A refused proposal is refused for the phase when its
# cavity is intact; when the cavity was changed earlier in the pass the edge
# or face is deferred to the next pass instead.
function apply_swap!(
    memory::PhaseMemory,
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions,
    created::CreatedEntities,
    proposal::Union{CavityProposal,Nothing},
    key,
    refused::Set,
    cavity_intact::Bool,
)
    if proposal === nothing
        cavity_intact && push!(refused, key)
        return false
    end
    apply_accepted!(memory, model, topology, options, created, proposal)
    return true
end

# Apply an accepted proposal within a pass: the topology, the record of the
# entities the operation created, and, when it moved nodes, the record of
# the moved nodes and the memory of the phase.  The elements around a moved
# node that the operation keeps have changed, so their stored values are
# evaluated again and the refusals of their edges and faces are forgotten;
# the elements the operation added are evaluated at the compaction.
function apply_accepted!(
    memory::PhaseMemory,
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions,
    created::CreatedEntities,
    proposal::CavityProposal,
)
    apply!(topology, proposal; metric=model.metric_field)
    record!(created, proposal.new_connectivity)
    isempty(proposal.moved_nodes) && return nothing
    union!(created.moved, proposal.moved_nodes)
    changed = Int[]
    for n in proposal.moved_nodes, e in node_elements(topology, n)
        (topology.element_alive[e] && !(e in changed)) && push!(changed, e)
    end
    forget_elements!(memory, topology, changed)
    isempty(memory.density) || energy_densities!(memory.density, model, topology, changed)
    options.shape_by_quality && !isempty(memory.quality) && element_qualities!(memory.quality, model, topology, changed)
    return nothing
end

# The topology is compacted after every operator, so that each one starts
# from a current adjacency: the operations of one operator are kept
# independent by the dead elements of their cavities, but the next operator
# would otherwise not see the elements they made.  An operator that accepted
# nothing left the topology unchanged (every accepted operation adds
# elements from index `first_new` on), so there is nothing to compact.
function compact_phase!(
    memory::PhaseMemory, model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions, first_new::Int
)
    first_new > length(topology.element_alive) && return nothing
    forget_changed!(memory, topology, first_new)
    node_map, element_map = compact!(topology)
    compact_metric!(model.metric_field, findall(>(0), node_map))
    remap!(memory, node_map)
    remap_element_values!(memory, model, topology, element_map, first_new, options)
    return nothing
end

# The result of an operator that did not run.
const NO_OPERATOR_RESULT = (accepted=0, decrease=0.0, time=0.0)

# Run one operator of a topology pass, `operator_pass!`, a function without
# arguments that returns the number of accepted operations and the decrease
# of the energy, and compact the topology after it.  Returns the number, the
# decrease, and the wall time of the operator with the compaction.
function run_operator!(
    operator_pass!::Function,
    memory::PhaseMemory,
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions,
)
    start = time()
    first_new = length(topology.element_alive) + 1
    accepted, decrease = operator_pass!()
    compact_phase!(memory, model, topology, options, first_new)
    return (accepted=accepted, decrease=decrease, time=time() - start)
end

function swap_operator!(memory::PhaseMemory, model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions)
    options.swaps || return NO_OPERATOR_RESULT
    return run_operator!(() -> swap_pass!(model, topology, options; memory), memory, model, topology, options)
end

# The collapses and then the splits of a pass, each when enabled.
function size_operators!(
    memory::PhaseMemory, model::SolidMechanics, topology::MeshTopology, options::AdaptivityOptions, shape::Bool
)
    collapses = splits = NO_OPERATOR_RESULT
    if options.collapses
        collapses = run_operator!(memory, model, topology, options) do
            return collapse_pass!(model, topology, options, shape; memory)
        end
    end
    if options.splits
        splits = run_operator!(memory, model, topology, options) do
            return split_pass!(model, topology, options, shape; memory)
        end
    end
    return collapses, splits
end

# The topology phase: passes of operations until none is accepted or the
# cap is reached.  Compacts the topology after every operator, so the
# adjacency is current at the start of the next.  Returns the number of
# accepted operations and the total decrease of the energy.
# The time of each pass and of each operator in it is logged after the pass
# and, when times are given, added to them.
function topology_phase!(
    model::SolidMechanics,
    topology::MeshTopology,
    options::AdaptivityOptions;
    times::Union{AdaptivityTimes,Nothing}=nothing,
)
    total_accepted = 0
    total_decrease = 0.0
    memory = PhaseMemory()
    previous_swaps = 1
    for pass in 1:options.maximum_passes
        pass_start = time()
        # The edges outside the length band of a prescribed target are
        # collapsed or split in every pass.  Collapses and splits driven by
        # the shape alone remove or add resolution, so they are tried only in
        # a pass where no swap was accepted (in the pass before, when the size
        # operations come first).  The decrease of the pass is summed in the
        # order in which the operators run.
        if options.size_first
            shape = options.shape_every_pass || previous_swaps == 0
            collapses, splits = size_operators!(memory, model, topology, options, shape)
            swaps = swap_operator!(memory, model, topology, options)
            decrease = collapses.decrease + splits.decrease + swaps.decrease
        else
            swaps = swap_operator!(memory, model, topology, options)
            shape = options.shape_every_pass || swaps.accepted == 0
            collapses, splits = size_operators!(memory, model, topology, options, shape)
            decrease = swaps.decrease + collapses.decrease + splits.decrease
        end
        previous_swaps = swaps.accepted
        accepted = swaps.accepted + collapses.accepted + splits.accepted
        norma_logf(
            0,
            :info,
            "Topology pass %d: %d operations accepted (%d swaps, %d collapses, %d splits), energy decrease %.6e",
            pass,
            accepted,
            swaps.accepted,
            collapses.accepted,
            splits.accepted,
            decrease,
        )
        time_pass = time() - pass_start
        norma_log(
            0,
            :time,
            "Topology pass $pass time = $(format_time(time_pass)) (swaps $(format_time(swaps.time)), " *
            "collapses $(format_time(collapses.time)), splits $(format_time(splits.time)))",
        )
        if times !== nothing
            times.passes += time_pass
            times.swaps += swaps.time
            times.collapses += collapses.time
            times.splits += splits.time
        end
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

# The name of the file at the given position of a sequence that ParaView
# opens as one: the first file keeps its name, and the others append -s_
# and the position in four or more digits.
function sequence_file_name(file::AbstractString, position::Integer)
    return position == 1 ? String(file) : file * "-s_" * lpad(position, 4, '0')
end

# The coupled loop: smoothing on the current mesh, a topology phase, a new
# mesh written to disk, and smoothing again on it, until a topology phase
# accepts nothing or the outer iterations are exhausted.  Each adapted mesh
# is a separate Exodus file, numbered after the output name, and so is the
# output of each smoothing phase: the first under the output name, the
# following ones with the suffixes -s_0002, -s_0003, and so on, which
# ParaView opens as one sequence.  The times written in a file continue
# from the file before it, one time step after its last, since ParaView
# orders the sequence by time and shows only one file at a repeated time.
function run_adaptive(params::Parameters)
    options = AdaptivityOptions(get(params, "adaptivity", Parameters()))
    output_file = params["output mesh file"]
    name = stripped_name(output_file)
    run_start = time()
    times = AdaptivityTimes()
    times.setup += @elapsed sim = create_simulation(params)
    times.smoothing += @elapsed run(sim)
    model = sim.model
    model isa SolidMechanics && model.mesh_smoothing || norma_abort("Adaptivity requires a mesh smoothing model")
    for iteration in 1:options.outer_iterations
        build_start = time()
        topology = build_topology(model)
        energy_before = sum(
            energy_densities(model, topology) .*
            ideal_element_volumes(model, 1, topology.connectivity, topology.positions),
        )
        time_build = time() - build_start
        times.build += time_build
        passes_before = times.passes
        accepted, decrease = topology_phase!(model, topology, options; times)
        time_passes = times.passes - passes_before
        observer = TOPOLOGY_PHASE_OBSERVER[]
        observer === nothing || observer(model, topology, accepted, decrease)
        norma_logf(
            0,
            :info,
            "Adaptivity iteration %d: %d operations accepted, energy %.6e -> %.6e",
            iteration,
            accepted,
            energy_before,
            energy_before - decrease,
        )
        time_write = 0.0
        if accepted > 0
            mesh_file = "$name-adapted-$iteration.g"
            time_write = @elapsed write_topology(
                topology, mesh_file; nodal_variables=metric_nodal_variables(model.metric_field)
            )
            times.write += time_write
        end
        norma_log(
            0,
            :time,
            "Topology phase time = $(format_time(time_build + time_passes + time_write)) " *
            "(build $(format_time(time_build)), passes $(format_time(time_passes)), " *
            "mesh write $(format_time(time_write)))",
        )
        accepted == 0 && break
        next_params = deepcopy(params)
        next_params["input mesh file"] = mesh_file
        next_params["output mesh file"] = sequence_file_name(output_file, iteration + 1)
        controller = sim.controller
        time_offset = get(sim.params, "exodus_time_offset", 0.0)
        elapsed = controller.time - controller.initial_time
        next_params["exodus_time_offset"] = time_offset + elapsed + controller.time_step
        times.setup += @elapsed sim = create_simulation(next_params)
        times.smoothing += @elapsed run(sim)
        model = sim.model
    end
    norma_log(0, :time, adaptivity_time_report(time() - run_start, times))
    return sim
end

function adaptivity_time_report(wall::Float64, times::AdaptivityTimes)::String
    return "Adaptivity Time = $(format_time(wall)) (smoothing $(format_time(times.smoothing)), " *
           "topology $(format_time(times.build + times.passes)), mesh write $(format_time(times.write)), " *
           "setup $(format_time(times.setup)))"
end
