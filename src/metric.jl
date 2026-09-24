# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Metric tensor targets for energetic mesh smoothing
# (docs/notes/ems-anisotropic): the sources of the field, their evaluation on
# an element, the factorization of a tensor into sizes and a rotation with a
# frame that follows from node to node, and the nodal representation that the
# adaptivity loop carries through its operations and writes to the adapted
# mesh.

# The unit regular tetrahedron of the smoothing rules, with edges of unit
# length; the ideal element is this shape scaled by the target size.
const UNIT_TETRAHEDRON = SMatrix{3,4,Float64,12}(
    0.5 / sqrt(2.0) * [
        1 -1 -1 1
        1 -1 1 -1
        1 1 -1 -1
    ],
)

# Order of the six components of a symmetric tensor in the input and the
# output: M_xx, M_yy, M_zz, M_xy, M_yz, M_zx.
const METRIC_COMPONENT_SUFFIXES = ("xx", "yy", "zz", "xy", "yz", "zx")

# Relative separation below which two eigenvalues count as repeated, so the
# in-plane phase of their eigenvectors is fixed by the reference frame.
const METRIC_REPEATED_EIGENVALUE_TOLERANCE = 1.0e-6

# Rotation angle between the frames of a node and its parent above which the
# principal frame is reported as not followed.
const METRIC_FRAME_JUMP_ANGLE = π / 4

function compile_field_expression(expression)
    return eval(build_function(eval(Meta.parse(string(expression))), [t, x, y, z]; expression=Val(false)))
end

function symmetric_from_components(c::AbstractVector{Float64})
    return SMatrix{3,3,Float64,9}(c[1], c[4], c[6], c[4], c[2], c[5], c[6], c[5], c[3])
end

function components_from_symmetric(M::SMatrix{3,3,Float64,9})
    return SVector{6,Float64}(M[1, 1], M[2, 2], M[3, 3], M[1, 2], M[2, 3], M[3, 1])
end

# Matrix logarithm of a symmetric positive definite tensor, or nothing when
# the tensor is not positive definite.
function log_symmetric(M::SMatrix{3,3,Float64,9})
    e = eigen(Symmetric(M))
    all(λ -> isfinite(λ) && λ > 0.0, e.values) || return nothing
    V = e.vectors
    return SMatrix{3,3,Float64,9}(V * Diagonal(log.(e.values)) * V')
end

function exp_symmetric(L::SMatrix{3,3,Float64,9})
    e = eigen(Symmetric(L))
    V = e.vectors
    return SMatrix{3,3,Float64,9}(V * Diagonal(exp.(e.values)) * V')
end

# Sizes and rotation of a metric tensor M = R diag(1/h²) R', with the frame
# chosen closest to `reference` among the equivalent ones: the columns of R
# are the eigenvectors of M in the order and with the signs that minimize the
# rotation from the reference, and for repeated eigenvalues the free in-plane
# phase is taken from the projection of the reference frame onto the
# eigenspace.  The energy of the smoother does not depend on this choice; the
# interpolation of sizes and rotation vectors between nodes does.  Returns
# nothing when M is not positive definite.
function principal_of_tensor(M::SMatrix{3,3,Float64,9}, reference::SMatrix{3,3,Float64,9}=SMatrix{3,3,Float64,9}(I))
    e = eigen(Symmetric(M))
    λ = e.values
    all(v -> isfinite(v) && v > 0.0, λ) || return nothing
    V = e.vectors
    tol = METRIC_REPEATED_EIGENVALUE_TOLERANCE * λ[3]
    if λ[3] - λ[1] ≤ tol
        # Isotropic: any frame is principal; keep the reference.
        h = 1.0 / sqrt((λ[1] + λ[2] + λ[3]) / 3.0)
        return SVector{3,Float64}(h, h, h), reference
    end
    if λ[2] - λ[1] ≤ tol || λ[3] - λ[2] ≤ tol
        # Transversely isotropic: the odd eigenvector is the axis, the plane
        # normal to it takes its in-plane phase from the reference.
        axis_index = λ[2] - λ[1] ≤ tol ? 3 : 1
        pair = axis_index == 3 ? (1, 2) : (2, 3)
        n = SVector{3,Float64}(V[:, axis_index])
        h_axis = 1.0 / sqrt(λ[axis_index])
        h_pair = 1.0 / sqrt((λ[pair[1]] + λ[pair[2]]) / 2.0)
        # The reference column most aligned with the axis takes the axis.
        dots = SVector{3,Float64}(dot(reference[:, 1], n), dot(reference[:, 2], n), dot(reference[:, 3], n))
        j = argmax(abs.(dots))
        n = dots[j] < 0.0 ? -n : n
        j1, j2 = j == 1 ? (2, 3) : (j == 2 ? (3, 1) : (1, 2))
        p = SVector{3,Float64}(reference[:, j1])
        p = p - dot(p, n) * n
        norm_p = norm(p)
        if norm_p ≤ eps(Float64)
            # The reference column lies along the axis; take the pair's own eigenvector.
            p = SVector{3,Float64}(V[:, pair[1]])
        else
            p = p / norm_p
        end
        q = cross(n, p)
        columns = (SVector{3,Float64}(zeros(3)), SVector{3,Float64}(zeros(3)), SVector{3,Float64}(zeros(3)))
        columns = Base.setindex(columns, n, j)
        columns = Base.setindex(columns, p, j1)
        columns = Base.setindex(columns, q, j2)
        R = SMatrix{3,3,Float64,9}(hcat(columns[1], columns[2], columns[3]))
        sizes = (0.0, 0.0, 0.0)
        sizes = Base.setindex(sizes, h_axis, j)
        sizes = Base.setindex(sizes, h_pair, j1)
        sizes = Base.setindex(sizes, h_pair, j2)
        return SVector{3,Float64}(sizes), R
    end
    # Distinct eigenvalues: the closest of the 24 proper signed permutations of
    # the eigenvectors, measured by the trace of reference' R (the cosine of
    # the rotation angle between the frames).
    best_trace = -Inf
    best_R = SMatrix{3,3,Float64,9}(I)
    best_h = SVector{3,Float64}(1.0, 1.0, 1.0)
    h = SVector{3,Float64}(1.0 / sqrt(λ[1]), 1.0 / sqrt(λ[2]), 1.0 / sqrt(λ[3]))
    for perm in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)), s1 in (1.0, -1.0),
        s2 in (1.0, -1.0), s3 in (1.0, -1.0)

        R = SMatrix{3,3,Float64,9}(hcat(s1 * V[:, perm[1]], s2 * V[:, perm[2]], s3 * V[:, perm[3]]))
        det(R) > 0.0 || continue
        trace = tr(reference' * R)
        if trace > best_trace
            best_trace = trace
            best_R = R
            best_h = SVector{3,Float64}(h[perm[1]], h[perm[2]], h[perm[3]])
        end
    end
    return best_h, best_R
end

function metric_tensor(h::SVector{3,Float64}, R::SMatrix{3,3,Float64,9})
    return SMatrix{3,3,Float64,9}(R * Diagonal(SVector{3,Float64}(1.0 / h[1]^2, 1.0 / h[2]^2, 1.0 / h[3]^2)) * R')
end

# Nodes adjacent to every node through the elements of a model.
function node_neighbors(model::SolidMechanics)
    num_nodes = size(model.reference, 2)
    neighbors = [Int[] for _ in 1:num_nodes]
    for block in model.blocks
        connectivity = block.connectivity
        for e in 1:size(connectivity, 2), i in 1:size(connectivity, 1), j in 1:size(connectivity, 1)
            i == j && continue
            a, b = connectivity[i, e], connectivity[j, e]
            b in neighbors[a] || push!(neighbors[a], b)
        end
    end
    return neighbors
end

# Sizes (3 × n) and rotation vectors (3 × n) of tensors given at the nodes,
# with the frame followed from node to node: the nodes are visited breadth
# first through the adjacency, each frame is chosen closest to the frame of
# the node it was reached from, and the first node of every component takes
# the frame closest to the global axes.  Frames that turn by more than
# METRIC_FRAME_JUMP_ANGLE from their parent are counted and reported, since
# the rotation vectors then cannot be interpolated between those nodes.
function principal_of_nodal_tensors(tensors::AbstractVector{SMatrix{3,3,Float64,9}}, neighbors::Vector{Vector{Int}})
    num_nodes = length(tensors)
    sizes = zeros(3, num_nodes)
    rotations = zeros(3, num_nodes)
    frames = Vector{SMatrix{3,3,Float64,9}}(undef, num_nodes)
    visited = falses(num_nodes)
    jumps = 0
    queue = Int[]
    for start in 1:num_nodes
        visited[start] && continue
        visited[start] = true
        result = principal_of_tensor(tensors[start])
        result === nothing && norma_abort("The nodal metric tensor at node $start is not positive definite")
        sizes[:, start], frames[start] = result
        rotations[:, start] = rv_of_rt(frames[start])
        empty!(queue)
        push!(queue, start)
        while !isempty(queue)
            parent = popfirst!(queue)
            for n in neighbors[parent]
                visited[n] && continue
                visited[n] = true
                result = principal_of_tensor(tensors[n], frames[parent])
                result === nothing && norma_abort("The nodal metric tensor at node $n is not positive definite")
                sizes[:, n], frames[n] = result
                # The rotation vector is continued past a half turn from the
                # parent's, so that adjacent vectors can be interpolated.
                rotations[:, n] = rv_continue(rv_of_rt(frames[n]), SVector{3,Float64}(rotations[:, parent]))
                cosine = clamp((tr(frames[parent]' * frames[n]) - 1.0) / 2.0, -1.0, 1.0)
                acos(cosine) > METRIC_FRAME_JUMP_ANGLE && (jumps += 1)
                push!(queue, n)
            end
        end
    end
    return sizes, rotations, jumps
end

function read_metric_variables(mesh::ExodusDatabase, names::Vector{String}, time_index::Int)
    available = read_exodus_names(mesh, NodalVariable)
    values = Vector{Vector{Float64}}(undef, length(names))
    for (i, name) in enumerate(names)
        name in available || norma_abort(
            "Nodal variable \"$name\" of the metric field is not in the input mesh '$(mesh.file_name)'; " *
            "available: $(join(available, ", "))",
        )
        values[i] = Float64.(Exodus.read_values(mesh, NodalVariable, time_index, name))
    end
    return values
end

function metric_names(params, key::String, count::Int)
    names = get(params, key, nothing)
    if !(names isa AbstractVector) || length(names) != count
        norma_abort("\"$key\" in \"metric field\" must list $count nodal variable names")
    end
    return String.(names)
end

function metric_time_index(params, mesh::ExodusDatabase)
    num_steps = Exodus.read_number_of_time_steps(mesh)
    num_steps ≥ 1 || norma_abort("The input mesh '$(mesh.file_name)' has no time steps to read nodal metric data from")
    time_index = get(params, "time index", num_steps)
    if !(time_index isa Integer) || time_index < 1 || time_index > num_steps
        norma_abort("\"time index\" in \"metric field\" must be an integer between 1 and $num_steps")
    end
    return Int(time_index)
end

# Compile the anisotropic target of `smooth reference: metric field` (or
# `metric field unrestricted`) from its block: `sizes` with an optional
# `rotation vector` (expressions in t, x, y, z), `tensor` (six expressions),
# `nodal sizes` with an optional `nodal rotation vector` (nodal variables of
# the input mesh), or `nodal tensor` (six nodal variables), the last two read
# at `time index` (default: the last step).  Nodal tensors are factored into
# sizes and rotation vectors with frames followed through the mesh
# (`interpolation: principal`, the default) or kept as logarithms
# (`interpolation: log-Euclidean`).  Returns `nothing` for the other
# smoothing modes.
function create_metric_field(smooth_reference::String, params, mesh::Union{ExodusDatabase,Nothing}=nothing)
    restricted = smooth_reference == "metric field"
    if !restricted && smooth_reference != "metric field unrestricted"
        return nothing
    end
    if !(params isa AbstractDict)
        norma_abort(
            "smooth reference = \"$smooth_reference\" requires a \"metric field\" block with \"sizes\", " *
            "\"tensor\", \"nodal sizes\", or \"nodal tensor\" under the model parameters",
        )
    end
    forms = [key for key in ("sizes", "tensor", "nodal sizes", "nodal tensor") if haskey(params, key)]
    if length(forms) != 1
        norma_abort(
            "\"metric field\" requires exactly one of \"sizes\", \"tensor\", \"nodal sizes\", \"nodal tensor\"; " *
            "got $(isempty(forms) ? "none" : join(forms, ", "))",
        )
    end
    form = forms[1]
    if haskey(params, "rotation vector") && form != "sizes"
        norma_abort("\"rotation vector\" in \"metric field\" goes with \"sizes\"")
    end
    if haskey(params, "nodal rotation vector") && form != "nodal sizes"
        norma_abort("\"nodal rotation vector\" in \"metric field\" goes with \"nodal sizes\"")
    end
    if haskey(params, "interpolation") && form != "nodal tensor"
        norma_abort("\"interpolation\" in \"metric field\" goes with \"nodal tensor\"")
    end
    if haskey(params, "time index") && !startswith(form, "nodal")
        norma_abort("\"time index\" in \"metric field\" goes with \"nodal sizes\" or \"nodal tensor\"")
    end
    if form == "sizes"
        sizes = params["sizes"]
        if !(sizes isa AbstractVector) || length(sizes) != 3
            norma_abort("\"metric field\" requires \"sizes\": three expressions in t, x, y, z (the principal sizes)")
        end
        size_funs = (
            compile_field_expression(sizes[1]), compile_field_expression(sizes[2]), compile_field_expression(sizes[3])
        )
        rotation = get(params, "rotation vector", nothing)
        if rotation === nothing
            rotation_funs = nothing
        elseif rotation isa AbstractVector && length(rotation) == 3
            rotation_funs = (
                compile_field_expression(rotation[1]),
                compile_field_expression(rotation[2]),
                compile_field_expression(rotation[3]),
            )
        else
            norma_abort("\"rotation vector\" in \"metric field\" must have three expressions in t, x, y, z")
        end
        return MetricField(PrincipalMetricFunctions(size_funs, rotation_funs), restricted)
    elseif form == "tensor"
        components = params["tensor"]
        if !(components isa AbstractVector) || length(components) != 6
            norma_abort(
                "\"tensor\" in \"metric field\" must have six expressions in t, x, y, z: the components " *
                join(("M_" * s for s in METRIC_COMPONENT_SUFFIXES), ", "),
            )
        end
        funs = ntuple(i -> compile_field_expression(components[i]), 6)
        return MetricField(TensorMetricFunctions(funs), restricted)
    end
    mesh === nothing && norma_abort("\"$form\" in \"metric field\" needs the input mesh to read nodal variables from")
    time_index = metric_time_index(params, mesh)
    if form == "nodal sizes"
        size_names = metric_names(params, "nodal sizes", 3)
        values = read_metric_variables(mesh, size_names, time_index)
        num_nodes = length(values[1])
        log_sizes = zeros(3, num_nodes)
        for i in 1:3, n in 1:num_nodes
            h = values[i][n]
            if !isfinite(h) || h ≤ 0.0
                norma_abort("Nodal metric size \"$(size_names[i])\" must be positive and finite; got $h at node $n")
            end
            log_sizes[i, n] = log(h)
        end
        rotation = zeros(3, num_nodes)
        rotation_names = String[]
        if haskey(params, "nodal rotation vector")
            rotation_names = metric_names(params, "nodal rotation vector", 3)
            values = read_metric_variables(mesh, rotation_names, time_index)
            for i in 1:3
                rotation[i, :] = values[i]
            end
            all(isfinite, rotation) || norma_abort("Nodal metric rotation vectors must be finite")
        end
        return MetricField(NodalPrincipalMetric(log_sizes, rotation, size_names, rotation_names, String[]), restricted)
    end
    names = metric_names(params, "nodal tensor", 6)
    values = read_metric_variables(mesh, names, time_index)
    num_nodes = length(values[1])
    tensors = Vector{SMatrix{3,3,Float64,9}}(undef, num_nodes)
    for n in 1:num_nodes
        tensors[n] = symmetric_from_components(SVector{6,Float64}(values[1][n], values[2][n], values[3][n],
            values[4][n], values[5][n], values[6][n]))
    end
    interpolation = get(params, "interpolation", "principal")
    if interpolation == "principal"
        neighbors = node_neighbors_from_mesh(mesh, num_nodes)
        sizes, rotations, jumps = principal_of_nodal_tensors(tensors, neighbors)
        if jumps > 0
            norma_logf(
                0,
                :warning,
                "The principal frame of the nodal metric tensor turns by more than %.0f degrees between %d " *
                "pairs of adjacent nodes; the rotation vectors are interpolated across those turns. " *
                "Consider \"interpolation: log-Euclidean\".",
                rad2deg(METRIC_FRAME_JUMP_ANGLE),
                jumps,
            )
        end
        return MetricField(NodalPrincipalMetric(log.(sizes), rotations, String[], String[], names), restricted)
    elseif interpolation == "log-Euclidean"
        log_tensor = zeros(6, num_nodes)
        for n in 1:num_nodes
            L = log_symmetric(tensors[n])
            L === nothing && norma_abort("The nodal metric tensor at node $n is not positive definite")
            log_tensor[:, n] = components_from_symmetric(L)
        end
        return MetricField(NodalTensorMetric(log_tensor, names), restricted)
    end
    return norma_abort(
        "\"interpolation\" in \"metric field\" must be \"principal\" or \"log-Euclidean\"; got \"$interpolation\""
    )
end

# Node adjacency of the input mesh, from the connectivity of its blocks.
function node_neighbors_from_mesh(mesh::ExodusDatabase, num_nodes::Int)
    neighbors = [Int[] for _ in 1:num_nodes]
    for block in Exodus.read_sets(mesh, Block)
        raw = get_block_connectivity(mesh, block.id)
        connectivity = reshape(Int.(vec(raw)), (size(raw, 2), size(raw, 1)))
        for e in 1:size(connectivity, 2), i in 1:size(connectivity, 1), j in 1:size(connectivity, 1)
            i == j && continue
            a, b = Int(connectivity[i, e]), Int(connectivity[j, e])
            b in neighbors[a] || push!(neighbors[a], b)
        end
    end
    return neighbors
end

# Sizes and rotation of the metric on one element: for the function sources
# at the point given (the centroid of the element in the mesh that samples the
# target), for the nodal sources the mean over the element's nodes.  Returns
# (h, R) with h positive and finite, or aborts with the location.
function principal_metric(source::PrincipalMetricFunctions, ::Any, point::SVector{3,Float64}, time::Float64)
    args = (time, point[1], point[2], point[3])
    h = SVector{3,Float64}(source.sizes[1](args), source.sizes[2](args), source.sizes[3](args))
    for i in 1:3
        if !isfinite(h[i]) || h[i] ≤ 0.0
            norma_abort(
                "Metric field sizes must be strictly positive and finite; got $(h[i]) for size $i at centroid " *
                "($(point[1]), $(point[2]), $(point[3])) and time $time",
            )
        end
    end
    if source.rotation === nothing
        return h, SMatrix{3,3,Float64,9}(I)
    end
    v = SVector{3,Float64}(source.rotation[1](args), source.rotation[2](args), source.rotation[3](args))
    return h, rt_of_rv(v)
end

function principal_metric(source::TensorMetricFunctions, ::Any, point::SVector{3,Float64}, time::Float64)
    args = (time, point[1], point[2], point[3])
    M = symmetric_from_components(SVector{6,Float64}(ntuple(i -> source.components[i](args), 6)))
    result = principal_of_tensor(M)
    result === nothing && norma_abort(
        "The metric tensor must be positive definite; it is not at centroid ($(point[1]), $(point[2]), " *
        "$(point[3])) and time $time",
    )
    return result
end

function nodal_indices_required(node_indices)
    node_indices === nothing && norma_abort("A nodal metric needs the node indices of the element")
    return node_indices
end

function principal_metric(source::NodalPrincipalMetric, node_indices, ::SVector{3,Float64}, ::Float64)
    nodes = nodal_indices_required(node_indices)
    weight = 1.0 / length(nodes)
    log_h = zero(SVector{3,Float64})
    v = zero(SVector{3,Float64})
    for n in nodes
        log_h += weight * SVector{3,Float64}(source.log_sizes[1, n], source.log_sizes[2, n], source.log_sizes[3, n])
        v += weight * SVector{3,Float64}(source.rotation[1, n], source.rotation[2, n], source.rotation[3, n])
    end
    return exp.(log_h), rt_of_rv(v)
end

function principal_metric(source::NodalTensorMetric, node_indices, ::SVector{3,Float64}, ::Float64)
    nodes = nodal_indices_required(node_indices)
    weight = 1.0 / length(nodes)
    L = zero(SVector{6,Float64})
    for n in nodes
        L += weight * SVector{6,Float64}(ntuple(i -> source.log_tensor[i, n], 6))
    end
    result = principal_of_tensor(exp_symmetric(symmetric_from_components(L)))
    result === nothing && norma_abort("The interpolated nodal metric tensor is not positive definite at nodes $nodes")
    return result
end

# Ideal reference element and metric factor for a TETRA4 element under a metric
# field.  The metric is evaluated once per element, at the centroid of the
# element in the mesh that samples the target for the function sources and
# from the element's nodes for the nodal sources, so the target is fixed
# during a solve and the assembled force is the exact gradient of the energy.
# With principal sizes h_i and rotation R, the metric factor is
# F_M = diag(1/h_i) R' (so F_M' F_M = M) and the ideal element is the unit
# regular tetrahedron Y mapped by F_M⁻¹ = R diag(h_i): scaled by h_i along the
# global axes, then rotated onto the principal directions.  Returns
# (X, F_M, F_M⁻¹).  The restricted rule scales the three sizes uniformly so
# that the ideal volume is never smaller than the volume of the original
# element, the anisotropic form of the `size field` floor: the smoother
# cannot create elements and should not ask one to shrink below what the mesh
# topology can accommodate.
function create_metric_reference(
    metric_field::MetricField, element_ref_pos::AbstractMatrix{Float64}, time::Float64; node_indices=nothing
)::Tuple{SMatrix{3,4,Float64,12},SMatrix{3,3,Float64,9},SMatrix{3,3,Float64,9}}
    X = SMatrix{3,4,Float64,12}(element_ref_pos)
    centroid = (X[:, 1] + X[:, 2] + X[:, 3] + X[:, 4]) / 4.0
    # The source is dispatched on dynamically; the declared result keeps the
    # caller typed.
    h, R = principal_metric(
        metric_field.source, node_indices, centroid, time
    )::Tuple{SVector{3,Float64},SMatrix{3,3,Float64,9}}
    if metric_field.restricted
        u = X[:, 2] - X[:, 1]
        v = X[:, 3] - X[:, 1]
        w = X[:, 4] - X[:, 1]
        element_volume = dot(u, cross(v, w)) / 6.0
        ideal_volume = h[1] * h[2] * h[3] / (6.0 * sqrt(2.0))
        if ideal_volume < element_volume
            h = h * cbrt(element_volume / ideal_volume)
        end
    end
    F_M = SMatrix{3,3,Float64,9}(Diagonal(SVector{3,Float64}(1.0 / h[1], 1.0 / h[2], 1.0 / h[3]))) * R'
    F_M_inv = R * SMatrix{3,3,Float64,9}(Diagonal(h))
    return F_M_inv * UNIT_TETRAHEDRON, F_M, F_M_inv
end

# Nodal metric data for the output: sizes (3 × n) and rotation vectors
# (3 × n, or nothing when the source has no rotation), at the given nodal
# positions and time for the function sources and from the carried values
# for the nodal sources.  Tensor sources are factored with frames followed
# through the mesh.
function metric_output_has_rotation(source::MetricSource)
    return !(source isa PrincipalMetricFunctions && source.rotation === nothing)
end

function nodal_metric_output(model::SolidMechanics, positions::AbstractMatrix{Float64}, time::Float64)
    source = model.metric_field.source
    num_nodes = size(positions, 2)
    if source isa PrincipalMetricFunctions
        sizes = zeros(3, num_nodes)
        rotation = source.rotation === nothing ? nothing : zeros(3, num_nodes)
        for n in 1:num_nodes
            args = (time, positions[1, n], positions[2, n], positions[3, n])
            for i in 1:3
                sizes[i, n] = source.sizes[i](args)
                rotation === nothing || (rotation[i, n] = source.rotation[i](args))
            end
        end
        return sizes, rotation
    elseif source isa TensorMetricFunctions
        tensors = Vector{SMatrix{3,3,Float64,9}}(undef, num_nodes)
        for n in 1:num_nodes
            args = (time, positions[1, n], positions[2, n], positions[3, n])
            tensors[n] = symmetric_from_components(SVector{6,Float64}(ntuple(i -> source.components[i](args), 6)))
        end
        sizes, rotation, _ = principal_of_nodal_tensors(tensors, node_neighbors(model))
        return sizes, rotation
    elseif source isa NodalPrincipalMetric
        return exp.(Matrix(source.log_sizes)), Matrix(source.rotation)
    end
    tensors = [exp_symmetric(symmetric_from_components(SVector{6,Float64}(ntuple(i -> source.log_tensor[i, n], 6))))
               for n in 1:num_nodes]
    sizes, rotation, _ = principal_of_nodal_tensors(tensors, node_neighbors(model))
    return sizes, rotation
end

# Names of the nodal variables of the two tensor forms of the metric at the
# nodes: the metric tensor M = R diag(1/h^2) R', whose quadratic form gives
# the squared length of a vector in the target, and the target tensor
# T = R diag(h) R', whose eigenvectors are the principal directions and whose
# eigenvalues are the target edge lengths, which is the form a tool that
# reads "scaled eigenvector" tensors expects.
const METRIC_TENSOR_NAMES = ["metric_$s" for s in METRIC_COMPONENT_SUFFIXES]
const TARGET_TENSOR_NAMES = ["target_$s" for s in METRIC_COMPONENT_SUFFIXES]

# The two tensor forms at every node, as 6 x n matrices in the component
# order of METRIC_COMPONENT_SUFFIXES.
function nodal_metric_tensors(model::SolidMechanics, positions::AbstractMatrix{Float64}, time::Float64)
    sizes, rotation = nodal_metric_output(model, positions, time)
    num_nodes = size(sizes, 2)
    metric = zeros(6, num_nodes)
    target = zeros(6, num_nodes)
    for n in 1:num_nodes
        h = SVector{3,Float64}(sizes[1, n], sizes[2, n], sizes[3, n])
        R = rotation === nothing ? SMatrix{3,3,Float64,9}(I) :
            rt_of_rv(SVector{3,Float64}(rotation[1, n], rotation[2, n], rotation[3, n]))
        metric[:, n] = components_from_symmetric(metric_tensor(h, R))
        target[:, n] = components_from_symmetric(SMatrix{3,3,Float64,9}(R * Diagonal(h) * R'))
    end
    return metric, target
end

# Nodal metric data carried through the adaptivity loop.  The data of the
# node that splits the edge (a, b) is the mean of the ends, in the space the
# source interpolates in; nothing for the function sources, which are
# sampled where the node lands.
function metric_node_data(field::Union{MetricField,Nothing}, a::Int, b::Int)
    field === nothing && return nothing
    source = field.source
    if source isa NodalPrincipalMetric
        log_sizes = 0.5 * (source.log_sizes[:, a] + source.log_sizes[:, b])
        rotation = 0.5 * (source.rotation[:, a] + source.rotation[:, b])
        return vcat(log_sizes, rotation)
    elseif source isa NodalTensorMetric
        return 0.5 * (source.log_tensor[:, a] + source.log_tensor[:, b])
    end
    return nothing
end

# A view of the field with one trial node appended, without copying the
# carried data (MatrixWithColumn), for the energies of a proposed split.
function metric_with_node(field::Union{MetricField,Nothing}, data::Union{Nothing,Vector{Float64}})
    (field === nothing || data === nothing) && return field
    source = field.source
    if source isa NodalPrincipalMetric
        log_sizes = MatrixWithColumn(source.log_sizes, SVector{3,Float64}(data[1], data[2], data[3]))
        rotation = MatrixWithColumn(source.rotation, SVector{3,Float64}(data[4], data[5], data[6]))
        return MetricField(
            NodalPrincipalMetric(log_sizes, rotation, source.size_names, source.rotation_names, source.tensor_names),
            field.restricted,
        )
    elseif source isa NodalTensorMetric
        log_tensor = MatrixWithColumn(source.log_tensor, SVector{6,Float64}(data))
        return MetricField(NodalTensorMetric(log_tensor, source.names), field.restricted)
    end
    return field
end

# Append the data of a node that was added to the mesh.
function add_metric_node!(field::Union{MetricField,Nothing}, data::Union{Nothing,Vector{Float64}})
    (field === nothing || data === nothing) && return field
    source = field.source
    if source isa NodalPrincipalMetric
        source.log_sizes = hcat(source.log_sizes, data[1:3])
        source.rotation = hcat(source.rotation, data[4:6])
    elseif source isa NodalTensorMetric
        source.log_tensor = hcat(source.log_tensor, data)
    end
    return field
end

# Keep the data of the alive nodes, in their compacted order.
function compact_metric!(field::Union{MetricField,Nothing}, alive_nodes::AbstractVector{Int})
    field === nothing && return field
    source = field.source
    if source isa NodalPrincipalMetric
        source.log_sizes = source.log_sizes[:, alive_nodes]
        source.rotation = source.rotation[:, alive_nodes]
    elseif source isa NodalTensorMetric
        source.log_tensor = source.log_tensor[:, alive_nodes]
    end
    return field
end

# Nodal variables that carry the metric on an adapted mesh, under the names
# the input named them by, so the input file applies unchanged to the
# adapted mesh; empty for the function sources.
function metric_nodal_variables(field::Union{MetricField,Nothing})
    variables = Dict{String,Vector{Float64}}()
    field === nothing && return variables
    source = field.source
    if source isa NodalPrincipalMetric
        num_nodes = size(source.log_sizes, 2)
        if isempty(source.tensor_names)
            for i in 1:3
                variables[source.size_names[i]] = exp.(source.log_sizes[i, :])
            end
            for i in 1:length(source.rotation_names)
                variables[source.rotation_names[i]] = source.rotation[i, :]
            end
        else
            components = zeros(6, num_nodes)
            for n in 1:num_nodes
                h = exp.(SVector{3,Float64}(source.log_sizes[1, n], source.log_sizes[2, n], source.log_sizes[3, n]))
                R = rt_of_rv(SVector{3,Float64}(source.rotation[1, n], source.rotation[2, n], source.rotation[3, n]))
                components[:, n] = components_from_symmetric(metric_tensor(h, R))
            end
            for i in 1:6
                variables[source.tensor_names[i]] = components[i, :]
            end
        end
    elseif source isa NodalTensorMetric
        num_nodes = size(source.log_tensor, 2)
        components = zeros(6, num_nodes)
        for n in 1:num_nodes
            L = symmetric_from_components(SVector{6,Float64}(ntuple(i -> source.log_tensor[i, n], 6)))
            components[:, n] = components_from_symmetric(exp_symmetric(L))
        end
        for i in 1:6
            variables[source.names[i]] = components[i, :]
        end
    end
    return variables
end
