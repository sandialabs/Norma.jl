# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

include("constitutive.jl")
include("interpolation.jl")
include("ics_bcs.jl")

using Base.Threads: @threads, threadid, nthreads
using NPZ

function LinearOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return LinearOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end

function QuadraticOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return QuadraticOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end

function CubicOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return CubicOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end 

function SolidMechanics(params::Parameters)
    input_mesh = params["input_mesh"]
    model_params = params["model"]
    coords = read_coordinates(input_mesh)
    num_nodes = Exodus.num_nodes(input_mesh.init)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    current = Matrix{Float64}(undef, 3, num_nodes)
    velocity = Matrix{Float64}(undef, 3, num_nodes)
    acceleration = Matrix{Float64}(undef, 3, num_nodes)
    for node in 1:num_nodes
        reference[:, node] = coords[:, node]
        current[:, node] = coords[:, node]
        velocity[:, node] = [0.0, 0.0, 0.0]
        acceleration[:, node] = [0.0, 0.0, 0.0]
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    if num_blks_params ≠ num_blks
        norma_abortf(
            "Number of blocks in mesh %s (%d) must be equal to number of blocks in materials %s (%d).",
            model_params["mesh"],
            num_blks,
            model_params["material"],
            num_blks_params,
        )
    end
    elem_block_names = Exodus.read_names(input_mesh, Block)
    materials = Vector{Solid}(undef, 0)
    kinematics = Undefined
    for elem_block_name in elem_block_names
        material_name = material_blocks[elem_block_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        if kinematics == Undefined
            kinematics = get_kinematics(material_model)
        else
            if kinematics ≠ get_kinematics(material_model)
                norma_abortf(
                    "Material of type %s has inconsistent kinematics %s compared to previous materials of type %s.",
                    string(typeof(material_model)),
                    string(get_kinematics(material_model)),
                    string(kinematics),
                )
            end
        end
        push!(materials, material_model)
    end
    time = 0.0
    failed = false
    internal_force = zeros(3 * num_nodes)
    boundary_force = zeros(3 * num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(3 * num_nodes)
    stress = Vector{Vector{Vector{Vector{Float64}}}}()
    stored_energy = Vector{Vector{Float64}}()
    for block in blocks
        block_id = block.id
        element_type_string, num_block_elems, _, _, _, _ = Exodus.read_block_parameters(input_mesh, block_id)
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        block_stress = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ in 1:num_block_elems
            element_stress = Vector{Vector{Float64}}()
            for _ in 1:num_points
                push!(element_stress, zeros(6))
            end
            push!(block_stress, element_stress)
            element_stored_energy = 0.0
            push!(block_stored_energy, element_stored_energy)
        end
        push!(stress, block_stress)
        push!(stored_energy, block_stored_energy)
    end
    strain_energy = 0.0
    stiffness = spzeros(0, 0)
    diag_stiffness = Float64[]
    mass = spzeros(0, 0)
    lumped_mass = Float64[]
    body_force = Float64[]
    compute_stiffness = true
    compute_diag_stiffness = true
    compute_mass = true
    compute_lumped_mass = true
    mesh_smoothing = get(params, "mesh smoothing", false)
    smooth_reference = get(model_params, "smooth reference", "")
    inclined_support = false
    num_dofs = 3 * num_nodes
    global_transform = sparse(I, num_dofs, num_dofs)

    return SolidMechanics(
        input_mesh,
        materials,
        reference,
        current,
        velocity,
        acceleration,
        internal_force,
        boundary_force,
        boundary_conditions,
        stress,
        stored_energy,
        strain_energy,
        stiffness,
        diag_stiffness,
        mass,
        lumped_mass,
        body_force,
        free_dofs,
        time,
        compute_stiffness,
        compute_diag_stiffness,
        compute_mass,
        compute_lumped_mass,
        failed,
        mesh_smoothing,
        smooth_reference,
        inclined_support,
        global_transform,
        kinematics,
    )
end

function create_model(params::Parameters)
    model_params = params["model"]
    model_name = model_params["type"]
    if model_name == "solid mechanics"
        return SolidMechanics(params)
    elseif model_name == "mesh smoothing"
        params["mesh smoothing"] = true
        return SolidMechanics(params)
    elseif model_name == "linear opinf rom"
        return LinearOpInfRom(params)
    elseif model_name == "quadratic opinf rom"
        return QuadraticOpInfRom(params)
    elseif model_name == "cubic opinf rom"
        return CubicOpInfRom(params)

    else
        norma_abort("Unknown type of model : $model_name")
    end
end

function create_smooth_reference(smooth_reference::String, element_type::ElementType, elem_ref_pos::Matrix{Float64})
    if element_type == TETRA4
        u = elem_ref_pos[:, 2] - elem_ref_pos[:, 1]
        v = elem_ref_pos[:, 3] - elem_ref_pos[:, 1]
        w = elem_ref_pos[:, 4] - elem_ref_pos[:, 1]

        if smooth_reference == "equal volume"
            h = equal_volume_tet_h(u, v, w)
        elseif smooth_reference == "average edge length"
            h = avg_edge_length_tet_h(u, v, w)
        elseif smooth_reference == "max"
            h = max(equal_volume_tet_h(u, v, w), avg_edge_length_tet_h(u, v, w))
        else
            norma_abort("Unknown type of mesh smoothing reference : $smooth_reference")
        end

        c = h * 0.5 / sqrt(2.0)
        A = [
            1 -1 -1 1
            1 -1 1 -1
            1 1 -1 -1
        ]
        return c * A
    else
        norma_abort("Unknown element type")
    end
end

function equal_volume_tet_h(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    h = cbrt(sqrt(2.0) * dot(u, cross(v, w)))
    return h
end

function avg_edge_length_tet_h(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    h = (norm(u) + norm(v) + norm(w) + norm(u - v) + norm(u - w) + norm(v - w)) / 6.0
    return h
end

function characteristic_element_length_centroid(nodal_coordinates::Matrix{Float64})::Float64
    centroid = sum(nodal_coordinates, dims=2) / size(nodal_coordinates, 2)
    total = 0.0
    @inbounds for i in 1:size(nodal_coordinates, 2)
        δ = nodal_coordinates[:, i] - centroid
        total += norm(δ)
    end
    return 2 * total / size(nodal_coordinates, 2)  # Approximate diameter
end

function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    stable_time_step = Inf
    for block_index in 1:num_blks
        material = materials[block_index]
        ρ = material.ρ
        M = get_p_wave_modulus(material)
        wave_speed = sqrt(M / ρ)
        minimum_block_edge_length = Inf
        block = blocks[block_index]
        block_id = block.id
        elem_block_conn = get_block_connectivity(input_mesh, block_id)
        num_block_elems, num_elem_nodes = size(elem_block_conn)
        for block_elem_index in 1:num_block_elems
            conn_indices = ((block_elem_index - 1) * num_elem_nodes + 1):(block_elem_index * num_elem_nodes)
            node_indices = elem_block_conn[conn_indices]
            elem_cur_pos = model.current[:, node_indices]
            minimum_elem_edge_length = characteristic_element_length_centroid(elem_cur_pos)
            minimum_block_edge_length = min(minimum_block_edge_length, minimum_elem_edge_length)
        end
        block_stable_time_step = integrator.CFL * minimum_block_edge_length / wave_speed
        stable_time_step = min(stable_time_step, block_stable_time_step)
    end
    if stable_time_step < integrator.time_step
        norma_logf(
            0,
            :warning,
            "Δt = %.3e exceeds stable Δt = %.3e — using stable step.",
            integrator.time_step,
            stable_time_step,
        )
    end
    integrator.time_step = min(stable_time_step, integrator.time_step)
    return nothing
end

function voigt_cauchy_from_stress(_::Solid, P::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9}, J::Float64)
    σ = F * P' ./ J
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function voigt_cauchy_from_stress(_::Linear_Elastic, σ::SMatrix{3,3,Float64,9}, _::SMatrix{3,3,Float64,9}, _::Float64)
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function dense(indices::Vector{Int64}, values::Vector{Float64}, vector_size::Int64)
    dense_vector = zeros(vector_size)
    @inbounds for i in 1:length(indices)
        dense_vector[indices[i]] += values[i]
    end
    return dense_vector
end

@generated function create_element_matrix(::Type{T}, ::Val{N}) where {T,N}
    dof_per_node = 3
    total_dofs = dof_per_node * N
    quote
        MMatrix{$total_dofs,$total_dofs,$T}(undef)
    end
end

@generated function create_element_vector(::Type{T}, ::Val{N}) where {T,N}
    dof_per_node = 3
    total_dofs = dof_per_node * N
    quote
        MVector{$total_dofs,$T}(undef)
    end
end

function create_gradient_operator(dNdX::SMatrix{3,N,T})::SMatrix{9,3N,T} where {N,T}
    B = MMatrix{9,3N,T}(undef)
    fill!(B, zero(T))
    @inbounds for i in 1:3         # i = direction of derivative
        for a in 1:N              # a = local node index
            # Place dNdX[:, a] into the appropriate 3×1 column
            B[(3 * (i - 1) + 1):(3 * i), (3 * (a - 1) + i)] = dNdX[:, a]
        end
    end
    return SMatrix{9,3N,T}(B)
end

function create_coo_vector(capacity::Int64)
    index = Vector{Int64}(undef, capacity)
    vals = Vector{Float64}(undef, capacity)
    return COOVector(index, vals, 0)
end

function create_coo_matrix(capacity::Int64)
    rows = Vector{Int64}(undef, capacity)
    cols = Vector{Int64}(undef, capacity)
    vals = Vector{Float64}(undef, capacity)
    return COOMatrix(rows, cols, vals, 0)
end

function ensure_capacity!(vector::COOVector, needed::Int64)
    current = length(vector.index)
    required = vector.len + needed
    if required > current
        newcap = max(required, ceil(Int64, 1.5 * current))
        resize!(vector.index, newcap)
        resize!(vector.vals, newcap)
    end
    return nothing
end

function ensure_capacity!(matrix::COOMatrix, needed::Int64)
    current = length(matrix.rows)
    required = matrix.len + needed
    if required > current
        newcap = max(required, ceil(Int64, 1.5 * current))
        resize!(matrix.rows, newcap)
        resize!(matrix.cols, newcap)
        resize!(matrix.vals, newcap)
    end
    return nothing
end

function assemble!(global_vector::COOVector, element_vector::AbstractVector{Float64}, dofs::AbstractVector{Int64})
    ndofs = length(dofs)
    ensure_capacity!(global_vector, ndofs)
    idx = global_vector.len + 1
    @inbounds for i in 1:ndofs
        global_vector.index[idx] = dofs[i]
        global_vector.vals[idx] = element_vector[i]
        idx += 1
    end
    global_vector.len += ndofs
    return nothing
end

function assemble!(global_matrix::COOMatrix, element_matrix::AbstractMatrix{Float64}, dofs::AbstractVector{Int64})
    ndofs = length(dofs)
    n2 = ndofs * ndofs
    ensure_capacity!(global_matrix, n2)
    idx = global_matrix.len + 1
    @inbounds for i in 1:ndofs
        I = dofs[i]
        @inbounds for j in 1:ndofs
            global_matrix.rows[idx] = I
            global_matrix.cols[idx] = dofs[j]
            global_matrix.vals[idx] = element_matrix[i, j]
            idx += 1
        end
    end
    global_matrix.len += n2
    return nothing
end

function count_coo_matrix_nnz(model::SolidMechanics)
    mesh = model.mesh
    blocks = Exodus.read_sets(mesh, Block)
    num_blocks = length(blocks)
    total = 0
    for block_index in 1:num_blocks
        block = blocks[block_index]
        block_id = block.id
        element_block_conn = get_block_connectivity(mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_conn)
        total += num_block_elements * num_element_nodes * num_element_nodes * 9
    end
    return total
end

function merge_threadlocal_coo_vectors(coo_vectors::Vector{COOVector}, num_dof::Int64)
    # Trimmed slices
    index = vcat((v.index[1:(v.len)] for v in coo_vectors)...)
    vals = vcat((v.vals[1:(v.len)] for v in coo_vectors)...)
    return dense(index, vals, num_dof)
end

function merge_threadlocal_coo_matrices(coo_matrices::Vector{COOMatrix}, num_dof::Int64)
    # Trimmed slices
    rows = vcat((m.rows[1:(m.len)] for m in coo_matrices)...)
    cols = vcat((m.cols[1:(m.len)] for m in coo_matrices)...)
    vals = vcat((m.vals[1:(m.len)] for m in coo_matrices)...)
    return sparse(rows, cols, vals, num_dof, num_dof)
end

function create_threadlocal_element_matrices(::Type{T}, ::Val{N}) where {T,N}
    local_mats = Vector{typeof(create_element_matrix(T, Val(N)))}(undef, nthreads())
    for i in 1:nthreads()
        local_mats[i] = create_element_matrix(T, Val(N))
    end
    return local_mats
end

function create_threadlocal_element_vectors(::Type{T}, ::Val{N}) where {T,N}
    local_vecs = Vector{typeof(create_element_vector(T, Val(N)))}(undef, nthreads())
    for i in 1:nthreads()
        local_vecs[i] = create_element_vector(T, Val(N))
    end
    return local_vecs
end

function create_threadlocal_coo_vectors(num_dofs::Int64)
    nthreads = Threads.nthreads()
    capacity = ceil(Int64, num_dofs / nthreads)
    coo_vectors = Vector{COOVector}(undef, nthreads)
    for i in 1:nthreads
        coo_vectors[i] = create_coo_vector(capacity)
    end
    return coo_vectors
end

function create_threadlocal_coo_matrices(coo_matrix_nnz::Int64)
    nthreads = Threads.nthreads()
    capacity = ceil(Int64, coo_matrix_nnz / nthreads)
    coo_matrices = Vector{COOMatrix}(undef, nthreads)
    for i in 1:nthreads
        coo_matrices[i] = create_coo_matrix(capacity)
    end
    return coo_matrices
end

function add_internal_force!(Fi::MVector{M,T}, grad_op::SMatrix{9,M,T}, stress::SVector{9,T}, dV::T) where {M,T}
    @einsum Fi[i] += grad_op[j, i] * stress[j] * dV
end

function add_diag_stiff!(K::MVector{M,T}, grad_op::SMatrix{9,M,T}, C::SMatrix{9,9,T}, dV::T) where {M,T}
    @inbounds for i in 1:M
        g = grad_op[:, i]
        K[i] += dV * dot(g, C * g)
    end
    return nothing
end

function add_lumped_mass!(M::MVector{R,T}, Nξ::SVector{N,T}, density::T, dV::T) where {R,N,T}
    @assert R == 3N
    s = sum(Nξ)
    w = (density * dV) .* Nξ .* s

    @inbounds for a in 1:N
        idx = 3 * (a - 1) + 1
        m = w[a]
        M[idx] += m
        M[idx + 1] += m
        M[idx + 2] += m
    end
    return nothing
end

function compute_flags(model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    is_implicit_dynamic = integrator isa Newmark
    is_explicit_dynamic = integrator isa CentralDifference
    is_implicit_static = integrator isa QuasiStatic
    is_dynamic = is_implicit_dynamic || is_explicit_dynamic
    is_implicit = is_implicit_dynamic || is_implicit_static
    is_hessian_opt = solver isa HessianMinimizer
    is_matrix_free = solver isa SteepestDescent
    need_diag_stiffness = is_implicit && is_matrix_free
    need_lumped_mass = is_explicit_dynamic || (is_implicit_dynamic && is_matrix_free)
    need_stiffness = is_implicit && is_hessian_opt
    need_mass = is_dynamic && is_hessian_opt
    compute_diag_stiffness = need_diag_stiffness && model.compute_diag_stiffness
    compute_lumped_mass = need_lumped_mass && model.compute_lumped_mass
    compute_stiffness = need_stiffness && model.compute_stiffness
    compute_mass = need_mass && model.compute_mass
    mesh_smoothing = model.mesh_smoothing

    return EvaluationFlags(
        is_dynamic,
        is_implicit,
        is_hessian_opt,
        is_matrix_free,
        need_diag_stiffness,
        need_lumped_mass,
        need_stiffness,
        need_mass,
        compute_diag_stiffness,
        compute_lumped_mass,
        compute_stiffness,
        compute_mass,
        mesh_smoothing,
    )
end

function create_threadlocal_arrays(model::SolidMechanics, flags::EvaluationFlags)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    energy = zeros(nthreads())
    internal_force = create_threadlocal_coo_vectors(num_dofs)

    diag_stiffness = create_threadlocal_coo_vectors(flags.compute_diag_stiffness ? num_dofs : 0)
    if flags.compute_diag_stiffness && model.kinematics == Infinitesimal
        model.compute_diag_stiffness = false
    end

    lumped_mass = create_threadlocal_coo_vectors(flags.compute_lumped_mass ? num_dofs : 0)
    if flags.compute_lumped_mass
        model.compute_lumped_mass = false
    end

    if flags.compute_stiffness || flags.compute_mass
        coo_matrix_nnz = count_coo_matrix_nnz(model)
    else
        coo_matrix_nnz = 0
    end

    stiffness = create_threadlocal_coo_matrices(flags.compute_stiffness ? coo_matrix_nnz : 0)
    if flags.compute_stiffness && model.kinematics == Infinitesimal
        model.compute_stiffness = false
    end

    mass = create_threadlocal_coo_matrices(flags.compute_mass ? coo_matrix_nnz : 0)
    if flags.compute_mass
        model.compute_mass = false
    end
    return SMThreadLocalArrays(energy, internal_force, diag_stiffness, lumped_mass, stiffness, mass)
end

function create_element_threadlocal_arrays(num_element_nodes::Int64, flags::EvaluationFlags)
    valN = Val(num_element_nodes)
    energy = zeros(nthreads())
    dofs = create_threadlocal_element_vectors(Int64, valN)
    internal_force = create_threadlocal_element_vectors(Float64, valN)
    diag_stiffness = create_threadlocal_element_vectors(Float64, flags.compute_diag_stiffness ? valN : Val(0))
    lumped_mass = create_threadlocal_element_vectors(Float64, flags.compute_lumped_mass ? valN : Val(0))
    stiffness = create_threadlocal_element_matrices(Float64, flags.compute_stiffness ? valN : Val(0))
    mass = create_threadlocal_element_matrices(Float64, flags.compute_mass ? valN : Val(0))
    return SMElementThreadLocalArrays(energy, dofs, internal_force, diag_stiffness, lumped_mass, stiffness, mass)
end

function reset_element_threadlocal_arrays!(
    element_arrays_tl::SMElementThreadLocalArrays,
    element_block_conn::Matrix{<:Integer},
    block_element_index::Integer,
    flags::EvaluationFlags,
)
    t = threadid()
    num_element_nodes = size(element_block_conn, 2)
    conn_indices = ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
    node_indices = element_block_conn[conn_indices]
    element_arrays_tl.dofs[t] = reshape(3 .* node_indices' .- [2, 1, 0], :)
    element_arrays_tl.energy[t] = 0.0
    fill!(element_arrays_tl.internal_force[t], 0.0)
    if flags.compute_diag_stiffness == true
        fill!(element_arrays_tl.diag_stiffness[t], 0.0)
    end
    if flags.compute_lumped_mass == true
        fill!(element_arrays_tl.lumped_mass[t], 0.0)
    end
    if flags.compute_stiffness == true
        fill!(element_arrays_tl.stiffness[t], 0.0)
    end
    if flags.compute_mass == true
        fill!(element_arrays_tl.mass[t], 0.0)
    end
    return node_indices
end

function compute_element_threadlocal_arrays!(
    element_arrays_tl::SMElementThreadLocalArrays,
    Np::SVector{N,T},
    dNdX::SMatrix{3,N,T},
    W::T,
    P::SMatrix{3,3,T,9},
    AA::SArray{Tuple{3,3,3,3},T},
    density::T,
    dvol::T,
    flags::EvaluationFlags,
) where {T,N}
    t = threadid()
    grad_op = create_gradient_operator(dNdX)
    stress = SVector{9,Float64}(P)
    element_arrays_tl.energy[t] += W * dvol
    add_internal_force!(element_arrays_tl.internal_force[t], grad_op, stress, dvol)
    if flags.compute_diag_stiffness == true
        moduli = second_from_fourth(AA)
        add_diag_stiff!(element_arrays_tl.diag_stiffness[t], grad_op, moduli, dvol)
    end
    if flags.compute_lumped_mass == true
        add_lumped_mass!(element_arrays_tl.lumped_mass[t], Np, density, dvol)
    end
    if flags.compute_stiffness == true
        moduli = second_from_fourth(AA)
        element_arrays_tl.stiffness[t] += grad_op' * moduli * grad_op * dvol
    end
    if flags.compute_mass == true
        reduced_mass = Np * Np' * density * dvol
        mass = element_arrays_tl.mass[t]
        for i in 1:3
            mass[i:3:end, i:3:end] .+= reduced_mass
        end
    end
    return nothing
end

function assemble_element_threadlocal_arrays!(
    arrays_tl::SMThreadLocalArrays, element_arrays_tl::SMElementThreadLocalArrays, flags::EvaluationFlags
)
    t = threadid()
    arrays_tl.energy[t] += element_arrays_tl.energy[t]
    assemble!(arrays_tl.internal_force[t], element_arrays_tl.internal_force[t], element_arrays_tl.dofs[t])
    if flags.compute_diag_stiffness == true
        assemble!(arrays_tl.diag_stiffness[t], element_arrays_tl.diag_stiffness[t], element_arrays_tl.dofs[t])
    end
    if flags.compute_lumped_mass == true
        assemble!(arrays_tl.lumped_mass[t], element_arrays_tl.lumped_mass[t], element_arrays_tl.dofs[t])
    end
    if flags.compute_stiffness == true
        assemble!(arrays_tl.stiffness[t], element_arrays_tl.stiffness[t], element_arrays_tl.dofs[t])
    end
    if flags.compute_mass == true
        assemble!(arrays_tl.mass[t], element_arrays_tl.mass[t], element_arrays_tl.dofs[t])
    end
    return nothing
end

function merge_threadlocal_arrays(
    model::SolidMechanics, arrays_tl::SMThreadLocalArrays, num_dofs::Int64, flags::EvaluationFlags
)
    model.strain_energy = sum(arrays_tl.energy)
    model.internal_force = merge_threadlocal_coo_vectors(arrays_tl.internal_force, num_dofs)
    if flags.compute_diag_stiffness == true
        model.diag_stiffness = merge_threadlocal_coo_vectors(arrays_tl.diag_stiffness, num_dofs)
    end
    if flags.compute_lumped_mass == true
        model.lumped_mass = merge_threadlocal_coo_vectors(arrays_tl.lumped_mass, num_dofs)
    end
    if flags.compute_stiffness == true
        model.stiffness = merge_threadlocal_coo_matrices(arrays_tl.stiffness, num_dofs)
    end
    if flags.compute_mass == true
        model.mass = merge_threadlocal_coo_matrices(arrays_tl.mass, num_dofs)
    end
    return nothing
end

function evaluate(model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    flags = compute_flags(model, integrator, solver)
    arrays_tl = create_threadlocal_arrays(model, flags)
    materials = model.materials
    input_mesh = model.mesh
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    body_force_vector = zeros(num_dofs)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blocks = length(blocks)
    for block_index in 1:num_blocks
        material = materials[block_index]
        density = material.ρ
        block = blocks[block_index]
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_conn = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_conn)
        element_arrays_tl = create_element_threadlocal_arrays(num_element_nodes, flags)
        @threads for block_element_index in 1:num_block_elements
            node_indices = reset_element_threadlocal_arrays!(
                element_arrays_tl, element_block_conn, block_element_index, flags
            )
            if flags.mesh_smoothing == true
                element_reference_position = create_smooth_reference(
                    model.smooth_reference, element_type, model.reference[:, node_indices]
                )
            else
                element_reference_position = model.reference[:, node_indices]
            end
            element_current_position = model.current[:, node_indices]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dNdX = dXdξ \ dNdξ
                F = SMatrix{3,3,Float64,9}(dNdX * element_current_position')
                J = det(F)
                if J ≤ 0.0 || isfinite(J) == false
                    model.failed = true
                    model.compute_stiffness = model.compute_diag_stiffness = true
                    model.compute_mass = model.compute_lumped_mass = true
                    norma_log(0, :error, "Non-positive Jacobian detected!")
                    norma_log(0, :error, "This may indicate element distortion.")
                    return nothing
                end
                W, P, AA = constitutive(material, F)
                ip_weight = ip_weights[point]
                det_dXdξ = det(dXdξ)
                dvol = det_dXdξ * ip_weight
                compute_element_threadlocal_arrays!(element_arrays_tl, Np, dNdX, W, P, AA, density, dvol, flags)
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[block_index][block_element_index][point] = voigt_cauchy
            end
            t = threadid()
            model.stored_energy[block_index][block_element_index] = element_arrays_tl.energy[t]
            assemble_element_threadlocal_arrays!(arrays_tl, element_arrays_tl, flags)
        end
    end
    merge_threadlocal_arrays(model, arrays_tl, num_dofs, flags)
    model.body_force = body_force_vector
    if flags.mesh_smoothing == true
        model.internal_force -= integrator.velocity
    end
    return nothing
end

function get_block_connectivity(mesh::ExodusDatabase, block_id::Integer)
    _, num_elems, num_nodes, _, _, _ = Exodus.read_block_parameters(mesh, Int32(block_id))
    conn = Exodus.read_block_connectivity(mesh, Int32(block_id), num_elems * num_nodes)
    return reshape(conn, (num_elems, num_nodes))
end
