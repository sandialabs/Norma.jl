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
        error(
            "number of blocks in mesh ",
            model_params["mesh"],
            " (",
            num_blks,
            ") must be equal to number of blocks in materials file ",
            model_params["material"],
            " (",
            num_blks_params,
            ")",
        )
    end
    elem_blk_names = Exodus.read_names(input_mesh, Block)
    materials = Vector{Solid}(undef, 0)
    kinematics = Undefined
    for elem_blk_name in elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        if kinematics == Undefined
            kinematics = get_kinematics(material_model)
        else
            if kinematics ≠ get_kinematics(material_model)
                error(
                    "Material ",
                    typeof(material_model),
                    " has inconsistent kinematics ",
                    get_kinematics(material_model),
                    " than previous materials of type ",
                    kinematics,
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
        blk_id = block.id
        element_type_string, num_blk_elems, _, _, _, _ = Exodus.read_block_parameters(input_mesh, blk_id)
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        block_stress = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ in 1:num_blk_elems
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
    mesh_smoothing = params["mesh smoothing"]
    if mesh_smoothing == true
        smooth_reference = model_params["smooth reference"]
    else
        smooth_reference = ""
    end
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

function HeatConduction(params::Parameters)
    input_mesh = params["input_mesh"]
    model_params = params["model"]
    coords = read_coordinates(input_mesh)
    num_nodes = Exodus.num_nodes(input_mesh.init)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    temperature = Vector{Float64}(undef, num_nodes)
    rate = Vector{Float64}(undef, num_nodes)
    for node in 1:num_nodes
        reference[:, node] = coords[:, node]
        temperature[node] = 0.0
        rate[node] = 0.0
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    if num_blks_params ≠ num_blks
        error(
            "number of blocks in mesh ",
            model_params["mesh"],
            " (",
            num_blks,
            ") must be equal to number of blocks in materials file ",
            model_params["material"],
            " (",
            num_blks_params,
            ")",
        )
    end
    elem_blk_names = Exodus.read_names(input_mesh, Block)
    materials = Vector{Thermal}(undef, 0)
    for elem_blk_name in elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    failed = false
    internal_heat_flux = zeros(num_nodes)
    boundary_heat_flux = zeros(num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(num_nodes)
    flux = Vector{Vector{Vector{Vector{Float64}}}}()
    stored_energy = Vector{Vector{Float64}}()
    for block in blocks
        blk_id = block.id
        element_type_string, num_blk_elems, _, _, _, _ = Exodus.read_block_parameters(input_mesh, blk_id)
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        block_flux = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ in 1:num_blk_elems
            element_flux = Vector{Vector{Float64}}()
            for _ in 1:num_points
                push!(element_flux, zeros(3))
            end
            push!(block_flux, element_flux)
            element_stored_energy = 0.0
            push!(block_stored_energy, element_stored_energy)
        end
        push!(flux, block_flux)
        push!(stored_energy, block_stored_energy)
    end
    return HeatConduction(
        input_mesh,
        materials,
        reference,
        temperature,
        rate,
        internal_heat_flux,
        boundary_heat_flux,
        boundary_conditions,
        flux,
        stored_energy,
        free_dofs,
        time,
        failed,
    )
end

function create_model(params::Parameters)
    model_params = params["model"]
    model_name = model_params["type"]
    if model_name == "solid mechanics"
        params["mesh smoothing"] = false
        return SolidMechanics(params)
    elseif model_name == "mesh smoothing"
        params["mesh smoothing"] = true
        return SolidMechanics(params)
    elseif model_name == "heat conduction"
        return HeatConduction(params)
    elseif model_name == "linear opinf rom"
        return LinearOpInfRom(params)
    elseif model_name == "quadratic opinf rom"
        return QuadraticOpInfRom(params)

    else
        error("Unknown type of model : ", model_name)
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
            error("Unknown type of mesh smoothing reference : ", smooth_reference)
        end

        c = h * 0.5 / sqrt(2.0)
        A = [
            1 -1 -1 1
            1 -1 1 -1
            1 1 -1 -1
        ]
        return c * A
    else
        error("Unknown element type")
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

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64}, edges::Vector{Tuple{Int64,Int64}})
    minimum_edge_length = Inf
    for edge in edges
        node_a = edge[1]
        node_b = edge[2]
        edge_vector = nodal_coordinates[:, node_a] - nodal_coordinates[:, node_b]
        distance = norm(edge_vector)
        minimum_edge_length = min(minimum_edge_length, distance)
    end
    return minimum_edge_length
end

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64}, element_type::ElementType)
    if element_type == TETRA4
        edges = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4), (2, 4)]
        return get_minimum_edge_length(nodal_coordinates, edges)
    elseif element_type == HEX8
        edges = [(1, 4), (1, 5), (4, 8), (5, 8), (2, 3), (2, 6), (3, 7), (6, 7), (1, 2), (3, 4), (5, 6), (7, 8)]
        return get_minimum_edge_length(nodal_coordinates, edges)
    else
        error("Invalid element type: ", element_type)
    end
end

function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    stable_time_step = Inf
    for blk_index in 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        M = get_p_wave_modulus(material)
        wave_speed = sqrt(M / ρ)
        minimum_blk_edge_length = Inf
        block = blocks[blk_index]
        blk_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, blk_id)[1]
        element_type = element_type_from_string(element_type_string)
        elem_blk_conn = get_block_connectivity(input_mesh, blk_id)
        num_blk_elems, num_elem_nodes = size(elem_blk_conn)
        for blk_elem_index in 1:num_blk_elems
            conn_indices = ((blk_elem_index - 1) * num_elem_nodes + 1):(blk_elem_index * num_elem_nodes)
            node_indices = elem_blk_conn[conn_indices]
            elem_cur_pos = model.current[:, node_indices]
            minimum_elem_edge_length = get_minimum_edge_length(elem_cur_pos, element_type)
            minimum_blk_edge_length = min(minimum_blk_edge_length, minimum_elem_edge_length)
        end
        blk_stable_time_step = integrator.CFL * minimum_blk_edge_length / wave_speed
        stable_time_step = min(stable_time_step, blk_stable_time_step)
    end
    integrator.stable_time_step = stable_time_step
    if stable_time_step < integrator.user_time_step
        @printf(
            "❗ Δt = %.3e exceeds stable Δt = %.3e — using stable step.\n", integrator.user_time_step, stable_time_step
        )
    end
    integrator.time_step = min(stable_time_step, integrator.user_time_step)
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

@generated function create_reduced_element_matrix(::Type{T}, ::Val{N}) where {T,N}
    quote
        MMatrix{$N,$N,$T}(undef)
    end
end

@generated function create_element_vector(::Type{T}, ::Val{N}) where {T,N}
    dof_per_node = 3
    total_dofs = dof_per_node * N
    quote
        MVector{$total_dofs,$T}(undef)
    end
end

@generated function create_reduced_element_vector(::Type{T}, ::Val{N}) where {T,N}
    quote
        MVector{$N,$T}(undef)
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

function row_sum_lump(A::SMatrix{N,N,T}) where {N,T}
    return SVector{N,T}(sum(A; dims=2)[:, 1])
end

function evaluate(model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    is_implicit_dynamic = integrator isa Newmark
    is_explicit_dynamic = integrator isa CentralDifference
    is_implicit_static = integrator isa QuasiStatic
    is_dynamic = is_implicit_dynamic || is_explicit_dynamic
    is_implicit = is_implicit_dynamic || is_implicit_static
    is_hessian_opt = solver isa HessianMinimizer
    is_matrix_free = solver isa SteepestDescent
    need_diag_stiffness = is_implicit == true && is_matrix_free == true
    need_lumped_mass = is_explicit_dynamic == true || (is_implicit_dynamic == true && is_matrix_free == true)
    need_stiffness = is_implicit == true && is_hessian_opt == true
    need_mass = is_dynamic == true && is_hessian_opt == true
    compute_diag_stiffness = need_diag_stiffness == true && model.compute_diag_stiffness == true
    compute_lumped_mass = need_lumped_mass == true && model.compute_lumped_mass == true
    compute_stiffness = need_stiffness == true && model.compute_stiffness == true
    compute_mass = need_mass == true && model.compute_mass == true
    materials = model.materials
    input_mesh = model.mesh
    mesh_smoothing = model.mesh_smoothing
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    energy_tl = zeros(nthreads())
    internal_force_tl = create_threadlocal_coo_vectors(num_dofs)
    if compute_diag_stiffness == true
        diag_stiffness_tl = create_threadlocal_coo_vectors(num_dofs)
        if model.kinematics == Infinitesimal
            model.compute_diag_stiffness = false
        end
    end
    if compute_lumped_mass == true
        lumped_mass_tl = create_threadlocal_coo_vectors(num_dofs)
        model.compute_lumped_mass = false
    end
    if compute_stiffness == true || compute_mass == true
        coo_matrix_nnz = count_coo_matrix_nnz(model)
    end
    if compute_stiffness == true
        stiffness_tl = create_threadlocal_coo_matrices(coo_matrix_nnz)
        if model.kinematics == Infinitesimal
            model.compute_stiffness = false
        end
    end
    if compute_mass == true
        mass_tl = create_threadlocal_coo_matrices(coo_matrix_nnz)
        model.compute_mass = false
    end
    body_force_vector = zeros(num_dofs)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blocks = length(blocks)
    for block_index in 1:num_blocks
        material = materials[block_index]
        if is_dynamic == true
            density = material.ρ
        end
        block = blocks[block_index]
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_conn = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_conn)
        element_dofs_tl = create_threadlocal_element_vectors(Int64, Val(num_element_nodes))
        element_internal_force_tl = create_threadlocal_element_vectors(Float64, Val(num_element_nodes))
        if compute_diag_stiffness == true
            element_diag_stiffness_tl = create_threadlocal_element_vectors(Float64, Val(num_element_nodes))
        end
        if compute_lumped_mass == true
            element_lumped_mass_tl = create_threadlocal_element_vectors(Float64, Val(num_element_nodes))
        end
        if compute_stiffness == true
            element_stiffness_tl = create_threadlocal_element_matrices(Float64, Val(num_element_nodes))
        end
        if compute_mass == true
            element_mass_tl = create_threadlocal_element_matrices(Float64, Val(num_element_nodes))
        end
        @threads for block_element_index in 1:num_block_elements
            t = threadid()
            element_energy = 0.0
            element_dofs = element_dofs_tl[t]
            element_internal_force = element_internal_force_tl[t]
            fill!(element_internal_force, 0.0)
            if compute_diag_stiffness == true
                element_diag_stiffness = element_diag_stiffness_tl[t]
                fill!(element_diag_stiffness, 0.0)
            end
            if compute_lumped_mass == true
                element_lumped_mass = element_lumped_mass_tl[t]
                fill!(element_lumped_mass, 0.0)
            end
            if compute_stiffness == true
                element_stiffness = element_stiffness_tl[t]
                fill!(element_stiffness, 0.0)
            end
            if compute_mass == true
                element_mass = element_mass_tl[t]
                fill!(element_mass, 0.0)
            end
            conn_indices = ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_conn[conn_indices]
            if mesh_smoothing == true
                element_reference_position = create_smooth_reference(
                    model.smooth_reference, element_type, model.reference[:, node_indices]
                )
            else
                element_reference_position = model.reference[:, node_indices]
            end
            element_current_position = model.current[:, node_indices]
            element_dofs = reshape(3 .* node_indices' .- [2, 1, 0], :)
            for point in 1:num_points
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dNdX = dXdξ \ dNdξ
                F = SMatrix{3,3,Float64,9}(dNdX * element_current_position')
                J = det(F)
                if J ≤ 0.0 || isfinite(J) == false
                    model.failed = true
                    model.compute_stiffness = model.compute_diag_stiffness = true
                    model.compute_mass = model.compute_lumped_mass = true
                    println("⛔️ Non-positive Jacobian detected!")
                    println("⛔️ This may indicate element distortion.")
                    println("⏮️  Attempting to recover...")
                    return nothing
                end
                W, P, A = constitutive(material, F)
                stress = SVector{9,Float64}(P)
                ip_weight = ip_weights[point]
                det_dXdξ = det(dXdξ)
                dvol = det_dXdξ * ip_weight
                element_energy += W * dvol
                grad_op = create_gradient_operator(dNdX)
                @einsum element_internal_force[i] += grad_op[j, i] * stress[j] * dvol
                if compute_diag_stiffness == true
                    moduli = second_from_fourth(A)
                    @einsum element_diag_stiffness[i] += grad_op[m, i] * moduli[m, n] * grad_op[n, i] * dvol
                end
                if compute_lumped_mass == true
                    Nξ = N[:, point]
                    reduced_mass = Nξ * Nξ' * density * dvol
                    element_lumped_mass[1:3:end] += row_sum_lump(reduced_mass)
                end
                if compute_stiffness == true
                    moduli = second_from_fourth(A)
                    element_stiffness += grad_op' * moduli * grad_op * dvol
                end
                if compute_mass == true
                    Nξ = N[:, point]
                    reduced_mass = Nξ * Nξ' * density * dvol
                    element_mass[1:3:end, 1:3:end] += reduced_mass
                end
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[block_index][block_element_index][point] = voigt_cauchy
            end
            if compute_lumped_mass == true
                element_lumped_mass[3:3:end] .= element_lumped_mass[2:3:end] .= element_lumped_mass[1:3:end]
            end
            if compute_mass == true
                element_mass[3:3:end, 3:3:end] .= element_mass[2:3:end, 2:3:end] .= element_mass[1:3:end, 1:3:end]
            end
            energy_tl[t] += element_energy
            model.stored_energy[block_index][block_element_index] = element_energy
            assemble!(internal_force_tl[t], element_internal_force, element_dofs)
            if compute_diag_stiffness == true
                assemble!(diag_stiffness_tl[t], element_diag_stiffness, element_dofs)
            end
            if compute_lumped_mass == true
                assemble!(lumped_mass_tl[t], element_lumped_mass, element_dofs)
            end
            if compute_stiffness == true
                assemble!(stiffness_tl[t], element_stiffness, element_dofs)
            end
            if compute_mass == true
                assemble!(mass_tl[t], element_mass, element_dofs)
            end
        end
    end
    model.strain_energy = sum(energy_tl)
    model.body_force = body_force_vector
    model.internal_force = merge_threadlocal_coo_vectors(internal_force_tl, num_dofs)
    if compute_diag_stiffness == true
        model.diag_stiffness = merge_threadlocal_coo_vectors(diag_stiffness_tl, num_dofs)
    end
    if compute_lumped_mass == true
        model.lumped_mass = merge_threadlocal_coo_vectors(lumped_mass_tl, num_dofs)
    end
    if compute_stiffness == true
        model.stiffness = merge_threadlocal_coo_matrices(stiffness_tl, num_dofs)
    end
    if compute_mass == true
        model.mass = merge_threadlocal_coo_matrices(mass_tl, num_dofs)
    end
    if mesh_smoothing == true
        model.internal_force -= integrator.velocity
    end
    return nothing
end
