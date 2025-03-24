# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

include("constitutive.jl")
include("interpolation.jl")
include("ics_bcs.jl")

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
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return LinearOpInfRom(
        opinf_model,
        basis,
        reduced_state,
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
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return QuadraticOpInfRom(
        opinf_model,
        basis,
        reduced_state,
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
        element_type, num_blk_elems, _, _, _, _ = Exodus.read_block_parameters(
            input_mesh, blk_id
        )
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
        free_dofs,
        time,
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
        element_type, num_blk_elems, _, _, _, _ = Exodus.read_block_parameters(
            input_mesh, blk_id
        )
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

function create_smooth_reference(
    smooth_reference::String, element_type::String, elem_ref_pos::Matrix{Float64}
)
    if element_type == "TETRA4"
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

function get_minimum_edge_length(
    nodal_coordinates::Matrix{Float64}, edges::Vector{Tuple{Int64,Int64}}
)
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

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64}, element_type::String)
    if element_type == "TETRA4"
        edges = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4), (2, 4)]
        return get_minimum_edge_length(nodal_coordinates, edges)
    elseif element_type == "HEX8"
        edges = [
            (1, 4),
            (1, 5),
            (4, 8),
            (5, 8),
            (2, 3),
            (2, 6),
            (3, 7),
            (6, 7),
            (1, 2),
            (3, 4),
            (5, 6),
            (7, 8),
        ]
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
        element_type = Exodus.read_block_parameters(input_mesh, blk_id)[1]
        elem_blk_conn = get_block_connectivity(input_mesh, blk_id)
        num_blk_elems, num_elem_nodes = size(elem_blk_conn)
        for blk_elem_index in 1:num_blk_elems
            conn_indices =
                ((blk_elem_index - 1) * num_elem_nodes + 1):(blk_elem_index * num_elem_nodes)
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
        println(
            "Warning: Estimated stable time step: ",
            stable_time_step,
            " < provided time step: ",
            integrator.user_time_step,
        )
    end
    return integrator.time_step = min(stable_time_step, integrator.user_time_step)
end

function voigt_cauchy_from_stress(
    _::Solid, P::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9}, J::Float64
)
    # Compute the Cauchy stress tensor
    σ = F * P' ./ J

    # Return as an SVector for efficient indexing and stack-allocation
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function voigt_cauchy_from_stress(
    _::Linear_Elastic, σ::SMatrix{3,3,Float64,9}, _::SMatrix{3,3,Float64,9}, _::Float64
)
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function assemble!(
    rows::Vector{Int64},
    global_vector::Vector{Float64},
    elem_vector::AbstractVector{Float64},
    dofs::AbstractVector{Int64},
)
    ndofs = length(dofs)

    # old length and new length
    old_len = length(rows)
    new_len = old_len + ndofs

    # Resize each vector in one step
    resize!(rows, new_len)
    resize!(global_vector, new_len)

    idx = old_len + 1
    @inbounds for i in 1:ndofs
        I = dofs[i]
        rows[idx] = I
        global_vector[idx] = elem_vector[i]
        idx += 1
    end
    return nothing
end

function assemble!(
    rows::Vector{Int64},
    cols::Vector{Int64},
    global_matrix::Vector{Float64},
    elem_matrix::AbstractMatrix{Float64},
    dofs::AbstractVector{Int64},
)
    ndofs = length(dofs)
    n2 = ndofs * ndofs

    # old length and new length
    old_len = length(rows)
    new_len = old_len + n2

    # Resize each vector in one step
    resize!(rows, new_len)
    resize!(cols, new_len)
    resize!(global_matrix, new_len)

    idx = old_len + 1
    @inbounds for i in 1:ndofs
        I = dofs[i]
        @inbounds for j in 1:ndofs
            rows[idx] = I
            cols[idx] = dofs[j]
            global_matrix[idx] = elem_matrix[i, j]
            idx += 1
        end
    end
    return nothing
end

function dense(indices::Vector{Int64}, values::Vector{Float64}, vector_size::Int64)
    dense_vector = zeros(vector_size)
    @inbounds for i in 1:length(indices)
        dense_vector[indices[i]] += values[i]
    end
    return dense_vector
end

function create_element_matrix(::Type{T}, ::Val{4}) where {T}
    return MMatrix{12,12,T}(undef)
end

function create_element_matrix(::Type{T}, ::Val{8}) where {T}
    return MMatrix{24,24,T}(undef)
end

function create_element_matrix(::Type{T}, ::Val{10}) where {T}
    return MMatrix{30,30,T}(undef)
end

function create_element_vector(::Type{T}, ::Val{4}) where {T}
    return MVector{12,T}(undef)
end

function create_element_vector(::Type{T}, ::Val{8}) where {T}
    return MVector{24,T}(undef)
end

function create_element_vector(::Type{T}, ::Val{10}) where {T}
    return MVector{30,T}(undef)
end

function create_gradient_operator(::Type{T}, ::Val{4}) where {T}
    return MMatrix{9,12,T}(undef)
end

function create_gradient_operator(::Type{T}, ::Val{8}) where {T}
    return MMatrix{9,24,T}(undef)
end

function create_gradient_operator(::Type{T}, ::Val{10}) where {T}
    return MMatrix{9,30,T}(undef)
end

function create_coo_matrix()
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    vals = Vector{Float64}()
    return rows, cols, vals
end

function create_coo_vector()
    index = Vector{Int64}()
    vals = Vector{Float64}()
    return index, vals
end

function merge_threadlocal_coo_vectors(
    index_tl::Vector{Vector{Int64}}, vals_tl::Vector{Vector{Float64}}, num_dof::Int64
)
    index = vcat(index_tl...)
    vals = vcat(vals_tl...)
    return dense(index, vals, num_dof)
end

function merge_threadlocal_coo_matrices(
    rows_tl::Vector{Vector{Int64}},
    cols_tl::Vector{Vector{Int64}},
    vals_tl::Vector{Vector{Float64}},
    num_dof::Int64,
)
    rows = vcat(rows_tl...)
    cols = vcat(cols_tl...)
    vals = vcat(vals_tl...)
    return sparse(rows, cols, vals, num_dof, num_dof)
end

using Base.Threads: @threads, threadid, nthreads

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

function create_threadlocal_gradient_operators(::Type{T}, ::Val{N}) where {T,N}
    local_ops = Vector{typeof(create_gradient_operator(T, Val(N)))}(undef, nthreads())
    for i in 1:nthreads()
        local_ops[i] = create_gradient_operator(T, Val(N))
    end
    return local_ops
end

function create_threadlocal_coo_matrices()
    rows_tl = Vector{Vector{Int64}}(undef, nthreads())
    cols_tl = Vector{Vector{Int64}}(undef, nthreads())
    vals_tl = Vector{Vector{Float64}}(undef, nthreads())
    for i in 1:nthreads()
        rows_tl[i], cols_tl[i], vals_tl[i] = create_coo_matrix()
    end
    return rows_tl, cols_tl, vals_tl
end

function create_threadlocal_coo_vectors()
    index_tl = Vector{Vector{Int64}}(undef, nthreads())
    vals_tl = Vector{Vector{Float64}}(undef, nthreads())
    for i in 1:nthreads()
        index_tl[i], vals_tl[i] = create_coo_vector()
    end
    return index_tl, vals_tl
end

function evaluate(integrator::TimeIntegrator, model::SolidMechanics)
    is_implicit_dynamic = integrator isa Newmark
    is_explicit_dynamic = integrator isa CentralDifference
    is_implicit_static = integrator isa QuasiStatic
    is_dynamic = is_implicit_dynamic || is_explicit_dynamic
    is_implicit = is_implicit_dynamic || is_implicit_static
    materials = model.materials
    input_mesh = model.mesh
    mesh_smoothing = model.mesh_smoothing
    num_nodes = size(model.reference, 2)
    num_dof = 3 * num_nodes
    energy_tl = zeros(nthreads())
    index_int_force_tl, internal_force_tl = create_threadlocal_coo_vectors()
    if is_explicit_dynamic == true
        index_lumped_mass_tl, lumped_mass_tl = create_threadlocal_coo_vectors()
    end
    if is_implicit == true
        rows_stiff_tl, cols_stiff_tl, stiffness_tl = create_threadlocal_coo_matrices()
    end
    if is_implicit_dynamic == true
        rows_mass_tl, cols_mass_tl, mass_tl = create_threadlocal_coo_matrices()
    end
    body_force_vector = zeros(num_dof)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    for blk_index in 1:num_blks
        material = materials[blk_index]
        if is_dynamic == true
            ρ = material.ρ
        end
        block = blocks[blk_index]
        blk_id = block.id
        element_type = Exodus.read_block_parameters(input_mesh, blk_id)[1]
        num_points = default_num_int_pts(element_type)
        N, dN, elem_weights = isoparametric(element_type, num_points)
        elem_blk_conn = get_block_connectivity(input_mesh, blk_id)
        num_blk_elems, num_elem_nodes = size(elem_blk_conn)
        elem_dofs_tl = create_threadlocal_element_vectors(Int64, Val(num_elem_nodes))
        element_internal_force_tl = create_threadlocal_element_vectors(
            Float64, Val(num_elem_nodes)
        )
        if is_explicit_dynamic == true
            element_lumped_mass_tl = create_threadlocal_element_vectors(
                Float64, Val(num_elem_nodes)
            )
        end
        if is_implicit == true
            element_stiffness_tl = create_threadlocal_element_matrices(
                Float64, Val(num_elem_nodes)
            )
        end
        if is_implicit_dynamic == true
            element_mass_tl = create_threadlocal_element_matrices(
                Float64, Val(num_elem_nodes)
            )
        end
        B_tl = create_threadlocal_gradient_operators(Float64, Val(num_elem_nodes))
        @threads for blk_elem_index in 1:num_blk_elems
            t = threadid()
            elem_dofs = elem_dofs_tl[t]
            element_internal_force = element_internal_force_tl[t]
            if is_explicit_dynamic == true
                element_lumped_mass = element_lumped_mass_tl[t]
            end
            if is_implicit == true
                element_stiffness = element_stiffness_tl[t]
            end
            if is_implicit_dynamic == true
                element_mass = element_mass_tl[t]
            end
            B = B_tl[t]
            conn_indices =
                ((blk_elem_index - 1) * num_elem_nodes + 1):(blk_elem_index * num_elem_nodes)
            node_indices = elem_blk_conn[conn_indices]
            if mesh_smoothing == true
                elem_ref_pos = create_smooth_reference(
                    model.smooth_reference, element_type, model.reference[:, node_indices]
                )
            else
                elem_ref_pos = model.reference[:, node_indices]
            end
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            fill!(element_internal_force, 0.0)
            if is_explicit_dynamic == true
                fill!(element_lumped_mass, 0.0)
            end
            if is_implicit == true
                fill!(element_stiffness, 0.0)
            end
            if is_implicit_dynamic == true
                fill!(element_mass, 0.0)
            end
            elem_dofs = reshape(3 .* node_indices' .- [2, 1, 0], :)
            for point in 1:num_points
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * elem_ref_pos')
                dNdX = dXdξ \ dNdξ
                F = SMatrix{3,3,Float64,9}(dNdX * elem_cur_pos')
                J = det(F)
                gradient_operator!(B, dNdX)
                j = det(dXdξ)
                if J ≤ 0.0
                    model.failed = true
                    @info "Non-positive Jacobian detected! This may indicate element distortion. Attempting to recover by adjusting time step size..."
                    if is_implicit_static == true
                        return 0.0,
                        zeros(num_dof), zeros(num_dof),
                        spzeros(num_dof, num_dof)
                    elseif is_implicit_dynamic == true
                        return 0.0,
                        zeros(num_dof), zeros(num_dof), spzeros(num_dof, num_dof),
                        spzeros(num_dof, num_dof)
                    elseif is_explicit_dynamic == true
                        return 0.0, zeros(num_dof), zeros(num_dof), zeros(num_dof)
                    else
                        error("Unknown type of time integrator", typeof(integrator))
                    end
                end
                W, P, A = constitutive(material, F)
                stress = P[1:9]
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                if is_implicit == true
                    moduli = second_from_fourth(A)
                    element_stiffness += B' * moduli * B * j * w
                end
                if is_dynamic == true
                    Nξ = N[:, point]
                    reduced_mass = Nξ * Nξ' * ρ * j * w
                end
                if is_implicit_dynamic == true
                    element_mass[1:3:end, 1:3:end] += reduced_mass
                end
                if is_explicit_dynamic == true
                    reduced_lumped_mass = sum(reduced_mass; dims=2)
                    element_lumped_mass[1:3:end] += reduced_lumped_mass
                end
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[blk_index][blk_elem_index][point] = voigt_cauchy
            end
            if is_implicit_dynamic == true
                element_mass[3:3:end, 3:3:end] .=
                    element_mass[2:3:end, 2:3:end] .= element_mass[1:3:end, 1:3:end]
            end
            if is_explicit_dynamic == true
                element_lumped_mass[3:3:end] .=
                    element_lumped_mass[2:3:end] .= element_lumped_mass[1:3:end]
            end
            energy_tl[t] += element_energy
            model.stored_energy[blk_index][blk_elem_index] = element_energy
            assemble!(
                index_int_force_tl[t],
                internal_force_tl[t],
                element_internal_force,
                elem_dofs,
            )
            if is_explicit_dynamic == true
                assemble!(
                    index_lumped_mass_tl[t],
                    lumped_mass_tl[t],
                    element_lumped_mass,
                    elem_dofs,
                )
            end
            if is_implicit == true
                assemble!(
                    rows_stiff_tl[t],
                    cols_stiff_tl[t],
                    stiffness_tl[t],
                    element_stiffness,
                    elem_dofs,
                )
            end
            if is_implicit_dynamic == true
                assemble!(
                    rows_mass_tl[t], cols_mass_tl[t], mass_tl[t], element_mass, elem_dofs
                )
            end
        end
    end
    energy = sum(energy_tl)
    internal_force_vector = merge_threadlocal_coo_vectors(
        index_int_force_tl, internal_force_tl, num_dof
    )
    if is_explicit_dynamic == true
        lumped_mass_vector = merge_threadlocal_coo_vectors(
            index_lumped_mass_tl, lumped_mass_tl, num_dof
        )
    end
    if is_implicit == true
        stiffness_matrix = merge_threadlocal_coo_matrices(
            rows_stiff_tl, cols_stiff_tl, stiffness_tl, num_dof
        )
    end
    if is_implicit_dynamic == true
        mass_matrix = merge_threadlocal_coo_matrices(
            rows_mass_tl, cols_mass_tl, mass_tl, num_dof
        )
    end
    if mesh_smoothing == true
        internal_force_vector -= integrator.velocity
    end
    model.internal_force = internal_force_vector
    if is_implicit_static == true
        return energy, internal_force_vector, body_force_vector, stiffness_matrix
    elseif is_implicit_dynamic == true
        return energy,
        internal_force_vector, body_force_vector, stiffness_matrix,
        mass_matrix
    elseif is_explicit_dynamic == true
        return energy, internal_force_vector, body_force_vector, lumped_mass_vector
    else
        error("Unknown type of time integrator", typeof(integrator))
    end
end
