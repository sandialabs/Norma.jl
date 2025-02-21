# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
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
    LinearOpInfRom(
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
    for node ∈ 1:num_nodes
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
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        if kinematics == Undefined
            kinematics = get_kinematics(material_model)
        else
            if kinematics ≠ get_kinematics(material_model)
                error("Material ", typeof(material_model), " has inconsistent kinematics ",
                    get_kinematics(material_model), " than previous materials of type ", kinematics)
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
    for block ∈ blocks
        blk_id = block.id
        element_type, num_blk_elems, _, _, _, _ =
            Exodus.read_block_parameters(input_mesh, blk_id)
        num_points = default_num_int_pts(element_type)
        block_stress = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ ∈ 1:num_blk_elems
            element_stress = Vector{Vector{Float64}}()
            for _ ∈ 1:num_points
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

    # BRP: define a global transform for inclined support
    inclined_support = false
    global_transform = sparse(Diagonal(ones(3 * num_nodes)))

    SolidMechanics(
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
    for node ∈ 1:num_nodes
        reference[:, node] = coords[:, node]
        temperature[node] = 0.0
        rate[node] = 0.0
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(elem_blk_ids)
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
    for elem_blk_name ∈ elem_blk_names
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
    for block ∈ blocks
        blk_id = block.id
        element_type, num_blk_elems, _, _, _, _ =
            Exodus.read_block_parameters(input_mesh, blk_id)
        num_points = default_num_int_pts(element_type)
        block_flux = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ ∈ 1:num_blk_elems
            element_flux = Vector{Vector{Float64}}()
            for _ ∈ 1:num_points
                push!(element_flux, zeros(3))
            end
            push!(block_flux, element_flux)
            element_stored_energy = 0.0
            push!(block_stored_energy, element_stored_energy)
        end
        push!(flux, block_flux)
        push!(stored_energy, block_stored_energy)
    end
    HeatConduction(
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

    else
        error("Unknown type of model : ", model_name)
    end
end

function create_smooth_reference(smooth_reference::String, element_type::String, elem_ref_pos::Matrix{Float64})
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
    nodal_coordinates::Matrix{Float64},
    edges::Vector{Tuple{Int64, Int64}},
)
    minimum_edge_length = Inf
    for edge ∈ edges
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
    for blk_index ∈ 1:num_blks
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
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
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
    integrator.time_step = min(stable_time_step, integrator.user_time_step)
end

function voigt_cauchy_from_stress(
    _::Solid,
    P::Matrix{Float64},
    F::Matrix{Float64},
    J::Float64,
)
    σ = F * P' ./ J
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function voigt_cauchy_from_stress(
    _::Linear_Elastic,
    σ::Matrix{Float64},
    _::Matrix{Float64},
    _::Float64,
)
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function assemble(
    rows::Vector{Int64},
    cols::Vector{Int64},
    global_stiffness::Vector{Float64},
    element_stiffness::Matrix{Float64},
    dofs::Vector{Int64},
)
    num_dofs = length(dofs)
    for i ∈ 1:num_dofs
        I = dofs[i]
        for j ∈ 1:num_dofs
            J = dofs[j]
            push!(rows, I)
            push!(cols, J)
            push!(global_stiffness, element_stiffness[i, j])
        end
    end
end

function assemble(
    rows::Vector{Int64},
    cols::Vector{Int64},
    global_stiffness::Vector{Float64},
    global_mass::Vector{Float64},
    element_stiffness::Matrix{Float64},
    element_mass::Matrix{Float64},
    dofs::Vector{Int64},
)
    num_dofs = length(dofs)
    for i ∈ 1:num_dofs
        I = dofs[i]
        for j ∈ 1:num_dofs
            J = dofs[j]
            push!(rows, I)
            push!(cols, J)
            push!(global_mass, element_mass[i, j])
            push!(global_stiffness, element_stiffness[i, j])
        end
    end
end

function evaluate(integrator::TimeIntegrator, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    mesh_smoothing = model.mesh_smoothing
    num_nodes = size(model.reference)[2]
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    body_force = zeros(num_dof)
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    stiffness = Vector{Float64}()
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blks = length(blocks)
    if typeof(integrator) == Newmark
        mass = Vector{Float64}()
    end
    if typeof(integrator) == CentralDifference
        lumped_mass = zeros(num_dof)
    end

    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        if typeof(integrator) == Newmark || typeof(integrator) == CentralDifference
            ρ = material.ρ
        end
        block = blocks[blk_index]
        blk_id = block.id
        element_type = Exodus.read_block_parameters(input_mesh, blk_id)[1]
        num_points = default_num_int_pts(element_type)
        N, dNdξ, elem_weights = isoparametric(element_type, num_points)
        elem_blk_conn = get_block_connectivity(input_mesh, blk_id)
        num_blk_elems, num_elem_nodes = size(elem_blk_conn)
        num_elem_dofs = 3 * num_elem_nodes
        elem_dofs = zeros(Int64, num_elem_dofs)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            if mesh_smoothing == true
                elem_ref_pos =
                    create_smooth_reference(model.smooth_reference, element_type, model.reference[:, node_indices])
            else
                elem_ref_pos = model.reference[:, node_indices]
            end
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = zeros(num_elem_dofs)
            element_stiffness = zeros(num_elem_dofs, num_elem_dofs)
            if typeof(integrator) == Newmark
                element_mass = zeros(num_elem_dofs, num_elem_dofs)
            end
            if typeof(integrator) == CentralDifference
                element_lumped_mass = zeros(num_elem_dofs)
            end
            index_x = 1:3:num_elem_dofs.-2
            index_y = index_x .+ 1
            index_z = index_x .+ 2
            elem_dofs[index_x] = 3 .* node_indices .- 2
            elem_dofs[index_y] = 3 .* node_indices .- 1
            elem_dofs[index_z] = 3 .* node_indices
            for point ∈ 1:num_points
                dNdξₚ = dNdξ[:, :, point]
                dXdξ = dNdξₚ * elem_ref_pos'
                dxdξ = dNdξₚ * elem_cur_pos'
                if det(dxdξ) ≤ 0.0
                    model.failed = true
                    @info "Non-positive Jacobian detected! This may indicate element distortion. Attempting to recover by adjusting time step size..."
                    if typeof(integrator) == QuasiStatic
                        return 0.0, zeros(num_dof), zeros(num_dof), spzeros(num_dof, num_dof)
                    elseif typeof(integrator) == Newmark
                        return 0.0, zeros(num_dof), zeros(num_dof), spzeros(num_dof, num_dof), spzeros(num_dof, num_dof)
                    elseif typeof(integrator) == CentralDifference
                        return 0.0, zeros(num_dof), zeros(num_dof), zeros(num_dof)
                    else
                        error("Unknown type of time integrator", typeof(integrator))
                    end
                end
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξₚ
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                J = det(dxdX)
                F = dxdX
                W, P, A = constitutive(material, F)
                stress = P[1:9]
                moduli = second_from_fourth(A)
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                element_stiffness += B' * moduli * B * j * w
                if typeof(integrator) == Newmark
                    reduced_mass = N[:, point] * N[:, point]' * ρ * j * w
                    element_mass[index_x, index_x] += reduced_mass
                    element_mass[index_y, index_y] += reduced_mass
                    element_mass[index_z, index_z] += reduced_mass
                end
                if typeof(integrator) == CentralDifference
                    reduced_mass = N[:, point] * N[:, point]' * ρ * j * w
                    reduced_lumped_mass = sum(reduced_mass, dims = 2)
                    element_lumped_mass[index_x] += reduced_lumped_mass
                    element_lumped_mass[index_y] += reduced_lumped_mass
                    element_lumped_mass[index_z] += reduced_lumped_mass
                end
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[blk_index][blk_elem_index][point] = voigt_cauchy
            end
            energy += element_energy
            model.stored_energy[blk_index][blk_elem_index] = element_energy
            internal_force[elem_dofs] += element_internal_force
            if typeof(integrator) == QuasiStatic
                assemble(rows, cols, stiffness, element_stiffness, elem_dofs)
            end
            if typeof(integrator) == Newmark
                assemble(rows, cols, stiffness, mass, element_stiffness, element_mass, elem_dofs)
            end
            if typeof(integrator) == CentralDifference
                lumped_mass[elem_dofs] += element_lumped_mass
            end
        end
    end

    if typeof(integrator) == QuasiStatic || typeof(integrator) == Newmark
        stiffness_matrix = sparse(rows, cols, stiffness)
    end
    if typeof(integrator) == Newmark
        mass_matrix = sparse(rows, cols, mass)
    end
    if mesh_smoothing == true
        internal_force -= integrator.velocity
    end
    model.internal_force = internal_force
    if typeof(integrator) == QuasiStatic
        return energy, internal_force, body_force, stiffness_matrix
    elseif typeof(integrator) == Newmark
        return energy, internal_force, body_force, stiffness_matrix, mass_matrix
    elseif typeof(integrator) == CentralDifference
        return energy, internal_force, body_force, lumped_mass
    else
        error("Unknown type of time integrator", typeof(integrator))
    end
end
