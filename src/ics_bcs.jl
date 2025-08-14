# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

include("opinf/opinf_ics_bcs.jl")
using Exodus

@variables t x y z
D = Differential(t)

function SolidMechanicsDirichletBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)

    # Build symbolic expressions
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))

    # Compile them into functions
    disp_fun = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
    velo_fun = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
    acce_fun = eval(build_function(acce_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsDirichletBoundaryCondition(
        node_set_name, offset, node_set_id, node_set_node_indices, disp_fun, velo_fun, acce_fun
    )
end

function SolidMechanicsInclinedDirichletBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)

    # Build symbolic expressions
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))

    # Compile them into functions
    disp_fun = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
    velo_fun = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
    acce_fun = eval(build_function(acce_num, [t, x, y, z]; expression=Val(false)))

    reference_normal_expression_vector = bc_params["normal vector"]
    reference_normal_nums = [eval(Meta.parse(string(ref_exp))) for ref_exp in reference_normal_expression_vector]

    reference_funs = [
        eval(build_function(norm_num, [t, x, y, z]; expression=Val(false))) for norm_num in reference_normal_nums
    ]
    return SolidMechanicsInclinedDirichletBoundaryCondition(
        node_set_name, node_set_id, node_set_node_indices, disp_fun, velo_fun, acce_fun, reference_funs
    )
end

function SolidMechanicsNeumannBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    traction_num = eval(Meta.parse(expression))

    # Compile them into functions
    traction_fun = eval(build_function(traction_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsNeumannBoundaryCondition(
        side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, traction_fun
    )
end

function SolidMechanicsNeumannPressureBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    pressure_num = eval(Meta.parse(expression))

    # Compile them into functions
    pressure_fun = eval(build_function(pressure_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsNeumannPressureBoundaryCondition(
        side_set_name, side_set_id, num_nodes_per_side, side_set_node_indices, pressure_fun
    )
end

function SolidMechanicsOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    side_set_name::String,
    side_set_node_indices::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    variational::Bool,
)
    if coupled_subsim.model isa RomModel
        coupled_mesh = coupled_subsim.model.fom_model.mesh
    else
        coupled_mesh = coupled_subsim.model.mesh
    end
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type_string = Exodus.read_block_parameters(coupled_mesh, coupled_block_id)[1]
    element_type = element_type_from_string(element_type_string)
    coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    unique_node_indices = unique(side_set_node_indices)
    for node_index in unique_node_indices
        point = subsim.model.reference[:, node_index]
        node_indices, ξ, found = find_point_in_mesh(point, coupled_subsim.model, coupled_block_id, tol)
        if found == false
            norma_abortf(
                "Could not find subdomain %s point (%.4e, %.4e, %.4e) in subdomain %s",
                subsim.name,
                point[1],
                point[2],
                point[3],
                coupled_subsim.name,
            )
        end
        N = interpolate(element_type, ξ)[1]
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    return SolidMechanicsOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_node_indices,
        coupled_nodes_indices,
        interpolation_function_values,
        coupled_subsim,
        subsim,
        variational,
    )
end

function SolidMechanicsContactSchwarzBoundaryCondition(
    coupled_subsim::SingleDomainSimulation, input_mesh::ExodusDatabase, bc_params::Parameters
)
    side_set_name = bc_params["side set"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)
    coupled_side_set_name = bc_params["source side set"]
    is_dirichlet = true
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(input_mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(input_mesh, side_set_id)
    coupled_bc_index = 0
    rotation_matrix = I(3)
    active_contact = false
    swap_bcs = get(bc_params, "swap BC types", false)
    friction_type_string = bc_params["friction type"]
    if friction_type_string == "frictionless"
        friction_type = 0
    elseif friction_type_string == "tied"
        friction_type = 1
    else
        norma_abort("Unknown or not implemented friction type : $friction_type_string")
    end
    variational = get(bc_params, "variational", false)
    return SolidMechanicsContactSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_subsim,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        is_dirichlet,
        swap_bcs,
        rotation_matrix,
        active_contact,
        friction_type,
        variational,
    )
end

function SolidMechanicsNonOverlapSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    is_dirichlet::Bool,
    swap_bcs::Bool,
    variational::Bool,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsNonOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_subsim,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        is_dirichlet,
        swap_bcs,
        variational,
    )
end

function SMCouplingSchwarzBC(
    subsim::SingleDomainSimulation,
    coupled_subsim::SingleDomainSimulation,
    input_mesh::ExodusDatabase,
    bc_type::String,
    bc_params::Parameters,
)
    side_set_name = bc_params["side set"]
    coupled_block_name = get(bc_params, "source block", "")
    tol = get(bc_params, "search tolerance", 1.0e-06)
    coupled_side_set_name = bc_params["source side set"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    num_nodes_sides = Int64.(num_nodes_sides)
    side_set_node_indices = Int64.(side_set_node_indices)
    variational = get(bc_params, "variational", false)
    if bc_type == "Schwarz overlap"
        SolidMechanicsOverlapSchwarzBoundaryCondition(
            coupled_block_name, tol, side_set_name, side_set_node_indices, coupled_subsim, subsim, variational
        )
    elseif bc_type == "Schwarz nonoverlap"
        default_bc_type = get(bc_params, "default BC type", "Dirichlet")
        if default_bc_type == "Dirichlet"
            is_dirichlet = true
        elseif default_bc_type == "Neumann"
            is_dirichlet = false
        else
            norma_abort("Invalid string for 'default BC type'!  Valid options are 'Dirichlet' and 'Neumann'")
        end
        swap_bcs = get(bc_params, "swap BC types", false)
        SolidMechanicsNonOverlapSchwarzBoundaryCondition(
            input_mesh,
            side_set_name,
            coupled_side_set_name,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            coupled_subsim,
            is_dirichlet,
            swap_bcs,
            variational,
        )
    else
        norma_abort("Unknown boundary condition type : $bc_type")
    end
end


function apply_bc(model::SolidMechanics, bc::SolidMechanicsDirichletBoundaryCondition)
    for node_index in bc.node_set_node_indices
        txzy = (
            model.time, model.reference[1, node_index], model.reference[2, node_index], model.reference[3, node_index]
        )

        disp_val = bc.disp_fun(txzy)
        velo_val = bc.velo_fun(txzy)
        acce_val = bc.acce_fun(txzy)

        dof_index = 3 * (node_index - 1) + bc.offset
        model.current[bc.offset, node_index] = model.reference[bc.offset, node_index] + disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsNeumannBoundaryCondition)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc.traction_fun, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index in side_nodes
            bc_val = nodal_force_component[side_node_index]
            side_node_index += 1
            dof_index = 3 * (node_index - 1) + bc.offset
            model.boundary_force[dof_index] += bc_val
        end
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsNeumannPressureBoundaryCondition)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_pressure(side_coordinates, bc.pressure_fun, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index in side_nodes
            dof_indices = [3 * node_index - 2, 3 * node_index - 1, 3 * node_index]
            model.boundary_force[dof_indices] += nodal_force_component[:, side_node_index]
            side_node_index += 1
        end
    end
end

function compute_rotation_matrix(axis::SVector{3,Float64})::SMatrix{3,3,Float64}
    e1 = @SVector [1.0, 0.0, 0.0]
    angle_btwn = acos(dot(axis, e1))
    w = cross(axis, e1)
    s = clamp(norm(w), -1.0, 1.0)
    if isapprox(angle_btwn, 0.0; atol=1e-12)
        return @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    elseif isapprox(angle_btwn, π; atol=1e-12)
        return @SMatrix [
            -1.0 0.0 0.0
            0.0 -1.0 0.0
            0.0 0.0 1.0
        ]
    else
        θ = angle_btwn > π / 2 ? π - asin(s) : asin(s)
        m = normalize(w)
        rv = θ * m
        return rt_of_rv(rv)
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsInclinedDirichletBoundaryCondition)
    for node_index in bc.node_set_node_indices
        txzy = (
            model.time, model.reference[1, node_index], model.reference[2, node_index], model.reference[3, node_index]
        )

        disp_val = bc.disp_fun(txzy)
        velo_val = bc.velo_fun(txzy)
        acce_val = bc.acce_fun(txzy)

        # Local basis from reference normal
        normal_vals = [norm_fun(txzy) for norm_fun in bc.reference_funs]
        axis = normalize(SVector{3,Float64}(normal_vals))
        rotation_matrix = compute_rotation_matrix(axis)

        original_disp = model.current[:, node_index] - model.reference[:, node_index]
        original_velocity = model.velocity[:, node_index]
        original_acceleration = model.acceleration[:, node_index]

        local_original_displacement = MVector(rotation_matrix * original_disp)
        local_original_velocity = MVector(rotation_matrix * original_velocity)
        local_original_acceleration = MVector(rotation_matrix * original_acceleration)

        model.free_dofs[3 * node_index - 2] = false
        model.free_dofs[3 * node_index - 1] = true
        model.free_dofs[3 * node_index] = true

        local_original_displacement[1] = disp_val
        local_original_velocity[1] = velo_val
        local_original_acceleration[1] = acce_val

        model.velocity[:, node_index] = rotation_matrix' * SVector(local_original_velocity)
        model.acceleration[:, node_index] = rotation_matrix' * SVector(local_original_acceleration)
        model.current[:, node_index] =
            model.reference[:, node_index] + rotation_matrix' * SVector(local_original_displacement)

        global_base = 3 * (node_index - 1)
        model.global_transform[(global_base + 1):(global_base + 3), (global_base + 1):(global_base + 3)] =
            rotation_matrix
    end
end

function find_point_in_mesh(point::Vector{Float64}, model::SolidMechanics, block_id::Int, tol::Float64)
    mesh = model.mesh
    element_type_string = Exodus.read_block_parameters(mesh, Int32(block_id))[1]
    element_type = element_type_from_string(element_type_string)
    element_block_connectivity = get_block_connectivity(mesh, block_id)
    num_block_elements, num_element_nodes = size(element_block_connectivity)
    node_indices = Vector{Int64}()
    found = false
    ξ = zeros(length(point))
    for block_element_index in 1:num_block_elements
        connectivity_indices =
            ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
        node_indices = element_block_connectivity[connectivity_indices]
        element_ref_pos = model.reference[:, node_indices]
        ξ, found = is_inside(element_type, element_ref_pos, point, tol)
        if found == true
            break
        end
    end
    return node_indices, ξ, found
end


function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    if bc.is_dirichlet == true
        contact_variational_dbc(model, bc)
    else
        contact_variational_nbc(model, bc)
    end
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    coupling_pointwise_dbc(model, bc)
    return nothing
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    if bc.is_dirichlet == true
        coupling_variational_dbc(model, bc)
    else
        coupling_variational_nbc(model, bc)
    end
end


function coupling_pointwise_dbc(model::SolidMechanics, bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    get_coupled_field = if bc.coupled_subsim.model isa SolidMechanics
        (field -> getfield(bc.coupled_subsim.model, field))
    else
        (field -> getfield(bc.coupled_subsim.model.fom_model, field))
    end

    current = get_coupled_field(:current)
    velocity = get_coupled_field(:velocity)
    acceleration = get_coupled_field(:acceleration)

    unique_node_indices = unique(bc.side_set_node_indices)

    for i in eachindex(unique_node_indices)
        node_index = unique_node_indices[i]
        coupled_node_indices = bc.coupled_nodes_indices[i]
        N = bc.interpolation_function_values[i]

        model.current[:, node_index] = current[:, coupled_node_indices] * N
        model.velocity[:, node_index] = velocity[:, coupled_node_indices] * N
        model.acceleration[:, node_index] = acceleration[:, coupled_node_indices] * N

        dof_index = (3 * node_index - 2):(3 * node_index)
        model.free_dofs[dof_index] .= false
    end
end

function coupling_variational_dbc(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    nodal_curr, nodal_velo, nodal_acce = get_dst_curr_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        @inbounds model.current[:, i_global] = nodal_curr[:, i_local]
        @inbounds model.velocity[:, i_global] = nodal_velo[:, i_local]
        @inbounds model.acceleration[:, i_global] = nodal_acce[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        model.free_dofs[global_range] .= false
    end
end

function coupling_variational_nbc(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    nodal_force = get_dst_force(bc)
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        local_range = (3 * (i_local - 1) + 1):(3 * i_local)
        @inbounds model.boundary_force[global_range] += nodal_force[local_range]
    end
end

function get_internal_force(model::SolidMechanics)
    return model.internal_force
end


function set_internal_force!(model::SolidMechanics, force)
    return model.internal_force = force
end

function apply_bc(model::Model, bc::SolidMechanicsSchwarzBoundaryCondition)
    parent_sim = bc.coupled_subsim.params["parent_simulation"]
    controller = parent_sim.controller

    # Skip application if contact is inactive
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact
        return nothing
    end

    coupled_subsim = bc.coupled_subsim
    integrator = coupled_subsim.integrator
    coupled_model = coupled_subsim.model

    # Save current state
    saved_disp = integrator.displacement
    saved_velo = integrator.velocity
    saved_acce = integrator.acceleration
    saved_∂Ω_f = get_internal_force(coupled_model)

    # Fetch interpolation inputs
    time = model.time
    coupled_index = parent_sim.subsim_name_index_map[coupled_subsim.name]
    num_dofs = length(coupled_subsim.model.free_dofs)

    time_hist = controller.time_hist[coupled_index]
    disp_hist = controller.disp_hist[coupled_index]
    velo_hist = controller.velo_hist[coupled_index]
    acce_hist = controller.acce_hist[coupled_index]
    ∂Ω_f_hist = controller.∂Ω_f_hist[coupled_index]

    # Interpolate or use fallback
    if !isempty(time_hist)
        interp_disp = interpolate(time_hist, disp_hist, time)
        interp_velo = interpolate(time_hist, velo_hist, time)
        interp_acce = interpolate(time_hist, acce_hist, time)
        interp_∂Ω_f = interpolate(time_hist, ∂Ω_f_hist, time)
    elseif isempty(time_hist) && !isempty(controller.stop_disp[coupled_index])
        interp_disp = controller.stop_disp[coupled_index]
        interp_velo = controller.stop_velo[coupled_index]
        interp_acce = controller.stop_acce[coupled_index]
        interp_∂Ω_f = controller.stop_∂Ω_f[coupled_index]
    else
        interp_disp = zeros(num_dofs)
        interp_velo = zeros(num_dofs)
        interp_acce = zeros(num_dofs)
        interp_∂Ω_f = zeros(num_dofs)
    end

    # Assign interpolated force
    set_internal_force!(coupled_model, interp_∂Ω_f)

    # Apply relaxed update if needed
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition || bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
        θ = controller.relaxation_parameter
        iter = controller.iteration_number

        λ_u_prev = iter < 2 ? interp_disp : controller.lambda_disp[coupled_index]
        λ_v_prev = iter < 2 ? interp_velo : controller.lambda_velo[coupled_index]
        λ_a_prev = iter < 2 ? interp_acce : controller.lambda_acce[coupled_index]

        controller.lambda_disp[coupled_index] = θ * interp_disp + (1 - θ) * λ_u_prev
        controller.lambda_velo[coupled_index] = θ * interp_velo + (1 - θ) * λ_v_prev
        controller.lambda_acce[coupled_index] = θ * interp_acce + (1 - θ) * λ_a_prev

        integrator.displacement = controller.lambda_disp[coupled_index]
        integrator.velocity = controller.lambda_velo[coupled_index]
        integrator.acceleration = controller.lambda_acce[coupled_index]
    else
        integrator.displacement = interp_disp
        integrator.velocity = interp_velo
        integrator.acceleration = interp_acce
    end

    # Apply boundary condition detail
    copy_solution_source_targets(integrator, coupled_subsim.solver, coupled_subsim.model)
    apply_bc_detail(model, bc)

    # Restore previous state
    integrator.displacement = saved_disp
    integrator.velocity = saved_velo
    integrator.acceleration = saved_acce
    set_internal_force!(coupled_model, saved_∂Ω_f)

    copy_solution_source_targets(integrator, coupled_subsim.solver, coupled_subsim.model)
    return nothing
end

function transfer_normal_component(source::Vector{Float64}, target::Vector{Float64}, normal::Vector{Float64})
    normal_projection = normal * normal'
    tangent_projection = I(length(normal)) - normal_projection
    return tangent_projection * target + normal_projection * source
end

function contact_variational_dbc(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    nodal_curr, nodal_velo, nodal_acce = get_dst_curr_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    for (i_local, i_global) in enumerate(global_from_local_map)
        normal = normals[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        if bc.friction_type == 0
            @inbounds model.current[:, i_global] = transfer_normal_component(
                nodal_curr[:, i_local], model.current[:, i_global], normal
            )
            @inbounds model.velocity[:, i_global] = transfer_normal_component(
                nodal_velo[:, i_local], model.velocity[:, i_global], normal
            )
            @inbounds model.acceleration[:, i_global] = transfer_normal_component(
                nodal_acce[:, i_local], model.acceleration[:, i_global], normal
            )
            model.free_dofs[[3 * i_global - 2]] .= false
            model.free_dofs[[3 * i_global - 1]] .= true
            model.free_dofs[[3 * i_global]] .= true
        elseif bc.friction_type == 1
            @inbounds model.current[:, i_global] = nodal_curr[:, i_local]
            @inbounds model.velocity[:, i_global] = nodal_velo[:, i_local]
            @inbounds model.acceleration[:, i_global] = nodal_velo[:, i_local]
            model.free_dofs[global_range] .= false
        else
            norma_abort("Unknown or not implemented friction type.")
        end
        # Update the rotation matrix
        axis = SVector{3,Float64}(-normalize(normal))
        bc.rotation_matrix = compute_rotation_matrix(axis)
        model.global_transform[global_range, global_range] = bc.rotation_matrix
    end
end

function apply_naive_stabilized_bcs(subsim::SingleDomainSimulation)
    bcs = subsim.model.boundary_conditions
    for bc in bcs
        if bc isa SolidMechanicsContactSchwarzBoundaryCondition
            unique_node_indices = unique(bc.side_set_node_indices)
            for node_index in unique_node_indices
                subsim.model.acceleration[:, node_index] .= 0.0
            end
        end
    end
    copy_solution_source_targets(subsim.model, subsim.integrator, subsim.solver)
    return nothing
end

function contact_variational_nbc(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    friction_type = bc.friction_type
    nodal_force = get_dst_force(bc)
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        local_range = (3 * (i_local - 1) + 1):(3 * i_local)
        normal = normals[:, i_local]
        node_force = nodal_force[local_range]
        if friction_type == 0
            target = model.boundary_force[global_range]
            eff_node_force = transfer_normal_component(node_force, target, normal)
        else
            eff_node_force = node_force
        end
        @inbounds model.boundary_force[global_range] += eff_node_force
    end
end

function extract_local_vector(global_vector::Vector{Float64}, global_from_local_map::Vector{Int64}, dim::Int64)
    num_local_nodes = length(global_from_local_map)
    local_vector = Vector{Float64}(undef, dim * num_local_nodes)
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (dim * (i_global - 1) + 1):(dim * i_global)
        local_range = (dim * (i_local - 1) + 1):(dim * i_local)
        @inbounds local_vector[local_range] = global_vector[global_range]
    end
    return local_vector
end

function extract_local_vector(bc::SolidMechanicsSchwarzBoundaryCondition, global_vector::Vector{Float64}, dim::Int64)
    global_from_local_map = bc.global_from_local_map
    return extract_local_vector(global_vector, global_from_local_map, dim)
end

function compute_neumann_projector(dst_model::SolidMechanics, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_model = dst_bc.coupled_subsim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    H = get_square_projection_matrix(src_model, src_bc)
    L = get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_bc.neumann_projector = L * (H \ I)
    return nothing
end

function compute_dirichlet_projector(dst_model::SolidMechanics, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_model = dst_bc.coupled_subsim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    W = get_square_projection_matrix(dst_model, dst_bc)
    L = get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_bc.dirichlet_projector = (W \ I) * L
    return nothing
end

function get_dst_force(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = dst_bc.coupled_subsim
    src_model = src_sim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_force = src_model.internal_force
    src_force = -extract_local_vector(src_bc, src_global_force, 3)
    neumann_projector = dst_bc.neumann_projector
    num_dst_nodes = size(neumann_projector, 1)
    dst_force = zeros(3 * num_dst_nodes)
    dst_force[1:3:end] = neumann_projector * src_force[1:3:end]
    dst_force[2:3:end] = neumann_projector * src_force[2:3:end]
    dst_force[3:3:end] = neumann_projector * src_force[3:3:end]
    return dst_force
end

function get_dst_curr_velo_acce(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = dst_bc.coupled_subsim
    src_model = src_sim.model isa RomModel ? src_sim.model.fom_model : src_sim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
    num_src_nodes = length(src_global_from_local_map)
    src_curr = zeros(3, num_src_nodes)
    src_velo = zeros(3, num_src_nodes)
    src_acce = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_curr[:, i_local] = src_model.current[:, i_global]
        src_velo[:, i_local] = src_model.velocity[:, i_global]
        src_acce[:, i_local] = src_model.acceleration[:, i_global]
    end
    dirichlet_projector = dst_bc.dirichlet_projector
    num_dst_nodes = size(dirichlet_projector, 1)
    dst_curr = zeros(3, num_dst_nodes)
    dst_velo = zeros(3, num_dst_nodes)
    dst_acce = zeros(3, num_dst_nodes)
    for i in 1:3
        dst_curr[i, :] = dirichlet_projector * src_curr[i, :]
        dst_velo[i, :] = dirichlet_projector * src_velo[i, :]
        dst_acce[i, :] = dirichlet_projector * src_acce[i, :]
    end
    return dst_curr, dst_velo, dst_acce
end

function node_set_id_from_name(node_set_name::String, mesh::ExodusDatabase)
    node_set_names = Exodus.read_names(mesh, NodeSet)
    num_names = length(node_set_names)
    node_set_index = 0
    for index in 1:num_names
        if node_set_name == node_set_names[index]
            node_set_index = index
            break
        end
    end
    if node_set_index == 0
        norma_abort("node set $node_set_name cannot be found in mesh")
    end
    node_set_ids = Exodus.read_ids(mesh, NodeSet)
    node_set_id = node_set_ids[node_set_index]
    return Int64(node_set_id)
end

function side_set_id_from_name(side_set_name::String, mesh::ExodusDatabase)
    side_set_names = Exodus.read_names(mesh, SideSet)
    num_names = length(side_set_names)
    side_set_index = 0
    for index in 1:num_names
        if side_set_name == side_set_names[index]
            side_set_index = index
            break
        end
    end
    if side_set_index == 0
        norma_abort("side set $side_set_name cannot be found in mesh")
    end
    side_set_ids = Exodus.read_ids(mesh, SideSet)
    side_set_id = side_set_ids[side_set_index]
    return Int64(side_set_id)
end

function block_id_from_name(block_name::String, mesh::ExodusDatabase)
    block_names = Exodus.read_names(mesh, Block)
    num_names = length(block_names)
    block_index = 0
    for index in 1:num_names
        if block_name == block_names[index]
            block_index = index
            break
        end
    end
    if block_index == 0
        norma_abort("block $block_name cannot be found in mesh")
    end
    block_ids = Exodus.read_ids(mesh, Block)
    block_id = block_ids[block_index]
    return Int64(block_id)
end

function component_offset_from_string(name::String)
    offset = 0
    if name == "x"
        offset = 1
    elseif name == "y"
        offset = 2
    elseif name == "z"
        offset = 3
    else
        norma_abort("invalid component name $name")
    end
    return offset
end

function extract_value(value::Real)
    return Float64(value)
end

function extract_value(symbol::Num)
    return Float64(symbol.val)
end

function create_bcs(params::Parameters)
    boundary_conditions = Vector{BoundaryCondition}()
    if haskey(params, "boundary conditions") == false
        return boundary_conditions
    end
    input_mesh = params["input_mesh"]
    bc_params = params["boundary conditions"]
    inclined_support_nodes = Vector{Int64}()
    for (bc_type, bc_type_params) in bc_params
        for bc_setting_params in bc_type_params
            if bc_type == "Dirichlet"
                boundary_condition = SolidMechanicsDirichletBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann"
                boundary_condition = SolidMechanicsNeumannBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann pressure"
                boundary_condition = SolidMechanicsNeumannPressureBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz contact"
                sim = params["parent_simulation"]
                sim.controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SolidMechanicsContactSchwarzBoundaryCondition(
                    coupled_subsim, input_mesh, bc_setting_params
                )
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Inclined Dirichlet"
                boundary_condition = SolidMechanicsInclinedDirichletBoundaryCondition(input_mesh, bc_setting_params)
                append!(inclined_support_nodes, boundary_condition.node_set_node_indices)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz overlap" || bc_type == "Schwarz nonoverlap"
                sim = params["parent_simulation"]
                subsim_name = params["name"]
                subdomain_index = sim.subsim_name_index_map[subsim_name]
                subsim = sim.subsims[subdomain_index]
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            else
                norma_abort("Unknown boundary condition type : $bc_type")
            end
        end
    end
    # BRP: do not support applying multiple inclined support BCs to a single node
    duplicate_inclined_support_conditions = length(unique(inclined_support_nodes)) < length(inclined_support_nodes)
    if duplicate_inclined_support_conditions
        norma_abort("Cannot apply multiple inclined BCs to a single node.")
    end
    return boundary_conditions
end

function apply_bcs(model::SolidMechanics)
    model.boundary_force .= 0.0
    model.free_dofs .= true
    for boundary_condition in model.boundary_conditions
        apply_bc(model, boundary_condition)
    end
end


function assign_velocity!(
    velocity::Matrix{Float64}, offset::Int64, node_index::Int32, velo_val::Float64, context::String
)
    current_val = velocity[offset, node_index]
    velocity_already_defined = !(current_val ≈ 0.0)
    dissimilar_velocities = !(current_val ≈ velo_val)
    if velocity_already_defined && dissimilar_velocities
        norma_abortf(
            "Inconsistent velocity initial conditions for node %d: " *
            "attempted to assign velocity %s (v = (%.4e, %.4e, %.4e)), " *
            "which conflicts with an already assigned value (v = (%.4e, %.4e, %.4e)).",
            node_index,
            context,
            velo_val[1],
            velo_val[2],
            velo_val[3],
            current_val[1],
            current_val[2],
            current_val[3],
        )
    else
        velocity[offset, node_index] = velo_val
    end
    return nothing
end

function apply_ics(params::Parameters, model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    if haskey(params, "initial conditions") == false
        return nothing
    end
    input_mesh = params["input_mesh"]
    ic_params = params["initial conditions"]
    for (ic_type, ic_type_params) in ic_params
        for ic in ic_type_params
            node_set_name = ic["node set"]
            expression = ic["function"]
            component = ic["component"]
            offset = component_offset_from_string(component)
            node_set_id = node_set_id_from_name(node_set_name, input_mesh)
            node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
            # expression is an arbitrary function of t, x, y, z in the input file
            # Parse and expand derivatives only once
            if ic_type == "displacement"
                disp_num = eval(Meta.parse(expression))
                velo_num = expand_derivatives(D(disp_num))
            elseif ic_type == "velocity"
                disp_num = nothing  # Not needed
                velo_num = eval(Meta.parse(expression))
            else
                norma_abort(
                    "Invalid initial condition type: '$ic_type'. Supported types are: displacement or velocity."
                )
            end
            for node_index in node_set_node_indices
                values = Dict(
                    t => model.time,
                    x => model.reference[1, node_index],
                    y => model.reference[2, node_index],
                    z => model.reference[3, node_index],
                )
                disp_val = disp_num === nothing ? 0.0 : extract_value(substitute(disp_num, values))
                velo_val = velo_num === nothing ? 0.0 : extract_value(substitute(velo_num, values))
                if ic_type == "displacement"
                    model.current[offset, node_index] = model.reference[offset, node_index] + disp_val
                    non_zero_velocity = !(velo_val ≈ 0.0)
                    if non_zero_velocity
                        assign_velocity!(model.velocity, offset, node_index, velo_val, "derived from displacement")
                    end
                end
                if ic_type == "velocity"
                    assign_velocity!(model.velocity, offset, node_index, velo_val, "directly from velocity IC")
                end
            end
        end
    end
    copy_solution_source_targets(model, integrator, solver)
    return nothing
end


function pair_schwarz_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        model = subsim.model
        bcs = model.boundary_conditions
        for (bc_index, bc) in enumerate(bcs)
            pair_bc(bc, bc_index)
        end
    end
end

function pair_bc(_::SolidMechanicsRegularBoundaryCondition, _::Int64) end

function pair_bc(bc::SolidMechanicsSchwarzBoundaryCondition, bc_index::Int64)
    if bc isa SolidMechanicsOverlapSchwarzBoundaryCondition
        return nothing
    end
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = bc.coupled_subsim.model
    coupled_bcs = coupled_model.boundary_conditions
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name
            coupled_bc.is_dirichlet = !bc.is_dirichlet
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
            if coupled_bc.is_dirichlet == false
                coupled_model.inclined_support = false
            end
        end
    end
    return nothing
end
