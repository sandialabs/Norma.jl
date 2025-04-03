# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

@variables t, x, y, z
D = Differential(t)

function SMDirichletBC(input_mesh::ExodusDatabase, bc_params::Parameters)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))
    return SMDirichletBC(node_set_name, offset, node_set_id, node_set_node_indices, disp_num, velo_num, acce_num)
end

function SMDirichletInclined(input_mesh::ExodusDatabase, bc_params::Parameters)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    if isa(expression, AbstractVector)
        if (length(expression) != 3)
            error("Vectorized function must have 3 elements.")
        end
        if all(x -> x isa String, expression) == false
            error("All functions must be strings (including zeros).")
        end
        disp_expression = [
            eval(Meta.parse(expression[1])), eval(Meta.parse(expression[2])), eval(Meta.parse(expression[3]))
        ]
        off_axis_free = false
    else
        disp_expression = [eval(Meta.parse(expression)), eval(Meta.parse("0.0")), eval(Meta.parse("0.0"))]
        off_axis_free = true
    end

    velo_expression = [
        expand_derivatives(D(disp_expression[1])),
        expand_derivatives(D(disp_expression[2])),
        expand_derivatives(D(disp_expression[3])),
    ]
    acce_expression = [
        expand_derivatives(D(velo_expression[1])),
        expand_derivatives(D(velo_expression[2])),
        expand_derivatives(D(velo_expression[3])),
    ]

    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file

    reference_normal = bc_params["normal vector"]

    if all(x -> x isa String, reference_normal) == false
        # We'll cast the normal into a string expression
        reference_normal = [string(vc) for vc in reference_normal]
    end

    reference_normal = [
        eval(Meta.parse(reference_normal[1])),
        eval(Meta.parse(reference_normal[2])),
        eval(Meta.parse(reference_normal[3])),
    ]
    return SMDirichletInclined(
        node_set_name,
        node_set_id,
        node_set_node_indices,
        disp_expression,
        velo_expression,
        acce_expression,
        reference_normal,
        off_axis_free,
    )
end

function SMNeumannBC(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    traction_num = eval(Meta.parse(expression))
    return SMNeumannBC(side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, traction_num)
end

function SMContactSchwarzBC(coupled_subsim::SingleDomainSimulation, input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    _, num_nodes_per_side, side_set_node_indices = get_side_set_local_from_global_map(input_mesh, side_set_id)
    coupled_block_name = bc_params["source block"]
    coupled_bc_index = 0
    coupled_mesh = coupled_subsim.model.mesh
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    coupled_side_set_name = bc_params["source side set"]
    coupled_side_set_id = side_set_id_from_name(coupled_side_set_name, coupled_mesh)
    is_dirichlet = true
    transfer_operator = Matrix{Float64}(undef, 0, 0)
    rotation_matrix = I(3)
    active_contact = false
    swap_bcs = false
    if haskey(bc_params, "swap BC types") == true
        swap_bcs = bc_params["swap BC types"]
    end

    friction_type_string = bc_params["friction type"]
    if friction_type_string == "frictionless"
        friction_type = 0
    elseif friction_type_string == "tied"
        friction_type = 1
    else
        error("Unknown or not implemented friction type : ", friction_type_string)
    end

    return SMContactSchwarzBC(
        side_set_name,
        side_set_id,
        num_nodes_per_side,
        side_set_node_indices,
        coupled_subsim,
        coupled_bc_index,
        coupled_block_id,
        coupled_side_set_id,
        is_dirichlet,
        transfer_operator,
        rotation_matrix,
        active_contact,
        swap_bcs,
        friction_type,
    )
end

function SMNonOverlapSchwarzBC(
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    coupled_nodes_indices::Vector{Vector{Int64}},
    interpolation_function_values::Vector{Vector{Float64}},
    coupled_subsim::Simulation,
    subsim::Simulation,
    coupled_side_set_id::Int64,
    is_dirichlet::Bool,
    swap_bcs::Bool,
)
    transfer_operator = Matrix{Float64}(undef, 0, 0)
    return SMNonOverlapSchwarzBC(
        side_set_id,
        side_set_node_indices,
        coupled_nodes_indices,
        interpolation_function_values,
        coupled_subsim,
        subsim,
        coupled_side_set_id,
        is_dirichlet,
        swap_bcs,
        transfer_operator,
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
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    _, _, side_set_node_indices = get_side_set_local_from_global_map(input_mesh, side_set_id)
    coupled_block_name = bc_params["source block"]
    if typeof(coupled_subsim.model) <: RomModel
        coupled_mesh = coupled_subsim.model.fom_model.mesh
    else
        coupled_mesh = coupled_subsim.model.mesh
    end
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type = Exodus.read_block_parameters(coupled_mesh, coupled_block_id)[1]
    coupled_side_set_name = bc_params["source side set"]
    coupled_side_set_id = side_set_id_from_name(coupled_side_set_name, coupled_mesh)
    coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    tol = 1.0e-06
    if haskey(bc_params, "search tolerance") == true
        tol = bc_params["search tolerance"]
    end
    side_set_node_indices = unique(side_set_node_indices)
    for node_index in side_set_node_indices
        point = subsim.model.reference[:, node_index]
        node_indices, ξ, found = find_point_in_mesh(point, coupled_subsim.model, coupled_block_id, tol)
        if found == false
            error("Could not find subdomain ", subsim.name, " point ", point, " in subdomain ", coupled_subsim.name)
        end
        N = interpolate(element_type, ξ)[1]
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    is_dirichlet = true
    swap_bcs = false
    if bc_type == "Schwarz overlap"
        SMOverlapSchwarzBC(
            side_set_name,
            side_set_node_indices,
            coupled_nodes_indices,
            interpolation_function_values,
            coupled_subsim,
            subsim,
            is_dirichlet,
            swap_bcs,
        )
    elseif bc_type == "Schwarz nonoverlap"
        if haskey(bc_params, "default BC type") == true
            default_bc_type = bc_params["default BC type"]
            if default_bc_type == "Dirichlet"
                is_dirichlet = true
            elseif default_bc_type == "Neumann"
                is_dirichlet = false
            else
                error("Invalid string for 'default BC type'!  Valid options are 'Dirichlet' and 'Neumann'")
            end
        end
        if haskey(bc_params, "swap BC types") == true
            swap_bcs = bc_params["swap BC types"]
        end
        SMNonOverlapSchwarzBC(
            side_set_id,
            side_set_node_indices,
            coupled_nodes_indices,
            interpolation_function_values,
            coupled_subsim,
            subsim,
            coupled_side_set_id,
            is_dirichlet,
            swap_bcs,
        )
    else
        error("Unknown boundary condition type : ", bc_type)
    end
end

function apply_bc(model::OpInfModel, bc::SMDirichletBC)
    model.fom_model.time = model.time
    apply_bc(model.fom_model, bc)
    bc_vector = zeros(0)
    for node_index in bc.node_set_node_indices
        disp_val = model.fom_model.current[bc.offset, node_index] - model.fom_model.reference[bc.offset, node_index]
        push!(bc_vector, disp_val)
    end

    offset = bc.offset
    if offset == 1
        offset_name = "x"
    end
    if offset == 2
        offset_name = "y"
    end
    if offset == 3
        offset_name = "z"
    end

    op_name = "B_" * bc.node_set_name * "-" * offset_name
    bc_operator = model.opinf_rom[op_name]
    # SM Dirichlet BC are only defined on a single x,y,z
    return model.reduced_boundary_forcing[:] += bc_operator[1, :, :] * bc_vector
end

function apply_bc(model::SolidMechanics, bc::SMDirichletBC)
    for node_index in bc.node_set_node_indices
        values = Dict(
            t => model.time,
            x => model.reference[1, node_index],
            y => model.reference[2, node_index],
            z => model.reference[3, node_index],
        )
        disp_sym = substitute(bc.disp_num, values)
        velo_sym = substitute(bc.velo_num, values)
        acce_sym = substitute(bc.acce_num, values)
        disp_val = extract_value(disp_sym)
        velo_val = extract_value(velo_sym)
        acce_val = extract_value(acce_sym)
        dof_index = 3 * (node_index - 1) + bc.offset
        model.current[bc.offset, node_index] = model.reference[bc.offset, node_index] + disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SMNeumannBC)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc.traction_num, model.time)
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

function apply_bc(model::SolidMechanics, bc::SMDirichletInclined)
    e1 = [1.0, 0.0, 0.0]

    for node_index in bc.node_set_node_indices
        values = Dict(
            t => model.time,
            x => model.reference[1, node_index],
            y => model.reference[2, node_index],
            z => model.reference[3, node_index],
        )
        # The local basis is determined from a normal vector
        axis = bc.reference_normal
        axis = [extract_value(substitute(av, values)) for av in axis]

        axis = axis / norm(axis)
        w = cross(axis, e1)
        s = norm(w)
        angle_btwn = acos(dot(axis, e1))

        if (angle_btwn ≈ 0.0)
            rotation_matrix = I(3) * 1.0
        elseif (angle_btwn ≈ π)
            rotation_matrix = I(3) * 1.0
            rotation_matrix[1, 1] = -1.0
            rotation_matrix[2, 2] = -1.0
        else
            if angle_btwn > π / 2
                θ = π - asin(s)
            else
                θ = asin(s)
            end

            m = w / s
            rv = θ * m
            # Rotation is converted via the pseudo vector to rotation matrix
            rotation_matrix = MiniTensor.rt_from_rv(rv)
        end

        disp_val_loc = [extract_value(substitute(exp, values)) for exp in bc.disp_expression]
        velo_val_loc = [extract_value(substitute(exp, values)) for exp in bc.velo_expression]
        acce_val_loc = [extract_value(substitute(exp, values)) for exp in bc.acce_expression]

        # Inherit free values from model (required for explicit dynamics)
        original_disp = model.current[:, node_index] - model.reference[:, node_index]
        original_acceleration = model.acceleration[:, node_index]
        original_velocity = model.velocity[:, node_index]
        local_original_displacement = rotation_matrix * original_disp
        local_original_acceleration = rotation_matrix * original_acceleration
        local_original_velocity = rotation_matrix * original_velocity

        # Inclined support is only applied in local X
        if bc.off_axis_free == true
            model.free_dofs[3 * node_index - 2] = false
            model.free_dofs[3 * node_index - 1] = true
            model.free_dofs[3 * node_index] = true
            # Set local X direction to inclined support
            local_original_displacement[1] = disp_val_loc[1]
            local_original_acceleration[1] = acce_val_loc[1]
            local_original_velocity[1] = velo_val_loc[1]
        else
            model.free_dofs[3 * node_index - 2] = false
            model.free_dofs[3 * node_index - 1] = false
            model.free_dofs[3 * node_index] = false
            # Set local X direction to inclined support
            local_original_displacement = disp_val_loc
            local_original_acceleration = acce_val_loc
            local_original_velocity = velo_val_loc
        end

        original_velocity = rotation_matrix' * local_original_velocity
        original_acceleration = rotation_matrix' * local_original_acceleration
        disp_vector_glob = rotation_matrix' * local_original_displacement

        model.current[:, node_index] = model.reference[:, node_index] + disp_vector_glob
        model.velocity[:, node_index] = original_velocity
        model.acceleration[:, node_index] = original_acceleration
        global_base = 3 * (node_index - 1) # Block index in global stiffness
        model.global_transform[(global_base + 1):(global_base + 3), (global_base + 1):(global_base + 3)] =
            rotation_matrix
    end
end

function find_point_in_mesh(point::Vector{Float64}, model::SolidMechanics, blk_id::Int, tol::Float64)
    mesh = model.mesh
    element_type = Exodus.read_block_parameters(mesh, Int32(blk_id))[1]
    elem_blk_conn = get_block_connectivity(mesh, blk_id)
    num_blk_elems, num_elem_nodes = size(elem_blk_conn)
    node_indices = Vector{Int64}()
    found = false
    ξ = zeros(length(point))
    for blk_elem_index in 1:num_blk_elems
        conn_indices = ((blk_elem_index - 1) * num_elem_nodes + 1):(blk_elem_index * num_elem_nodes)
        node_indices = elem_blk_conn[conn_indices]
        elem_ref_pos = model.reference[:, node_indices]
        ξ, found = is_inside(element_type, elem_ref_pos, point, tol)
        if found == true
            break
        end
    end
    return node_indices, ξ, found
end

function find_point_in_mesh(point::Vector{Float64}, model::RomModel, blk_id::Int, tol::Float64)
    node_indices, ξ, found = find_point_in_mesh(point, model.fom_model, blk_id, tol)
    return node_indices, ξ, found
end

function apply_bc_detail(model::SolidMechanics, bc::SMContactSchwarzBC)
    if bc.is_dirichlet == true
        apply_sm_schwarz_contact_dirichlet(model, bc)
    else
        apply_sm_schwarz_contact_neumann(model, bc)
    end
end

function apply_bc_detail(model::SolidMechanics, bc::CouplingSchwarzBoundaryCondition)
    if bc.is_dirichlet == true
        apply_sm_schwarz_coupling_dirichlet(model, bc)
    else
        apply_sm_schwarz_coupling_neumann(model, bc)
    end
end

function apply_bc_detail(model::OpInfModel, bc::CouplingSchwarzBoundaryCondition)
    if (typeof(bc.coupled_subsim.model) == SolidMechanics)
        ## Apply BC to the FOM vector
        apply_bc_detail(model.fom_model, bc)

        # populate our own BC vector
        bc_vector = zeros(3, length(bc.side_set_node_indices))
        for i in 1:length(bc.side_set_node_indices)
            node_index = bc.side_set_node_indices[i]
            bc_vector[:, i] = model.fom_model.current[:, node_index] - model.fom_model.reference[:, node_index]
        end
        op_name = "B_" * bc.side_set_name
        bc_operator = model.opinf_rom[op_name]
        for i in 1:3
            model.reduced_boundary_forcing[:] += bc_operator[i, :, :] * bc_vector[i, :]
        end
    else
        throw("ROM-ROM coupling not supported yet")
    end
end

function apply_sm_schwarz_coupling_dirichlet(model::SolidMechanics, bc::CouplingSchwarzBoundaryCondition)
    if (typeof(bc.coupled_subsim.model) == SolidMechanics)
        for i in 1:length(bc.side_set_node_indices)
            node_index = bc.side_set_node_indices[i]
            coupled_node_indices = bc.coupled_nodes_indices[i]
            N = bc.interpolation_function_values[i]
            elem_posn = bc.coupled_subsim.model.current[:, coupled_node_indices]
            elem_velo = bc.coupled_subsim.model.velocity[:, coupled_node_indices]
            elem_acce = bc.coupled_subsim.model.acceleration[:, coupled_node_indices]
            point_posn = elem_posn * N
            point_velo = elem_velo * N
            point_acce = elem_acce * N
            model.current[:, node_index] = point_posn
            model.velocity[:, node_index] = point_velo
            model.acceleration[:, node_index] = point_acce
            dof_index = [3 * node_index - 2, 3 * node_index - 1, 3 * node_index]
            model.free_dofs[dof_index] .= false
        end
    elseif (typeof(bc.coupled_subsim.model) <: RomModel)
        for i in 1:length(bc.side_set_node_indices)
            node_index = bc.side_set_node_indices[i]
            coupled_node_indices = bc.coupled_nodes_indices[i]
            N = bc.interpolation_function_values[i]
            elem_posn = bc.coupled_subsim.model.fom_model.current[:, coupled_node_indices]
            elem_velo = bc.coupled_subsim.model.fom_model.velocity[:, coupled_node_indices]
            elem_acce = bc.coupled_subsim.model.fom_model.acceleration[:, coupled_node_indices]
            point_posn = elem_posn * N
            point_velo = elem_velo * N
            point_acce = elem_acce * N
            model.current[:, node_index] = point_posn
            model.velocity[:, node_index] = point_velo
            model.acceleration[:, node_index] = point_acce
            dof_index = [3 * node_index - 2, 3 * node_index - 1, 3 * node_index]
            model.free_dofs[dof_index] .= false
        end
    end
end

function apply_sm_schwarz_coupling_neumann(model::SolidMechanics, bc::CouplingSchwarzBoundaryCondition)
    schwarz_tractions = get_dst_traction(bc)
    global_from_local_map = get_side_set_global_from_local_map(model.mesh, bc.side_set_id)
    num_local_nodes = length(global_from_local_map)
    for local_node in 1:num_local_nodes
        global_node = global_from_local_map[local_node]
        node_tractions = schwarz_tractions[:, local_node]
        model.boundary_force[(3 * global_node - 2):(3 * global_node)] += node_tractions
    end
end

function apply_bc(model::SolidMechanics, bc::SchwarzBoundaryCondition)
    global_sim = bc.coupled_subsim.params["global_simulation"]
    schwarz_controller = global_sim.schwarz_controller
    if typeof(bc) == SMContactSchwarzBC && schwarz_controller.active_contact == false
        return nothing
    end
    empty_history = length(global_sim.schwarz_controller.time_hist) == 0
    same_step = schwarz_controller.same_step == true
    if empty_history == true
        apply_bc_detail(model, bc)
        return nothing
    end
    # Save solution of coupled simulation
    saved_disp = bc.coupled_subsim.integrator.displacement
    saved_velo = bc.coupled_subsim.integrator.velocity
    saved_acce = bc.coupled_subsim.integrator.acceleration
    saved_∂Ω_f = bc.coupled_subsim.model.internal_force
    time = model.time
    coupled_name = bc.coupled_subsim.name
    coupled_index = global_sim.subsim_name_index_map[coupled_name]
    time_hist = global_sim.schwarz_controller.time_hist[coupled_index]
    disp_hist = global_sim.schwarz_controller.disp_hist[coupled_index]
    velo_hist = global_sim.schwarz_controller.velo_hist[coupled_index]
    acce_hist = global_sim.schwarz_controller.acce_hist[coupled_index]
    ∂Ω_f_hist = global_sim.schwarz_controller.∂Ω_f_hist[coupled_index]
    interp_disp = same_step == true ? disp_hist[end] : interpolate(time_hist, disp_hist, time)
    interp_velo = same_step == true ? velo_hist[end] : interpolate(time_hist, velo_hist, time)
    interp_acce = same_step == true ? acce_hist[end] : interpolate(time_hist, acce_hist, time)
    interp_∂Ω_f = same_step == true ? ∂Ω_f_hist[end] : interpolate(time_hist, ∂Ω_f_hist, time)
    bc.coupled_subsim.model.internal_force = interp_∂Ω_f
    if typeof(bc) == SMContactSchwarzBC || typeof(bc) == SMNonOverlapSchwarzBC
        relaxation_parameter = global_sim.schwarz_controller.relaxation_parameter
        schwarz_iteration = global_sim.schwarz_controller.iteration_number
        if schwarz_iteration == 1
            lambda_dispᵖʳᵉᵛ = zeros(length(interp_disp))
            lambda_veloᵖʳᵉᵛ = zeros(length(interp_velo))
            lambda_acceᵖʳᵉᵛ = zeros(length(interp_acce))
        else
            lambda_dispᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_disp[coupled_index]
            lambda_veloᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_velo[coupled_index]
            lambda_acceᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_acce[coupled_index]
        end
        bc.coupled_subsim.integrator.displacement =
            global_sim.schwarz_controller.lambda_disp[coupled_index] =
                relaxation_parameter * interp_disp + (1 - relaxation_parameter) * lambda_dispᵖʳᵉᵛ
        bc.coupled_subsim.integrator.velocity =
            global_sim.schwarz_controller.lambda_velo[coupled_index] =
                relaxation_parameter * interp_velo + (1 - relaxation_parameter) * lambda_veloᵖʳᵉᵛ
        bc.coupled_subsim.integrator.acceleration =
            global_sim.schwarz_controller.lambda_acce[coupled_index] =
                relaxation_parameter * interp_acce + (1 - relaxation_parameter) * lambda_acceᵖʳᵉᵛ
    else
        bc.coupled_subsim.integrator.displacement = interp_disp
        bc.coupled_subsim.integrator.velocity = interp_velo
        bc.coupled_subsim.integrator.acceleration = interp_acce
    end
    copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
    apply_bc_detail(model, bc)
    bc.coupled_subsim.integrator.displacement = saved_disp
    bc.coupled_subsim.integrator.velocity = saved_velo
    bc.coupled_subsim.integrator.acceleration = saved_acce
    bc.coupled_subsim.model.internal_force = saved_∂Ω_f
    return copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
end

function apply_bc(model::RomModel, bc::SchwarzBoundaryCondition)
    global_sim = bc.coupled_subsim.params["global_simulation"]
    schwarz_controller = global_sim.schwarz_controller
    if typeof(bc) == SMContactSchwarzBC && schwarz_controller.active_contact == false
        return nothing
    end
    empty_history = length(global_sim.schwarz_controller.time_hist) == 0
    same_step = schwarz_controller.same_step == true
    if empty_history == true
        apply_bc_detail(model, bc)
        return nothing
    end
    # Save solution of coupled simulation
    saved_disp = bc.coupled_subsim.integrator.displacement
    saved_velo = bc.coupled_subsim.integrator.velocity
    saved_acce = bc.coupled_subsim.integrator.acceleration
    saved_∂Ω_f = bc.coupled_subsim.model.internal_force
    time = model.time
    coupled_name = bc.coupled_subsim.name
    coupled_index = global_sim.subsim_name_index_map[coupled_name]
    time_hist = global_sim.schwarz_controller.time_hist[coupled_index]
    disp_hist = global_sim.schwarz_controller.disp_hist[coupled_index]
    velo_hist = global_sim.schwarz_controller.velo_hist[coupled_index]
    acce_hist = global_sim.schwarz_controller.acce_hist[coupled_index]
    ∂Ω_f_hist = global_sim.schwarz_controller.∂Ω_f_hist[coupled_index]
    interp_disp = same_step == true ? disp_hist[end] : interpolate(time_hist, disp_hist, time)
    interp_velo = same_step == true ? velo_hist[end] : interpolate(time_hist, velo_hist, time)
    interp_acce = same_step == true ? acce_hist[end] : interpolate(time_hist, acce_hist, time)
    interp_∂Ω_f = same_step == true ? ∂Ω_f_hist[end] : interpolate(time_hist, ∂Ω_f_hist, time)
    if typeof(bc.coupled_subsim.model) == SolidMechanics
        bc.coupled_subsim.model.internal_force = interp_∂Ω_f
    elseif typeof(bc.coupled_subsim.model) <: RomModel
        bc.coupled_subsim.model.fom_model.internal_force = interp_∂Ω_f
    end

    if typeof(bc) == SMContactSchwarzBC || typeof(bc) == SMNonOverlapSchwarzBC
        relaxation_parameter = global_sim.schwarz_controller.relaxation_parameter
        schwarz_iteration = global_sim.schwarz_controller.iteration_number
        if schwarz_iteration == 1
            lambda_dispᵖʳᵉᵛ = zeros(length(interp_disp))
            lambda_veloᵖʳᵉᵛ = zeros(length(interp_velo))
            lambda_acceᵖʳᵉᵛ = zeros(length(interp_acce))
        else
            lambda_dispᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_disp[coupled_index]
            lambda_veloᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_velo[coupled_index]
            lambda_acceᵖʳᵉᵛ = global_sim.schwarz_controller.lambda_acce[coupled_index]
        end
        bc.coupled_subsim.integrator.displacement =
            global_sim.schwarz_controller.lambda_disp[coupled_index] =
                relaxation_parameter * interp_disp + (1 - relaxation_parameter) * lambda_dispᵖʳᵉᵛ
        bc.coupled_subsim.integrator.velocity =
            global_sim.schwarz_controller.lambda_velo[coupled_index] =
                relaxation_parameter * interp_velo + (1 - relaxation_parameter) * lambda_veloᵖʳᵉᵛ
        bc.coupled_subsim.integrator.acceleration =
            global_sim.schwarz_controller.lambda_acce[coupled_index] =
                relaxation_parameter * interp_acce + (1 - relaxation_parameter) * lambda_acceᵖʳᵉᵛ
    else
        bc.coupled_subsim.integrator.displacement = interp_disp
        bc.coupled_subsim.integrator.velocity = interp_velo
        bc.coupled_subsim.integrator.acceleration = interp_acce
    end
    # Copies from integrator to model
    copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
    apply_bc_detail(model, bc)
    bc.coupled_subsim.integrator.displacement = saved_disp
    bc.coupled_subsim.integrator.velocity = saved_velo
    bc.coupled_subsim.integrator.acceleration = saved_acce
    bc.coupled_subsim.model.internal_force = saved_∂Ω_f
    # Copy from integrator to model
    return copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
end

function transfer_normal_component(source::Vector{Float64}, target::Vector{Float64}, normal::Vector{Float64})
    normal_projection = normal * normal'
    tangent_projection = I(length(normal)) - normal_projection
    return tangent_projection * target + normal_projection * source
end

function apply_sm_schwarz_contact_dirichlet(model::SolidMechanics, bc::SMContactSchwarzBC)
    side_set_node_indices = unique(bc.side_set_node_indices)
    for node_index in side_set_node_indices
        point = model.current[:, node_index]
        new_point, ξ, _, closest_face_node_indices, closest_normal, _ = project_point_to_side_set(
            point, bc.coupled_subsim.model, bc.coupled_side_set_id
        )

        # The local basis comes from the closest_normal:
        axis = -closest_normal / norm(-closest_normal)

        e1 = [1.0, 0.0, 0.0]
        angle_btwn = acos(dot(axis, e1))
        w = cross(axis, e1)
        s = norm(w)
        if (angle_btwn ≈ 0.0)
            rotation_matrix = I(3) * 1.0
        elseif (angle_btwn ≈ π)
            rotation_matrix = I(3) * 1.0
            rotation_matrix[1, 1] = -1.0
            rotation_matrix[2, 2] = -1.0
        else
            if angle_btwn > π / 2
                θ = π - asin(s)
            else
                θ = asin(s)
            end
            m = w / s
            rv = θ * m
            # Rotation is converted via the pseudo vector to rotation matrix
            bc.rotation_matrix = MiniTensor.rt_from_rv(rv)
        end

        model.current[:, node_index] = new_point
        num_nodes = length(closest_face_node_indices)
        element_type = get_element_type(2, num_nodes)
        N, _, _ = interpolate(element_type, ξ)
        source_velo = bc.coupled_subsim.model.velocity[:, closest_face_node_indices] * N
        source_acce = bc.coupled_subsim.model.acceleration[:, closest_face_node_indices] * N
        model.free_dofs[[3 * node_index - 2]] .= false
        model.free_dofs[[3 * node_index - 1]] .= true
        model.free_dofs[[3 * node_index]] .= true
        if bc.friction_type == 0
            model.velocity[:, node_index] = transfer_normal_component(
                source_velo, model.velocity[:, node_index], closest_normal
            )
            model.acceleration[:, node_index] = transfer_normal_component(
                source_acce, model.acceleration[:, node_index], closest_normal
            )
        elseif bc.friction_type == 1
            model.velocity[:, node_index] = source_velo
            model.acceleration[:, node_index] = source_acce
        else
            error("Unknown or not implemented friction type.")
        end
        global_base = 3 * (node_index - 1) # Block index in global stiffness
        model.global_transform[(global_base + 1):(global_base + 3), (global_base + 1):(global_base + 3)] =
            bc.rotation_matrix
    end
end

function apply_naive_stabilized_bcs(subsim::SingleDomainSimulation)
    bcs = subsim.model.boundary_conditions
    for bc in bcs
        if typeof(bc) == SMContactSchwarzBC
            side_set_node_indices = unique(bc.side_set_node_indices)
            for node_index in side_set_node_indices
                subsim.model.acceleration[:, node_index] = zeros(3)
            end
        end
    end
    return copy_solution_source_targets(subsim.model, subsim.integrator, subsim.solver)
end

function apply_sm_schwarz_contact_neumann(model::SolidMechanics, bc::SMContactSchwarzBC)
    schwarz_tractions = get_dst_traction(bc)
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    global_from_local_map = get_side_set_global_from_local_map(model.mesh, bc.side_set_id)
    num_local_nodes = length(global_from_local_map)
    for local_node in 1:num_local_nodes
        global_node = global_from_local_map[local_node]
        node_tractions = schwarz_tractions[:, local_node]
        normal = normals[:, local_node]
        if bc.friction_type == 0
            model.boundary_force[(3 * global_node - 2):(3 * global_node)] += transfer_normal_component(
                node_tractions, model.boundary_force[(3 * global_node - 2):(3 * global_node)], normal
            )
        elseif bc.friction_type == 1
            model.boundary_force[(3 * global_node - 2):(3 * global_node)] += node_tractions
        else
            error("Unknown or not implemented friction type.")
        end
    end
end

function local_traction_from_global_force(mesh::ExodusDatabase, side_set_id::Integer, global_force::Vector{Float64})
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    num_local_nodes = length(global_from_local_map)
    local_traction = zeros(3, num_local_nodes)
    for local_node in 1:num_local_nodes
        global_node = global_from_local_map[local_node]
        local_traction[:, local_node] = global_force[(3 * global_node - 2):(3 * global_node)]
    end
    return local_traction
end

function update_transfer_operator(dst_model::SolidMechanics, dst_bc::SchwarzBoundaryCondition)
    src_side_set_id = dst_bc.coupled_side_set_id
    src_model = dst_bc.coupled_subsim.model
    dst_side_set_id = dst_bc.side_set_id
    square_projection_matrix = get_square_projection_matrix(src_model, src_side_set_id)
    rectangular_projection_matrix = get_rectangular_projection_matrix(
        src_model, src_side_set_id, dst_model, dst_side_set_id
    )
    return dst_bc.transfer_operator = rectangular_projection_matrix * (square_projection_matrix \ I)
end

function get_dst_traction(dst_bc::SchwarzBoundaryCondition)
    src_mesh = dst_bc.coupled_subsim.model.mesh
    src_side_set_id = dst_bc.coupled_side_set_id
    src_global_force = -dst_bc.coupled_subsim.model.internal_force
    src_local_traction = local_traction_from_global_force(src_mesh, src_side_set_id, src_global_force)
    num_dst_nodes = size(dst_bc.transfer_operator, 1)
    dst_traction = zeros(3, num_dst_nodes)
    dst_traction[1, :] = dst_bc.transfer_operator * src_local_traction[1, :]
    dst_traction[2, :] = dst_bc.transfer_operator * src_local_traction[2, :]
    dst_traction[3, :] = dst_bc.transfer_operator * src_local_traction[3, :]
    return dst_traction
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
        error("node set ", node_set_name, " cannot be found in mesh")
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
        error("side set ", side_set_name, " cannot be found in mesh")
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
        error("block ", block_name, " cannot be found in mesh")
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
        error("invalid component name ", name)
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
                boundary_condition = SMDirichletBC(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann"
                boundary_condition = SMNeumannBC(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz contact"
                sim = params["global_simulation"]
                sim.schwarz_controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMContactSchwarzBC(coupled_subsim, input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Inclined Dirichlet"
                boundary_condition = SMDirichletInclined(input_mesh, bc_setting_params)
                append!(inclined_support_nodes, boundary_condition.node_set_node_indices)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz overlap" || bc_type == "Schwarz nonoverlap"
                sim = params["global_simulation"]
                subsim_name = params["name"]
                subdomain_index = sim.subsim_name_index_map[subsim_name]
                subsim = sim.subsims[subdomain_index]
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            else
                error("Unknown boundary condition type : ", bc_type)
            end
        end
    end
    # BRP: do not support applying multiple inclined support BCs to a single node
    duplicate_inclined_support_conditions = length(unique(inclined_support_nodes)) < length(inclined_support_nodes)
    if duplicate_inclined_support_conditions
        throw(error("Cannot apply multiple inclined BCs to a single node."))
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

function apply_bcs(model::RomModel)
    model.reduced_boundary_forcing[:] .= 0.0
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
        error(
            "Inconsistent velocity initial conditions (ICs) for node ",
            node_index,
            ": attempted to assign velocity ",
            context,
            " (v = ",
            velo_val,
            ")",
            " which conflicts with an already assigned value (v = ",
            current_val,
            ").",
        )
    else
        velocity[offset, node_index] = velo_val
    end
    return nothing
end

function apply_ics(params::Parameters, model::SolidMechanics)
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
                error("Invalid initial condition type: '$ic_type'. Supported types are: displacement or velocity.")
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
end

function apply_ics(params::Parameters, model::RomModel)
    apply_ics(params, model.fom_model)

    if haskey(params, "initial conditions") == false
        return nothing
    end
    n_var, n_node, n_mode = model.basis.size
    n_var_fom, n_node_fom = size(model.fom_model.current)

    # Make sure basis is the right size
    if n_var != n_var_fom || n_node != n_node_fom
        throw("Basis is wrong size")
    end

    # project onto basis
    for k in 1:n_mode
        model.reduced_state[k] = 0.0
        model.reduced_velocity[k] = 0.0
        for j in 1:n_node
            for n in 1:n_var
                model.reduced_state[k] +=
                    model.basis[n, j, k] *
                    (model.fom_model.current[n, j] - model.fom_model.reference[n, j])
                model.reduced_velocity[k] +=
                    model.basis[n, j, k] *
                    (model.fom_model.velocity[n, j])
            end
        end
    end
end

function pair_schwarz_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        model = subsim.model
        name = subsim.name
        bcs = model.boundary_conditions
        for bc in bcs
            pair_bc(name, bc)
        end
    end
end

function pair_bc(_::String, _::RegularBoundaryCondition) end

function pair_bc(name::String, bc::ContactSchwarzBoundaryCondition)
    coupled_model = bc.coupled_subsim.model
    coupled_bcs = coupled_model.boundary_conditions
    for coupled_bc in coupled_bcs
        if is_coupled_to_current(name, coupled_bc) == true
            coupled_bc.is_dirichlet = !bc.is_dirichlet
            if coupled_bc.is_dirichlet == false
                coupled_model.inclined_support = false
            end
        end
    end
end

function pair_bc(name::String, bc::CouplingSchwarzBoundaryCondition)
    if typeof(bc) == SMNonOverlapSchwarzBC
        coupled_model = bc.coupled_subsim.model
        coupled_bcs = coupled_model.boundary_conditions
        for coupled_bc in coupled_bcs
            if is_coupled_to_current(name, coupled_bc) == true
                coupled_bc.is_dirichlet = !bc.is_dirichlet
            end
        end
    end
end

function is_coupled_to_current(_::String, _::RegularBoundaryCondition)
    return false
end

function is_coupled_to_current(name::String, coupled_bc::ContactSchwarzBoundaryCondition)
    return name == coupled_bc.coupled_subsim.name
end

function is_coupled_to_current(name::String, coupled_bc::CouplingSchwarzBoundaryCondition)
    return name == coupled_bc.coupled_subsim.name
end
