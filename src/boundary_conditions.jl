# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.


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

function SolidMechanicsRobinBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    robin_parameter = bc_params["robin parameter"]
    tol = 1.0e-16
    if (abs(robin_parameter) < tol) 
        norma_abort(
                "The robin parameter is close to zero.  Robin BC is equivalent " *
                "to Neumann BC. Use Neumann BC in input file.")
    end 
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    rhs_num = eval(Meta.parse(expression))

    # Compile them into functions
    traction_fun = eval(build_function(rhs_num, [t, x, y, z]; expression=Val(false)))
 
    #We want to set traction + robin_parameter * disp = traction_fun
    return SolidMechanicsRobinBoundaryCondition(
        side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, 
        traction_fun, robin_parameter
    )
end


function SolidMechanicsOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    use_weak::Bool,
)
    mesh = get_fom_model(subsim).mesh
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_mesh = get_fom_model(coupled_subsim).mesh
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
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    return SolidMechanicsOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_nodes_indices,
        interpolation_function_values,
        coupled_block_name,
        tol,
        dirichlet_projector,
        use_weak,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    mesh::ExodusDatabase,
    side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    impedance::Float64,
    robin_parameter::Float64,
    impedance_scale::Vector{Float64},
)
    # Pointwise interpolation infrastructure (same as regular overlap)
    coupled_mesh = get_fom_model(coupled_subsim).mesh
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
                subsim.name, point[1], point[2], point[3], coupled_subsim.name,
            )
        end
        N = interpolate(element_type, ξ)[1]
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    # Surface projector infrastructure (for weak force application)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    square_projector = Matrix{Float64}(undef, 0, 0)
    return SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        coupled_nodes_indices,
        interpolation_function_values,
        local_from_global_map,
        global_from_local_map,
        square_projector,
        impedance,
        robin_parameter,
        impedance_scale,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsContactSchwarzBoundaryCondition(
    subsim::SingleDomainSimulation,
    coupled_subsim::SingleDomainSimulation,
    input_mesh::ExodusDatabase,
    bc_params::Parameters,
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
    return SolidMechanicsContactSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        is_dirichlet,
        swap_bcs,
        rotation_matrix,
        active_contact,
        friction_type,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsRobinSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation, 
    robin_parameter::Float64,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsRobinSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        square_projector,
        robin_parameter,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsImpedanceSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    impedance::Float64,
    robin_parameter::Float64,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsImpedanceSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        square_projector,
        impedance,
        robin_parameter,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsNonOverlapSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    subsim::Simulation,
    coupled_subsim::Simulation,
    is_dirichlet::Bool,
    swap_bcs::Bool,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
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
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        square_projector,
        is_dirichlet,
        swap_bcs,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
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
    if bc_type == "Schwarz overlap"
        use_weak = get(bc_params, "weak", false)
        SolidMechanicsOverlapSchwarzBoundaryCondition(
            coupled_block_name, tol, side_set_name, side_set_id, side_set_node_indices,
            num_nodes_sides, coupled_subsim, subsim, use_weak
        )
    elseif bc_type == "Schwarz DN nonoverlap"
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
            subsim,
            coupled_subsim,
            is_dirichlet,
            swap_bcs,
        )
    elseif bc_type == "Schwarz RR nonoverlap"
        robin_parameter = Float64(bc_params["robin parameter"])
        SolidMechanicsRobinSchwarzBoundaryCondition(
            input_mesh,
            side_set_name,
            coupled_side_set_name,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            coupled_subsim,
            subsim,
            robin_parameter,
        )
    elseif bc_type == "Schwarz impedance nonoverlap" || bc_type == "Schwarz impedance overlap"
        # Impedance Z = √(ρ(λ + 2μ)), computed from material properties
        model_params = subsim.params["model"]
        mat_params = model_params["material"]
        mat_blocks = mat_params["blocks"]
        mat_name = first(values(mat_blocks))
        mat_props = mat_params[mat_name]
        E = Float64(mat_props["elastic modulus"])
        ν = Float64(mat_props["Poisson's ratio"])
        ρ = Float64(mat_props["density"])
        λ_lame = E * ν / ((1 + ν) * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        impedance = sqrt(ρ * (λ_lame + 2μ))
        robin_parameter = Float64(get(bc_params, "robin parameter", 0.0))
        raw_scale = get(bc_params, "impedance scale", 1.0)
        if raw_scale isa AbstractVector
            impedance_scale = Float64.(raw_scale)
        else
            impedance_scale = [Float64(raw_scale)]
        end
        if bc_type == "Schwarz impedance nonoverlap"
            SolidMechanicsImpedanceSchwarzBoundaryCondition(
                input_mesh,
                side_set_name,
                coupled_side_set_name,
                side_set_id,
                side_set_node_indices,
                num_nodes_sides,
                coupled_subsim,
                subsim,
                impedance,
                robin_parameter,
            )
        else
            coupled_block_name = bc_params["source block"]
            tol = Float64(get(bc_params, "search tolerance", 1.0e-06))
            SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
                coupled_block_name,
                tol,
                input_mesh,
                side_set_name,
                side_set_id,
                side_set_node_indices,
                num_nodes_sides,
                coupled_subsim,
                subsim,
                impedance,
                robin_parameter,
                impedance_scale,
            )
        end
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
        model.displacement[bc.offset, node_index] = disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsNeumannRobinBoundaryCondition)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc.traction_fun, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index in side_nodes
            dof_index = 3 * (node_index - 1) + bc.offset
            model.boundary_force[dof_index] += nodal_force_component[side_node_index]
            side_node_index += 1
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

function apply_robin_bcs_internal_force!(model::SolidMechanics)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinBoundaryCondition || continue
        ss_node_index = 1
        for side in bc.num_nodes_per_side
            side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
            side_coordinates = model.reference[:, side_nodes]
            K_side = bc.robin_parameter * get_side_set_nodal_stiffness(side_coordinates, model.time)
            ss_node_index += side
            for (i, node_i) in enumerate(side_nodes)
                dof_i = 3 * (node_i - 1) + bc.offset
                for (j, node_j) in enumerate(side_nodes)
                    u_j = model.displacement[bc.offset, node_j]
                    model.internal_force[dof_i] += K_side[i, j] * u_j
                end
            end
        end
    end
end

function build_robin_stiffness(model::SolidMechanics)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_robin = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinBoundaryCondition || continue
        ss_node_index = 1
        for side in bc.num_nodes_per_side
            side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
            side_coordinates = model.reference[:, side_nodes]
            K_side = bc.robin_parameter * get_side_set_nodal_stiffness(side_coordinates, model.time)
            ss_node_index += side
            for (i, node_i) in enumerate(side_nodes)
                dof_i = 3 * (node_i - 1) + bc.offset
                for (j, node_j) in enumerate(side_nodes)
                    dof_j = 3 * (node_j - 1) + bc.offset
                    K_robin[dof_i, dof_j] += K_side[i, j]
                end
            end
        end
    end
    return K_robin
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


