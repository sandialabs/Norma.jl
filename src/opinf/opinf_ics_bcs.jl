# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using PyCall
import NPZ

@variables t, x, y, z
D = Differential(t)


function SMOpInfDirichletBC(input_mesh::ExodusDatabase, bc_params::Dict{String,Any})
    fom_bc = SMDirichletBC(input_mesh,bc_params)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))

    opinf_model_file = bc_params["model-file"]
    py""" 
    import torch
    def get_model(model_file):
      return torch.load(model_file)
    """
    model = py"get_model"(opinf_model_file)

    opinf_model_file = bc_params["model-file"]
    basis_file = bc_params["basis-file"]
    basis = NPZ.npzread(basis_file)
    basis = basis["basis"]

    SMOpInfDirichletBC(
        node_set_name,
        offset,
        node_set_id,
        node_set_node_indices,
        disp_num,
        velo_num,
        acce_num,
        fom_bc,
        model,
        basis,
    )
end


function SMOpInfCouplingSchwarzBC(
    subsim::SingleDomainSimulation,
    coupled_subsim::SingleDomainSimulation,
    input_mesh::ExodusDatabase,
    bc_type::String,
    bc_params::Parameters,
)
    fom_bc = SMCouplingSchwarzBC(subsim,coupled_subsim,input_mesh,"Schwarz overlap",bc_params) 
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

    opinf_model_directory = bc_params["model-directory"]
    py""" 
    import torch
    def get_model(model_file):
      return torch.load(model_file)
    """
    ensemble_size = bc_params["ensemble-size"]
    model = []
    for i in 1:ensemble_size
      tmp =  py"get_model"(opinf_model_directory * "/BC-" * side_set_name * "-" * string(i-1) * ".pt")
      push!(model,tmp)
    end

    basis_file = bc_params["model-directory"] * "/nn-opinf-basis-" * side_set_name * ".npz" 
    basis = NPZ.npzread(basis_file)
    basis = basis["basis"]

    if bc_type == "OpInf Schwarz overlap"
        SMOpInfOverlapSchwarzBC(
            side_set_name,
            side_set_node_indices,
            coupled_nodes_indices,
            interpolation_function_values,
            coupled_subsim,
            subsim,
            is_dirichlet,
            swap_bcs,
            fom_bc,
            model,
            basis
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

function apply_bc(model::NeuralNetworkOpInfRom, bc::SMOpInfDirichletBC)
    model.fom_model.time = model.time
    apply_bc(model.fom_model,bc.fom_bc)
    bc_vector = zeros(0)
    for node_index ∈ bc.fom_bc.node_set_node_indices
        dof_index = 3 * (node_index - 1) + bc.fom_bc.offset
        disp_val = model.fom_model.current[bc.fom_bc.offset,node_index] - model.fom_model.reference[bc.fom_bc.offset, node_index]
        push!(bc_vector,disp_val)
    end

    py"""
    import numpy as np
    def setup_inputs(x):
        xi = np.zeros((1,x.size))
        xi[0] = x
        inputs = torch.tensor(xi)
        return inputs
    """

    reduced_bc_vector = bc.basis[1,:,:]' * bc_vector 
    print(size(reduced_bc_vector),size(bc.basis),size(bc_vector))
    model_inputs = py"setup_inputs"(reduced_bc_vector)
    ensemble_size = size(bc.nn_model)[1]
    for i in 1:ensemble_size
      reduced_forcing = bc.nn_model[i].forward(model_inputs)
      reduced_forcing = reduced_forcing.detach().numpy()[1,:]
      model.reduced_boundary_forcing[:] += reduced_forcing 
    end
    model.reduced_boundary_forcing[:] = model.reduced_boundary_forcing[:] ./ ensemble_size
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


function apply_bc_detail(model::NeuralNetworkOpInfRom, bc::CouplingSchwarzBoundaryCondition)
    if (typeof(bc.coupled_subsim.model) == SolidMechanics)
        ## Apply BC to the FOM vector
        apply_bc_detail(model.fom_model, bc.fom_bc)

        # populate our own BC vector
        bc_vector = zeros(3, length(bc.fom_bc.side_set_node_indices))
        for i in 1:length(bc.fom_bc.side_set_node_indices)
            node_index = bc.fom_bc.side_set_node_indices[i]
            bc_vector[:, i] = model.fom_model.current[:, node_index] - model.fom_model.reference[:, node_index]
        end  

        py"""
        import numpy as np
        def setup_inputs(x):
            xi = np.zeros((1,x.size))
            xi[0] = x
            inputs = torch.tensor(xi)
            return inputs
        """

        reduced_bc_vector = zeros(bc.basis.size[3])
        for i in 1:3
            reduced_bc_vector[:] += bc.basis[i,:,:]' * bc_vector[i,:] 
        end  
        model_inputs = py"setup_inputs"(reduced_bc_vector)
        ensemble_size = size(bc.nn_model)[1]
        local_reduced_forcing = zeros(model.reduced_boundary_forcing.size[1])
        for i in 1:ensemble_size
          reduced_forcing = bc.nn_model[i].forward(model_inputs)
          reduced_forcing = reduced_forcing.detach().numpy()[1,:]
          local_reduced_forcing[:] += reduced_forcing 
        end
        local_reduced_forcing[:] = local_reduced_forcing[:] ./ ensemble_size
        model.reduced_boundary_forcing[:] += local_reduced_forcing

    else 
        throw("ROM-ROM coupling not supported yet")
    end  
end

function apply_bc(model::RomModel, bc::SchwarzBoundaryCondition)
    global_sim = bc.coupled_subsim.params["global_simulation"]
    controller = global_sim.controller
    if bc isa SMContactSchwarzBC && controller.active_contact == false
        return nothing
    end
    empty_history = length(controller.time_hist) == 0
    same_step = controller.same_step
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
    time_hist = controller.time_hist[coupled_index]
    disp_hist = controller.disp_hist[coupled_index]
    velo_hist = controller.velo_hist[coupled_index]
    acce_hist = controller.acce_hist[coupled_index]
    ∂Ω_f_hist = controller.∂Ω_f_hist[coupled_index]
    interp_disp = same_step == true ? disp_hist[end] : interpolate(time_hist, disp_hist, time)
    interp_velo = same_step == true ? velo_hist[end] : interpolate(time_hist, velo_hist, time)
    interp_acce = same_step == true ? acce_hist[end] : interpolate(time_hist, acce_hist, time)
    interp_∂Ω_f = same_step == true ? ∂Ω_f_hist[end] : interpolate(time_hist, ∂Ω_f_hist, time)
    if bc.coupled_subsim.model isa SolidMechanics
        bc.coupled_subsim.model.internal_force = interp_∂Ω_f
    elseif bc.coupled_subsim.model isa RomModel
        bc.coupled_subsim.model.fom_model.internal_force = interp_∂Ω_f
    end

    if bc isa SMContactSchwarzBC || bc isa SMNonOverlapSchwarzBC
        relaxation_parameter = controller.relaxation_parameter
        schwarz_iteration = controller.iteration_number
        if schwarz_iteration == 1
            lambda_dispᵖʳᵉᵛ = zeros(length(interp_disp))
            lambda_veloᵖʳᵉᵛ = zeros(length(interp_velo))
            lambda_acceᵖʳᵉᵛ = zeros(length(interp_acce))
        else
            lambda_dispᵖʳᵉᵛ = controller.lambda_disp[coupled_index]
            lambda_veloᵖʳᵉᵛ = controller.lambda_velo[coupled_index]
            lambda_acceᵖʳᵉᵛ = controller.lambda_acce[coupled_index]
        end
        bc.coupled_subsim.integrator.displacement =
            controller.lambda_disp[coupled_index] =
                relaxation_parameter * interp_disp + (1 - relaxation_parameter) * lambda_dispᵖʳᵉᵛ
        bc.coupled_subsim.integrator.velocity =
            controller.lambda_velo[coupled_index] =
                relaxation_parameter * interp_velo + (1 - relaxation_parameter) * lambda_veloᵖʳᵉᵛ
        bc.coupled_subsim.integrator.acceleration =
            controller.lambda_acce[coupled_index] =
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


function apply_bcs(model::RomModel)
    model.reduced_boundary_forcing[:] .= 0.0
    for boundary_condition in model.boundary_conditions
        apply_bc(model, boundary_condition)
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


function find_point_in_mesh(point::Vector{Float64}, model::RomModel, blk_id::Int, tol::Float64)
    node_indices, ξ, found = find_point_in_mesh(point, model.fom_model, blk_id, tol)
    return node_indices, ξ, found
end

