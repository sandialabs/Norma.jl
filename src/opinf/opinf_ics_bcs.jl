# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

        
function SolidMechanicsOpInfDirichletBC(input_mesh::ExodusDatabase, bc_params::Dict{String,Any})
    fom_bc = SolidMechanicsDirichletBoundaryCondition(input_mesh,bc_params)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))
    
    opinf_model_directory = bc_params["model-directory"]
    py""" 
    import torch
    def get_model(model_file):
      import os
      assert os.path.isfile(model_file) ,  print(model_file + " cannot be found" )
      return torch.load(model_file,weights_only=False)
    """
    ensemble_size = bc_params["ensemble-size"]
    
    if offset == 1
        offset_name = "x"
    end
    if offset == 2
        offset_name = "y"
    end
    if offset == 3
        offset_name = "z"
    end


    model = [] 
    for i in 1:ensemble_size
      tmp =  py"get_model"(opinf_model_directory * "/BC-" * node_set_name * "-" * offset_name * "-" * string(i-1) * ".pt")
      push!(model,tmp)
    end     
            
    basis_file = bc_params["model-directory"] * "/nn-opinf-basis-" * node_set_name * "-" * offset_name * ".npz"
    basis = NPZ.npzread(basis_file)
    basis = basis["basis"]
            
            
    SolidMechanicsOpInfDirichletBC(
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
    
function SolidMechanicsOpInfOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    side_set_name::String,
    side_set_node_indices::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    variational::Bool,
    bc_params::Parameters,
)
    fom_bc = SolidMechanicsOverlapSchwarzBoundaryCondition(coupled_block_name,tol,side_set_name,side_set_node_indices,coupled_subsim,subsim,variational)
    opinf_model_directory = bc_params["model-directory"]
    py""" 
    import torch
    def get_model(model_file):
      return torch.load(model_file,weights_only=False)
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

    return SolidMechanicsOpInfOverlapSchwarzBoundaryCondition(
        fom_bc.name,
        fom_bc.side_set_node_indices,
        fom_bc.coupled_nodes_indices,
        fom_bc.interpolation_function_values,
        coupled_subsim,
        subsim,
        variational,
        fom_bc,
        model,
        basis
    )
end

function SMOpInfCouplingSchwarzBC(
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
    SolidMechanicsOpInfOverlapSchwarzBoundaryCondition(
        coupled_block_name, tol, side_set_name, side_set_node_indices, coupled_subsim, subsim, variational, bc_params
        )
end

function apply_bc(model::NeuralNetworkOpInfRom, bc::SolidMechanicsOpInfDirichletBC)
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
    model_inputs = py"setup_inputs"(reduced_bc_vector)
    ensemble_size = size(bc.nn_model)[1]
    for i in 1:ensemble_size
      reduced_forcing = bc.nn_model[i].forward(model_inputs)
      reduced_forcing = reduced_forcing.detach().numpy()[1,:]
      model.reduced_boundary_forcing[:] += reduced_forcing
    end 
    model.reduced_boundary_forcing[:] = model.reduced_boundary_forcing[:] ./ ensemble_size
end     

   
function apply_bc_detail(model::NeuralNetworkOpInfRom, bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
    if (typeof(bc.coupled_subsim.model) == SolidMechanics)
        ## Apply BC to the FOM vector
        apply_bc_detail(model.fom_model, bc.fom_bc)

        # populate our own BC vector
        unique_node_indices = unique(bc.fom_bc.side_set_node_indices)
        bc_vector = zeros(3, length(unique_node_indices))
        for i in 1:length(unique_node_indices)
            node_index = unique_node_indices[i]
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

        reduced_bc_vector = zeros(size(bc.basis)[3])
        for i in 1:3
            reduced_bc_vector[:] += bc.basis[i,:,:]' * bc_vector[i,:]
        end
        model_inputs = py"setup_inputs"(reduced_bc_vector)
        ensemble_size = size(bc.nn_model)[1]
        local_reduced_forcing = zeros(size(model.reduced_boundary_forcing)[1])
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


function apply_ics(params::Parameters, model::RomModel, integrator::TimeIntegrator, solver::Solver)
    ## Need to create a fake time integrator and solver for the FOM IC routine
    apply_ics(params, model.fom_model, integrator.fom_integrator, solver.fom_solver)

    if haskey(params, "initial conditions") == false
        return nothing
    end
    n_var, n_node, n_mode = model.basis.size
    n_var_fom, n_node_fom = size(model.fom_model.current)

    # Make sure basis is the right size
    if n_var != n_var_fom || n_node != n_node_fom
        norma_abort("Basis is wrong size")
    end

    # project onto basis
    for k in 1:n_mode
        model.reduced_state[k] = 0.0
        model.reduced_velocity[k] = 0.0
        for j in 1:n_node
            for n in 1:n_var
                model.reduced_state[k] +=
                    model.basis[n, j, k] * (model.fom_model.current[n, j] - model.fom_model.reference[n, j])
                model.reduced_velocity[k] += model.basis[n, j, k] * (model.fom_model.velocity[n, j])
            end
        end
    end
end

function apply_bcs(model::RomModel)
    model.reduced_boundary_forcing[:] .= 0.0
    for boundary_condition in model.boundary_conditions
        apply_bc(model, boundary_condition)
    end
end

function get_internal_force(model::RomModel)
    return model.fom_model.internal_force
end

function set_internal_force!(model::RomModel, force)
    return model.fom_model.internal_force = force
end

function apply_bc(model::OpInfModel, bc::SolidMechanicsDirichletBoundaryCondition)
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
    op_name = "B_" * bc.name * "-" * offset_name
    bc_operator = model.opinf_rom[op_name]
    # SM Dirichlet BC are only defined on a single x,y,z
    return model.reduced_boundary_forcing[:] += bc_operator[1, :, :] * bc_vector
end

function apply_bc_detail(model::OpInfModel, bc::SolidMechanicsCouplingSchwarzBoundaryCondition)
    if bc.coupled_subsim.model isa SolidMechanics || bc.coupled_subsim.model isa OpInfModel
        ## Apply BC to the FOM vector
        apply_bc_detail(model.fom_model, bc)

        unique_node_indices = unique(bc.side_set_node_indices)

        # populate our own BC vector
        bc_vector = zeros(3, length(unique_node_indices))
        for i in eachindex(unique_node_indices)
            node_index = unique_node_indices[i]
            bc_vector[:, i] = model.fom_model.current[:, node_index] - model.fom_model.reference[:, node_index]
        end
        op_name = "B_" * bc.name
        bc_operator = model.opinf_rom[op_name]
        for i in 1:3
            model.reduced_boundary_forcing[:] += bc_operator[i, :, :] * bc_vector[i, :]
        end
    else
        throw("ROM-ROM coupling not supported yet")
    end
end

function find_point_in_mesh(point::Vector{Float64}, model::RomModel, block_id::Int, tol::Float64)
    node_indices, ξ, found = find_point_in_mesh(point, model.fom_model, block_id, tol)
    return node_indices, ξ, found
end
