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
  model = load_torch_models([
    opinf_model_directory * "/BC-" * node_set_name * "-" * offset_name * "-" * string(i - 1) * ".pt" for i in 1:ensemble_size
  ])
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
  side_set_id::Int64,
  side_set_node_indices::Vector{Int64},
  num_nodes_sides::Vector{Int64},
  coupled_subsim::Simulation,
  subsim::Simulation,
  bc_params::Parameters,
)
  # FIX 1: default must be "" (String), not false (Bool), to match the
  # SolidMechanicsOverlapSchwarzBoundaryCondition constructor signature
  # which changed compute_overlap_l2_error from Bool to String on main.
  compute_overlap_l2_error = get(bc_params, "compute overlap L2 relative error", "")
  fom_bc = SolidMechanicsOverlapSchwarzBoundaryCondition(
    coupled_block_name, tol, side_set_name, side_set_id, side_set_node_indices,
    num_nodes_sides, coupled_subsim, subsim, false, compute_overlap_l2_error
  )
  opinf_model_directory = bc_params["model-directory"]
  ensemble_size = bc_params["ensemble-size"]
  model = load_torch_models([
    opinf_model_directory * "/BC-" * side_set_name * "-" * string(i - 1) * ".pt" for i in 1:ensemble_size
  ])
  basis_file = bc_params["model-directory"] * "/nn-opinf-basis-" * side_set_name * ".npz"
  basis = NPZ.npzread(basis_file)
  basis = basis["basis"]
  return SolidMechanicsOpInfOverlapSchwarzBoundaryCondition(
    fom_bc.name,
    fom_bc.side_set_node_indices,
    fom_bc.coupled_nodes_indices,
    fom_bc.interpolation_function_values,
    fom_bc,
    model,
    basis,
    subsim.parent,
    subsim.handle,
    coupled_subsim.handle,
  )
end

function compute_overlap_l2_error!(bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
  return compute_overlap_l2_error!(bc.fom_bc)
end

function stores_overlap_l2_error(bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
  # FIX 2: fom_bc.compute_overlap_l2_error is now a String, not a Bool;
  # use !isempty() to return a Bool as callers expect.
  return !isempty(bc.fom_bc.compute_overlap_l2_error)
end

function overlap_l2_error_field(bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
  return bc.fom_bc.compute_overlap_l2_error
end

function get_overlap_l2_error(bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
  return bc.fom_bc.overlap_l2_error
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
  # Overlap coupling uses a source block, not a source side set; the latter is
  # not needed here (see SMCouplingSchwarzBC).
  side_set_id = side_set_id_from_name(side_set_name, input_mesh)
  num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
  num_nodes_sides = Int64.(num_nodes_sides)
  side_set_node_indices = Int64.(side_set_node_indices)
  SolidMechanicsOpInfOverlapSchwarzBoundaryCondition(
    coupled_block_name, tol, side_set_name, side_set_id, side_set_node_indices,
    num_nodes_sides, coupled_subsim, subsim, bc_params
  )
end

## apply_bc(::NeuralNetworkOpInfRom, ::SolidMechanicsOpInfDirichletBC) and
## apply_bc_detail(::NeuralNetworkOpInfRom, ::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
## run PyTorch forward passes and are defined in ext/NormaPyTorchExt.jl. They are
## only reachable once a NeuralNetworkOpInfRom has been constructed, which already
## requires the PyCall extension.

function apply_ics(params::Parameters, model::RomModel, integrator::TimeIntegrator, solver::Solver)
  ## Need to create a fake time integrator and solver for the FOM IC routine
  apply_ics(params, model.fom_model, integrator.fom_integrator, solver.fom_solver)
  # The FOM displacement/velocity arrays are nonzero (and must be projected
  # onto the reduced basis below) in two mutually-exclusive cases: an explicit
  # `initial conditions:` block (just applied above), or a `restart:` snapshot
  # (applied directly to model.fom_model.displacement/velocity inside
  # SolidMechanics() construction, before this function ever runs — see
  # model.jl and process_restart!() in simulation.jl). Without this check,
  # reduced_state/reduced_velocity would silently stay at their
  # zero-initialized default on restart, discarding the restored FOM state.
  has_ics = haskey(params, "initial conditions")
  is_restart = get(params, "restart_info", nothing) !== nothing
  if !has_ics && !is_restart
    return nothing
  end
  n_var, n_node, n_mode = size(model.basis)
  n_var_fom, n_node_fom = size(model.fom_model.displacement)
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
          model.basis[n, j, k] * model.fom_model.displacement[n, j]
        model.reduced_velocity[k] += model.basis[n, j, k] * (model.fom_model.velocity[n, j])
      end
    end
  end
  if is_restart
    norma_log(0, :restart, "Projected restarted FOM displacement/velocity onto the reduced basis.")
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
    disp_val = model.fom_model.displacement[bc.offset, node_index]
    push!(bc_vector, disp_val)
  end
  # SM Dirichlet BC are only defined on a single x,y,z
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
  # Compute contribution from BC
  op_name = "B_" * bc.name * "-" * offset_name
  bc_operator = model.opinf_rom[op_name]
  scale = get(model.opinf_rom, "input-scale-" * op_name, 1.0)
  bc_vector .= scale .* bc_vector
  bc_contrib = bc_operator[1, :, :] * bc_vector
  model.reduced_boundary_forcing[:] += bc_contrib
end

function apply_bc_detail(model::OpInfModel, bc::SolidMechanicsCouplingSchwarzBoundaryCondition)
  coupled_model = coupled_subsim_of(bc).model
  if coupled_model isa SolidMechanics || coupled_model isa OpInfModel
    ## Apply BC to the FOM vector
    apply_bc_detail(model.fom_model, bc)
    if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition && !bc.is_dirichlet
      nodal_force = get_dst_force(bc)
      num_nodes = length(bc.global_from_local_map)
      force = reshape(nodal_force, (3, num_nodes))
      # apply scaling
      scale = get(model.opinf_rom, "input-scale-B_N_" * bc.name, 1.0)
      force .= scale .* force
      # apply operator to BC internal force vector
      op_name = "B_N_" * bc.name
      bc_operator = model.opinf_rom[op_name]
      for i in 1:3
        bc_contrib = bc_operator[i, :, :] * force[i, :]
        model.reduced_boundary_forcing[:] += bc_contrib
      end
    else
      unique_node_indices = unique(bc.side_set_node_indices)
      # populate our own BC vector
      bc_vector = zeros(3, length(unique_node_indices))
      for i in eachindex(unique_node_indices)
        node_index = unique_node_indices[i]
        bc_vector[:, i] = model.fom_model.displacement[:, node_index]
      end
      # apply pre-scaling to input
      scale = get(model.opinf_rom, "input-scale-B_" * bc.name, 1.0)
      bc_vector .= scale .* bc_vector
      op_name = "B_" * bc.name
      bc_operator = model.opinf_rom[op_name]
      for i in 1:3
        model.reduced_boundary_forcing[:] += bc_operator[i, :, :] * bc_vector[i, :]
      end
    end
  else
    throw("ROM-ROM coupling not supported yet")
  end
end

function find_point_in_mesh(point::Vector{Float64}, model::RomModel, block_id::Int, tol::Float64)
  node_indices, ξ, found = find_point_in_mesh(point, model.fom_model, block_id, tol)
  return node_indices, ξ, found
end
