# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

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

