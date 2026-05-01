# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function write_stop_csv(sim::SingleDomainSimulation, model::RomModel)
    stop = sim.controller.stop
    index_string = "-" * string(stop; pad=4)
    prefix = sim.name * "-"
    reduced_states_filename = prefix * "reduced_states" * index_string * ".csv"
    writedlm(reduced_states_filename, model.reduced_state)
    write_stop_csv(sim, model.fom_model)
    return nothing
end

function write_stop_exodus(sim::SingleDomainSimulation, model::RomModel)
    integrator = sim.integrator
    solver = sim.solver
    #Re-construct full state
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration

    for i in 1:size(model.fom_model.displacement)[2]
        x_dof_index = 3 * (i - 1) + 1
        y_dof_index = 3 * (i - 1) + 2
        z_dof_index = 3 * (i - 1) + 3
        if model.fom_model.free_dofs[x_dof_index]
            model.fom_model.displacement[1, i] = model.basis[1, i, :]'displacement
            model.fom_model.velocity[1, i] = model.basis[1, i, :]'velocity
            model.fom_model.acceleration[1, i] = model.basis[1, i, :]'acceleration
        end

        if model.fom_model.free_dofs[y_dof_index]
            model.fom_model.displacement[2, i] = model.basis[2, i, :]'displacement
            model.fom_model.velocity[2, i] = model.basis[2, i, :]'velocity
            model.fom_model.acceleration[2, i] = model.basis[2, i, :]'acceleration
        end

        if model.fom_model.free_dofs[z_dof_index]
            model.fom_model.displacement[3, i] = model.basis[3, i, :]'displacement
            model.fom_model.velocity[3, i] = model.basis[3, i, :]'velocity
            model.fom_model.acceleration[3, i] = model.basis[3, i, :]'acceleration
        end
    end
    evaluate(model.fom_model, integrator.fom_integrator, solver.fom_solver)
    write_stop_exodus(sim, model.fom_model)
    return nothing
end

function write_sideset_stop_csv(sim::SingleDomainSimulation, model::RomModel)
    integrator = sim.integrator
    solver = sim.solver

    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration

    has_nonoverlap = any(
        bc -> bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition,
        model.boundary_conditions,
    )

    # Re-construct full state only when we need internal force for nonoverlap output.
    if has_nonoverlap
        for i in 1:size(model.fom_model.current)[2]
            x_dof_index = 3 * (i - 1) + 1
            y_dof_index = 3 * (i - 1) + 2
            z_dof_index = 3 * (i - 1) + 3
            if model.fom_model.free_dofs[x_dof_index]
                model.fom_model.current[1, i] = model.basis[1, i, :]'displacement + model.fom_model.reference[1, i]
                model.fom_model.velocity[1, i] = model.basis[1, i, :]'velocity
                model.fom_model.acceleration[1, i] = model.basis[1, i, :]'acceleration
            end

            if model.fom_model.free_dofs[y_dof_index]
                model.fom_model.current[2, i] = model.basis[2, i, :]'displacement + model.fom_model.reference[2, i]
                model.fom_model.velocity[2, i] = model.basis[2, i, :]'velocity
                model.fom_model.acceleration[2, i] = model.basis[2, i, :]'acceleration
            end

            if model.fom_model.free_dofs[z_dof_index]
                model.fom_model.current[3, i] = model.basis[3, i, :]'displacement + model.fom_model.reference[3, i]
                model.fom_model.velocity[3, i] = model.basis[3, i, :]'velocity
                model.fom_model.acceleration[3, i] = model.basis[3, i, :]'acceleration
            end
        end

        # Ensure internal force is consistent for force sideset output
        evaluate(model.fom_model, integrator.fom_integrator, solver.fom_solver)
    else
        # Only reconstruct nodes needed for sideset output
        nodes = Set{Int64}()
        for bc in model.boundary_conditions
            if bc isa SolidMechanicsDirichletBoundaryCondition
                for n in bc.node_set_node_indices
                    push!(nodes, n)
                end
            elseif bc isa SolidMechanicsOverlapSchwarzBoundaryCondition ||
                bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                for n in bc.side_set_node_indices
                    push!(nodes, n)
                end
            end
        end

        for i in nodes
            x_dof_index = 3 * (i - 1) + 1
            y_dof_index = 3 * (i - 1) + 2
            z_dof_index = 3 * (i - 1) + 3
            if model.fom_model.free_dofs[x_dof_index]
                model.fom_model.current[1, i] = model.basis[1, i, :]'displacement + model.fom_model.reference[1, i]
                model.fom_model.velocity[1, i] = model.basis[1, i, :]'velocity
                model.fom_model.acceleration[1, i] = model.basis[1, i, :]'acceleration
            end

            if model.fom_model.free_dofs[y_dof_index]
                model.fom_model.current[2, i] = model.basis[2, i, :]'displacement + model.fom_model.reference[2, i]
                model.fom_model.velocity[2, i] = model.basis[2, i, :]'velocity
                model.fom_model.acceleration[2, i] = model.basis[2, i, :]'acceleration
            end

            if model.fom_model.free_dofs[z_dof_index]
                model.fom_model.current[3, i] = model.basis[3, i, :]'displacement + model.fom_model.reference[3, i]
                model.fom_model.velocity[3, i] = model.basis[3, i, :]'velocity
                model.fom_model.acceleration[3, i] = model.basis[3, i, :]'acceleration
            end
        end
    end

    stop = sim.controller.stop
    index_string = "-" * string(stop; pad=4)
    prefix = sim.name * "-"
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsDirichletBoundaryCondition
            node_set_name = bc.name
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
            curr_filename = prefix * node_set_name * "-" * offset_name * "-curr" * index_string * ".csv"
            disp_filename = prefix * node_set_name * "-" * offset_name * "-disp" * index_string * ".csv"
            velo_filename = prefix * node_set_name * "-" * offset_name * "-velo" * index_string * ".csv"
            acce_filename = prefix * node_set_name * "-" * offset_name * "-acce" * index_string * ".csv"
            writedlm(curr_filename, model.fom_model.current[bc.offset, bc.node_set_node_indices])
            writedlm(velo_filename, model.fom_model.velocity[bc.offset, bc.node_set_node_indices])
            writedlm(acce_filename, model.fom_model.acceleration[bc.offset, bc.node_set_node_indices])
            writedlm(
                disp_filename,
                model.fom_model.current[bc.offset, bc.node_set_node_indices] -
                model.fom_model.reference[bc.offset, bc.node_set_node_indices],
            )
        elseif bc isa SolidMechanicsOverlapSchwarzBoundaryCondition ||
            bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
            side_set_name = bc.name
            curr_filename = prefix * side_set_name * "-curr" * index_string * ".csv"
            disp_filename = prefix * side_set_name * "-disp" * index_string * ".csv"
            velo_filename = prefix * side_set_name * "-velo" * index_string * ".csv"
            acce_filename = prefix * side_set_name * "-acce" * index_string * ".csv"
            unique_indices = unique(bc.side_set_node_indices)
            writedlm_nodal_array(curr_filename, model.fom_model.current[:, unique_indices])
            writedlm_nodal_array(velo_filename, model.fom_model.velocity[:, unique_indices])
            writedlm_nodal_array(acce_filename, model.fom_model.acceleration[:, unique_indices])
            writedlm_nodal_array(
                disp_filename,
                model.fom_model.current[:, unique_indices] - model.fom_model.reference[:, unique_indices],
            )

            if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                force_filename = prefix * side_set_name * "-force" * index_string * ".csv"
                # See ics_bcs.get_dst_force() for this
                # Will get projected with Neumann projector onto destination simulation boundary
                force_global = model.fom_model.internal_force
                force = extract_local_vector(bc, force_global, 3)
                num_nodes = length(bc.global_from_local_map)
                force_out = reshape(force, (3, num_nodes))
                writedlm_nodal_array(force_filename, force_out)
            end
        end
    end
    return nothing
end
