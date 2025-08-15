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
    #Re-construct full state
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration

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
    write_stop_exodus(sim, model.fom_model)
    return nothing
end

