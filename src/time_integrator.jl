# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using DelimitedFiles
using Format

function adaptive_stepping_parameters(integrator_params::Parameters)
    has_minimum = haskey(integrator_params, "minimum time step")
    has_decrease = haskey(integrator_params, "decrease factor")
    has_maximum = haskey(integrator_params, "maximum time step")
    has_increase = haskey(integrator_params, "increase factor")
    has_any = has_minimum || has_decrease || has_maximum || has_increase
    has_all = has_minimum && has_decrease && has_maximum && has_increase
    if has_any == true && has_all == false
        error(
            "Adaptive time stepping requires 4 parameters: minimum and maximum time steps and decrease and increase factors",
        )
    elseif has_any == true && has_all == true
        minimum_time_step = integrator_params["minimum time step"]
        decrease_factor = integrator_params["decrease factor"]
        maximum_time_step = integrator_params["maximum time step"]
        increase_factor = integrator_params["increase factor"]
    else
        minimum_time_step = maximum_time_step = integrator_params["time step"]
        decrease_factor = increase_factor = 1.0
    end
    if minimum_time_step > maximum_time_step
        error("Minimum time step ", minimum_time_step, " must be greater or equal than maximum ", maximum_time_step)
    end
    if decrease_factor > 1.0
        error("Decrease factor ", decrease_factor, " must be less or equal to one.")
    end
    if increase_factor < 1.0
        error("Increase factor ", increase_factor, " must be greater or equal to one.")
    end
    return minimum_time_step, decrease_factor, maximum_time_step, increase_factor
end

function QuasiStatic(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    runtime_step = 0.0
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = initial_time
    stop = 0
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    stored_energy = 0.0
    initial_equilibrium = false
    if haskey(integrator_params, "initial equilibrium") == true
        initial_equilibrium = integrator_params["initial equilibrium"]
    end
    return QuasiStatic(
        initial_time,
        final_time,
        time_step,
        runtime_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        prev_time,
        time,
        stop,
        displacement,
        velocity,
        acceleration,
        stored_energy,
        initial_equilibrium,
    )
end

function Newmark(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    runtime_step = 0.0
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = initial_time
    stop = 0
    β = integrator_params["β"]
    γ = integrator_params["γ"]
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    return Newmark(
        initial_time,
        final_time,
        time_step,
        runtime_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        prev_time,
        time,
        stop,
        β,
        γ,
        displacement,
        velocity,
        acceleration,
        disp_pre,
        velo_pre,
        stored_energy,
        kinetic_energy,
    )
end

function CentralDifference(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = 0.0
    runtime_step = 0.0
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    stable_time_step = 0.0
    user_time_step = integrator_params["time step"]
    time = prev_time = initial_time
    stop = 0
    CFL = integrator_params["CFL"]
    γ = integrator_params["γ"]
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    return CentralDifference(
        initial_time,
        final_time,
        time_step,
        runtime_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        user_time_step,
        stable_time_step,
        prev_time,
        time,
        stop,
        CFL,
        γ,
        displacement,
        velocity,
        acceleration,
        stored_energy,
        kinetic_energy,
    )
end

function create_time_integrator(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    integrator_name = integrator_params["type"]
    if integrator_name == "quasi static"
        return QuasiStatic(params, model)
    elseif integrator_name == "Newmark"
        return Newmark(params, model)
    elseif integrator_name == "central difference"
        return CentralDifference(params, model)
    else
        error("Unknown type of time integrator : ", integrator_name)
    end
end

function is_static(integrator::TimeIntegrator)
    return integrator isa QuasiStatic
end

function is_dynamic(integrator::TimeIntegrator)
    return is_static(integrator) == false
end

function initialize(integrator::Newmark, solver::HessianMinimizer, model::RomModel)
    # Compute initial accelerations
    evaluate(model.fom_model, integrator, solver)
    free = model.fom_model.free_dofs
    internal_force = model.fom_model.internal_force
    external_force = model.fom_model.body_force + model.fom_model.boundary_force
    if model.fom_model.inclined_support == true
        external_force = model.fom_model.global_transform * external_force
        internal_force = model.fom_model.global_transform * internal_force
    end
    inertial_force = external_force - internal_force

    # project onto basis
    n_var, n_node, n_mode = size(model.basis)
    num_nodes = n_node
    #_, num_nodes = size(model.fom_model.reference)
    acceleration = zeros(3 * num_nodes)

    acceleration[free] = model.fom_model.mass[free, free] \ inertial_force[free]

    integrator.displacement[:] = model.reduced_state[:]
    integrator.velocity[:] = model.reduced_velocity[:]
    solver.solution[:] = model.reduced_state[:]

    for k in 1:n_mode
        integrator.acceleration[k] = 0.0
        for j in 1:n_node
            for n in 1:n_var
                integrator.acceleration[k] += model.basis[n, j, k] * acceleration[3 * (j - 1) + n]
            end
        end
    end
    return nothing
end

function predict(integrator::Newmark, solver::Solver, model::RomModel)
    dt = integrator.time_step
    beta = integrator.β
    gamma = integrator.γ
    integrator.disp_pre[:] =
        integrator.displacement[:] =
            integrator.displacement[:] +
            dt * integrator.velocity +
            1.0 / 2.0 * dt * dt * (1.0 - 2.0 * beta) * integrator.acceleration
    integrator.velo_pre[:] = integrator.velocity[:] += dt * (1.0 - gamma) * integrator.acceleration
    solver.solution[:] = integrator.displacement[:]
    model.reduced_state[:] = integrator.displacement[:]
    return nothing
end

function correct(integrator::Newmark, solver::Solver, model::RomModel)
    dt = integrator.time_step
    beta = integrator.β
    gamma = integrator.γ
    integrator.displacement[:] = solver.solution[:]
    integrator.acceleration[:] = (solver.solution - integrator.disp_pre) / (beta * dt * dt)
    integrator.velocity[:] = integrator.velo_pre + dt * gamma * integrator.acceleration[:]
    model.reduced_state[:] = solver.solution[:]
    return nothing
end

function initialize(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    if integrator.initial_equilibrium == true
        norma_log(0, :equilibrium, "Establishing Initial Equilibrium...")
        solve(integrator, solver, model)
        if model.failed == true
            error("Finite element model failed to establish initial equlibrium")
        end
    end
    return nothing
end

function predict(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    return nothing
end

function correct(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    copy_solution_source_targets(solver, model, integrator)
    return nothing
end

function initialize(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    evaluate(model, integrator, solver)
    if model.failed == true
        error("Finite element model failed to initialize")
    end
    internal_force = model.internal_force
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    inertial_force = external_force - internal_force
    kinetic_energy = 0.5 * dot(integrator.velocity, model.mass, integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = model.strain_energy
    integrator.acceleration[free] = solve_linear(model.mass[free, free], inertial_force[free])
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function initialize(integrator::Newmark, solver::MatrixFree, model::SolidMechanics)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    evaluate(model, integrator, solver)
    if model.failed == true
        error("Finite element model failed to initialize")
    end
    internal_force = model.internal_force
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    inertial_force = external_force - internal_force
    kinetic_energy = 0.5 * model.lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = model.strain_energy
    integrator.acceleration[free] = inertial_force[free] ./ model.lumped_mass[free]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function predict(integrator::Newmark, solver::Solver, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    fixed = .!free
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    uᵖʳᵉ[free] = u[free] += Δt * v[free] + (0.5 - β) * Δt * Δt * a[free]
    vᵖʳᵉ[free] = v[free] += (1.0 - γ) * Δt * a[free]
    uᵖʳᵉ[fixed] = u[fixed]
    vᵖʳᵉ[fixed] = v[fixed]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function correct(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement = solver.solution
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    integrator.acceleration[free] = (u[free] - uᵖʳᵉ[free]) / β / Δt / Δt
    integrator.velocity[free] = vᵖʳᵉ[free] + γ * Δt * integrator.acceleration[free]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function initialize(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    set_time_step(integrator, model)
    evaluate(model, integrator, solver)
    if model.failed == true
        error("The finite element model has failed to initialize")
    end
    internal_force = model.internal_force
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    kinetic_energy = 0.5 * model.lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = model.strain_energy
    inertial_force = external_force - internal_force
    integrator.acceleration[free] = inertial_force[free] ./ model.lumped_mass[free]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function predict(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    set_time_step(integrator, model)
    Δt = integrator.time_step
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    u[free] += Δt * v[free] + 0.5 * Δt * Δt * a[free]
    v[free] += (1.0 - γ) * Δt * a[free]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function correct(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    Δt = integrator.time_step
    γ = integrator.γ
    a = integrator.acceleration = solver.solution
    free = model.free_dofs
    integrator.velocity[free] += γ * Δt * a[free]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function initialize_writing(params::Parameters, integrator::DynamicTimeIntegrator, model::RomModel)
    initialize_writing(params, integrator, model.fom_model)
    return nothing
end

function initialize_writing(params::Parameters, integrator::TimeIntegrator, _::SolidMechanics)
    output_mesh = params["output_mesh"]

    # global variables
    num_global_vars = Exodus.read_number_of_variables(output_mesh, GlobalVariable)
    num_global_vars += 1  # adjust as necessary
    Exodus.write_number_of_variables(output_mesh, GlobalVariable, num_global_vars)

    runtime_step_index = 1
    Exodus.write_name(output_mesh, GlobalVariable, Int32(runtime_step_index), "runtime_step")

    # setup nodal variables
    num_node_vars = 6
    node_var_names = ["refe_x", "refe_y", "refe_z", "disp_x", "disp_y", "disp_z"]
    if is_dynamic(integrator) == true
        num_node_vars += 6
        append!(node_var_names, ["velo_x", "velo_y", "velo_z", "acce_x", "acce_y", "acce_z"])
    end
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_names(output_mesh, NodalVariable, node_var_names)

    # get maximum number of quadrature points
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for block in blocks
        blk_id = block.id
        element_type_string = Exodus.read_block_parameters(output_mesh, blk_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end

    num_element_vars = 6 * max_num_int_points + 1
    Exodus.write_number_of_variables(output_mesh, ElementVariable, num_element_vars)

    el_var_names = String[]
    for point in 1:max_num_int_points
        ip_str = "_" * string(point)
        push!(el_var_names, "stress_xx" * ip_str)
        push!(el_var_names, "stress_yy" * ip_str)
        push!(el_var_names, "stress_zz" * ip_str)
        push!(el_var_names, "stress_yz" * ip_str)
        push!(el_var_names, "stress_xz" * ip_str)
        push!(el_var_names, "stress_xy" * ip_str)
    end
    push!(el_var_names, "stored_energy")
    Exodus.write_names(output_mesh, ElementVariable, el_var_names)
    return nothing
end

function finalize_writing(params::Parameters)
    input_mesh = params["input_mesh"]
    Exodus.close(input_mesh)
    output_mesh = params["output_mesh"]
    Exodus.close(output_mesh)
    return nothing
end

function writedlm_nodal_array(filename::String, nodal_array::Matrix{Float64})
    open(filename, "w") do io
        for col in 1:size(nodal_array, 2)
            # Write each column as a comma-separated line
            println(io, join(nodal_array[:, col], ","))
        end
    end
    return nothing
end

function write_step(params::Parameters, integrator::TimeIntegrator, model::Model)
    stop = integrator.stop
    exodus_interval = get(params, "Exodus output interval", 1)
    if exodus_interval > 0 && stop % exodus_interval == 0
        write_step_exodus(params, integrator, model)
    end
    csv_interval = get(params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        sim_id = 1
        if haskey(params, "global_simulation") == true
            sim_id = params["global_simulation"].subsim_name_index_map[params["name"]]
        end
        write_step_csv(integrator, model, sim_id)
        if haskey(params, "CSV write sidesets") == true
            write_sideset_step_csv(integrator, model, sim_id)
        end
    end
    return nothing
end

function write_step_csv(integrator::TimeIntegrator, model::SolidMechanics, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop; pad=4)
    sim_id_string = string(sim_id; pad=2) * "-"
    if stop == 0
        refe_filename = sim_id_string * "refe" * ".csv"
        writedlm_nodal_array(refe_filename, model.reference)
    end
    free_dofs_filename = sim_id_string * "free_dofs" * index_string * ".csv"
    writedlm(free_dofs_filename, model.free_dofs)
    curr_filename = sim_id_string * "curr" * index_string * ".csv"
    writedlm_nodal_array(curr_filename, model.current)
    disp_filename = sim_id_string * "disp" * index_string * ".csv"
    writedlm_nodal_array(disp_filename, model.current - model.reference)
    time_filename = sim_id_string * "time" * index_string * ".csv"
    writedlm(time_filename, integrator.time, '\n')
    potential_filename = sim_id_string * "potential" * index_string * ".csv"
    writedlm(potential_filename, integrator.stored_energy, '\n')
    if is_dynamic(integrator) == true
        velo_filename = sim_id_string * "velo" * index_string * ".csv"
        writedlm_nodal_array(velo_filename, model.velocity)
        acce_filename = sim_id_string * "acce" * index_string * ".csv"
        writedlm_nodal_array(acce_filename, model.acceleration)
        kinetic_filename = sim_id_string * "kinetic" * index_string * ".csv"
        writedlm(kinetic_filename, integrator.kinetic_energy, '\n')
    end
    return nothing
end

function write_sideset_step_csv(integrator::DynamicTimeIntegrator, model::SolidMechanics, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop; pad=4)
    sim_id_string = string(sim_id; pad=2) * "-"
    for bc in model.boundary_conditions
        if bc isa SMDirichletBC
            node_set_name = bc.node_set_name
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
            curr_filename = sim_id_string * node_set_name * "-" * offset_name * "-curr" * index_string * ".csv"
            disp_filename = sim_id_string * node_set_name * "-" * offset_name * "-disp" * index_string * ".csv"
            velo_filename = sim_id_string * node_set_name * "-" * offset_name * "-velo" * index_string * ".csv"
            acce_filename = sim_id_string * node_set_name * "-" * offset_name * "-acce" * index_string * ".csv"
            writedlm(curr_filename, model.current[bc.offset, bc.node_set_node_indices])
            writedlm(velo_filename, model.velocity[bc.offset, bc.node_set_node_indices])
            writedlm(acce_filename, model.acceleration[bc.offset, bc.node_set_node_indices])
            writedlm(
                disp_filename,
                model.current[bc.offset, bc.node_set_node_indices] -
                model.reference[bc.offset, bc.node_set_node_indices],
            )
        elseif bc isa SMOverlapSchwarzBC
            side_set_name = bc.side_set_name
            curr_filename = sim_id_string * side_set_name * "-curr" * index_string * ".csv"
            disp_filename = sim_id_string * side_set_name * "-disp" * index_string * ".csv"
            velo_filename = sim_id_string * side_set_name * "-velo" * index_string * ".csv"
            acce_filename = sim_id_string * side_set_name * "-acce" * index_string * ".csv"
            writedlm_nodal_array(curr_filename, model.current[:, bc.side_set_node_indices])
            writedlm_nodal_array(velo_filename, model.velocity[:, bc.side_set_node_indices])
            writedlm_nodal_array(acce_filename, model.acceleration[:, bc.side_set_node_indices])
            writedlm_nodal_array(
                disp_filename, model.current[:, bc.side_set_node_indices] - model.reference[:, bc.side_set_node_indices]
            )
        end
    end
    return nothing
end

function write_step_csv(integrator::DynamicTimeIntegrator, model::RomModel, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop; pad=4)
    sim_id_string = string(sim_id; pad=2) * "-"
    reduced_states_filename = sim_id_string * "reduced_states" * index_string * ".csv"
    writedlm(reduced_states_filename, model.reduced_state)
    write_step_csv(integrator, model.fom_model, sim_id)
    return nothing
end

function write_step_exodus(params::Parameters, integrator::DynamicTimeIntegrator, model::RomModel)
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
    write_step_exodus(params, integrator, model.fom_model)
    return nothing
end

function write_step_exodus(params::Parameters, integrator::TimeIntegrator, model::SolidMechanics)
    time = integrator.time
    stop = integrator.stop
    runtime_step = integrator.runtime_step
    time_index = stop + 1
    output_mesh = params["output_mesh"]
    Exodus.write_time(output_mesh, time_index, time)

    # global variables
    Exodus.write_values(output_mesh, GlobalVariable, time_index, [runtime_step])

    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_x", refe_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_y", refe_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_z", refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_x", disp_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_y", disp_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_z", disp_z)
    if is_dynamic(integrator) == true
        velocity = model.velocity
        velo_x = velocity[1, :]
        velo_y = velocity[2, :]
        velo_z = velocity[3, :]
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_x", velo_x)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_y", velo_y)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_z", velo_z)
        acceleration = model.acceleration
        acce_x = acceleration[1, :]
        acce_y = acceleration[2, :]
        acce_z = acceleration[3, :]
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_x", acce_x)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_y", acce_y)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_z", acce_z)
    end
    stress = model.stress
    stored_energy = model.stored_energy
    blocks = Exodus.read_sets(output_mesh, Block)
    for (block, block_stress, block_stored_energy) in zip(blocks, stress, stored_energy)
        blk_id = block.id
        element_type_string, num_blk_elems, _, _, _, _ = Exodus.read_block_parameters(output_mesh, blk_id)
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        stress_xx = zeros(num_blk_elems, num_points)
        stress_yy = zeros(num_blk_elems, num_points)
        stress_zz = zeros(num_blk_elems, num_points)
        stress_yz = zeros(num_blk_elems, num_points)
        stress_xz = zeros(num_blk_elems, num_points)
        stress_xy = zeros(num_blk_elems, num_points)
        for blk_elem_index in 1:num_blk_elems
            element_stress = block_stress[blk_elem_index]
            for point in 1:num_points
                point_stress = element_stress[point]
                stress_xx[blk_elem_index, point] = point_stress[1]
                stress_yy[blk_elem_index, point] = point_stress[2]
                stress_zz[blk_elem_index, point] = point_stress[3]
                stress_yz[blk_elem_index, point] = point_stress[4]
                stress_xz[blk_elem_index, point] = point_stress[5]
                stress_xy[blk_elem_index, point] = point_stress[6]
            end
        end
        for point in 1:num_points
            ip_str = "_" * string(point)
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_xx" * ip_str, stress_xx[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_yy" * ip_str, stress_yy[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_zz" * ip_str, stress_zz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_yz" * ip_str, stress_yz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_xz" * ip_str, stress_xz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(blk_id), "stress_xy" * ip_str, stress_xy[:, point]
            )
        end
        Exodus.write_values(
            output_mesh, ElementVariable, time_index, Int64(blk_id), "stored_energy", block_stored_energy
        )
    end
    return nothing
end
