# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Printf
using YAML

include("interpolation_types.jl")
include("simulation_types.jl")
include("minitensor.jl")
include("model.jl")
include("time_integrator.jl")
include("solver.jl")
include("io.jl")

function create_simulation(input_file::String)
    norma_log(0, :setup, "Reading from " * input_file)
    params = YAML.load_file(input_file; dicttype=Parameters)
    basename = stripped_name(input_file)
    params["name"] = basename
    return create_simulation(params)
end

function create_simulation(params::Parameters)
    sim_type = params["type"]
    if sim_type == "single"
        sim = SingleDomainSimulation(params)
        create_bcs(sim)
        return sim
    elseif sim_type == "multi"
        sim = MultiDomainSimulation(params)
        create_bcs(sim)
        return sim
    else
        norma_abort("Unknown type of simulation: $sim_type")
    end
end

function SingleDomainSimulation(params::Parameters)
    basename = params["name"]
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    rm(output_mesh_file; force=true)
    input_mesh = Exodus.ExodusDatabase(input_mesh_file, "r")
    Exodus.copy(input_mesh, output_mesh_file)
    output_mesh = Exodus.ExodusDatabase(output_mesh_file, "rw")
    params["output_mesh"] = output_mesh
    params["input_mesh"] = input_mesh
    controller = create_controller(params)
    model = create_model(params)
    integrator = create_time_integrator(params, model)
    solver = create_solver(params, model)
    failed = false
    return SingleDomainSimulation(basename, params, controller, integrator, solver, model, failed)
end

function MultiDomainSimulation(params::Parameters)
    basename = params["name"]
    domain_paths = params["domains"]
    subsims = Vector{SingleDomainSimulation}()
    controller = create_controller(params)
    initial_time = controller.initial_time
    final_time = controller.final_time
    time_step = controller.time_step
    exodus_interval = get(params, "Exodus output interval", 1)
    csv_interval = get(params, "CSV output interval", 0)
    subsim_name_index_map = Dict{String,Int64}()
    subsim_index = 1
    for domain_path in domain_paths
        domain_name = stripped_name(domain_path)
        norma_log(4, :domain, domain_name)
        subparams = YAML.load_file(domain_path; dicttype=Parameters)
        subparams["name"] = domain_name
        integrator_params = subparams["time integrator"]
        subsim_time_step = get(integrator_params, "time step", time_step)
        subsim_time_step = min(subsim_time_step, time_step)
        integrator_params["initial time"] = initial_time
        integrator_params["final time"] = final_time
        integrator_params["time step"] = subsim_time_step
        subparams["Exodus output interval"] = exodus_interval
        subparams["CSV output interval"] = csv_interval
        subsim = SingleDomainSimulation(subparams)
        params[domain_name] = subsim.params
        push!(subsims, subsim)
        subsim_name_index_map[domain_name] = subsim_index
        subsim_index += 1
    end
    num_domains = length(subsims)
    failed = false
    sim = MultiDomainSimulation(basename, params, controller, num_domains, subsims, subsim_name_index_map, failed)
    for subsim in sim.subsims
        subsim.params["parent_simulation"] = sim
    end
    return sim
end

function SolidMultiDomainTimeController(params::Parameters)
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["absolute tolerance"]
    relative_tolerance = params["relative tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
    absolute_error = relative_error = 0.0
    time = prev_time = initial_time
    stop = 0
    converged = false
    iteration_number = 0
    stop_disp = [Vector{Float64}[] for _ in 1:num_domains]
    stop_velo = [Vector{Float64}[] for _ in 1:num_domains]
    stop_acce = [Vector{Float64}[] for _ in 1:num_domains]
    stop_∂Ω_f = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_disp = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_velo = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_acce = [Vector{Float64}[] for _ in 1:num_domains]
    time_hist = [Float64[] for _ in 1:num_domains]
    disp_hist = [Vector{Float64}[] for _ in 1:num_domains]
    velo_hist = [Vector{Float64}[] for _ in 1:num_domains]
    acce_hist = [Vector{Float64}[] for _ in 1:num_domains]
    ∂Ω_f_hist = [Vector{Float64}[] for _ in 1:num_domains]
    relaxation_parameter = get(params, "relaxation parameter", 1.0)
    naive_stabilized = get(params, "naive stabilized", false)
    lambda_disp = [Vector{Float64}[] for _ in 1:num_domains]
    lambda_velo = [Vector{Float64}[] for _ in 1:num_domains]
    lambda_acce = [Vector{Float64}[] for _ in 1:num_domains]
    is_schwarz = true
    schwarz_contact = false
    active_contact = false
    contact_hist = Vector{Bool}[]
    schwarz_iters = zeros(Int64, num_stops - 1)

    return SolidMultiDomainTimeController(
        minimum_iterations,
        maximum_iterations,
        absolute_tolerance,
        relative_tolerance,
        absolute_error,
        relative_error,
        initial_time,
        final_time,
        time_step,
        time,
        prev_time,
        num_stops,
        stop,
        converged,
        iteration_number,
        stop_disp,
        stop_velo,
        stop_acce,
        stop_∂Ω_f,
        schwarz_disp,
        schwarz_velo,
        schwarz_acce,
        time_hist,
        disp_hist,
        velo_hist,
        acce_hist,
        ∂Ω_f_hist,
        relaxation_parameter,
        naive_stabilized,
        lambda_disp,
        lambda_velo,
        lambda_acce,
        is_schwarz,
        schwarz_contact,
        active_contact,
        contact_hist,
        schwarz_iters,
    )
end

function SolidSingleDomainTimeController(params::Parameters)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
    time = prev_time = initial_time
    stop = 0
    return SolidSingleDomainTimeController(initial_time, final_time, time_step, time, prev_time, num_stops, stop)
end

function create_controller(params::Parameters)
    sim_type = params["type"]
    if sim_type == "single"
        return SolidSingleDomainTimeController(params)
    elseif sim_type == "multi"
        return SolidMultiDomainTimeController(params)
    else
        norma_abort("Unknown type of simulation: $sim_type")
    end
end

function create_bcs(sim::SingleDomainSimulation)
    boundary_conditions = create_bcs(sim.params)
    for bc in boundary_conditions
        if bc isa SolidMechanicsInclinedDirichletBoundaryCondition ||
            bc isa SolidMechanicsContactSchwarzBoundaryCondition
            sim.model.inclined_support = true
            break
        end
    end
    return sim.model.boundary_conditions = boundary_conditions
end

function create_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        create_bcs(subsim)
    end
    pair_schwarz_bcs(sim)
    return nothing
end

function evolve(sim::Simulation)
    sync_control_time(sim)
    initialize(sim)
    initialize_writing(sim)
    write_stop(sim)
    while true
        advance_control_time(sim)
        sync_control_time(sim)
        advance_control(sim)
        write_stop(sim)
        if stop_evolve(sim) == true
            break
        end
    end
    finalize_writing(sim)
    return nothing
end

function stop_evolve(sim::Simulation)
    time = sim.controller.time
    final_time = sim.controller.final_time
    return isapprox(time, final_time; rtol=1e-9, atol=1e-12) || time > final_time
end

function solve_contact(sim::MultiDomainSimulation)
    if sim.controller.active_contact == true
        schwarz(sim)
    else
        advance_independent(sim)
    end
end

function decrease_time_step(sim::SingleDomainSimulation)
    norma_log(0, :recover, "Attempting to recover by decreasing time step...")
    time_step = sim.integrator.time_step
    decrease_factor = sim.integrator.decrease_factor
    if decrease_factor == 1.0
        norma_abortf(
            "Cannot adapt time step %.4e because decrease factor is %.4e. Enable adaptive time stepping.",
            time_step,
            decrease_factor,
        )
    end
    new_time_step = decrease_factor * time_step
    minimum_time_step = sim.integrator.minimum_time_step
    if new_time_step < minimum_time_step
        norma_abortf("Cannot adapt time step to %.4e because minimum is %.4e.", new_time_step, minimum_time_step)
    end
    sim.integrator.time_step = new_time_step
    return nothing
end

function increase_time_step(sim::SingleDomainSimulation)
    increase_factor = sim.integrator.increase_factor
    if increase_factor > 1.0
        time_step = sim.integrator.time_step
        maximum_time_step = sim.integrator.maximum_time_step
        new_time_step = min(increase_factor * time_step, maximum_time_step)
        if new_time_step > time_step
            sim.integrator.time_step = new_time_step
        end
    end
    return nothing
end

function advance_one_step(sim::SingleDomainSimulation)
    while true
        prev_time = sim.integrator.prev_time
        time = sim.integrator.time
        time_step = sim.integrator.time_step
        norma_logf(4, :advance, "Time = [%.4e, %.4e] : Δt = %.4e", prev_time, time, time_step)
        apply_bcs(sim)
        solve(sim)
        if sim.failed == false
            increase_time_step(sim)
            save_curr_state(sim)
            break
        end
        restore_prev_state(sim)
        decrease_time_step(sim)
        advance_time(sim)
        sim.failed = sim.model.failed = sim.solver.failed = false
    end
    return nothing
end

function advance_control(sim::SingleDomainSimulation)
    subcycle(sim)
    return nothing
end

function advance_control(sim::MultiDomainSimulation)
    if sim.controller.schwarz_contact == false
        schwarz(sim)
        return nothing
    end
    save_stop_state(sim)
    solve_contact(sim)
    was_in_contact = sim.controller.active_contact
    detect_contact(sim)
    if sim.controller.active_contact ≠ was_in_contact
        if was_in_contact == true
            norma_log(0, :contact, "Released — reattempting control step.")
        else
            norma_log(0, :contact, "Initiated — reattempting control step.")
        end
        restore_stop_state(sim)
        solve_contact(sim)
    end
    return nothing
end

function apply_ics(sim::SingleDomainSimulation)
    apply_ics(sim.params, sim.model, sim.integrator, sim.solver)
    return nothing
end

function apply_ics(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_ics(subsim)
    end
end

function apply_bcs(sim::SingleDomainSimulation)
    apply_bcs(sim.model)
    return nothing
end

function apply_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_bcs(subsim)
    end
end

function initialize(sim::SingleDomainSimulation)
    apply_ics(sim)
    apply_bcs(sim)
    initialize(sim.integrator, sim.solver, sim.model)
    save_curr_state(sim)
    return nothing
end

function initialize(sim::MultiDomainSimulation)
    initialize_bc_projectors(sim)
    apply_ics(sim)
    apply_bcs(sim)
    for subsim in sim.subsims
        initialize(subsim.integrator, subsim.solver, subsim.model)
        save_curr_state(subsim)
    end
    detect_contact(sim)
    return nothing
end

function solve(sim::SingleDomainSimulation)
    solve(sim.integrator, sim.solver, sim.model)
    return sim.failed = sim.model.failed
end

function initialize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        initialize_writing(subsim)
    end
end

function finalize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        finalize_writing(subsim)
    end
end

function sync_control_time(sim::SingleDomainSimulation)
    set_initial_subcycle_time(sim)
    return nothing
end

function sync_control_time(sim::MultiDomainSimulation)
    set_initial_subcycle_time(sim)
    controller = sim.controller
    initial_time = controller.initial_time
    final_time = controller.final_time
    time_step = controller.time_step
    time = controller.time
    prev_time = controller.prev_time
    num_stops = controller.num_stops
    stop = controller.stop
    for subsim in sim.subsims
        subsim.controller.initial_time = initial_time
        subsim.controller.final_time = final_time
        subsim.controller.time_step = time_step
        subsim.controller.time = time
        subsim.controller.prev_time = prev_time
        subsim.controller.num_stops = num_stops
        subsim.controller.stop = stop
    end
    return nothing
end

function set_initial_subcycle_time(sim::SingleDomainSimulation)
    time = sim.controller.prev_time
    sim.model.time = sim.integrator.time = time
    return nothing
end

function set_initial_subcycle_time(sim::MultiDomainSimulation)
    time = sim.controller.prev_time
    for subsim in sim.subsims
        subsim.model.time = subsim.integrator.time = time
    end
    return nothing
end

function get_adjusted_timestep(t::Float64, dt::Float64, t_stop::Float64, eps::Float64=0.01)::Float64
    t_next = t + dt
    gap = t_stop - t
    tol = eps * dt
    return if isapprox(t_next, t_stop; rtol=1e-6, atol=1e-12) || t_next > t_stop
        gap
    elseif t_next ≥ t_stop - tol
        gap / 2
    else
        dt
    end
end

function advance_time(sim::SingleDomainSimulation)
    controller = sim.controller
    integrator = sim.integrator
    time = integrator.time
    time_stop = controller.time
    time_step = get_adjusted_timestep(time, integrator.time_step, time_stop)
    next_time = time + time_step
    integrator.prev_time = time
    integrator.time = sim.model.time = next_time
    integrator.time_step = time_step
    return nothing
end

function advance_control_time(sim::Simulation)
    controller = sim.controller
    controller.prev_time = controller.time
    stop = controller.stop + 1
    initial_time = controller.initial_time
    time_step = controller.time_step
    next_time = stop * time_step + initial_time
    controller.time = next_time
    controller.stop = stop
    return nothing
end

function advance_independent(sim::MultiDomainSimulation)
    sim.controller.iteration_number = 0
    sim.controller.is_schwarz = false
    save_stop_state(sim)
    set_initial_subcycle_time(sim)
    subcycle(sim)
    return nothing
end

function schwarz(sim::MultiDomainSimulation)
    iteration_number = 1
    sim.controller.is_schwarz = true
    save_stop_state(sim)
    save_schwarz_state(sim)
    reset_histories(sim)
    swap_swappable_bcs(sim)

    while true
        norma_log(0, :schwarz, "Iteration [$iteration_number]")
        sim.controller.iteration_number = iteration_number
        set_initial_subcycle_time(sim)
        subcycle(sim)
        ΔU, Δu = update_schwarz_convergence_criterion(sim)
        raw_status = sim.controller.converged ? "[CONVERGED]" : "[CONVERGING]"
        status = colored_status(raw_status)
        norma_logf(
            0,
            :schwarz,
            "Criterion [%d] %s = %.3e : %s = %.3e : %s",
            iteration_number,
            "|ΔU|",
            ΔU,
            "|ΔU|/|U|",
            Δu,
            status,
        )
        if stop_schwarz(sim, iteration_number + 1) == true
            plural = iteration_number == 1 ? "" : "s"
            norma_log(0, :schwarz, "Performed $iteration_number Schwarz Iteration" * plural)
            sim.controller.schwarz_iters[sim.controller.stop] = iteration_number
            break
        end
        iteration_number += 1
        save_schwarz_state(sim)
        restore_stop_state(sim)
    end
end

function save_curr_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    if sim.model.inclined_support == true
        global_transform_T = sim.model.global_transform'
        integrator.prev_disp = global_transform_T * integrator.displacement
        integrator.prev_velo = global_transform_T * integrator.velocity
        integrator.prev_acce = global_transform_T * integrator.acceleration
    else
        integrator.prev_disp = copy(integrator.displacement)
        integrator.prev_velo = copy(integrator.velocity)
        integrator.prev_acce = copy(integrator.acceleration)
    end
    return integrator.prev_∂Ω_f = copy(sim.model.internal_force)
end

function restore_prev_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    integrator.time = integrator.prev_time
    if sim.model.inclined_support == true
        global_transform = sim.model.global_transform
        integrator.displacement = global_transform * integrator.prev_disp
        integrator.velocity = global_transform * integrator.prev_velo
        integrator.acceleration = global_transform * integrator.prev_acce
    else
        integrator.displacement = copy(integrator.prev_disp)
        integrator.velocity = copy(integrator.prev_velo)
        integrator.acceleration = copy(integrator.prev_acce)
    end
    sim.model.internal_force = copy(integrator.prev_∂Ω_f)
    copy_solution_source_targets(sim.integrator, sim.solver, sim.model)
    return nothing
end

function save_stop_state(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    subsims = sim.subsims
    for i in 1:num_domains
        subsim = subsims[i]
        # If this model has inclined support on, we need to rotate the integrator values
        if subsim.model.inclined_support == true
            global_transform_T = subsim.model.global_transform'
            controller.stop_disp[i] = global_transform_T * subsim.integrator.displacement
            controller.stop_velo[i] = global_transform_T * subsim.integrator.velocity
            controller.stop_acce[i] = global_transform_T * subsim.integrator.acceleration
        else
            controller.stop_disp[i] = copy(subsim.integrator.displacement)
            controller.stop_velo[i] = copy(subsim.integrator.velocity)
            controller.stop_acce[i] = copy(subsim.integrator.acceleration)
        end
        controller.stop_∂Ω_f[i] = copy(subsim.model.internal_force)
    end
end

function restore_stop_state(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    subsims = sim.subsims
    for i in 1:num_domains
        subsim = subsims[i]
        subsim.integrator.time = controller.prev_time
        # If this model has inclined support on, we need to rotate the integrator values
        if subsim.model.inclined_support == true
            global_transform = subsim.model.global_transform
            subsim.integrator.displacement = global_transform * controller.stop_disp[i]
            subsim.integrator.velocity = global_transform * controller.stop_velo[i]
            subsim.integrator.acceleration = global_transform * controller.stop_acce[i]
        else
            subsim.integrator.displacement = copy(controller.stop_disp[i])
            subsim.integrator.velocity = copy(controller.stop_velo[i])
            subsim.integrator.acceleration = copy(controller.stop_acce[i])
        end
        subsim.model.internal_force = copy(controller.stop_∂Ω_f[i])
        copy_solution_source_targets(subsim.integrator, subsim.solver, subsim.model)
    end
end

function save_schwarz_state(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = sim.num_domains
    for i in 1:num_domains
        subsim = subsims[i]
        # Note: Integrator values are not rotated due to inclined support as the
        # schwarz variables are only used for Schwarz convergence which are compared
        # to simulation integrator values.
        controller.schwarz_disp[i] = copy(subsim.integrator.displacement)
        controller.schwarz_velo[i] = copy(subsim.integrator.velocity)
        controller.schwarz_acce[i] = copy(subsim.integrator.acceleration)
    end
end

function swap_swappable_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        swap_swappable_bcs(subsim)
    end
end

function swap_swappable_bcs(sim::SingleDomainSimulation)
    for bc in sim.model.boundary_conditions
        if bc isa SolidMechanicsContactSchwarzBoundaryCondition ||
            bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
            if (bc.swap_bcs == true)
                bc.is_dirichlet = !bc.is_dirichlet
            end
        end
    end
end

function subcycle(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    for subsim_index in 1:num_domains
        subsim = sim.subsims[subsim_index]
        norma_log(4, :domain, subsim.name)
        reset_history(controller, subsim_index)
        while true
            advance_time(subsim)
            advance_one_step(subsim)
            if controller.active_contact == true && controller.naive_stabilized == true
                apply_naive_stabilized_bcs(subsim)
            end
            save_history_snapshot(controller, subsim, subsim_index)
            if stop_subcyle(subsim) == true
                break
            end
        end
    end
    return nothing
end

function subcycle(sim::SingleDomainSimulation)
    while true
        advance_time(sim)
        advance_one_step(sim)
        if stop_subcyle(sim) == true
            break
        end
    end
    return nothing
end

function stop_subcyle(sim::SingleDomainSimulation)
    return isapprox(sim.integrator.time, sim.controller.time; rtol=1.0e-06, atol=1.0e-12)
end

function reset_histories(sim::MultiDomainSimulation)
    controller = sim.controller
    if controller.is_schwarz == false
        return nothing
    end
    num_domains = sim.num_domains
    for subsim_index in 1:num_domains
        reset_history(controller, subsim_index)
    end
    return nothing
end

function reset_history(controller::MultiDomainTimeController, subsim_index::Int64)
    if controller.is_schwarz == false
        return nothing
    end
    empty!(controller.time_hist[subsim_index])
    empty!(controller.disp_hist[subsim_index])
    empty!(controller.velo_hist[subsim_index])
    empty!(controller.acce_hist[subsim_index])
    empty!(controller.∂Ω_f_hist[subsim_index])
    return nothing
end

function save_history_snapshot(controller::MultiDomainTimeController, sim::SingleDomainSimulation, subsim_index::Int64)
    if controller.is_schwarz == false
        return nothing
    end
    push!(controller.time_hist[subsim_index], sim.integrator.time)
    push!(controller.disp_hist[subsim_index], copy(sim.integrator.displacement))
    push!(controller.velo_hist[subsim_index], copy(sim.integrator.velocity))
    push!(controller.acce_hist[subsim_index], copy(sim.integrator.acceleration))
    push!(controller.∂Ω_f_hist[subsim_index], copy(sim.model.internal_force))
    return nothing
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = sim.num_domains
    norms_disp = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i in 1:num_domains
        Δt = controller.time_step
        x_prev = controller.schwarz_disp[i] + Δt * controller.schwarz_velo[i]
        x_curr = subsims[i].integrator.displacement + Δt * subsims[i].integrator.velocity
        norms_disp[i] = norm(x_curr)
        norms_diff[i] = norm(x_curr - x_prev)
    end
    norm_disp = norm(norms_disp)
    norm_diff = norm(norms_diff)
    controller.absolute_error = norm_diff
    controller.relative_error = norm_disp > 0.0 ? norm_diff / norm_disp : norm_diff
    conv_abs = controller.absolute_error ≤ controller.absolute_tolerance
    conv_rel = controller.relative_error ≤ controller.relative_tolerance
    controller.converged = conv_abs || conv_rel
    return controller.absolute_error, controller.relative_error
end

function stop_schwarz(sim::MultiDomainSimulation, iteration_number::Int64)
    if sim.controller.absolute_error == 0.0
        return true
    end
    exceeds_minimum_iterations = iteration_number > sim.controller.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations = iteration_number > sim.controller.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return sim.controller.converged
end

function check_overlap(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    distance_tol = 0.0
    parametric_tol = 1.0e-06
    overlap = false
    unique_node_indices = unique(bc.side_set_node_indices)
    coupled_model = bc.coupled_subsim.model
    coupled_bc = coupled_model.boundary_conditions[bc.coupled_bc_index]
    coupled_side_set_id = coupled_bc.side_set_id

    # Tolerance value to calculate closest point projection, to avoid projection failure
    # Set to 10 times the minimum characteristic edge length of the side set nodes
    # TODO (BRP): perhaps calculate this "characteristic size" a priori, and only recalculate for finite kinematics
    nodal_points = model.current[:, unique_node_indices]
    # Compute minimum distances between the nodes in the side set
    max_distance_tolerance =
        10 * minimum([
            minimum([norm(p1 - p2) for (j, p2) in enumerate(eachcol(nodal_points)) if i != j]) for
            (i, p1) in enumerate(eachcol(nodal_points))
        ])

    for node_index in unique_node_indices
        point = model.current[:, node_index]
        # Precompute the face node distance
        face_nodes, face_node_indices, min_distance = closest_face_to_point(point, coupled_model, coupled_side_set_id)
        if min_distance > max_distance_tolerance
            continue
        end
        _, ξ, _, coupled_face_node_indices, _, distance = project_point_to_side_set(
            point, face_nodes, face_node_indices
        )
        num_nodes_coupled_side = length(coupled_face_node_indices)
        parametric_dim = length(ξ)
        element_type = get_element_type(parametric_dim, num_nodes_coupled_side)
        overlap = distance ≤ distance_tol && is_inside_parametric(element_type, ξ, parametric_tol)
        if overlap == true
            break
        end
    end
    return overlap
end

function check_compression(
    mesh::ExodusDatabase, model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition
)
    compression_tol = 0.0
    compression = false
    nodal_forces = get_dst_force(bc)
    normals = compute_normal(mesh, bc.side_set_id, model)
    global_from_local_map = bc.global_from_local_map
    for local_node in eachindex(global_from_local_map)
        local_range = (3 * (local_node - 1) + 1):(3 * local_node)
        nodal_force = nodal_forces[local_range]
        normal = normals[:, local_node]
        normal_force = dot(nodal_force, normal)
        compressive_force = normal_force ≤ compression_tol
        if compressive_force == true
            compression = true
            break
        end
    end
    return compression
end

function initialize_bc_projectors(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if !(
                bc isa SolidMechanicsContactSchwarzBoundaryCondition ||
                bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
            )
                continue
            end
            compute_dirichlet_projector(subsim.model, bc)
            compute_neumann_projector(subsim.model, bc)
        end
    end
end

function detect_contact(sim::MultiDomainSimulation)
    if sim.controller.schwarz_contact == false
        return nothing
    end
    num_domains = sim.num_domains
    persistence = sim.controller.active_contact
    contact_domain = falses(num_domains)
    for domain in 1:num_domains
        subsim = sim.subsims[domain]
        mesh = subsim.params["input_mesh"]
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SolidMechanicsContactSchwarzBoundaryCondition
                if persistence == true
                    compression = check_compression(mesh, subsim.model, bc)
                    contact_domain[domain] = compression == true
                else
                    overlap = check_overlap(subsim.model, bc)
                    contact_domain[domain] = overlap == true
                end
            end
        end
    end

    sim.controller.active_contact = any(contact_domain)
    for domain in 1:num_domains
        subsim = sim.subsims[domain]
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SolidMechanicsContactSchwarzBoundaryCondition
                bc.active_contact = sim.controller.active_contact
            end
        end
    end
    if sim.controller.active_contact == true
        norma_log(0, :contact, "Detected")
    end
    resize!(sim.controller.contact_hist, sim.controller.stop + 1)
    sim.controller.contact_hist[sim.controller.stop + 1] = sim.controller.active_contact
    write_schwarz_params_csv(sim)
    return nothing
end

function write_schwarz_params_csv(sim::MultiDomainSimulation)
    stop = sim.controller.stop
    csv_interval = get(sim.params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        index_string = "-" * string(stop; pad=4)
        contact_filename = "contact" * index_string * ".csv"
        writedlm(contact_filename, sim.controller.active_contact, '\n')
        iters_filename = "iterations" * index_string * ".csv"
        writedlm(iters_filename, sim.controller.iteration_number, '\n')
    end
end
