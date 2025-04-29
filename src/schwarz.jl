# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function SolidMultiDomainController(params::Parameters)
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    input_time_step = params["time step"]
    num_stops = max(round(Int64, (final_time - initial_time) / input_time_step) + 1, 2)
    time_step = (final_time - initial_time) / (num_stops - 1)
    norma_logf(0, :time, "Time Step = %.4e : Adjusted = %.4e : Total Stops = %d", input_time_step, time_step, num_stops)
    absolute_error = relative_error = 0.0
    time = prev_time = initial_time
    same_step = true
    stop = 0
    converged = false
    iteration_number = 0
    stop_disp = Vector{Vector{Float64}}(undef, num_domains)
    stop_velo = Vector{Vector{Float64}}(undef, num_domains)
    stop_acce = Vector{Vector{Float64}}(undef, num_domains)
    stop_∂Ω_f = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_disp = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_velo = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_acce = Vector{Vector{Float64}}(undef, num_domains)
    time_hist = Vector{Vector{Float64}}()
    disp_hist = Vector{Vector{Vector{Float64}}}()
    velo_hist = Vector{Vector{Vector{Float64}}}()
    acce_hist = Vector{Vector{Vector{Float64}}}()
    ∂Ω_f_hist = Vector{Vector{Vector{Float64}}}()
    if haskey(params, "relaxation parameter") == true
        relaxation_parameter = params["relaxation parameter"]
    else
        relaxation_parameter = 1.0
    end
    if haskey(params, "naive stabilized") == true
        naive_stabilized = params["naive stabilized"]
    else
        naive_stabilized = false
    end
    lambda_disp = Vector{Vector{Float64}}(undef, num_domains)
    lambda_velo = Vector{Vector{Float64}}(undef, num_domains)
    lambda_acce = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_contact = false
    active_contact = false
    contact_hist = Vector{Bool}()

    csv_interval = get(params, "CSV output interval", 0)
    if csv_interval > 0
        iterations = params["maximum iterations"]
        convergence_hist = zeros(Float64, iterations, 2)
    else
        convergence_hist = Array{Float64}(undef, 0, 0)
    end
    return SolidMultiDomainController(
        num_domains,
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
        same_step,
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
        schwarz_contact,
        active_contact,
        contact_hist,
        convergence_hist,
    )
end

function create_controller(params::Parameters)
    return SolidMultiDomainController(params)
end

function advance_independent(sim::MultiDomainSimulation)
    sim.controller.iteration_number = 0
    save_stop_solutions(sim)
    set_subcycle_times(sim)
    synchronize(sim)
    is_schwarz = false
    return subcycle(sim, is_schwarz)
end

function schwarz(sim::MultiDomainSimulation)
    iteration_number = 1
    save_stop_solutions(sim)
    save_schwarz_solutions(sim)
    resize_histories(sim)
    set_subcycle_times(sim)
    swap_swappable_bcs(sim)
    is_schwarz = true

    csv_interval = get(sim.params, "CSV output interval", 0)
    if csv_interval > 0
        sim.controller.convergence_hist .= 0.0
    end

    while true
        norma_log(0, :schwarz, "Iteration [$iteration_number]")
        sim.controller.iteration_number = iteration_number
        synchronize(sim)
        subcycle(sim, is_schwarz)
        ΔU, Δu = update_schwarz_convergence_criterion(sim)
        if csv_interval > 0
            sim.controller.convergence_hist[iteration_number, 1] = ΔU
            sim.controller.convergence_hist[iteration_number, 2] = Δu
        end
        raw_status = sim.controller.converged ? "[DONE]" : "[WAIT]"
        status = colored_status(raw_status)
        norma_logf(0, :schwarz, "Convergence [%d] %s = %.3e : %s = %.3e : %s", iteration_number, "|ΔU|", ΔU, "|ΔU|/|U|", Δu, status)
                if stop_schwarz(sim, iteration_number + 1) == true
            plural = iteration_number == 1 ? "" : "s"
            norma_log(0, :schwarz, "Performed $iteration_number Schwarz Iteration" * plural)
            break
        end
        iteration_number += 1
        save_schwarz_solutions(sim)
        restore_stop_solutions(sim)
    end
end

function save_stop_solutions(sim::SingleDomainSimulation)
    sim.controller.prev_time = sim.integrator.prev_time
    sim.controller.time = sim.integrator.time
    sim.controller.stop = sim.integrator.stop
    if sim.model.inclined_support == true
        global_transform_T = sim.model.global_transform'
        sim.controller.stop_disp = global_transform_T * sim.integrator.displacement
        sim.controller.stop_velo = global_transform_T * sim.integrator.velocity
        sim.controller.stop_acce = global_transform_T * sim.integrator.acceleration
    else
        sim.controller.stop_disp = copy(sim.integrator.displacement)
        sim.controller.stop_velo = copy(sim.integrator.velocity)
        sim.controller.stop_acce = copy(sim.integrator.acceleration)
    end
    return sim.controller.stop_∂Ω_f = copy(sim.model.internal_force)
end

function save_stop_solutions(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    for i in 1:(controller.num_domains)
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

function restore_stop_solutions(sim::SingleDomainSimulation)
    sim.integrator.prev_time = sim.controller.prev_time
    sim.integrator.time = sim.controller.time
    sim.integrator.stop = sim.controller.stop
    if sim.model.inclined_support == true
        global_transform = sim.model.global_transform
        sim.integrator.displacement = global_transform * sim.controller.stop_disp
        sim.integrator.velocity = global_transform * sim.controller.stop_velo
        sim.integrator.acceleration = global_transform * sim.controller.stop_acce
    else
        sim.integrator.displacement = copy(sim.controller.stop_disp)
        sim.integrator.velocity = copy(sim.controller.stop_velo)
        sim.integrator.acceleration = copy(sim.controller.stop_acce)
    end
    sim.model.internal_force = copy(sim.controller.stop_∂Ω_f)
    return copy_solution_source_targets(sim.integrator, sim.solver, sim.model)
end

function restore_stop_solutions(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    for i in 1:(controller.num_domains)
        subsim = subsims[i]
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

function save_schwarz_solutions(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    for i in 1:(controller.num_domains)
        subsim = subsims[i]
        # Note: Integrator values are not rotated due to inclined support as the 
        # schwarz variables are only used for Schwarz convergence which are compared
        # to simulation integrator values.
        controller.schwarz_disp[i] = copy(subsim.integrator.displacement)
        controller.schwarz_velo[i] = copy(subsim.integrator.velocity)
        controller.schwarz_acce[i] = copy(subsim.integrator.acceleration)
    end
end

function set_subcycle_times(sim::MultiDomainSimulation)
    initial_time = sim.controller.prev_time
    final_time = sim.controller.time
    for subsim in sim.subsims
        subsim.integrator.initial_time = initial_time
        subsim.integrator.time = initial_time
        subsim.integrator.final_time = final_time
    end
end

function swap_swappable_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        swap_swappable_bcs(subsim)
    end
end

function swap_swappable_bcs(sim::SingleDomainSimulation)
    for bc in sim.model.boundary_conditions
        if bc isa ContactSchwarzBoundaryCondition || bc isa CouplingSchwarzBoundaryCondition
            if (bc.swap_bcs == true)
                bc.is_dirichlet = !bc.is_dirichlet
            end
        end
    end
end

function subcycle(sim::MultiDomainSimulation, is_schwarz::Bool)
    subsim_index = 1
    for subsim in sim.subsims
        norma_log(4, :domain, subsim.name)
        norma_logf(4, :initial, "Time = %.4e", subsim.integrator.time)
        stop_index = 1
        while true
            advance_time(subsim)
            if stop_evolve(subsim) == true
                break
            end
            subsim.model.time = subsim.integrator.time
            advance(subsim)
            norma_logf(4, :advance, "Time = %.4e : Δt = %.4e", subsim.integrator.time, subsim.integrator.time_step)
            if sim.controller.active_contact == true && sim.controller.naive_stabilized == true
                apply_naive_stabilized_bcs(subsim)
            end
            stop_index += 1
            if is_schwarz == true
                save_history_snapshot(sim.controller, sim.subsims, subsim_index, stop_index)
            end
        end
        subsim_index += 1
    end
end

function resize_histories(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = controller.num_domains
    resize!(controller.time_hist, num_domains)
    resize!(controller.disp_hist, num_domains)
    resize!(controller.velo_hist, num_domains)
    resize!(controller.acce_hist, num_domains)
    resize!(controller.∂Ω_f_hist, num_domains)
    for i in 1:num_domains
        num_steps = round(Int64, controller.time_step / subsims[i].integrator.time_step)
        Δt = controller.time_step / num_steps
        num_stops = num_steps + 1
        subsims[i].integrator.time_step = Δt
        controller.time_hist[i] = Vector{Float64}(undef, num_stops)
        controller.disp_hist[i] = Vector{Vector{Float64}}(undef, num_stops)
        controller.velo_hist[i] = Vector{Vector{Float64}}(undef, num_stops)
        controller.acce_hist[i] = Vector{Vector{Float64}}(undef, num_stops)
        controller.∂Ω_f_hist[i] = Vector{Vector{Float64}}(undef, num_stops)
        for stop in 1:num_stops
            controller.time_hist[i][stop] = controller.prev_time + (stop - 1) * Δt
            controller.disp_hist[i][stop] = copy(controller.stop_disp[i])
            controller.velo_hist[i][stop] = copy(controller.stop_velo[i])
            controller.acce_hist[i][stop] = copy(controller.stop_acce[i])
            controller.∂Ω_f_hist[i][stop] = copy(controller.stop_∂Ω_f[i])
        end
    end
end

function save_history_snapshot(
    controller::MultiDomainController, sims::Vector{SingleDomainSimulation}, subsim_index::Int64, stop_index::Int64
)
    controller.disp_hist[subsim_index][stop_index] = copy(sims[subsim_index].integrator.displacement)
    controller.velo_hist[subsim_index][stop_index] = copy(sims[subsim_index].integrator.velocity)
    controller.acce_hist[subsim_index][stop_index] = copy(sims[subsim_index].integrator.acceleration)
    controller.∂Ω_f_hist[subsim_index][stop_index] = copy(sims[subsim_index].model.internal_force)
    return nothing
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = controller.num_domains
    norms_disp = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i in 1:num_domains
        Δt = controller.time_step
        xᵖʳᵉᵛ = controller.schwarz_disp[i] + Δt * controller.schwarz_velo[i]
        xᶜᵘʳʳ = subsims[i].integrator.displacement + Δt * subsims[i].integrator.velocity
        norms_disp[i] = norm(xᶜᵘʳʳ)
        norms_diff[i] = norm(xᶜᵘʳʳ - xᵖʳᵉᵛ)
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

function check_overlap(model::SolidMechanics, bc::SMContactSchwarzBC)
    distance_tol = 0.0
    parametric_tol = 1.0e-06
    overlap = false
    for node_index in bc.side_set_node_indices
        point = model.current[:, node_index]
        _, ξ, _, coupled_face_node_indices, _, distance = project_point_to_side_set(
            point, bc.coupled_subsim.model, bc.coupled_side_set_id
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

function check_compression(mesh::ExodusDatabase, model::SolidMechanics, bc::SMContactSchwarzBC)
    compression_tol = 0.0
    compression = false
    nodal_reactions = get_dst_traction(bc)
    normals = compute_normal(mesh, bc.side_set_id, model)
    global_from_local_map = get_side_set_global_from_local_map(mesh, bc.side_set_id)
    num_local_nodes = length(global_from_local_map)
    for local_node in 1:num_local_nodes
        nodal_reaction = nodal_reactions[:, local_node]
        normal = normals[:, local_node]
        normal_traction = dot(nodal_reaction, normal)
        compressive_traction = normal_traction ≤ compression_tol
        if compressive_traction == true
            compression = true
            break
        end
    end
    return compression
end

function initialize_transfer_operators(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if !(bc isa SMContactSchwarzBC || bc isa SMNonOverlapSchwarzBC)
                continue
            end
            compute_transfer_operator(subsim.model, bc)
        end
    end
end

function update_transfer_operators(sim::MultiDomainSimulation)
    is_contact = sim.controller.schwarz_contact
    for subsim in sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if !(bc isa SMContactSchwarzBC || bc isa SMNonOverlapSchwarzBC)
                continue
            end
            if is_contact == true || subsim.model.kinematics == Infinitesimal
                compute_transfer_operator(subsim.model, bc)
            end
        end
    end
end

function detect_contact(sim::MultiDomainSimulation)
    if sim.controller.schwarz_contact == false
        return nothing
    end
    num_domains = sim.controller.num_domains
    persistence = sim.controller.active_contact
    contact_domain = falses(num_domains)
    for domain in 1:num_domains
        subsim = sim.subsims[domain]
        mesh = subsim.params["input_mesh"]
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SMContactSchwarzBC
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
            if bc isa SMContactSchwarzBC
                bc.active_contact = sim.controller.active_contact
            end
        end
    end
    if sim.controller.active_contact == true
        norma_log(0, :contact, "Detected")
    end
    resize!(sim.controller.contact_hist, sim.controller.stop + 1)
    sim.controller.contact_hist[sim.controller.stop + 1] = sim.controller.active_contact
    write_scharz_params_csv(sim)
    return sim.controller.active_contact
end

function write_scharz_params_csv(sim::MultiDomainSimulation)
    stop = sim.controller.stop
    csv_interval = get(sim.params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        index_string = "-" * string(stop; pad=4)
        contact_filename = "contact" * index_string * ".csv"
        writedlm(contact_filename, sim.controller.active_contact, '\n')
        iters_filename = "iterations" * index_string * ".csv"
        writedlm(iters_filename, sim.controller.iteration_number, '\n')
        conv_filename = "convergence_values" * index_string * ".csv"
        writedlm(conv_filename, sim.controller.convergence_hist, '\n')
    end
end
