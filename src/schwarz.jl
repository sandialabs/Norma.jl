function SolidSchwarzController(params::Dict{String,Any})
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    absolute_error = relative_error = 0.0
    time = prev_time = initial_time
    same_step = get(params, "same time step for domains", false)
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
    SolidSchwarzController(
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
    )
end

function create_schwarz_controller(params::Dict{String,Any})
    type = params["subdomains type"]
    if type == "static solid mechanics" || type == "dynamic solid mechanics" || type == "dynamic linear opinf rom"
        return SolidSchwarzController(params)
    else
        error("Unknown type of Schwarz controller : ", type)
    end
end

function advance_independent(sim::MultiDomainSimulation)
    sim.schwarz_controller.iteration_number = 0
    save_stop_solutions(sim)
    set_subcycle_times(sim)
    synchronize(sim)
    is_schwarz = false
    subcycle(sim, is_schwarz)
end

function schwarz(sim::MultiDomainSimulation)
    iteration_number = 1
    save_stop_solutions(sim)
    save_schwarz_solutions(sim)
    resize_histories(sim)
    set_subcycle_times(sim)
    is_schwarz = true
    while true
        println("Schwarz iteration=", iteration_number)
        sim.schwarz_controller.iteration_number = iteration_number
        synchronize(sim)
        subcycle(sim, is_schwarz)
        if sim.failed == true
            return
        end
        iteration_number += 1
        ΔX, Δx = update_schwarz_convergence_criterion(sim)
        println("Schwarz criterion |ΔX|=", ΔX, " |ΔX|/|X|=", Δx)
        if stop_schwarz(sim, iteration_number) == true
            println("Performed ", iteration_number - 1, " Schwarz iterations")
            break
        end
        save_schwarz_solutions(sim)
        restore_stop_solutions(sim)
    end
end

function save_stop_solutions(sim::MultiDomainSimulation)
    save_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function save_stop_solutions(
    schwarz_controller::SolidSchwarzController,
    sims::Vector{SingleDomainSimulation},
)

    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.stop_disp[i] = deepcopy(sims[i].integrator.displacement)
        schwarz_controller.stop_velo[i] = deepcopy(sims[i].integrator.velocity)
        schwarz_controller.stop_acce[i] = deepcopy(sims[i].integrator.acceleration)
        if typeof(sims[i].model) == SolidMechanics
          schwarz_controller.stop_∂Ω_f[i] = deepcopy(sims[i].model.internal_force)
        elseif typeof(sims[i].model) == LinearOpInfRom
          schwarz_controller.stop_∂Ω_f[i] = deepcopy(sims[i].model.fom_model.internal_force)
        end
    end
end


function restore_stop_solutions(sim::MultiDomainSimulation)
    restore_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function restore_stop_solutions(
    schwarz_controller::SolidSchwarzController,
    sims::Vector{SingleDomainSimulation},
)
    for i ∈ 1:schwarz_controller.num_domains
        sims[i].integrator.displacement = deepcopy(schwarz_controller.stop_disp[i])
        sims[i].integrator.velocity = deepcopy(schwarz_controller.stop_velo[i])
        sims[i].integrator.acceleration = deepcopy(schwarz_controller.stop_acce[i])
        if typeof(sims[i].model) == SolidMechanics
            sims[i].model.internal_force = deepcopy(schwarz_controller.stop_∂Ω_f[i])
        elseif typeof(sims[i].model) == LinearOpInfRom
            sims[i].model.fom_model.internal_force = deepcopy(schwarz_controller.stop_∂Ω_f[i])
        end
        copy_solution_source_targets(sims[i].integrator, sims[i].solver, sims[i].model)
    end
end

function save_schwarz_solutions(sim::MultiDomainSimulation)
    save_schwarz_solutions(sim.schwarz_controller, sim.subsims)
end

function save_schwarz_solutions(
    schwarz_controller::SolidSchwarzController,
    sims::Vector{SingleDomainSimulation},
)
    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.schwarz_disp[i] = deepcopy(sims[i].integrator.displacement)
        schwarz_controller.schwarz_velo[i] = deepcopy(sims[i].integrator.velocity)
        schwarz_controller.schwarz_acce[i] = deepcopy(sims[i].integrator.acceleration)
    end
end

function set_subcycle_times(sim::MultiDomainSimulation)
    initial_time = sim.schwarz_controller.prev_time
    final_time = sim.schwarz_controller.time
    for subsim ∈ sim.subsims
        subsim.integrator.initial_time = initial_time
        subsim.integrator.time = initial_time
        subsim.integrator.final_time = final_time
    end
end

function subcycle(sim::MultiDomainSimulation, is_schwarz::Bool)
    subsim_index = 1
    for subsim ∈ sim.subsims
        println("subcycle ", subsim.name)
        stop_index = 1
        while true
            advance_time(subsim)
            if stop_evolve(subsim) == true
                break
            end
            subsim.model.time = subsim.integrator.time
            @debug "Subdomain : $(subsim.name), Time : $(subsim.model.time)"
            @debug "Before applying BCs"
            apply_bcs(subsim)
            @debug "After applying BCs and before advance"
            advance(subsim)
            @debug "After advance"
            if subsim.failed == true
                sim.failed = true
                return
            end
            if sim.schwarz_controller.active_contact == true &&
               sim.schwarz_controller.naive_stabilized == true
                apply_naive_stabilized_bcs(subsim)
            end
            stop_index += 1
            if is_schwarz == true
                save_history_snapshot(
                    sim.schwarz_controller,
                    sim.subsims,
                    subsim_index,
                    stop_index,
                )
            end
        end
        subsim_index += 1
    end
end

function resize_histories(sim::MultiDomainSimulation)
    resize_histories(sim.schwarz_controller, sim.subsims)
end

function resize_histories(
    schwarz_controller::SolidSchwarzController,
    sims::Vector{SingleDomainSimulation},
)
    num_domains = schwarz_controller.num_domains
    resize!(schwarz_controller.time_hist, num_domains)
    resize!(schwarz_controller.disp_hist, num_domains)
    resize!(schwarz_controller.velo_hist, num_domains)
    resize!(schwarz_controller.acce_hist, num_domains)
    resize!(schwarz_controller.∂Ω_f_hist, num_domains)
    for subsim ∈ 1:num_domains
        num_steps =
            round(Int64, schwarz_controller.time_step / sims[subsim].integrator.time_step)
        Δt = schwarz_controller.time_step / num_steps
        num_stops = num_steps + 1
        sims[subsim].integrator.time_step = Δt
        schwarz_controller.time_hist[subsim] = Vector{Float64}(undef, num_stops)
        schwarz_controller.disp_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.velo_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.acce_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.∂Ω_f_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        for stop ∈ 1:num_stops
            schwarz_controller.time_hist[subsim][stop] =
                schwarz_controller.prev_time + (stop - 1) * Δt
            schwarz_controller.disp_hist[subsim][stop] =
                deepcopy(schwarz_controller.stop_disp[subsim])
            schwarz_controller.velo_hist[subsim][stop] =
                deepcopy(schwarz_controller.stop_velo[subsim])
            schwarz_controller.acce_hist[subsim][stop] =
                deepcopy(schwarz_controller.stop_acce[subsim])
            schwarz_controller.∂Ω_f_hist[subsim][stop] =
                deepcopy(schwarz_controller.stop_∂Ω_f[subsim])
        end
    end
end

function save_history_snapshot(
    schwarz_controller::SchwarzController,
    sims::Vector{SingleDomainSimulation},
    subsim_index::Int64,
    stop_index::Int64,
)
    schwarz_controller.disp_hist[subsim_index][stop_index] =
        deepcopy(sims[subsim_index].integrator.displacement)
    schwarz_controller.velo_hist[subsim_index][stop_index] =
        deepcopy(sims[subsim_index].integrator.velocity)
    schwarz_controller.acce_hist[subsim_index][stop_index] =
        deepcopy(sims[subsim_index].integrator.acceleration)
    if typeof(sims[subsim_index].model) == SolidMechanics
      schwarz_controller.∂Ω_f_hist[subsim_index][stop_index] =
          deepcopy(sims[subsim_index].model.internal_force)
    elseif typeof(sims[subsim_index].model) == LinearOpInfRom
      schwarz_controller.∂Ω_f_hist[subsim_index][stop_index] =
          deepcopy(sims[subsim_index].model.fom_model.internal_force)
    end
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    return update_schwarz_convergence_criterion(sim.schwarz_controller, sim.subsims)
end

function update_schwarz_convergence_criterion(
    schwarz_controller::SolidSchwarzController,
    sims::Vector{SingleDomainSimulation},
)
    num_domains = schwarz_controller.num_domains
    norms_disp = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i ∈ 1:num_domains
        Δt = schwarz_controller.time_step
        xᵖʳᵉᵛ = schwarz_controller.schwarz_disp[i] + Δt * schwarz_controller.schwarz_velo[i]
        xᶜᵘʳʳ = sims[i].integrator.displacement + Δt * sims[i].integrator.velocity
        norms_disp[i] = norm(xᶜᵘʳʳ)
        norms_diff[i] = norm(xᶜᵘʳʳ - xᵖʳᵉᵛ)
    end
    norm_disp = norm(norms_disp)
    norm_diff = norm(norms_diff)
    schwarz_controller.absolute_error = norm_diff
    schwarz_controller.relative_error = norm_disp > 0.0 ? norm_diff / norm_disp : norm_diff
    conv_abs = schwarz_controller.absolute_error ≤ schwarz_controller.absolute_tolerance
    conv_rel = schwarz_controller.relative_error ≤ schwarz_controller.relative_tolerance
    schwarz_controller.converged = conv_abs || conv_rel
    return schwarz_controller.absolute_error, schwarz_controller.relative_error
end

function stop_schwarz(sim::MultiDomainSimulation, iteration_number::Int64)
    for subsim ∈ sim.subsims
        if subsim.solver.failed == true
            return true
        end
    end
    if sim.schwarz_controller.absolute_error == 0.0
        return true
    end
    exceeds_minimum_iterations =
        iteration_number > sim.schwarz_controller.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations =
        iteration_number > sim.schwarz_controller.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return sim.schwarz_controller.converged
end

function check_overlap(model::SolidMechanics, bc::SMContactSchwarzBC)
    distance_tol = 0.0
    parametric_tol = 1.0e-06
    overlap = false
    for node_index ∈ bc.side_set_node_indices
        point = model.current[:, node_index]
        _, ξ, _, coupled_face_node_indices, _, distance =
            project_point_to_side_set(point, bc.coupled_subsim.model, bc.coupled_side_set_id)
        num_nodes_coupled_side = length(coupled_face_node_indices)
        parametric_dim = length(ξ)
        element_type = get_element_type(parametric_dim, num_nodes_coupled_side)
        overlap =
            distance ≤ distance_tol && is_inside_parametric(element_type, ξ, parametric_tol)
        if overlap == true
            break
        end
    end
    return overlap
end

function check_compression(
    mesh::ExodusDatabase,
    model::SolidMechanics,
    bc::SMContactSchwarzBC,
)
    compression_tol = 0.0
    compression = false
    nodal_reactions = get_dst_traction(bc)
    normals = compute_normal(mesh, bc.side_set_id, model)
    global_from_local_map = get_side_set_global_from_local_map(mesh, bc.side_set_id)
    num_local_nodes = length(global_from_local_map)
    for local_node ∈ 1:num_local_nodes
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
    is_contact = sim.schwarz_controller.schwarz_contact
    for subsim ∈ sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc ∈ bcs
            if typeof(bc) ≠ SMContactSchwarzBC && typeof(bc) ≠ SMNonOverlapSchwarzBC
                continue
            end
            update_transfer_operator(subsim.model, bc)
        end
    end
end

function update_transfer_operators(sim::MultiDomainSimulation)
    is_contact = sim.schwarz_controller.schwarz_contact
    for subsim ∈ sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc ∈ bcs
            if typeof(bc) ≠ SMContactSchwarzBC && typeof(bc) ≠ SMNonOverlapSchwarzBC
                continue
            end
            if is_contact == true || subsim.model.kinematics == Infinitesimal
                update_transfer_operator(subsim.model, bc)
            end
        end
    end
end

function detect_contact(sim::MultiDomainSimulation)
    if sim.schwarz_controller.schwarz_contact == false
        return
    end
    num_domains = sim.schwarz_controller.num_domains
    persistence = sim.schwarz_controller.active_contact
    contact_domain = falses(num_domains)
    for domain ∈ 1:num_domains
        subsim = sim.subsims[domain]
        mesh = subsim.params["input_mesh"]
        bcs = subsim.model.boundary_conditions
        for bc ∈ bcs
            if typeof(bc) == SMContactSchwarzBC
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
    sim.schwarz_controller.active_contact = any(contact_domain)
    println("contact ", sim.schwarz_controller.active_contact)
    resize!(sim.schwarz_controller.contact_hist, sim.schwarz_controller.stop + 1)
    sim.schwarz_controller.contact_hist[sim.schwarz_controller.stop+1] =
        sim.schwarz_controller.active_contact
    write_scharz_params_csv(sim)
    return sim.schwarz_controller.active_contact
end

function write_scharz_params_csv(sim::MultiDomainSimulation)
    stop = sim.schwarz_controller.stop
    csv_interval = get(sim.params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        index_string = "-" * string(stop, pad=4)
        contact_filename = "contact" * index_string * ".csv"
        writedlm(contact_filename, sim.schwarz_controller.active_contact, '\n')
        iters_filename = "iterations" * index_string * ".csv"
        writedlm(iters_filename, sim.schwarz_controller.iteration_number, '\n')
    end
end
