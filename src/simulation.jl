# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Printf
using YAML


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
        initialize_storage(sim)
        return sim
    elseif sim_type == "multi"
        sim = MultiDomainSimulation(params)
        create_bcs(sim)
        initialize_storage(sim)
        return sim
    else
        norma_abort("Unknown type of simulation: $sim_type")
    end
end

function initialize_storage(sim::SingleDomainSimulation)
    model = sim.model
    integrator = sim.integrator
    solver = sim.solver
    model isa SolidMechanics || return nothing
    num_dof = length(model.free_dofs)
    integrator.displacement = unsafe_wrap(Vector{Float64}, pointer(model.displacement), num_dof)
    integrator.velocity = unsafe_wrap(Vector{Float64}, pointer(model.velocity), num_dof)
    integrator.acceleration = unsafe_wrap(Vector{Float64}, pointer(model.acceleration), num_dof)
    if solver isa ExplicitSolver
        solver.solution = unsafe_wrap(Vector{Float64}, pointer(model.acceleration), num_dof)
    else
        solver.solution = unsafe_wrap(Vector{Float64}, pointer(model.displacement), num_dof)
    end
    return nothing
end

function initialize_storage(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        initialize_storage(subsim)
    end
    return nothing
end

function SingleDomainSimulation(params::Parameters)
    t_setup = time()
    basename = params["name"]
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    norma_log(0, :setup, "Input:  $input_mesh_file")
    norma_log(0, :setup, "Output: $output_mesh_file")
    rm(output_mesh_file; force=true)
    input_mesh = Exodus.ExodusDatabase(input_mesh_file, "r")
    Exodus.copy(input_mesh, output_mesh_file)
    output_mesh = Exodus.ExodusDatabase(output_mesh_file, "rw")
    params["output_mesh"] = output_mesh
    params["input_mesh"] = input_mesh
    n_nodes = Exodus.num_nodes(input_mesh.init)
    n_elems = Exodus.num_elements(input_mesh.init)
    norma_logf(0, :setup, "Mesh:   %d nodes, %d elements", n_nodes, n_elems)
    controller = create_controller(params)
    norma_logf(0, :setup, "Time:   [%.2e, %.2e], Δt = %.2e, %d steps",
               controller.initial_time, controller.final_time,
               controller.time_step, controller.num_stops - 1)
    model = create_model(params)
    _log_materials(model, input_mesh)
    integrator = create_time_integrator(params, model)
    solver = create_solver(params, model)
    norma_logf(0, :setup, "Solver: %s, %s", _integrator_name(integrator), _solver_name(solver))
    norma_log(0, :setup, "Setup complete ($(format_time(time() - t_setup)))")
    failed = false
    return SingleDomainSimulation(basename, params, controller, integrator, solver, model, failed, nothing, nothing)
end

# True when this subsim belongs to a MultiDomainSimulation (i.e., is coupled).
is_coupled(sim::SingleDomainSimulation) = sim.parent !== nothing

# Resolve a domain name to its stable handle within a MultiDomainSimulation.
handle(sim::MultiDomainSimulation, name::String) = sim.handle_by_name[name]

_integrator_name(::QuasiStatic) = "Quasi-static"
_integrator_name(::Newmark) = "Newmark"
_integrator_name(::CentralDifference) = "Central difference"
_integrator_name(ig) = replace(string(typeof(ig)), r"^.*\." => "")

_solver_name(::HessianMinimizer) = "Newton (Hessian minimizer)"
_solver_name(::SteepestDescent) = "Steepest descent"
_solver_name(::ExplicitSolver) = "Explicit"
_solver_name(s) = replace(string(typeof(s)), r"^.*\." => "")

function log_dof_counts(model::SolidMechanics)
    n_dofs = length(model.free_dofs)
    n_free = count(model.free_dofs)
    norma_logf(0, :setup, "DOFs:   %d total, %d free, %d constrained",
               n_dofs, n_free, n_dofs - n_free)
end

log_dof_counts(::Model) = nothing

_material_name(m) = replace(string(typeof(m)), r"^.*\." => "")

# Format a single material struct field as "name = value".  Floats get
# scientific notation for readability; integers (e.g. SethHill's m, n)
# print as-is.
function _format_material_field(name::Symbol, value)
    if value isa AbstractFloat
        return @sprintf("%s = %.3e", String(name), value)
    else
        return string(name) * " = " * string(value)
    end
end

function _log_materials(model::SolidMechanics, input_mesh)
    block_names = Exodus.read_names(input_mesh, Block)
    for (block_name, material) in zip(block_names, model.materials)
        props = join((_format_material_field(f, getfield(material, f))
                      for f in fieldnames(typeof(material))), ", ")
        norma_log(0, :setup, "Block \"$block_name\": $(_material_name(material)) : $props")
    end
end

_log_materials(::Model, _) = nothing

function MultiDomainSimulation(params::Parameters)
    basename = params["name"]
    domain_paths = params["domains"]
    subsims = Vector{SingleDomainSimulation}()
    controller = create_controller(params)
    initial_time = controller.initial_time
    final_time = controller.final_time
    integrator_dt = params["time step"]
    exodus_interval = Float64(get(params, "Exodus output interval", integrator_dt))
    csv_interval = Float64(get(params, "CSV output interval", 0.0))
    handle_by_name = Dict{String,DomainHandle}()
    name_by_handle = String[]
    subsim_index = 1
    for domain_path in domain_paths
        domain_name = stripped_name(domain_path)
        norma_log(4, :domain, domain_name)
        subparams = YAML.load_file(domain_path; dicttype=Parameters)
        subparams["name"] = domain_name
        integrator_params = subparams["time integrator"]
        subsim_time_step = get(integrator_params, "time step", integrator_dt)
        subsim_time_step = min(subsim_time_step, integrator_dt)
        integrator_params["initial time"] = initial_time
        integrator_params["final time"] = final_time
        integrator_params["time step"] = subsim_time_step
        subparams["Exodus output interval"] = exodus_interval
        subparams["CSV output interval"] = csv_interval
        subsim = SingleDomainSimulation(subparams)
        params[domain_name] = subsim.params
        push!(subsims, subsim)
        handle_by_name[domain_name] = DomainHandle(subsim_index)
        push!(name_by_handle, domain_name)
        subsim_index += 1
    end
    num_domains = length(subsims)
    swaps = parse_swap_plans(params)
    failed = false
    sim = MultiDomainSimulation(
        basename, params, controller, num_domains, subsims,
        handle_by_name, name_by_handle, swaps, failed,
    )
    for (i, subsim) in enumerate(sim.subsims)
        subsim.parent = sim
        subsim.handle = DomainHandle(i)
    end
    validate_swap_plans(sim)
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
    # Multi-domain: control step = integrator dt (Schwarz coupling each step).
    # Output timing is handled by _is_output_time in write_stop.
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
    has_relaxation_key = haskey(params, "relaxation")
    has_relaxation_parameter = haskey(params, "relaxation parameter")
    if has_relaxation_key && has_relaxation_parameter
        norma_abort(
            "Schwarz controller: specify either `relaxation: aitken` or `relaxation parameter: <float>`, not both.",
        )
    end
    relaxation_method = :fixed
    relaxation_parameter = 1.0
    if has_relaxation_key
        relaxation_value = params["relaxation"]
        if relaxation_value isa AbstractString && lowercase(relaxation_value) == "aitken"
            relaxation_method = :aitken
        else
            norma_abort("Schwarz controller: unsupported `relaxation: $(relaxation_value)` (only `aitken` is recognized).")
        end
    elseif has_relaxation_parameter
        relaxation_parameter = Float64(params["relaxation parameter"])
    end
    naive_stabilized = get(params, "naive stabilized", false)
    lambda_disp = [Vector{Float64}[] for _ in 1:num_domains]
    lambda_velo = [Vector{Float64}[] for _ in 1:num_domains]
    lambda_acce = [Vector{Float64}[] for _ in 1:num_domains]
    aitken_prev_residual_disp = [Vector{Float64}() for _ in 1:num_domains]
    aitken_prev_residual_velo = [Vector{Float64}() for _ in 1:num_domains]
    aitken_prev_residual_acce = [Vector{Float64}() for _ in 1:num_domains]
    aitken_theta_disp = ones(Float64, num_domains)
    is_schwarz = true
    schwarz_contact = false
    active_contact = false
    contact_hist = Vector{Bool}[]
    schwarz_iters = zeros(Int64, num_stops - 1)
    use_interface_predictor = get(params, "interface predictor", false)
    predictor_disp = [Vector{Float64}() for _ in 1:num_domains]
    predictor_velo = [Vector{Float64}() for _ in 1:num_domains]
    predictor_acce = [Vector{Float64}() for _ in 1:num_domains]
    predictor_∂Ω_f = [Vector{Float64}() for _ in 1:num_domains]
    prev_stop_disp = [Vector{Float64}() for _ in 1:num_domains]
    prev_stop_∂Ω_f = [Vector{Float64}() for _ in 1:num_domains]

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
        relaxation_method,
        naive_stabilized,
        lambda_disp,
        lambda_velo,
        lambda_acce,
        aitken_prev_residual_disp,
        aitken_prev_residual_velo,
        aitken_prev_residual_acce,
        aitken_theta_disp,
        is_schwarz,
        schwarz_contact,
        active_contact,
        contact_hist,
        schwarz_iters,
        use_interface_predictor,
        predictor_disp,
        predictor_velo,
        predictor_acce,
        predictor_∂Ω_f,
        prev_stop_disp,
        prev_stop_∂Ω_f,
    )
end

function SolidSingleDomainTimeController(params::Parameters)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    # Control step = integrator dt always. Output timing is handled
    # by _is_output_time in write_stop, not by controller stops.
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
    boundary_conditions = _create_bcs(sim)
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
    t_batch = time()
    while true
        advance_control_time(sim)
        maybe_apply_swaps!(sim)
        sync_control_time(sim)
        advance_control(sim)
        write_stop(sim; wall_time=time() - t_batch)
        t_batch = time()
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
            "Cannot adapt time step %.2e because decrease factor is %.2e. Enable adaptive time stepping.",
            time_step,
            decrease_factor,
        )
    end
    new_time_step = decrease_factor * time_step
    minimum_time_step = sim.integrator.minimum_time_step
    if new_time_step < minimum_time_step
        norma_abortf("Cannot adapt time step to %.2e because minimum is %.2e.", new_time_step, minimum_time_step)
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
    is_explicit = sim.integrator isa ExplicitDynamicTimeIntegrator
    while true
        prev_time = sim.integrator.prev_time
        time = sim.integrator.time
        time_step = sim.integrator.time_step
        # For implicit, print each step's time interval (Newton iters follow).
        # For explicit, suppress — only output stops are printed.
        if !is_explicit
            norma_logf(4, :advance, "Time = [%.2e, %.2e] : Δt = %.2e", prev_time, time, time_step)
        end
        apply_bcs(sim)
        # Defense-in-depth: constitutive model evaluate() catches math errors
        # and sets model.failed.  This outer catch handles the same errors if
        # they escape from any other point in the solve path.
        try
            solve(sim)
        catch e
            e isa _MATH_ERRORS || rethrow()
            norma_logf(4, :solve, "Caught %s during solve — treating as step failure", typeof(e))
            sim.failed = sim.model.failed = sim.solver.failed = true
        end
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
    log_dof_counts(sim.model)
    initialize(sim.integrator, sim.solver, sim.model)
    save_curr_state(sim)
    return nothing
end

function initialize(sim::MultiDomainSimulation)
    initialize_bc_projectors(sim)
    apply_ics(sim)
    for (subsim_index, subsim) in enumerate(sim.subsims)
        if subsim.model isa SolidMechanics
            save_history_snapshot(sim.controller, subsim, subsim_index)
        end
    end
    apply_bcs(sim)
    for subsim in sim.subsims
        log_dof_counts(subsim.model)
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

function get_adjusted_timestep(t::Float64, dt::Float64, t_stop::Float64)::Float64
    gap = t_stop - t
    gap <= 0.0 && return 0.0
    t_next = t + dt
    return t_stop - t_next <= 0.5 * dt ? gap : dt
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
    iteration_number = 0
    sim.controller.is_schwarz = true
    save_stop_state(sim)
    save_schwarz_state(sim)
    reset_histories(sim)
    swap_swappable_bcs(sim)

    if sim.controller.use_interface_predictor
        compute_interface_predictor!(sim)
    end

    while true
        norma_log(0, :schwarz, "Iteration [$iteration_number]")
        sim.controller.iteration_number = iteration_number
        set_initial_subcycle_time(sim)
        subcycle(sim)
        ΔU, Δu = update_schwarz_convergence_criterion(sim)
        if iteration_number == 0
            # Initial Schwarz pass: there is no prior iterate to compare against,
            # so the relative criterion is not yet a Schwarz convergence measure.
            # Report the absolute update for reference and skip the relative test.
            sim.controller.converged = false
            norma_logf(
                0,
                :schwarz,
                "Criterion [%d] |ΔU| = %.2e : |ΔU|/|U| = %8s : %s",
                iteration_number,
                ΔU,
                "—",
                colored_status("[INITIAL]"),
            )
            # Early exit: if the initial absolute update already meets the absolute
            # tolerance, no further Schwarz iterations are needed. The relative test
            # still cannot be applied on iteration 0 (no prior iterate).
            if ΔU ≤ sim.controller.absolute_tolerance
                norma_log(0, :schwarz, "Performed 0 Schwarz Iterations")
                sim.controller.schwarz_iters[sim.controller.stop] = 0
                break
            end
        else
            raw_status = sim.controller.converged ? "[CONVERGED]" : "[CONVERGING]"
            status = colored_status(raw_status)
            norma_logf(
                0,
                :schwarz,
                "Criterion [%d] |ΔU| = %.2e : |ΔU|/|U| = %.2e : %s",
                iteration_number,
                ΔU,
                Δu,
                status,
            )
            if stop_schwarz(sim, iteration_number + 1) == true
                plural = iteration_number == 1 ? "" : "s"
                norma_log(0, :schwarz, "Performed $iteration_number Schwarz Iteration" * plural)
                sim.controller.schwarz_iters[sim.controller.stop] = iteration_number
                break
            end
        end
        iteration_number += 1
        save_schwarz_state(sim)
        restore_stop_state(sim)
    end
end


function save_curr_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    integrator.prev_disp = copy(integrator.displacement)
    integrator.prev_velo = copy(integrator.velocity)
    integrator.prev_acce = copy(integrator.acceleration)
    return integrator.prev_∂Ω_f = copy(sim.model.internal_force)
end

function restore_prev_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    integrator.time = integrator.prev_time
    integrator.displacement .= integrator.prev_disp
    integrator.velocity .= integrator.prev_velo
    integrator.acceleration .= integrator.prev_acce
    sim.model.internal_force = copy(integrator.prev_∂Ω_f)
    return nothing
end

function save_stop_state(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    subsims = sim.subsims
    if controller.use_interface_predictor
        for i in 1:num_domains
            if !isempty(controller.stop_disp[i])
                controller.prev_stop_disp[i] = controller.stop_disp[i]
            end
            if !isempty(controller.stop_∂Ω_f[i])
                controller.prev_stop_∂Ω_f[i] = controller.stop_∂Ω_f[i]
            end
        end
    end
    for i in 1:num_domains
        subsim = subsims[i]
        controller.stop_disp[i] = copy(subsim.integrator.displacement)
        controller.stop_velo[i] = copy(subsim.integrator.velocity)
        controller.stop_acce[i] = copy(subsim.integrator.acceleration)
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
        subsim.integrator.displacement .= controller.stop_disp[i]
        subsim.integrator.velocity .= controller.stop_velo[i]
        subsim.integrator.acceleration .= controller.stop_acce[i]
        subsim.model.internal_force = copy(controller.stop_∂Ω_f[i])
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
    is_explicit = sim.integrator isa ExplicitDynamicTimeIntegrator
    t_last_log = time()
    while true
        advance_time(sim)
        advance_one_step(sim)
        if is_explicit && time() - t_last_log >= 1.0
            norma_logf(4, :progress, "Time = %.2e : Δt = %.2e",
                       sim.integrator.time, sim.integrator.time_step)
            t_last_log = time()
        end
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
    norms_pos = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i in 1:num_domains
        Δt = controller.time_step
        u_prev = controller.schwarz_disp[i] + Δt * controller.schwarz_velo[i]
        u_curr = subsims[i].integrator.displacement + Δt * subsims[i].integrator.velocity
        if subsims[i].model isa SolidMechanics
            X = vec(subsims[i].model.reference)
            norms_pos[i] = norm(X + u_curr)
        else
            norms_pos[i] = norm(u_curr)
        end
        norms_diff[i] = norm(u_curr - u_prev)
    end
    norm_pos = norm(norms_pos)
    norm_diff = norm(norms_diff)
    controller.absolute_error = norm_diff
    controller.relative_error = norm_pos > 0.0 ? norm_diff / norm_pos : norm_diff
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
    coupled_model = coupled_subsim_of(bc).model
    coupled_bc = coupled_model.boundary_conditions[bc.coupled_bc_index]
    coupled_side_set_id = coupled_bc.side_set_id

    # Tolerance value to calculate closest point projection, to avoid projection failure
    # Set to 10 times the minimum characteristic edge length of the side set nodes
    # TODO (BRP): perhaps calculate this "characteristic size" a priori, and only recalculate for finite kinematics
    nodal_points = model.reference[:, unique_node_indices] + model.displacement[:, unique_node_indices]
    # Compute minimum distances between the nodes in the side set
    max_distance_tolerance =
        10 * minimum([
            minimum([norm(p1 - p2) for (j, p2) in enumerate(eachcol(nodal_points)) if i != j]) for
            (i, p1) in enumerate(eachcol(nodal_points))
        ])

    for node_index in unique_node_indices
        point = model.reference[:, node_index] + model.displacement[:, node_index]
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
            if bc isa SolidMechanicsRobinSchwarzBoundaryCondition
                compute_robin_schwarz_projectors!(subsim.model, bc)
            elseif bc isa SolidMechanicsImpedanceSchwarzBoundaryCondition
                compute_impedance_schwarz_projectors!(subsim.model, bc)
            elseif bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
                compute_impedance_overlap_schwarz_projectors!(subsim.model, bc)
            elseif bc isa SolidMechanicsOverlapSchwarzBoundaryCondition && bc.use_weak
                coupled_model = get_fom_model(coupled_subsim_of(bc))
                W = get_square_projection_matrix(subsim.model, bc)
                L = get_overlap_rectangular_projection_matrix(
                    subsim.model, bc, coupled_model, bc.coupled_block_name, bc.search_tolerance
                )
                bc.dirichlet_projector = (W \ I) * L
            elseif bc isa SolidMechanicsContactSchwarzBoundaryCondition ||
                   bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                compute_dirichlet_projector(subsim.model, bc)
                compute_neumann_projector(subsim.model, bc)
                if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                    bc.square_projector = get_square_projection_matrix(subsim.model, bc)
                end
            end
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

function extract_interface_state(global_state::Vector{Float64}, bc::SolidMechanicsSchwarzBoundaryCondition)
    n_local = length(bc.global_from_local_map)
    local_mat = Matrix{Float64}(undef, 3, n_local)
    for (li, gi) in enumerate(bc.global_from_local_map)
        @inbounds local_mat[:, li] = global_state[(3gi - 2):(3gi)]
    end
    return local_mat
end

function write_interface_state!(global_state::Vector{Float64}, local_mat::Matrix{Float64}, bc::SolidMechanicsSchwarzBoundaryCondition)
    for (li, gi) in enumerate(bc.global_from_local_map)
        @inbounds global_state[(3gi - 2):(3gi)] = local_mat[:, li]
    end
    return nothing
end

function compute_interface_predictor!(sim::MultiDomainSimulation)
    controller = sim.controller
    processed = Set{Tuple{Int,Int}}()

    for (dom_k, subsim_k) in enumerate(sim.subsims)
        for bc_k in subsim_k.model.boundary_conditions
            bc_k isa SolidMechanicsNonOverlapSchwarzBoundaryCondition ||
            bc_k isa SolidMechanicsRobinSchwarzBoundaryCondition || continue

            dom_j = bc_k.coupled_handle.id
            pair = minmax(dom_k, dom_j)
            pair ∈ processed && continue
            push!(processed, pair)

            subsim_j = sim.subsims[dom_j]
            bc_j = subsim_j.model.boundary_conditions[bc_k.coupled_bc_index]

            # Projectors: P_kj maps from domain j's interface space to domain k's
            #             P_jk maps from domain k's interface space to domain j's
            P_kj = bc_k.dirichlet_projector
            P_jk = bc_j.dirichlet_projector

            # Start predictors from the current stop state (t_n)
            pred_disp_k = copy(controller.stop_disp[dom_k])
            pred_velo_k = copy(controller.stop_velo[dom_k])
            pred_acce_k = copy(controller.stop_acce[dom_k])
            pred_disp_j = copy(controller.stop_disp[dom_j])
            pred_velo_j = copy(controller.stop_velo[dom_j])
            pred_acce_j = copy(controller.stop_acce[dom_j])

            both_static = is_static(subsim_k.integrator) && is_static(subsim_j.integrator)

            if both_static
                # Quasi-static: velocities and accelerations are always zero.
                # Displacement predictor: linear extrapolation in time when prior step exists.
                isempty(controller.prev_stop_disp[dom_k]) && continue
                isempty(controller.prev_stop_disp[dom_j]) && continue

                u_k = extract_interface_state(controller.stop_disp[dom_k], bc_k)
                u_k_prev = extract_interface_state(controller.prev_stop_disp[dom_k], bc_k)
                u_j = extract_interface_state(controller.stop_disp[dom_j], bc_j)
                u_j_prev = extract_interface_state(controller.prev_stop_disp[dom_j], bc_j)

                # Each domain extrapolates its own interface displacement independently.
                # This correctly handles non-conforming meshes: domain k predicts its own
                # interface motion, and domain j predicts its own interface motion.
                u_k_pred = @. 2 * u_k - u_k_prev
                u_j_pred = @. 2 * u_j - u_j_prev

                write_interface_state!(pred_disp_k, u_k_pred, bc_k)
                write_interface_state!(pred_disp_j, u_j_pred, bc_j)
            else
                # Dynamic: velocity predictor via perfectly inelastic collision using
                # the interface mass matrices W (row-sum lumping).
                # W = ∫ N Nᵀ dΓ is purely geometric; it works for both explicit
                # (lumped_mass) and implicit (consistent mass) time integrators.
                W_k = bc_k.square_projector    # (n_k × n_k)
                W_j = bc_j.square_projector    # (n_j × n_j)

                # Lumped interface masses: row sums of W
                m_k = vec(sum(W_k; dims=2))    # (n_k,)
                m_j = vec(sum(W_j; dims=2))    # (n_j,)

                # Interface velocities at t_n from stop state
                v_k = extract_interface_state(controller.stop_velo[dom_k], bc_k)  # (3, n_k)
                v_j = extract_interface_state(controller.stop_velo[dom_j], bc_j)  # (3, n_j)

                # Project domain j's mass and velocity into domain k's interface space
                m_j_in_k = P_kj * m_j    # (n_k,)
                v_j_in_k = Matrix{Float64}(undef, 3, size(P_kj, 1))
                for d in 1:3
                    v_j_in_k[d, :] = P_kj * v_j[d, :]
                end

                # Perfectly inelastic collision: conserves linear momentum
                m_total = m_k .+ m_j_in_k
                v_pred_k = Matrix{Float64}(undef, 3, length(m_k))
                for d in 1:3
                    @. v_pred_k[d, :] = (m_k * v_k[d, :] + m_j_in_k * v_j_in_k[d, :]) / m_total
                end

                # Project predicted velocity from domain k's space back to domain j's space
                v_pred_j = Matrix{Float64}(undef, 3, size(P_jk, 1))
                for d in 1:3
                    v_pred_j[d, :] = P_jk * v_pred_k[d, :]
                end

                # Displacement predictor: velocity-consistent linear extrapolation
                # u_pred = u_n + Δt * v_pred  (each domain uses its own u_n as baseline)
                Δt = controller.time_step
                u_k = extract_interface_state(controller.stop_disp[dom_k], bc_k)
                u_j = extract_interface_state(controller.stop_disp[dom_j], bc_j)
                u_pred_k = @. u_k + Δt * v_pred_k
                u_pred_j = @. u_j + Δt * v_pred_j

                # Acceleration predictor: zero at the interface (standard Newmark prediction)
                # pred_acce_k and pred_acce_j already contain zero-valued stop_acce (from t_n
                # corrected acceleration); no modification needed here.

                write_interface_state!(pred_disp_k, u_pred_k, bc_k)
                write_interface_state!(pred_velo_k, v_pred_k, bc_k)
                write_interface_state!(pred_disp_j, u_pred_j, bc_j)
                write_interface_state!(pred_velo_j, v_pred_j, bc_j)
            end

            # Force predictor: linear temporal extrapolation f_pred = 2*f_n - f_{n-1}.
            # Which domain's force is read by the coupling:
            #   - Non-overlap Neumann (is_dirichlet == false) reads force FROM the coupled domain.
            #   - Robin reads force from BOTH coupled domains.
            # So set predictor_∂Ω_f[x] when domain x's force is read by the other domain.
            is_robin = bc_k isa SolidMechanicsRobinSchwarzBoundaryCondition

            pred_∂Ω_f_k = copy(controller.stop_∂Ω_f[dom_k])
            pred_∂Ω_f_j = copy(controller.stop_∂Ω_f[dom_j])

            if (is_robin || !bc_j.is_dirichlet) && !isempty(controller.prev_stop_∂Ω_f[dom_k])
                @. pred_∂Ω_f_k = 2 * controller.stop_∂Ω_f[dom_k] - controller.prev_stop_∂Ω_f[dom_k]
            end
            if (is_robin || !bc_k.is_dirichlet) && !isempty(controller.prev_stop_∂Ω_f[dom_j])
                @. pred_∂Ω_f_j = 2 * controller.stop_∂Ω_f[dom_j] - controller.prev_stop_∂Ω_f[dom_j]
            end

            controller.predictor_disp[dom_k] = pred_disp_k
            controller.predictor_velo[dom_k] = pred_velo_k
            controller.predictor_acce[dom_k] = pred_acce_k
            controller.predictor_∂Ω_f[dom_k] = pred_∂Ω_f_k
            controller.predictor_disp[dom_j] = pred_disp_j
            controller.predictor_velo[dom_j] = pred_velo_j
            controller.predictor_acce[dom_j] = pred_acce_j
            controller.predictor_∂Ω_f[dom_j] = pred_∂Ω_f_j
        end
    end
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
