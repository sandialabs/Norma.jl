# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function evolve(sim::Simulation)
    sync_time(sim)
    initialize(sim)
    initialize_writing(sim)
    write_step(sim)
    while true
        advance_time(sim)
        if stop_evolve(sim) == true
            break
        end
        sync_time(sim)
        advance(sim)
        write_step(sim)
    end
    return finalize_writing(sim)
end

function stop_evolve(sim::SingleDomainSimulation)
    return sim.integrator.time > sim.integrator.final_time
end

function stop_evolve(sim::MultiDomainSimulation)
    return sim.schwarz_controller.time > sim.schwarz_controller.final_time
end

function solve_contact(sim::MultiDomainSimulation)
    if sim.schwarz_controller.active_contact == true
        schwarz(sim)
    else
        advance_independent(sim)
    end
end

function decrease_time_step(sim::SingleDomainSimulation)
    time_step = sim.integrator.time_step
    minimum_time_step = sim.integrator.minimum_time_step
    decrease_factor = sim.integrator.decrease_factor
    if decrease_factor == 1.0
        error(
            "Cannot adapt time step ",
            time_step,
            " because decrease factor is ",
            decrease_factor,
            ". Enable adaptive time stepping.",
        )
    end
    new_time_step = decrease_factor * time_step
    if new_time_step < minimum_time_step
        error(
            "Cannot adapt time step to ",
            new_time_step,
            " because minimum is ",
            minimum_time_step,
        )
    end
    sim.integrator.time_step = new_time_step
    @info "Time step failure. Decreasing time step." time_step new_time_step
end

function increase_time_step(sim::SingleDomainSimulation)
    time_step = sim.integrator.time_step
    maximum_time_step = sim.integrator.maximum_time_step
    increase_factor = sim.integrator.increase_factor
    if increase_factor > 1.0
        new_time_step = min(increase_factor * time_step, maximum_time_step)
        if new_time_step > time_step
            sim.integrator.time_step = new_time_step
            @info "Time step success. Increasing time step." time_step new_time_step
        end
    end
end

function advance(sim::SingleDomainSimulation)
    while true
        apply_bcs(sim)
        solve(sim)
        if sim.failed == false
            increase_time_step(sim)
            break
        end
        restore_stop_solutions(sim)
        decrease_time_step(sim)
        advance_time(sim)
        sync_time(sim)
        sim.failed = sim.model.failed = sim.solver.failed = false
    end
end

function advance(sim::MultiDomainSimulation)
    start_runtimer(sim)
    update_transfer_operators(sim)
    if sim.schwarz_controller.schwarz_contact == false
        schwarz(sim)
        return nothing
    end
    save_stop_solutions(sim)
    solve_contact(sim)
    was_in_contact = sim.schwarz_controller.active_contact
    detect_contact(sim)
    if sim.schwarz_controller.active_contact ≠ was_in_contact
        if was_in_contact == true
            println("Contact released — reattempting control step.")
        else
            println("Contact initiated — reattempting control step.")
        end
        restore_stop_solutions(sim)
        solve_contact(sim)
    end
    return end_runtimer(sim)
end

function apply_ics(sim::SingleDomainSimulation)
    return apply_ics(sim.params, sim.model)
end

function apply_ics(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_ics(subsim)
    end
end

function apply_bcs(sim::SingleDomainSimulation)
    return apply_bcs(sim.model)
end

function apply_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_bcs(subsim)
    end
end

function initialize(sim::SingleDomainSimulation)
    apply_ics(sim)
    apply_bcs(sim)
    return initialize(sim.integrator, sim.solver, sim.model)
end

function initialize(sim::MultiDomainSimulation)
    initialize_transfer_operators(sim)
    apply_ics(sim)
    apply_bcs(sim)
    for subsim in sim.subsims
        initialize(subsim.integrator, subsim.solver, subsim.model)
    end
    return detect_contact(sim)
end

function solve(sim::SingleDomainSimulation)
    solve(sim.integrator, sim.solver, sim.model)
    return sim.failed = sim.model.failed
end

function initialize_writing(sim::SingleDomainSimulation)
    return initialize_writing(sim.params, sim.integrator, sim.model)
end

function initialize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        initialize_writing(subsim)
    end
end

function write_step(sim::SingleDomainSimulation)
    return write_step(sim.params, sim.integrator, sim.model)
end

function write_step(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.time
    stop = sim.schwarz_controller.stop
    for subsim in sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = stop
        write_step(subsim)
    end
end

function finalize_writing(sim::SingleDomainSimulation)
    return finalize_writing(sim.params)
end

function finalize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        finalize_writing(subsim)
    end
end

using Printf

function sync_time(sim::SingleDomainSimulation)
    synchronize(sim)
    stop = sim.integrator.stop
    initial_time = sim.integrator.prev_time
    final_time = sim.integrator.time
    if stop == 0
        @printf("Initializing: Stop 0 | Time = %.4e\n", final_time)
    else
        Δt = final_time - initial_time
        @printf(
            "Advancing to: Stop %d | Time = %.4e | Δt = %.4e\n",
            stop,
            final_time,
            Δt
        )
    end
end

function sync_time(sim::MultiDomainSimulation)
    synchronize(sim)
    stop = sim.schwarz_controller.stop
    initial_time = sim.schwarz_controller.prev_time
    final_time = sim.schwarz_controller.time
    if stop == 0
        @printf("Initializing: Stop 0 | Time = %.4e\n", final_time)
    else
        Δt = final_time - initial_time
        @printf(
            "Advancing to: Stop %d | Time = %.4e | Δt = %.4e\n",
            stop,
            final_time,
            Δt
        )
    end
end

function synchronize(sim::SingleDomainSimulation)
    return sim.model.time = sim.integrator.time
end

function synchronize(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.prev_time
    subsim_stop = 0
    for subsim in sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = subsim_stop
    end
end

function advance_time(sim::SingleDomainSimulation)
    if sim.failed == false
        save_stop_solutions(sim)
    end
    sim.integrator.prev_time = sim.integrator.time
    next_time = round(sim.integrator.time + sim.integrator.time_step; digits=12)
    sim.integrator.time = sim.model.time = next_time
    return sim.integrator.stop += 1
end

function advance_time(sim::MultiDomainSimulation)
    sim.schwarz_controller.prev_time = sim.schwarz_controller.time
    stop = sim.schwarz_controller.stop + 1
    final_time = sim.schwarz_controller.final_time
    initial_time = sim.schwarz_controller.initial_time
    num_stops = sim.schwarz_controller.num_stops
    next_time = round(
        (final_time - initial_time) * Float64(stop) / Float64(num_stops - 1) + initial_time;
        digits=12,
    )
    sim.schwarz_controller.time = next_time
    return sim.schwarz_controller.stop = stop
end

function start_runtimer(sim::SingleDomainSimulation)
    return sim.integrator.runtime_step = time()
end

function start_runtimer(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        subsim.integrator.runtime_step = time()
    end
end

function end_runtimer(sim::SingleDomainSimulation)
    return sim.integrator.runtime_step = time() - sim.integrator.runtime_step
end

function end_runtimer(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        subsim.integrator.runtime_step = time() - subsim.integrator.runtime_step
    end
end
