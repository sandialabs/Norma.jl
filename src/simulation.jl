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
include("schwarz.jl")
# For OpInf models
include("opinf/opinf_model.jl")


function create_simulation(input_file::String)
    norma_log(0, :setup, "Reading from " * input_file)
    params = YAML.load_file(input_file; dicttype=Parameters)
    params["name"] = input_file
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
        error("Unknown type of simulation: ", sim_type)
    end
end

function create_bcs(sim::SingleDomainSimulation)
    boundary_conditions = create_bcs(sim.params)
    for bc in boundary_conditions
        if bc isa SMDirichletInclined || bc isa SMContactSchwarzBC
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
    return pair_schwarz_bcs(sim)
end

function SingleDomainSimulation(params::Parameters)
    name = params["name"]
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
    return SingleDomainSimulation(name, params, controller, integrator, solver, model, failed)
end

function MultiDomainSimulation(params::Parameters)
    name = params["name"]
    domain_names = params["domains"]
    subsims = Vector{SingleDomainSimulation}()
    controller = create_controller(params)
    initial_time = controller.initial_time
    final_time = controller.final_time
    time_step = controller.time_step
    same_step = true
    exodus_interval = get(params, "Exodus output interval", 1)
    csv_interval = get(params, "CSV output interval", 0)
    subsim_name_index_map = Dict{String,Int64}()
    subsim_index = 1
    for domain_name in domain_names
        norma_log(4, :domain, domain_name)
        subparams = YAML.load_file(domain_name; dicttype=Parameters)
        subparams["name"] = domain_name
        ti_params = subparams["time integrator"]
        subsim_time_step = get(ti_params, "time step", time_step)
        subsim_time_step = min(subsim_time_step, time_step)
        ti_params["initial time"] = initial_time
        ti_params["final time"] = final_time
        ti_params["time step"] = subsim_time_step
        subparams["Exodus output interval"] = exodus_interval
        subparams["CSV output interval"] = csv_interval
        subsim = SingleDomainSimulation(subparams)
        if subsim.integrator.time_step ≠ time_step
            same_step = false
        end
        params[domain_name] = subsim.params
        push!(subsims, subsim)
        subsim_name_index_map[domain_name] = subsim_index
        subsim_index += 1
    end
    controller.same_step = same_step
    failed = false
    sim = MultiDomainSimulation(name, params, controller, subsims, subsim_name_index_map, failed)
    for subsim in sim.subsims
        subsim.params["global_simulation"] = sim
    end
    return sim
end

function SolidMultiDomainController(params::Parameters)
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
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
    relaxation_parameter = get(params, "relaxation parameter", 1.0)
    naive_stabilized = get(params, "naive stabilized", false)
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

function SolidSingleController(params::Parameters)
    ti_params = params["time integrator"]
    initial_time = ti_params["initial time"]
    final_time = ti_params["final time"]
    time_step = ti_params["time step"]
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
    time = prev_time = initial_time
    stop = 0
    stop_disp = Vector{Float64}()
    stop_velo = Vector{Float64}()
    stop_acce = Vector{Float64}()
    stop_∂Ω_f = Vector{Float64}()
    return SolidSingleController(
        initial_time,
        final_time,
        time_step,
        time,
        prev_time,
        num_stops,
        stop,
        stop_disp,
        stop_velo,
        stop_acce,
        stop_∂Ω_f,
    )
end

function create_controller(params::Parameters)
    sim_type = params["type"]
    if sim_type == "single"
        return SolidSingleController(params)
    elseif sim_type == "multi"
        return SolidMultiDomainController(params)
    else
        error("Unknown type of simulation: ", sim_type)
    end
end
