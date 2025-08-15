# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function RomNewmark(params::Parameters, model::RomModel)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    β = integrator_params["β"]
    γ = integrator_params["γ"]
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    fom_newmark = Newmark(params, model.fom_model)
    return RomNewmark(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        β,
        γ,
        displacement,
        velocity,
        acceleration,
        disp_pre,
        velo_pre,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        kinetic_energy,
        fom_newmark,
    )
end

function RomCentralDifference(params::Parameters, model::RomModel)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    CFL = integrator_params["CFL"]
    γ = integrator_params["γ"]
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    stored_energy = 0.0
    kinetic_energy = 0.0
    fom_integrator = CentralDifference(params, model.fom_model)
    return RomCentralDifference(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        CFL,
        γ,
        displacement,
        velocity,
        acceleration,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        kinetic_energy,
        fom_integrator,
    )
end

function create_time_integrator(params::Parameters, model::RomModel)
    integrator_params = params["time integrator"]
    integrator_name = integrator_params["type"]
    if integrator_name == "quasi static"
        norma_abort("This integrator is not supported for ROMs : $integrator_name")
    elseif integrator_name == "Newmark"
        return RomNewmark(params, model)
    elseif integrator_name == "central difference"
        return RomCentralDifference(params, model)
    else
        norma_abort("Unknown type of time integrator : $integrator_name")
    end
end


function initialize(integrator::RomNewmark, solver::RomHessianMinimizer, model::RomModel)
    # Compute initial accelerations
    initialize(integrator.fom_integrator, solver.fom_solver, model.fom_model)

    # project onto basis
    n_var, n_node, n_mode = size(model.basis)
    integrator.displacement[:] = model.reduced_state[:]
    integrator.velocity[:] = model.reduced_velocity[:]
    solver.solution[:] = model.reduced_state[:]

    for k in 1:n_mode
        integrator.acceleration[k] = 0.0
        for j in 1:n_node
            for n in 1:n_var
                integrator.acceleration[k] +=
                    model.basis[n, j, k] * integrator.fom_integrator.acceleration[3 * (j - 1) + n]
            end
        end
    end
    return nothing
end

function predict(integrator::RomNewmark, solver::Solver, model::RomModel)
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

function correct(integrator::RomNewmark, solver::Solver, model::RomModel)
    dt = integrator.time_step
    beta = integrator.β
    gamma = integrator.γ
    integrator.displacement[:] = solver.solution[:]
    integrator.acceleration[:] = (solver.solution - integrator.disp_pre) / (beta * dt * dt)
    integrator.velocity[:] = integrator.velo_pre + dt * gamma * integrator.acceleration[:]
    model.reduced_state[:] = solver.solution[:]
    return nothing
end

function initialize(integrator::RomCentralDifference, solver::RomExplicitSolver, model::RomModel)
    # Compute initial accelerations
    initialize(integrator.fom_integrator, solver.fom_solver, model.fom_model)
    # project onto basis
    n_var, n_node, n_mode = size(model.basis)
    integrator.displacement[:] = model.reduced_state[:]
    integrator.velocity[:] = model.reduced_velocity[:]
    solver.solution[:] = model.reduced_state[:]

    for k in 1:n_mode
        integrator.acceleration[k] = 0.0
        for j in 1:n_node
            for n in 1:n_var
                integrator.acceleration[k] +=
                    model.basis[n, j, k] * integrator.fom_integrator.acceleration[3 * (j - 1) + n]
            end
        end
    end
    return nothing
end

function predict(integrator::RomCentralDifference, solver::RomExplicitSolver, model::RomModel)
    dt = integrator.time_step
    gamma = integrator.γ
    integrator.displacement[:] =
        integrator.displacement[:] + dt * integrator.velocity + 1.0 / 2.0 * dt * dt * integrator.acceleration
    integrator.velocity[:] += dt * (1.0 - gamma) * integrator.acceleration
    model.reduced_state[:] = integrator.displacement[:]
    return nothing
end

function correct(integrator::RomCentralDifference, solver::RomExplicitSolver, model::RomModel)
    dt = integrator.time_step
    gamma = integrator.γ
    # Displacement doesn't change
    integrator.acceleration[:] = solver.solution[:]
    integrator.velocity[:] = integrator.velocity[:] + dt * gamma * integrator.acceleration[:]
    return nothing
end
