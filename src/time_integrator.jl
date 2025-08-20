# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function adaptive_stepping_parameters(integrator_params::Parameters)
    has_minimum = haskey(integrator_params, "minimum time step")
    has_decrease = haskey(integrator_params, "decrease factor")
    has_maximum = haskey(integrator_params, "maximum time step")
    has_increase = haskey(integrator_params, "increase factor")
    has_any = has_minimum || has_decrease || has_maximum || has_increase
    has_all = has_minimum && has_decrease && has_maximum && has_increase
    if has_any == true && has_all == false
        norma_abort(
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
        norma_abortf(
            "Minimum time step %.4e must be less than or equal to maximum %.4e.", minimum_time_step, maximum_time_step
        )
    end
    if decrease_factor > 1.0
        norma_abortf("Decrease factor %.4e must be less than or equal to one.", decrease_factor)
    end
    if increase_factor < 1.0
        norma_abortf("Increase factor %.4e must be greater than or equal to one.", increase_factor)
    end
    return minimum_time_step, decrease_factor, maximum_time_step, increase_factor
end

function QuasiStatic(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    stored_energy = 0.0
    initial_equilibrium = get(integrator_params, "initial equilibrium", false)
    return QuasiStatic(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        displacement,
        velocity,
        acceleration,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        initial_equilibrium,
    )
end


function Newmark(params::Parameters, model::Model)
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
    return Newmark(
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
    )
end


function CentralDifference(params::Parameters, model::Model)
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
    return CentralDifference(
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
    )
end

function create_time_integrator(params::Parameters, model::SolidMechanics)
    integrator_params = params["time integrator"]
    integrator_name = integrator_params["type"]
    if integrator_name == "quasi static"
        return QuasiStatic(params, model)
    elseif integrator_name == "Newmark"
        return Newmark(params, model)
    elseif integrator_name == "central difference"
        return CentralDifference(params, model)
    else
        norma_abort("Unknown type of time integrator : $integrator_name")
    end
end


function is_static(integrator::TimeIntegrator)
    return integrator isa QuasiStatic
end

function is_dynamic(integrator::TimeIntegrator)
    return is_static(integrator) == false
end


function initialize(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    if integrator.initial_equilibrium == true
        norma_log(0, :equilibrium, "Establishing Initial Equilibrium...")
        solve(integrator, solver, model)
        if model.failed == true
            norma_abort("Finite element model failed to establish initial equlibrium")
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
        norma_abort("Finite element model failed to initialize")
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
    atol = solver.linear_solver_absolute_tolerance
    rtol = solver.linear_solver_relative_tolerance
    integrator.acceleration[free] = solve_linear(model.mass[free, free], inertial_force[free], atol, rtol)
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
    u_pre = integrator.disp_pre
    v_pre = integrator.velo_pre
    u_pre[free] = u[free] += Δt * v[free] + (0.5 - β) * Δt * Δt * a[free]
    v_pre[free] = v[free] += (1.0 - γ) * Δt * a[free]
    u_pre[fixed] = u[fixed]
    v_pre[fixed] = v[fixed]
    copy_solution_source_targets(integrator, solver, model)
    return nothing
end

function correct(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement = solver.solution
    u_pre = integrator.disp_pre
    v_pre = integrator.velo_pre
    integrator.acceleration[free] = (u[free] - u_pre[free]) / β / Δt / Δt
    integrator.velocity[free] = v_pre[free] + γ * Δt * integrator.acceleration[free]
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
        norma_abort("The finite element model has failed to initialize")
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

include("opinf/opinf_time_integrator.jl")

