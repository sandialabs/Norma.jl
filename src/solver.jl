# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using LinearAlgebra

function create_step(solver_params::Parameters)
    step_name = solver_params["step"]
    if step_name == "full Newton"
        return NewtonStep(solver_params)
    elseif step_name == "explicit"
        return ExplicitStep(solver_params)
    elseif step_name == "steepest descent"
        return SteepestDescentStep(solver_params)
    else
        error("Unknown type of solver step: ", step_name)
    end
end

function HessianMinimizer(params::Parameters, model::Model)
    solver_params = params["solver"]
    num_dof = length(model.free_dofs)
    minimum_iterations = solver_params["minimum iterations"]
    maximum_iterations = solver_params["maximum iterations"]
    absolute_tolerance = solver_params["absolute tolerance"]
    relative_tolerance = solver_params["relative tolerance"]
    absolute_error = 0.0
    relative_error = 0.0
    value = 0.0
    gradient = zeros(num_dof)
    hessian = spzeros(num_dof, num_dof)
    solution = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    ls_backtrack_factor = 0.5
    ls_decrease_factor = 1.0e-04
    ls_max_iters = 16
    if haskey(solver_params, "line search backtrack factor")
        ls_backtrack_factor = solver_params["line search backtrack factor"]
    end
    if haskey(solver_params, "line search decrease factor")
        ls_decrease_factor = solver_params["line search decrease factor"]
    end
    if haskey(solver_params, "line search maximum iterations")
        ls_max_iters = solver_params["line search maximum iterations"]
    end
    line_search = BackTrackLineSearch(ls_backtrack_factor, ls_decrease_factor, ls_max_iters)
    HessianMinimizer(
        minimum_iterations,
        maximum_iterations,
        absolute_tolerance,
        relative_tolerance,
        absolute_error,
        relative_error,
        value,
        gradient,
        hessian,
        solution,
        initial_norm,
        converged,
        failed,
        step,
        line_search,
    )
end

function ExplicitSolver(params::Parameters, model::Model)
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    num_dof = length(model.free_dofs)
    value = 0.0
    gradient = zeros(num_dof)
    solution = zeros(num_dof)
    initial_guess = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    ExplicitSolver(
        value,
        gradient,
        solution,
        initial_guess,
        initial_norm,
        converged,
        failed,
        step,
    )
end

function SteepestDescent(params::Parameters, model::Model)
    solver_params = params["solver"]
    num_dof = length(model.free_dofs)
    minimum_iterations = solver_params["minimum iterations"]
    maximum_iterations = solver_params["maximum iterations"]
    absolute_tolerance = solver_params["absolute tolerance"]
    relative_tolerance = solver_params["relative tolerance"]
    absolute_error = 0.0
    relative_error = 0.0
    value = 0.0
    gradient = zeros(num_dof)
    solution = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    ls_backtrack_factor = 0.5
    ls_decrease_factor = 1.0e-04
    ls_max_iters = 16
    if haskey(solver_params, "line search backtrack factor")
        ls_backtrack_factor = solver_params["line search backtrack factor"]
    end
    if haskey(solver_params, "line search decrease factor")
        ls_decrease_factor = solver_params["line search decrease factor"]
    end
    if haskey(solver_params, "line search maximum iterations")
        ls_max_iters = solver_params["line search maximum iterations"]
    end
    line_search = BackTrackLineSearch(ls_backtrack_factor, ls_decrease_factor, ls_max_iters)
    SteepestDescent(
        minimum_iterations,
        maximum_iterations,
        absolute_tolerance,
        relative_tolerance,
        absolute_error,
        relative_error,
        value,
        gradient,
        solution,
        initial_norm,
        converged,
        failed,
        step,
        line_search,
    )
end

function NewtonStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    NewtonStep(step_length)
end

function ExplicitStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    ExplicitStep(step_length)
end

function SteepestDescentStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    SteepestDescentStep(step_length)
end

function create_solver(params::Parameters, model::Model)
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params, model)
    elseif solver_name == "explicit solver"
        return ExplicitSolver(params, model)
    elseif solver_name == "steepest descent"
        return SteepestDescent(params, model)
    else
        error("Unknown type of solver : ", solver_name)
    end
end

function copy_solution_source_targets(
    integrator::QuasiStatic,
    solver::Solver,
    model::SolidMechanics,
)
    displacement_local = integrator.displacement
    solver.solution = displacement_local
    # BRP: apply inclined support inverse transform
    if model.inclined_support == true
        displacement = model.global_transform' * displacement_local
    else
        displacement = displacement_local
    end

    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(
    solver::Solver,
    model::SolidMechanics,
    integrator::QuasiStatic,
)
    displacement_local = solver.solution
    integrator.displacement = displacement_local
    if model.inclined_support == true
        displacement = model.global_transform' * displacement_local
    else
        displacement = displacement_local
    end
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(
    integrator::Newmark,
    _::HessianMinimizer,
    model::LinearOpInfRom,
)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    # Clean this up; maybe make a free dofs 2d array or move to a basis in matrix format
    for i ∈ 1:size(model.fom_model.current)[2]
        x_dof_index = 3 * (i - 1) + 1
        y_dof_index = 3 * (i - 1) + 2
        z_dof_index = 3 * (i - 1) + 3
        if model.fom_model.free_dofs[x_dof_index]
            model.fom_model.current[1, i] =
                model.basis[1, i, :]'displacement + model.fom_model.reference[1, i]
            model.fom_model.velocity[1, i] = model.basis[1, i, :]'velocity
            model.fom_model.acceleration[1, i] = model.basis[1, i, :]'acceleration
        end

        if model.fom_model.free_dofs[y_dof_index]
            model.fom_model.current[2, i] =
                model.basis[2, i, :]'displacement + model.fom_model.reference[2, i]
            model.fom_model.velocity[2, i] = model.basis[2, i, :]'velocity
            model.fom_model.acceleration[2, i] = model.basis[2, i, :]'acceleration
        end

        if model.fom_model.free_dofs[z_dof_index]
            model.fom_model.current[3, i] =
                model.basis[3, i, :]'displacement + model.fom_model.reference[3, i]
            model.fom_model.velocity[3, i] = model.basis[3, i, :]'velocity
            model.fom_model.acceleration[3, i] = model.basis[3, i, :]'acceleration
        end
    end

end

function copy_solution_source_targets(
    model::SolidMechanics,
    integrator::QuasiStatic,
    solver::Solver,
)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        integrator.displacement[3*node-2:3*node] = nodal_displacement
    end
    # Convert integrator displacement from global to local
    if model.inclined_support == true
        integrator.displacement = model.global_transform * integrator.displacement
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(
    integrator::Newmark,
    solver::HessianMinimizer,
    model::SolidMechanics,
)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = displacement

    if model.inclined_support == true
        displacement = model.global_transform' * displacement
        velocity = model.global_transform' * velocity
        acceleration = model.global_transform' * integrator.acceleration
    end

    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(
    solver::HessianMinimizer,
    model::SolidMechanics,
    integrator::Newmark,
)
    displacement = solver.solution
    integrator.displacement = displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration

    if model.inclined_support == true
        displacement = model.global_transform' * displacement
        velocity = model.global_transform' * velocity
        acceleration = model.global_transform' * acceleration
    end

    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(
    model::SolidMechanics,
    integrator::Newmark,
    solver::HessianMinimizer,
)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]

        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            local_transform = model.global_transform[base+1:base+3, base+1:base+3]
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        integrator.displacement[3*node-2:3*node] = nodal_displacement
        integrator.velocity[3*node-2:3*node] = nodal_velocity
        integrator.acceleration[3*node-2:3*node] = nodal_acceleration
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(
    integrator::CentralDifference,
    solver::ExplicitSolver,
    model::SolidMechanics,
)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Local (integrator) to global (model), use transpose
            local_transform = model.global_transform[base+1:base+3, base+1:base+3]'
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(
    solver::ExplicitSolver,
    model::SolidMechanics,
    integrator::CentralDifference,
)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = solver.solution
    integrator.acceleration = acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Local (integrator) to global (model), use transpose
            local_transform = model.global_transform[base+1:base+3, base+1:base+3]'
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(
    model::SolidMechanics,
    integrator::CentralDifference,
    solver::ExplicitSolver,
)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Global (model) to local (integrator)
            local_transform = model.global_transform[base+1:base+3, base+1:base+3]
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        integrator.displacement[3*node-2:3*node] = nodal_displacement
        integrator.velocity[3*node-2:3*node] = nodal_velocity
        integrator.acceleration[3*node-2:3*node] = nodal_acceleration
    end
    solver.solution = integrator.acceleration
end

function evaluate(integrator::Newmark, solver::HessianMinimizer, model::LinearOpInfRom)
    beta = integrator.β
    gamma = integrator.γ
    dt = integrator.time_step

    #Ax* = b
    # Put in residual format
    #e = [x* - x] -> x* = x + e
    #Ax + Ae = b
    #Ax - b = -Ae
    #Ae = r, r = b - Ax 
    ##M uddot + Ku = f

    num_dof = length(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)
    LHS = I / (dt * dt * beta) + Matrix{Float64}(model.opinf_rom["K"])
    RHS =
        model.opinf_rom["f"] +
        model.reduced_boundary_forcing +
        1.0 / (dt * dt * beta) .* integrator.disp_pre

    residual = RHS - LHS * solver.solution
    solver.hessian[:, :] = LHS
    solver.gradient[:] = -residual
end

function evaluate(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    stored_energy, internal_force, body_force, stiffness_matrix =
        evaluate(integrator, model)
    if model.failed == true
        return
    end
    integrator.stored_energy = stored_energy
    solver.value = stored_energy
    external_force = body_force + model.boundary_force
    if model.inclined_support == true
        solver.gradient = model.global_transform * (internal_force - external_force)
        solver.hessian = model.global_transform * stiffness_matrix * model.global_transform'
    else
        solver.gradient = internal_force - external_force
        solver.hessian = stiffness_matrix
    end
end

function evaluate(integrator::QuasiStatic, solver::SteepestDescent, model::SolidMechanics)
    stored_energy, internal_force, body_force, _ = evaluate(integrator, model)
    if model.failed == true
        return
    end
    integrator.stored_energy = stored_energy
    solver.value = stored_energy
    external_force = body_force + model.boundary_force
    solver.gradient = internal_force - external_force
end

function evaluate(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    stored_energy, internal_force, body_force, stiffness_matrix, mass_matrix =
        evaluate(integrator, model)
    if model.failed == true
        return
    end
    integrator.stored_energy = stored_energy
    β = integrator.β
    Δt = integrator.time_step
    inertial_force = mass_matrix * integrator.acceleration
    kinetic_energy = 0.5 * dot(integrator.velocity, mass_matrix, integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    external_force = body_force + model.boundary_force
    if model.inclined_support == true
        stiffness_matrix =
            model.global_transform * stiffness_matrix * model.global_transform'
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    solver.hessian = stiffness_matrix + mass_matrix / β / Δt / Δt
    solver.value = stored_energy - external_force ⋅ integrator.displacement + kinetic_energy
    solver.gradient = internal_force - external_force + inertial_force
end

function evaluate(
    integrator::CentralDifference,
    solver::ExplicitSolver,
    model::SolidMechanics,
)
    stored_energy, internal_force, body_force, lumped_mass = evaluate(integrator, model)
    if model.failed == true
        return
    end

    integrator.stored_energy = stored_energy
    # Intertial force (local)
    inertial_force = lumped_mass .* integrator.acceleration
    kinetic_energy = 0.5 * lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    # Body force -> global
    # Model boundary force -> global
    external_force = body_force + model.boundary_force
    if model.inclined_support == true
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    # External and internal force in local
    solver.value = stored_energy - external_force ⋅ integrator.displacement + kinetic_energy
    # Graident -> local, local, local
    solver.gradient = internal_force - external_force + inertial_force
    solver.lumped_hessian = lumped_mass
end

# Taken from ELASTOPLASITICITY—PART II: GLOBALLY CONVERGENT SCHEMES, Perez-Foguet & Armero, 2002
function backtrack_line_search(
    integrator::TimeIntegrator,
    solver::Solver,
    model::SolidMechanics,
    direction::Vector{Float64},
)
    backtrack_factor = solver.line_search.backtrack_factor
    decrease_factor = solver.line_search.decrease_factor
    max_iters = solver.line_search.max_iters
    free = model.free_dofs
    resid = solver.gradient
    merit = 0.5 * dot(resid, resid)
    merit_prime = -2.0 * merit
    step_length = solver.step.step_length
    step = step_length * direction
    initial_solution = 1.0 * solver.solution
    for _ ∈ 1:max_iters
        merit_old = merit
        step = step_length * direction
        solver.solution[free] = initial_solution[free] + step[free]
        copy_solution_source_targets(solver, model, integrator)
        evaluate(integrator, solver, model)
        if model.failed == true
            return step
        end
        resid = solver.gradient
        merit = 0.5 * dot(resid, resid)
        if merit ≤ (1.0 - 2.0 * decrease_factor * step_length) * merit_old
            break
        end
        step_length =
            max(
                backtrack_factor * step_length,
                -0.5 * step_length * step_length * merit_prime,
            ) / (merit - merit_old - step_length * merit_prime)
    end
    return step
end

function compute_step(
    _::QuasiStatic,
    model::SolidMechanics,
    solver::HessianMinimizer,
    _::NewtonStep,
)
    free = model.free_dofs
    return -solver.hessian[free, free] \ solver.gradient[free]
end

function compute_step(
    _::Newmark,
    model::SolidMechanics,
    solver::HessianMinimizer,
    _::NewtonStep,
)
    free = model.free_dofs
    return -solver.hessian[free, free] \ solver.gradient[free]
end

function compute_step(
    _::DynamicTimeIntegrator,
    model::LinearOpInfRom,
    solver::HessianMinimizer,
    _::NewtonStep,
)
    return -solver.hessian \ solver.gradient
end

function compute_step(
    _::CentralDifference,
    model::SolidMechanics,
    solver::ExplicitSolver,
    _::ExplicitStep,
)
    free = model.free_dofs
    return -solver.gradient[free] ./ solver.lumped_hessian[free]
end

function compute_step(
    integrator::QuasiStatic,
    model::SolidMechanics,
    solver::SteepestDescent,
    _::SteepestDescentStep,
)
    free = model.free_dofs
    step = backtrack_line_search(integrator, solver, model, -solver.gradient)
    return step[free]
end

function update_solver_convergence_criterion(
    solver::HessianMinimizer,
    absolute_error::Float64,
)
    solver.absolute_error = absolute_error
    solver.relative_error =
        solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : absolute_error
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    solver.converged = converged_absolute || converged_relative
end

function update_solver_convergence_criterion(
    solver::SteepestDescent,
    absolute_error::Float64,
)
    solver.absolute_error = absolute_error
    solver.relative_error =
        solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : absolute_error
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    solver.converged = converged_absolute || converged_relative
end

function update_solver_convergence_criterion(solver::ExplicitSolver, _::Float64)
    solver.converged = true
end

function stop_solve(solver::HessianMinimizer, iteration_number::Int64)
    if solver.failed == true
        return true
    end
    zero_residual = solver.absolute_error == 0.0
    if zero_residual == true
        return true
    end
    exceeds_minimum_iterations = iteration_number > solver.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations = iteration_number > solver.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return solver.converged
end

function stop_solve(solver::SteepestDescent, iteration_number::Int64)
    if solver.failed == true
        return true
    end
    zero_residual = solver.absolute_error == 0.0
    if zero_residual == true
        return true
    end
    exceeds_minimum_iterations = iteration_number > solver.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations = iteration_number > solver.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return solver.converged
end

function stop_solve(_::ExplicitSolver, _::Int64)
    return true
end

function solve(integrator::TimeIntegrator, solver::Solver, model::Model)
    predict(integrator, solver, model)
    evaluate(integrator, solver, model)
    if model.failed == true
        return
    end
    residual = solver.gradient
    norm_residual = norm(residual[model.free_dofs])
    solver.initial_norm = norm_residual
    iteration_number = 0
    solver.failed = solver.failed || model.failed
    step_type = solver.step
    while true
        step = compute_step(integrator, model, solver, step_type)
        solver.solution[model.free_dofs] += step
        correct(integrator, solver, model)
        evaluate(integrator, solver, model)
        if model.failed == true
            return
        end
        residual = solver.gradient
        norm_residual = norm(residual[model.free_dofs])
        if iteration_number == 0
            println("|R|=", norm_residual)
        else
            println("|R|=", norm_residual, ", solver iteration=", iteration_number)
        end
        update_solver_convergence_criterion(solver, norm_residual)
        iteration_number += 1
        if stop_solve(solver, iteration_number) == true
            break
        end
    end
    if model.inclined_support == true
        solver.gradient = model.global_transform' * solver.gradient
    end
end
