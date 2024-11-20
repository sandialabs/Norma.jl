function create_step(solver_params::Dict{Any,Any})
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

function HessianMinimizer(params::Dict{Any,Any})
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
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
    )
end

function ExplicitSolver(params::Dict{Any,Any})
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
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

function SteepestDescent(params::Dict{Any,Any})
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
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
    )
end

function NewtonStep(params::Dict{Any,Any})
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    NewtonStep(step_length)
end

function ExplicitStep(params::Dict{Any,Any})
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    ExplicitStep(step_length)
end

function SteepestDescentStep(params::Dict{Any,Any})
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    if haskey(params, "max step length") == true
        max_step_length = params["max step length"]
    else
        max_step_length = step_length
    end
    if haskey(params, "min step length") == true
        min_step_length = params["min step length"]
    else
        min_step_length = step_length
    end
    if haskey(params, "step reduction factor") == true
        reduction_factor = params["reduction factor"]
    else
        reduction_factor = 0.5
    end
    if haskey(params, "step increase factor") == true
        increase_factor = params["increase factor"]
    else
        increase_factor = 1.1
    end
    SteepestDescentStep(
        AdaptiveStepLength(step_length),
        max_step_length,
        min_step_length,
        reduction_factor,
        increase_factor
   )
end

function create_solver(params::Dict{Any,Any})
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params)
    elseif solver_name == "explicit solver"
        return ExplicitSolver(params)
    elseif solver_name == "steepest descent"
        return SteepestDescent(params)
    else
        error("Unknown type of solver : ", solver_name)
    end
end

function copy_solution_source_targets(
    integrator::QuasiStatic,
    solver::Any,
    model::SolidMechanics,
)
    displacement = integrator.displacement
    solver.solution = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(
    solver::Any,
    model::SolidMechanics,
    integrator::QuasiStatic,
)
    displacement = solver.solution
    integrator.displacement = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(
    model::SolidMechanics,
    integrator::QuasiStatic,
    solver::Any,
)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        integrator.displacement[3*node-2:3*node] = nodal_displacement
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
        integrator.displacement[3*node-2:3*node] = nodal_displacement
        integrator.velocity[3*node-2:3*node] = nodal_velocity
        integrator.acceleration[3*node-2:3*node] = nodal_acceleration
    end
    solver.solution = integrator.acceleration
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
    solver.gradient = internal_force - external_force
    solver.hessian = stiffness_matrix
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
    inertial_force = lumped_mass .* integrator.acceleration
    kinetic_energy = 0.5 * lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    external_force = body_force + model.boundary_force
    solver.value = stored_energy - external_force ⋅ integrator.displacement + kinetic_energy
    solver.gradient = internal_force - external_force + inertial_force
    solver.lumped_hessian = lumped_mass
end

# Taken from ELASTOPLASITICITY—PART II: GLOBALLY CONVERGENT SCHEMES, Perez-Foguet & Armero, 2002
function backtrack_line_search(
    integrator::TimeIntegrator,
    solver::Any,
    model::SolidMechanics,
    direction::Vector{Float64},
)
    backtrack_factor = 0.1
    decrease_factor = 1.0e-04
    free = model.free_dofs
    resid = solver.gradient
    merit = 0.5 * dot(resid, resid)
    merit_prime = -2.0 * merit
    step_length = solver.step.step_length.value
    step = step_length * direction
    initial_solution = 1.0 * solver.solution
    max_ls_iters = 20
    for _ ∈ 1:max_ls_iters
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
    prev_solution = 1.0 * solver.solution
    while true
        step = compute_step(integrator, model, solver, step_type)
        if model.failed == true
            if solver.step.step_length isa AdaptiveStepLength
                if solver.step.reduction_factor < 1.0
                && solver.step.reduction_factor > 0.0
                && solver.step.step_length.value * solver.step.reduction_factor > solver.step.min_step_length
                    @info "Reducing step length from $(solver.step.step_length.value) to $(solver.step.step_length.value * solver.step.reduction_factor)"
                    solver.step.step_length.value *= solver.step.reduction_factor
                    model.failed = false
                    solver.solution = 1.0 * prev_solution
                    continue
                else
                    error("evaluation of model has failed and step length cannot be reduced")
                end
            else
                error("evaluation of model has failed")
            end
        end
        prev_solution = 1.0 * solver.solution
        solver.solution[model.free_dofs] += step
        correct(integrator, solver, model)
        evaluate(integrator, solver, model)
        if (solver.step.step_length isa AdaptiveStepLength
            && solver.step.increase_factor > 1.0
            && solver.step.step_length.value * solver.step.reduction_factor < solver.step.max_step_length)
            solver.step.step_length.value *= solver.step.increase_factor
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
end
