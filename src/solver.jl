# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

include("opinf/opinf_solver.jl")
using IterativeSolvers
using LinearAlgebra
using Printf


function HessianMinimizer(params::Parameters, model::Model)
    solver_params = params["solver"]
    num_dof = length(model.free_dofs)
    minimum_iterations = solver_params["minimum iterations"]
    maximum_iterations = solver_params["maximum iterations"]
    absolute_tolerance = solver_params["absolute tolerance"]
    relative_tolerance = solver_params["relative tolerance"]
    # See documentation for IterativeSolvers.cg and others for this
    default_abs_tol = 0.0
    default_rel_tol = sqrt(eps(Float64))
    linear_solver_absolute_tolerance = get(solver_params, "linear solver absolute tolerance", default_abs_tol)
    linear_solver_relative_tolerance = get(solver_params, "linear solver relative tolerance", default_rel_tol)
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
    use_line_search = get(solver_params, "use line search", false)
    ls_backtrack_factor = get(solver_params, "line search backtrack factor", 0.5)
    ls_decrease_factor = get(solver_params, "line search decrease factor", 1.0e-04)
    ls_max_iters = get(solver_params, "line search maximum iterations", 16)
    line_search = BackTrackLineSearch(ls_backtrack_factor, ls_decrease_factor, ls_max_iters)
    return HessianMinimizer(
        minimum_iterations,
        maximum_iterations,
        absolute_tolerance,
        relative_tolerance,
        absolute_error,
        relative_error,
        linear_solver_absolute_tolerance,
        linear_solver_relative_tolerance,
        value,
        gradient,
        hessian,
        solution,
        initial_norm,
        converged,
        failed,
        step,
        line_search,
        use_line_search,
    )
end

function ExplicitSolver(params::Parameters, model::Model)
    solver_params = params["solver"]
    num_dof = length(model.free_dofs)
    value = 0.0
    gradient = zeros(num_dof)
    lumped_hessian = zeros(num_dof)
    solution = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    return ExplicitSolver(value, gradient, lumped_hessian, solution, initial_norm, converged, failed, step)
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
    ls_backtrack_factor = get(solver_params, "line search backtrack factor", 0.5)
    ls_decrease_factor = get(solver_params, "line search decrease factor", 1.0e-04)
    ls_max_iters = get(solver_params, "line search maximum iterations", 16)
    line_search = BackTrackLineSearch(ls_backtrack_factor, ls_decrease_factor, ls_max_iters)
    line_search = BackTrackLineSearch(ls_backtrack_factor, ls_decrease_factor, ls_max_iters)
    return SteepestDescent(
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

function create_solver(params::Parameters, model::SolidMechanics)
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params, model)
    elseif solver_name == "explicit solver"
        return ExplicitSolver(params, model)
    elseif solver_name == "steepest descent"
        return SteepestDescent(params, model)
    else
        norma_abort("Unknown type of solver : $solver_name")
    end
end


function NewtonStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    return NewtonStep(step_length)
end

function ExplicitStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    return ExplicitStep(step_length)
end

function SteepestDescentStep(params::Parameters)
    if haskey(params, "step length") == true
        step_length = params["step length"]
    else
        step_length = 1.0
    end
    return SteepestDescentStep(step_length)
end

function create_step(solver_params::Parameters)
    step_name = solver_params["step"]
    if step_name == "full Newton"
        return NewtonStep(solver_params)
    elseif step_name == "explicit"
        return ExplicitStep(solver_params)
    elseif step_name == "steepest descent"
        return SteepestDescentStep(solver_params)
    else
        norma_abort("Unknown type of solver step: $step_name")
    end
end

function copy_solution_source_targets(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    displacement_local = integrator.displacement
    solver.solution = displacement_local
    if model.inclined_support == true
        displacement = model.global_transform' * displacement_local
    else
        displacement = displacement_local
    end
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
    return nothing
end

function copy_solution_source_targets(solver::Solver, model::SolidMechanics, integrator::QuasiStatic)
    displacement_local = solver.solution
    integrator.displacement = displacement_local
    if model.inclined_support == true
        displacement = model.global_transform' * displacement_local
    else
        displacement = displacement_local
    end
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
    return nothing
end


function copy_solution_source_targets(model::SolidMechanics, integrator::QuasiStatic, solver::Solver)
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        integrator.displacement[(3 * node - 2):(3 * node)] = nodal_displacement
    end
    if model.inclined_support == true
        integrator.displacement = model.global_transform * integrator.displacement
    end
    solver.solution = integrator.displacement
    return nothing
end

function copy_solution_source_targets(integrator::Newmark, solver::Solver, model::SolidMechanics)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = displacement

    if model.inclined_support == true
        displacement = model.global_transform' * displacement
        velocity = model.global_transform' * velocity
        acceleration = model.global_transform' * integrator.acceleration
    end

    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        nodal_velocity = velocity[(3 * node - 2):(3 * node)]
        nodal_acceleration = acceleration[(3 * node - 2):(3 * node)]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
    return nothing
end

function copy_solution_source_targets(solver::Solver, model::SolidMechanics, integrator::Newmark)
    displacement = solver.solution
    integrator.displacement = displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration

    if model.inclined_support == true
        displacement = model.global_transform' * displacement
        velocity = model.global_transform' * velocity
        acceleration = model.global_transform' * acceleration
    end

    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        nodal_velocity = velocity[(3 * node - 2):(3 * node)]
        nodal_acceleration = acceleration[(3 * node - 2):(3 * node)]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
    return nothing
end

function copy_solution_source_targets(model::SolidMechanics, integrator::Newmark, solver::Solver)
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]

        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            local_transform = model.global_transform[(base + 1):(base + 3), (base + 1):(base + 3)]
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        integrator.displacement[(3 * node - 2):(3 * node)] = nodal_displacement
        integrator.velocity[(3 * node - 2):(3 * node)] = nodal_velocity
        integrator.acceleration[(3 * node - 2):(3 * node)] = nodal_acceleration
    end
    solver.solution = integrator.displacement
    return nothing
end

function copy_solution_source_targets(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = acceleration
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        nodal_velocity = velocity[(3 * node - 2):(3 * node)]
        nodal_acceleration = acceleration[(3 * node - 2):(3 * node)]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Local (integrator) to global (model), use transpose
            local_transform = model.global_transform[(base + 1):(base + 3), (base + 1):(base + 3)]'
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
    return nothing
end

function copy_solution_source_targets(solver::ExplicitSolver, model::SolidMechanics, integrator::CentralDifference)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = solver.solution
    integrator.acceleration = acceleration
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = displacement[(3 * node - 2):(3 * node)]
        nodal_velocity = velocity[(3 * node - 2):(3 * node)]
        nodal_acceleration = acceleration[(3 * node - 2):(3 * node)]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Local (integrator) to global (model), use transpose
            local_transform = model.global_transform[(base + 1):(base + 3), (base + 1):(base + 3)]'
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
    return nothing
end

function copy_solution_source_targets(model::SolidMechanics, integrator::CentralDifference, solver::ExplicitSolver)
    num_nodes = size(model.reference, 2)
    for node in 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]
        if model.inclined_support == true
            base = 3 * (node - 1) # Block index in global stiffness
            # Global (model) to local (integrator)
            local_transform = model.global_transform[(base + 1):(base + 3), (base + 1):(base + 3)]
            nodal_displacement = local_transform * nodal_displacement
            nodal_velocity = local_transform * nodal_velocity
            nodal_acceleration = local_transform * nodal_acceleration
        end
        integrator.displacement[(3 * node - 2):(3 * node)] = nodal_displacement
        integrator.velocity[(3 * node - 2):(3 * node)] = nodal_velocity
        integrator.acceleration[(3 * node - 2):(3 * node)] = nodal_acceleration
    end
    solver.solution = integrator.acceleration
    return nothing
end


function evaluate(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    evaluate(model, integrator, solver)
    if model.failed == true
        return nothing
    end
    integrator.stored_energy = model.strain_energy
    solver.value = model.strain_energy
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        solver.gradient = model.global_transform * (model.internal_force - external_force)
        solver.hessian = model.global_transform * model.stiffness * model.global_transform'
    else
        solver.gradient = model.internal_force - external_force
        solver.hessian = model.stiffness
    end
    return nothing
end

function evaluate(integrator::QuasiStatic, solver::MatrixFree, model::SolidMechanics)
    evaluate(model, integrator, solver)
    if model.failed == true
        return nothing
    end
    integrator.stored_energy = model.strain_energy
    solver.value = model.strain_energy
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        solver.gradient = model.global_transform * (model.internal_force - external_force)
    else
        solver.gradient = model.internal_force - external_force
    end
    return nothing
end

function evaluate(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    evaluate(model, integrator, solver)
    if model.failed == true
        return nothing
    end
    integrator.stored_energy = model.strain_energy
    β = integrator.β
    Δt = integrator.time_step
    inertial_force = model.mass * integrator.acceleration
    kinetic_energy = 0.5 * dot(integrator.velocity, model.mass, integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    internal_force = model.internal_force
    external_force = model.body_force + model.boundary_force
    stiffness = model.stiffness
    if model.inclined_support == true
        global_transform = model.global_transform
        internal_force = global_transform * internal_force
        external_force = global_transform * external_force
        stiffness = global_transform * stiffness * global_transform'
    end
    solver.hessian = stiffness + model.mass / β / Δt / Δt
    solver.gradient = internal_force - external_force + inertial_force
    solver.value = model.strain_energy - external_force ⋅ integrator.displacement + kinetic_energy
    return nothing
end

function evaluate(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    evaluate(model, integrator, solver)
    if model.failed == true
        return nothing
    end
    integrator.stored_energy = model.strain_energy
    # Intertial force (local)
    inertial_force = model.lumped_mass .* integrator.acceleration
    kinetic_energy = 0.5 * model.lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    # Body force -> global
    # Model boundary force -> global
    internal_force = model.internal_force
    external_force = model.body_force + model.boundary_force
    if model.inclined_support == true
        external_force = model.global_transform * external_force
        internal_force = model.global_transform * internal_force
    end
    # External and internal force in local
    solver.value = model.strain_energy - external_force ⋅ integrator.displacement + kinetic_energy
    # Gradient -> local, local, local
    solver.gradient = internal_force - external_force + inertial_force
    solver.lumped_hessian = model.lumped_mass
    return nothing
end

function backtrack_line_search(
    integrator::TimeIntegrator, solver::Solver, model::SolidMechanics, direction::Vector{Float64}
)
    backtrack_factor = solver.line_search.backtrack_factor
    decrease_factor = solver.line_search.decrease_factor
    max_iters = solver.line_search.max_iters
    free = model.free_dofs
    resid = solver.gradient[free]
    merit_init = 0.5 * dot(resid, resid)
    step_length = solver.step.step_length
    step = step_length * direction
    initial_solution = 1.0 * solver.solution
    compute_stiffness = model.compute_stiffness
    compute_mass = model.compute_mass
    compute_lumped_mass = model.compute_lumped_mass
    model.compute_stiffness = false
    model.compute_mass = false
    model.compute_lumped_mass = false
    for iter in 1:max_iters
        norma_logf(8, :linesearch, "Line Search [%d] |ΔX| = %.3e", iter, step_length)
        step = step_length * direction
        solver.solution[free] = initial_solution[free] + step
        copy_solution_source_targets(solver, model, integrator)
        evaluate(integrator, solver, model)
        if model.failed == true
            step_length *= backtrack_factor
            continue
        end
        resid = solver.gradient[free]
        merit = 0.5 * dot(resid, resid)
        if merit ≤ merit_init + decrease_factor * step_length * dot(resid, direction)
            break
        end
        step_length *= backtrack_factor
    end
    solver.solution = initial_solution
    copy_solution_source_targets(solver, model, integrator)
    model.compute_stiffness = compute_stiffness
    model.compute_mass = compute_mass
    model.compute_lumped_mass = compute_lumped_mass
    return step
end

function solve_linear(A::SparseMatrixCSC{Float64}, b::Vector{Float64}, atol::Float64, rtol::Float64)
    return cg(A, b; abstol=atol, reltol=rtol)
end

function compute_step(integrator::TimeIntegrator, model::SolidMechanics, solver::HessianMinimizer, _::NewtonStep)
    free = model.free_dofs
    A = solver.hessian[free, free]
    b = solver.gradient[free]
    atol = solver.linear_solver_absolute_tolerance
    rtol = solver.linear_solver_relative_tolerance
    step = -solve_linear(A, b, atol, rtol)
    return solver.use_line_search ? backtrack_line_search(integrator, solver, model, step) : step
end

function compute_step(integrator::TimeIntegrator, model::SolidMechanics, solver::MatrixFree, _::SteepestDescentStep)
    free = model.free_dofs
    direction = LinearAlgebra.normalize(-solver.gradient[free])
    return backtrack_line_search(integrator, solver, model, direction)
end

function compute_step(_::CentralDifference, model::SolidMechanics, solver::ExplicitSolver, _::ExplicitStep)
    free = model.free_dofs
    return -solver.gradient[free] ./ solver.lumped_hessian[free]
end


function update_solver_convergence_criterion(solver::HessianMinimizer, absolute_error::Float64)
    solver.absolute_error = absolute_error
    solver.relative_error = solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : absolute_error
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    return solver.converged = converged_absolute || converged_relative
end


function update_solver_convergence_criterion(solver::SteepestDescent, absolute_error::Float64)
    solver.absolute_error = absolute_error
    solver.relative_error = solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : absolute_error
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    return solver.converged = converged_absolute || converged_relative
end

function update_solver_convergence_criterion(solver::ExplicitSolver, _::Float64)
    return solver.converged = true
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
    is_explicit_dynamic = integrator isa ExplicitDynamicTimeIntegrator
    predict(integrator, solver, model)
    evaluate(integrator, solver, model)
    if model.failed == true
        return nothing
    end
    residual = solver.gradient
    norm_residual = norm(residual[model.free_dofs])
    if is_explicit_dynamic == false
        raw_status = "[WAIT]"
        status = colored_status(raw_status)
        norma_logf(8, :solve, "Iteration [%d] %s = %.3e : %s = %.3e : %s", 0, "|R|", norm_residual, "|r|", 1.0, status)
    end
    solver.initial_norm = norm_residual
    iteration_number = 1
    solver.failed = solver.failed || model.failed
    step_type = solver.step
    while true
        step = compute_step(integrator, model, solver, step_type)
        solver.solution[model.free_dofs] += step
        correct(integrator, solver, model)
        evaluate(integrator, solver, model)
        if model.failed == true
            return nothing
        end
        residual = solver.gradient
        norm_residual = norm(residual[model.free_dofs])
        update_solver_convergence_criterion(solver, norm_residual)
        if is_explicit_dynamic == false
            raw_status = solver.converged ? "[DONE]" : "[WAIT]"
            status = colored_status(raw_status)
            norma_logf(
                8,
                :solve,
                "Iteration [%d] %s = %.3e : %s = %.3e : %s",
                iteration_number,
                "|R|",
                solver.absolute_error,
                "|r|",
                solver.relative_error,
                status,
            )
        end
        iteration_number += 1
        if stop_solve(solver, iteration_number) == true
            break
        end
    end
    if model.inclined_support == true
        solver.gradient = model.global_transform' * solver.gradient
    end
    return nothing
end
