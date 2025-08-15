# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function RomHessianMinimizer(params::Parameters, model::RomModel)
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
    fom_solver = HessianMinimizer(params, model.fom_model)
    return RomHessianMinimizer(
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
        fom_solver,
    )
end

function RomExplicitSolver(params::Parameters, model::RomModel)
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
    fom_solver = ExplicitSolver(params, model.fom_model)
    return RomExplicitSolver(
        value, gradient, lumped_hessian, solution, initial_norm, converged, failed, step, fom_solver
    )
end


function create_solver(params::Parameters, model::RomModel)
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return RomHessianMinimizer(params, model)
    elseif solver_name == "explicit solver"
        return RomExplicitSolver(params, model)
    elseif solver_name == "steepest descent"
        norma_abort("Not supported for ROM : $solver_name")
    else
        norma_abort("Unknown type of solver : $solver_name")
    end
end


function copy_solution_source_targets(integrator::DynamicTimeIntegrator, solver::Solver, model::RomModel)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    # Clean this up; maybe make a free dofs 2d array or move to a basis in matrix format
    for i in 1:size(model.fom_model.current)[2]
        x_dof_index = 3 * (i - 1) + 1
        y_dof_index = 3 * (i - 1) + 2
        z_dof_index = 3 * (i - 1) + 3
        if model.fom_model.free_dofs[x_dof_index]
            model.fom_model.current[1, i] = model.basis[1, i, :]'displacement + model.fom_model.reference[1, i]
            model.fom_model.velocity[1, i] = model.basis[1, i, :]'velocity
            model.fom_model.acceleration[1, i] = model.basis[1, i, :]'acceleration
        end

        if model.fom_model.free_dofs[y_dof_index]
            model.fom_model.current[2, i] = model.basis[2, i, :]'displacement + model.fom_model.reference[2, i]
            model.fom_model.velocity[2, i] = model.basis[2, i, :]'velocity
            model.fom_model.acceleration[2, i] = model.basis[2, i, :]'acceleration
        end

        if model.fom_model.free_dofs[z_dof_index]
            model.fom_model.current[3, i] = model.basis[3, i, :]'displacement + model.fom_model.reference[3, i]
            model.fom_model.velocity[3, i] = model.basis[3, i, :]'velocity
            model.fom_model.acceleration[3, i] = model.basis[3, i, :]'acceleration
        end
    end
    return nothing
end


function evaluate(integrator::RomNewmark, solver::RomHessianMinimizer, model::QuadraticOpInfRom)
    beta = integrator.β
    gamma = integrator.γ
    dt = integrator.time_step

    ##Quadratic OpInf
    num_dof = length(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)
    # Create tangent stiffness for quadratic operator
    H = model.opinf_rom["H"]
    x1 = kron(I, solver.solution)
    x2 = kron(solver.solution, I)
    Hprime = H * x1 + H * x2
    xsqr = kron(solver.solution, solver.solution)
    LHS_linear = I / (dt * dt * beta) + Matrix{Float64}(model.opinf_rom["K"])
    LHS_nonlinear = Hprime

    RHS = model.opinf_rom["f"] + model.reduced_boundary_forcing + 1.0 / (dt * dt * beta) .* integrator.disp_pre

    residual = RHS - LHS_linear * solver.solution - H * xsqr
    solver.hessian[:, :] = LHS_linear + LHS_nonlinear
    solver.gradient[:] = -residual
    return nothing
end

function evaluate(integrator::RomNewmark, solver::RomHessianMinimizer, model::CubicOpInfRom)
    beta = integrator.β
    gamma = integrator.γ
    dt = integrator.time_step

    ##Cubic OpInf
    num_dof = length(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)
    # Create tangent stiffness for quadratic operator
    H = model.opinf_rom["H"]
    x1 = kron(I, solver.solution)
    x2 = kron(solver.solution, I)
    Hprime = H * x1 + H * x2
    xsqr = kron(solver.solution, solver.solution)
    # Create tangent stiffness for cubic operator
    G = model.opinf_rom["G"]
    x3 = kron(I, solver.solution, solver.solution)
    x4 = kron(solver.solution, I, solver.solution)
    x5 = kron(solver.solution, solver.solution, I)
    Gprime = G * x3 + G * x4 + G * x5
    xcub = kron(solver.solution, solver.solution, solver.solution)
    # Create LHSs
    LHS_linear = I / (dt * dt * beta) + Matrix{Float64}(model.opinf_rom["K"])
    LHS_nonlinear = Hprime + Gprime

    RHS = model.opinf_rom["f"] + model.reduced_boundary_forcing + 1.0 / (dt * dt * beta) .* integrator.disp_pre

    residual = RHS - LHS_linear * solver.solution - H * xsqr - G * xcub
    solver.hessian[:, :] = LHS_linear + LHS_nonlinear
    solver.gradient[:] = -residual
    return nothing
end

function evaluate(integrator::RomCentralDifference, solver::RomExplicitSolver, model::LinearOpInfRom)
    gamma = integrator.γ
    dt = integrator.time_step
    ## Value for accelertaion - assumes bases are orthonormal to avoid solve
    solver.solution[:] =
        model.opinf_rom["f"] + model.reduced_boundary_forcing - model.opinf_rom["K"] * integrator.displacement[:]
    return nothing
end

function evaluate(integrator::RomNewmark, solver::RomHessianMinimizer, model::LinearOpInfRom)
    beta = integrator.β
    gamma = integrator.γ
    dt = integrator.time_step

    num_dof = length(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)
    LHS = I / (dt * dt * beta) + Matrix{Float64}(model.opinf_rom["K"])
    RHS = model.opinf_rom["f"] + model.reduced_boundary_forcing + 1.0 / (dt * dt * beta) .* integrator.disp_pre

    residual = RHS - LHS * solver.solution
    solver.hessian[:, :] = LHS
    solver.gradient[:] = -residual
    return nothing
end

function compute_step(_::DynamicTimeIntegrator, model::RomModel, solver::RomHessianMinimizer, _::NewtonStep)
    atol = solver.linear_solver_absolute_tolerance
    rtol = solver.linear_solver_relative_tolerance
    return -solve_linear(solver.hessian, solver.gradient, atol, rtol)
end

function compute_step(_::RomCentralDifference, model::RomModel, solver::RomExplicitSolver, _::ExplicitStep)
    free = model.free_dofs
    return zeros(size(model.free_dofs))
end

function update_solver_convergence_criterion(solver::RomHessianMinimizer, absolute_error::Float64)
    solver.absolute_error = absolute_error
    solver.relative_error = solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : absolute_error
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    return solver.converged = converged_absolute || converged_relative
end

function update_solver_convergence_criterion(solver::RomExplicitSolver, _::Float64)
    return solver.converged = true
end

function stop_solve(solver::RomHessianMinimizer, iteration_number::Int64)
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

function stop_solve(_::RomExplicitSolver, _::Int64)
    return true
end

