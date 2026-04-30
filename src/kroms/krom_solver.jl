function evaluate(integrator::RomNewmark, solver::RomHessianMinimizer, model::RBFKernelROM)
    beta = integrator.β
    gamma = integrator.γ
    dt = integrator.time_step

    ##RBF Kernel ROM
    num_dof = length(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)

    # Create tangent stiffness for quadratic operator
    val, LHS_nonlinear = evaluate_rbf_interpolant(model.interpolant, solver.solution)
    LHS_linear = I / (dt * dt * beta)

    RHS =
        model.forcing +
        model.reduced_boundary_forcing +
        1.0 / (dt * dt * beta) .* integrator.disp_pre

    residual = RHS - LHS_linear * solver.solution - val
    solver.hessian[:, :] = LHS_linear + LHS_nonlinear
    solver.gradient[:] = -residual
    return nothing
end
