# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Time of the components of one step on the initialized model of a case:
#   julia -t N --project=.. components.jl case.yaml
# Each component is called once to compile and then timed as the minimum of
# five calls, with the bytes it allocates. The explicit case reports the
# evaluation, the reductions, the stable step, and the kinematic updates; the
# implicit cases add the matrix bookkeeping and the linear solve with and
# without the Jacobi preconditioner.
using Norma, LinearAlgebra, SparseArrays
input = ARGS[1]
sim = Norma.create_simulation(input)
Norma.initialize(sim)
model, integrator, solver = sim.model, sim.integrator, sim.solver
ndofs = 3 * size(model.reference, 2)
free = model.free_dofs
function bench(name, f; n=5)
    f()
    t = minimum(@elapsed(f()) for _ in 1:n)
    a = @allocated f()
    println(
        "COMP threads=",
        Threads.nthreads(),
        " ",
        rpad(name, 28),
        " ms=",
        round(1000t; digits=3),
        " MB=",
        round(a / 1e6; digits=2),
    )
    return nothing
end
bench("apply_bcs", () -> Norma.apply_bcs(model))
bench("evaluate(model)", () -> Norma.evaluate(model, integrator, solver))
bench("evaluate(solver)", () -> Norma.evaluate(integrator, solver, model))
bench("reduce_element_vector", () -> Norma.reduce_element_vector(model.value_buffers, model.value_buffers.force, ndofs))
bench("norm(gradient[free])", () -> norm(solver.gradient[free]))
bench("predict", () -> Norma.predict(integrator, solver, model))
bench("correct", () -> Norma.correct(integrator, solver, model))
if integrator isa Norma.CentralDifference
    bench("stable_time_step", () -> Norma.stable_time_step(integrator, model))
    bench("compute_step", () -> Norma.compute_step(integrator, model, solver, solver.step))
else
    pattern = model.matrix_pattern
    bench("reduce_element_matrix", () -> Norma.reduce_element_matrix(pattern, pattern.stiffness, ndofs))
    bench(
        "newmark_hessian", () -> Norma.newmark_hessian(model.stiffness, model.mass, integrator.β, integrator.time_step)
    )
    bench("mass*acceleration", () -> model.mass * integrator.acceleration)
    bench("hessian[free,free]", () -> solver.hessian[free, free])
    bench("reduced_matrix!", () -> Norma.reduced_matrix!(solver, solver.hessian, free))
    A = solver.hessian[free, free]
    b = solver.gradient[free]
    bench("symmetric_to_tolerance", () -> Norma.symmetric_to_tolerance(A))
    bench("solve_linear", () -> Norma.solve_linear(solver, A, b))
    if solver.linear_solver == "cg"
        println("COMP cg iterations none=", solver.cg_workspace.iterations)
        solver.preconditioner = :jacobi
        bench("solve_linear jacobi", () -> Norma.solve_linear(solver, A, b))
        println("COMP cg iterations jacobi=", solver.cg_workspace.iterations)
        solver.preconditioner = :none
    end
    if solver.linear_solver == "cg"
        x, hist = Norma.IterativeSolvers.cg(
            Norma.SymmetricRowOperator(A),
            b;
            abstol=solver.linear_solver_absolute_tolerance,
            reltol=solver.linear_solver_relative_tolerance,
            log=true,
        )
        println("COMP cg iterations=", hist.iters, " n=", size(A, 1), " nnz=", nnz(A))
        bench("one mul! (row operator)", () -> mul!(x, Norma.SymmetricRowOperator(A), b))
        bench("one mul! (CSC)", () -> mul!(x, A, b))
    else
        println("COMP direct n=", size(A, 1), " nnz=", nnz(A), " factor_nnz=", nnz(solver.direct_solver_cache.factor))
        bench("cholesky (fresh)", () -> cholesky(Symmetric(A); check=false))
        F = cholesky(Symmetric(A); check=false)
        bench("cholesky! (reuse symbolic)", () -> cholesky!(F, Symmetric(A); check=false))
        bench("factor \\ b", () -> F \ b)
    end
end
