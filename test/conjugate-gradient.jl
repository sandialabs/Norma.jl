# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The conjugate gradient of the Hessian minimizer on its reused workspace,
# without and with the Jacobi preconditioner: it meets the stopping test of
# the library solver it replaced, agrees with the direct solution, repeats
# itself exactly, and gives the same displacements as the unpreconditioned
# solve on a Newmark run.
using LinearAlgebra
using SparseArrays
using Test

@testset "Conjugate Gradient" begin
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    sim = Norma.create_simulation("clamped.yaml")
    Norma.initialize(sim)
    model, integrator, solver = sim.model, sim.integrator, sim.solver
    Norma.predict(integrator, solver, model)
    Norma.evaluate(integrator, solver, model)
    free = model.free_dofs
    A = Norma.reduced_matrix!(solver, solver.hessian, free)
    b = copy(Norma.reduced_vector!(solver.reduced, solver.gradient))
    x_direct = cholesky(Symmetric(A)) \ b
    reltol = solver.linear_solver_relative_tolerance
    abstol = solver.linear_solver_absolute_tolerance
    iterations = Dict{Symbol,Int64}()
    for preconditioner in (:none, :jacobi, :cholesky)
        solver.preconditioner = preconditioner
        x = copy(Norma.conjugate_gradient!(solver, A, b))
        ws = solver.cg_workspace
        iterations[preconditioner] = ws.iterations
        @test 0 < ws.iterations < size(A, 1)
        @test norm(A * x - b) <= 2.0 * max(reltol * norm(b), abstol)
        @test norm(x - x_direct) / norm(x_direct) < 1.0e-6
        # The same solve again gives the same iterate, bit for bit.
        @test Norma.conjugate_gradient!(solver, A, b) == x
        @test (preconditioner == :jacobi) == !isempty(ws.inverse_diagonal)
    end
    @test iterations[:jacobi] <= iterations[:none]
    # A fresh factor is the exact inverse: the preconditioned iteration
    # converges in a few iterations.
    @test iterations[:cholesky] <= 3
    # The factor is kept across solves of the same matrix and pattern.
    factor = solver.direct_cache.factor
    Norma.conjugate_gradient!(solver, A, b)
    @test solver.direct_cache.factor === factor
    Norma.finalize_writing(sim)

    # The run with the Jacobi preconditioner against the run without.
    displacements = Dict{String,Vector{Float64}}()
    for preconditioner in ("none", "jacobi", "cholesky")
        input = read("clamped.yaml", String)
        write(
            "clamped-$preconditioner.yaml",
            replace(input, "solver:\n" => "solver:\n  preconditioner: $preconditioner\n"),
        )
        run = Norma.run("clamped-$preconditioner.yaml")
        @test run.solver.preconditioner == Symbol(preconditioner)
        displacements[preconditioner] = copy(run.integrator.displacement)
        rm("clamped-$preconditioner.yaml"; force=true)
    end
    for preconditioner in ("jacobi", "cholesky")
        difference = norm(displacements[preconditioner] - displacements["none"]) / norm(displacements["none"])
        @test difference < 1.0e-7
    end
    rm("clamped.yaml"; force=true)
    rm("clamped.g"; force=true)
    rm("clamped.e"; force=true)
end
