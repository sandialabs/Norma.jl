# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Sparse Cholesky of the first reduced tangent of a case under three orderings:
#   julia -t N --project=.. ordering.jl case.yaml
# For the default ordering of CHOLMOD, the nested dissection of METIS that
# Norma uses, and the natural ordering, one line with the nonzeros of the
# factor, the seconds of the full factorization (symbolic analysis and
# numeric phase), of the numeric phase alone with the analysis reused, and of
# one solve, each the minimum of three calls, and the relative residual of
# the solution. The BLAS thread count follows the Julia thread count.
using Norma, LinearAlgebra, SparseArrays
input = ARGS[1]
sim = Norma.create_simulation(input)
Norma.initialize(sim)
model, integrator, solver = sim.model, sim.integrator, sim.solver
Norma.predict(integrator, solver, model)
Norma.evaluate(integrator, solver, model)
A = Norma.reduced_matrix!(solver, solver.hessian, model.free_dofs)
b = copy(Norma.reduced_vector!(solver.reduced, solver.gradient))
n = size(A, 1)
function measure(name, perm)
    F = cholesky(Symmetric(A); perm=perm)
    t_factor = minimum(@elapsed(cholesky(Symmetric(A); perm=perm)) for _ in 1:3)
    t_numeric = minimum(@elapsed(cholesky!(F, Symmetric(A))) for _ in 1:3)
    t_solve = minimum(@elapsed(F \ b) for _ in 1:3)
    x = F \ b
    return println(
        "ORDERING threads=",
        Threads.nthreads(),
        " case=",
        input,
        " ordering=",
        name,
        " n=",
        n,
        " nnzA=",
        nnz(A),
        " nnzL=",
        nnz(sparse(F.L)),
        " factor_s=",
        round(t_factor; digits=3),
        " numeric_s=",
        round(t_numeric; digits=3),
        " solve_s=",
        round(t_solve; digits=4),
        " residual=",
        norm(A * x - b) / norm(b),
    )
end
measure("default", nothing)
t_ordering = minimum(@elapsed(Norma.fill_reducing_ordering(A)) for _ in 1:3)
println("ORDERING threads=", Threads.nthreads(), " case=", input, " metis_s=", round(t_ordering; digits=3))
measure("metis", Norma.fill_reducing_ordering(A))
measure("natural", collect(1:n))
