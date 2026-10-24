# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The sparse product of the first reduced tangent of a case with 64-bit and
# 32-bit column indices, and with one stored triangle on one thread:
#   julia -t N --project=.. product.jl case.yaml
# Each variant is the minimum of 50 products. The product with one triangle
# scatters into the rows of the transposed part, so it runs on one thread.
using Norma, LinearAlgebra, SparseArrays
input = ARGS[1]
sim = Norma.create_simulation(input)
Norma.initialize(sim)
model, integrator, solver = sim.model, sim.integrator, sim.solver
Norma.predict(integrator, solver, model)
Norma.evaluate(integrator, solver, model)
A = Norma.reduced_matrix!(solver, solver.hessian, model.free_dofs)
n = size(A, 1)
x = rand(n)
y = zeros(n)
@inline function rows!(y, rowptr::Vector{Ti}, colval::Vector{Ti}, nzval, x, rows) where {Ti}
    @inbounds for row in rows
        acc = 0.0
        for k in rowptr[row]:(rowptr[row + 1] - 1)
            acc += nzval[k] * x[colval[k]]
        end
        y[row] = acc
    end
    return nothing
end
function product!(y, rowptr, colval, nzval, x)
    Norma.chunk_results(Nothing, () -> nothing, length(y); chunk=Norma.vector_chunk(length(y))) do rows, _
        rows!(y, rowptr, colval, nzval, x, rows)
        return nothing
    end
    return y
end
function triangle!(y, colptr, rowval, nzval, x)
    fill!(y, 0.0)
    @inbounds for j in eachindex(y)
        xj = x[j]
        acc = 0.0
        for k in colptr[j]:(colptr[j + 1] - 1)
            i = rowval[k]
            v = nzval[k]
            acc += v * x[i]
            i != j && (y[i] += v * xj)
        end
        y[j] += acc
    end
    return y
end
function bench(name, f)
    f()
    t = minimum(@elapsed(f()) for _ in 1:50)
    return println(
        "PRODUCT threads=",
        Threads.nthreads(),
        " case=",
        input,
        " n=",
        n,
        " nnz=",
        nnz(A),
        " variant=",
        name,
        " ms=",
        round(1000t; digits=3),
    )
end
rowptr64, colval64, nzval = A.colptr, A.rowval, A.nzval
rowptr32, colval32 = Vector{Int32}(rowptr64), Vector{Int32}(colval64)
U = sparse(UpperTriangular(A))
bench("int64", () -> product!(y, rowptr64, colval64, nzval, x))
bench("int32", () -> product!(y, rowptr32, colval32, nzval, x))
bench("int64-one-thread", () -> rows!(y, rowptr64, colval64, nzval, x, 1:n))
bench("triangle-one-thread", () -> triangle!(y, U.colptr, U.rowval, U.nzval, x))
