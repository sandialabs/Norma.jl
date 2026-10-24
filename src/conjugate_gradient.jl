# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Preconditioned conjugate gradient on the reduced system of the Hessian
# minimizer. The library solver (IterativeSolvers.cg) allocated its state at
# every solve and ran its vector operations on one thread: of a 45-iteration
# solve on 45 000 unknowns on 12 threads, 7.2 ms were the threaded products
# and 6.7 ms the rest. Here the vectors live in a workspace reused across
# solves, and each iteration makes three passes, each in fixed chunks on all
# threads (chunk_results): the product with the matrix fused with the dot
# product u ⋅ A u, the update of the iterate and the residual fused with the
# dot products r ⋅ r and r ⋅ z, and the new direction. The chunks and the
# order of the partial sums are fixed, so the iterates are the same on every
# thread count. The stopping test is that of the library solver: the residual
# norm at or below max(reltol ‖r₀‖, abstol), at most n iterations.
using LinearAlgebra

# Below the threshold of the row product every pass runs inline, as a single
# chunk; above it, every pass uses the same chunks, so that each thread keeps
# working on the entries it wrote in the previous pass.
vector_chunk(n::Int64) = n < ROW_PRODUCT_THREAD_THRESHOLD ? n : VECTOR_CHUNK

function resize_workspace!(ws::ConjugateGradientWorkspace, n::Int64)
    for v in (ws.x, ws.r, ws.u, ws.c, ws.z)
        length(v) == n || resize!(v, n)
    end
    return nothing
end

# r ⋅ z by chunks, with the partial sums in chunk order.
function chunk_dot(r::Vector{Float64}, z::Vector{Float64})
    partials = chunk_results(Float64, () -> nothing, length(r); chunk=vector_chunk(length(r))) do range, _
        s = 0.0
        @inbounds for i in range
            s += r[i] * z[i]
        end
        return s
    end
    return sum(partials)
end

# c = A u and u ⋅ c, by chunks of rows of the symmetric matrix.
function product_and_dot!(c::Vector{Float64}, A::SymmetricRowOperator, u::Vector{Float64})
    rowptr, colval, nzval = A.rowptr, A.colval, A.nzval
    partials = chunk_results(Float64, () -> nothing, A.n; chunk=vector_chunk(A.n)) do rows, _
        row_product!(c, rowptr, colval, nzval, u, rows)
        s = 0.0
        @inbounds for i in rows
            s += u[i] * c[i]
        end
        return s
    end
    return sum(partials)
end

# x += α u and r -= α c, with r ⋅ r and r ⋅ z for z = D⁻¹ r (z = r without a
# preconditioner), in one pass.
function update_and_dots!(
    x::Vector{Float64},
    r::Vector{Float64},
    u::Vector{Float64},
    c::Vector{Float64},
    α::Float64,
    inverse_diagonal::Vector{Float64},
)
    jacobi = !isempty(inverse_diagonal)
    # The two forms are separate loops: a branch inside the loop kept it from
    # vectorizing (0.22 against 0.07 ms per pass on 45 000 entries).
    partials = chunk_results(NTuple{2,Float64}, () -> nothing, length(x); chunk=vector_chunk(length(x))) do range, _
        rr = 0.0
        rz = 0.0
        if jacobi
            @inbounds for i in range
                x[i] += α * u[i]
                r[i] -= α * c[i]
                ri = r[i]
                rr += ri * ri
                rz += ri * ri * inverse_diagonal[i]
            end
        else
            @inbounds for i in range
                x[i] += α * u[i]
                r[i] -= α * c[i]
                ri = r[i]
                rr += ri * ri
            end
            rz = rr
        end
        return (rr, rz)
    end
    return sum(first, partials), sum(last, partials)
end

# u = z + β u with z = D⁻¹ r.
function direction!(u::Vector{Float64}, r::Vector{Float64}, β::Float64, inverse_diagonal::Vector{Float64})
    jacobi = !isempty(inverse_diagonal)
    chunk_results(Nothing, () -> nothing, length(u); chunk=vector_chunk(length(u))) do range, _
        if jacobi
            @inbounds for i in range
                u[i] = r[i] * inverse_diagonal[i] + β * u[i]
            end
        else
            @inbounds for i in range
                u[i] = r[i] + β * u[i]
            end
        end
        return nothing
    end
    return nothing
end

function conjugate_gradient!(
    ws::ConjugateGradientWorkspace, A::SymmetricRowOperator, b::Vector{Float64}, abstol::Float64, reltol::Float64
)
    n = A.n
    length(b) == n || norma_abort("The right-hand side does not match the matrix.")
    resize_workspace!(ws, n)
    x, r, u, c = ws.x, ws.r, ws.u, ws.c
    inverse_diagonal = ws.inverse_diagonal
    fill!(x, 0.0)
    fill!(u, 0.0)
    fill!(c, 0.0)
    copyto!(r, b)
    rr, ρ = update_and_dots!(x, r, u, c, 0.0, inverse_diagonal)
    residual = sqrt(rr)
    tolerance = max(reltol * residual, abstol)
    β = 0.0
    iteration = 0
    while iteration < n && residual > tolerance
        direction!(u, r, β, inverse_diagonal)
        α = ρ / product_and_dot!(c, A, u)
        rr, ρ_new = update_and_dots!(x, r, u, c, α, inverse_diagonal)
        residual = sqrt(rr)
        β = ρ_new / ρ
        ρ = ρ_new
        iteration += 1
    end
    ws.iterations = iteration
    ws.residual = residual
    return x
end

# The conjugate gradient with a general preconditioner, applied by
# `apply!(z, r)` into the workspace vector z: the Cholesky factor of a
# previous matrix, whose application is a pair of triangular solves.
function conjugate_gradient!(
    ws::ConjugateGradientWorkspace,
    A::SymmetricRowOperator,
    b::Vector{Float64},
    abstol::Float64,
    reltol::Float64,
    apply!::Function,
)
    n = A.n
    length(b) == n || norma_abort("The right-hand side does not match the matrix.")
    resize_workspace!(ws, n)
    x, r, u, c, z = ws.x, ws.r, ws.u, ws.c, ws.z
    none = Float64[]
    fill!(x, 0.0)
    fill!(u, 0.0)
    fill!(c, 0.0)
    copyto!(r, b)
    rr, _ = update_and_dots!(x, r, u, c, 0.0, none)
    residual = sqrt(rr)
    tolerance = max(reltol * residual, abstol)
    ρ = 0.0
    iteration = 0
    while iteration < n && residual > tolerance
        apply!(z, r)
        ρ_new = chunk_dot(r, z)
        β = iteration == 0 ? 0.0 : ρ_new / ρ
        ρ = ρ_new
        direction!(u, z, β, none)
        α = ρ / product_and_dot!(c, A, u)
        rr, _ = update_and_dots!(x, r, u, c, α, none)
        residual = sqrt(rr)
        iteration += 1
    end
    ws.iterations = iteration
    ws.residual = residual
    return x
end

# The Cholesky factor of A as the preconditioner: kept in the direct solver
# cache and recomputed when the pattern changes, when the previous solve took
# more than `preconditioner_refresh` iterations (the factor has gone stale),
# or when there is none. With the symbolic analysis reused, a refactorization
# costs the numeric phase only. A matrix that is not positive definite falls
# back to the Jacobi preconditioner for that solve.
function cholesky_preconditioner!(solver::HessianMinimizer, A::SparseMatrixCSC{Float64,Int64})
    cache = solver.direct_cache
    ws = solver.cg_workspace
    factor = cache.factor
    same_pattern = factor !== nothing && cache.colptr === A.colptr && cache.rowval === A.rowval
    stale = ws.iterations > solver.preconditioner_refresh
    if factor === nothing || !same_pattern || stale
        As = Symmetric(A)
        factor = if same_pattern
            cholesky!(factor, As; check=false)
        else
            cholesky(As; perm=fill_reducing_ordering(A), check=false)
        end
        if issuccess(factor) == false
            norma_log(4, :solve, "Cholesky preconditioner: matrix not positive definite, using Jacobi")
            cache.factor = nothing
            return nothing
        end
        cache.factor = factor
        cache.colptr = A.colptr
        cache.rowval = A.rowval
    end
    return factor
end

# The inverse diagonal of the reduced matrix for the Jacobi preconditioner,
# from the diagonal positions of the reduced structure when A is the reduced
# matrix (the search of diag(A) and its sparse broadcast cost 3 ms, as much as
# a solve), else from diag(A).
function jacobi_preconditioner!(solver::HessianMinimizer, A::SparseMatrixCSC{Float64,Int64})
    ws = solver.cg_workspace
    n = size(A, 1)
    resize!(ws.inverse_diagonal, n)
    reduced = solver.reduced
    if reduced !== nothing && A.colptr === reduced.colptr
        nzval = A.nzval
        positions = reduced.diagonal
        @inbounds for i in 1:n
            ws.inverse_diagonal[i] = 1.0 / nzval[positions[i]]
        end
    else
        ws.inverse_diagonal .= 1.0 ./ Vector(diag(A))
    end
    return nothing
end

# The conjugate gradient of a Hessian minimizer on the reduced matrix A, with
# its preconditioner and tolerances. The iterate returned is the workspace
# vector, consumed by the caller before the next solve.
function conjugate_gradient!(solver::HessianMinimizer, A::SparseMatrixCSC{Float64,Int64}, b::Vector{Float64})
    ws = solver.cg_workspace
    abstol = solver.linear_solver_absolute_tolerance
    reltol = solver.linear_solver_relative_tolerance
    if solver.preconditioner == :cholesky
        factor = cholesky_preconditioner!(solver, A)
        if factor !== nothing
            empty!(ws.inverse_diagonal)
            apply!(z, r) = ldiv!(z, factor, r)
            return conjugate_gradient!(ws, SymmetricRowOperator(A), b, abstol, reltol, apply!)
        end
    end
    if solver.preconditioner == :none
        empty!(ws.inverse_diagonal)
    else
        jacobi_preconditioner!(solver, A)
    end
    return conjugate_gradient!(ws, SymmetricRowOperator(A), b, abstol, reltol)
end
