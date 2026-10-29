# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

abstract type Solver end
abstract type Minimizer <: Solver end
abstract type MatrixFree <: Solver end
abstract type Explicit <: Solver end
abstract type Step end
abstract type LineSearch end

using SparseArrays

struct BackTrackLineSearch <: LineSearch
    backtrack_factor::Float64
    decrease_factor::Float64
    max_iters::Int64
end

# Cached sparse factorization of the direct linear solver, a Cholesky factor
# of CHOLMOD or an LU factor of UMFPACK whose numeric part is freed after each
# solve. Its ordering and symbolic analysis are reused by a numeric
# refactorization while the pattern of the reduced system is unchanged, which
# is every Newton iteration once the free degrees of freedom are fixed.
mutable struct DirectSolverCache
    factor::Any
    colptr::Vector{Int64}
    rowval::Vector{Int64}
end

# The system on the free degrees of freedom of an assembled matrix. Extracting
# A[free, free] searched every column of the matrix at every Newton iteration
# (5.1 ms on 45 000 free degrees of freedom). The structure of the reduced
# matrix and the position in the full matrix of each of its nonzeros are
# computed once per full pattern (recognized by the identity of its structure
# arrays) and mask, and each iteration gathers the values. The reduced matrix
# shares these structure arrays across iterations, so the direct solver
# recognizes an unchanged pattern by identity and reuses its symbolic
# factorization.
mutable struct ReducedSystem
    full_colptr::Vector{Int64}
    full_rowval::Vector{Int64}
    free::BitVector
    free_index::Vector{Int64}
    colptr::Vector{Int64}
    rowval::Vector{Int64}
    nz_source::Vector{Int64}
    nzval::Vector{Float64}
    rhs::Vector{Float64}
    # Position of the diagonal entry of every column among the nonzeros of the
    # reduced matrix, for the Jacobi preconditioner.
    diagonal::Vector{Int64}
end

# Workspace of the conjugate gradient (conjugate_gradient.jl), reused across
# solves: the iterate, the residual, the search direction, its product with
# the matrix, and the inverse diagonal of the Jacobi preconditioner (empty
# without one). The last solve's iteration count and residual norm are kept
# for reports and tests.
mutable struct ConjugateGradientWorkspace
    x::Vector{Float64}
    r::Vector{Float64}
    u::Vector{Float64}
    c::Vector{Float64}
    z::Vector{Float64}
    inverse_diagonal::Vector{Float64}
    iterations::Int64
    residual::Float64
end

function ConjugateGradientWorkspace()
    return ConjugateGradientWorkspace(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], 0, 0.0)
end

mutable struct HessianMinimizer <: Minimizer
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    linear_solver::String
    linear_solver_absolute_tolerance::Float64
    linear_solver_relative_tolerance::Float64
    direct_cache::DirectSolverCache
    reduced::Union{Nothing,ReducedSystem}
    preconditioner::Symbol
    # Iterations of a solve preconditioned by the Cholesky factor above which
    # the factor is recomputed at the next solve.
    preconditioner_refresh::Int64
    cg_workspace::ConjugateGradientWorkspace
    unconverged_action::String
    value::Float64
    gradient::Vector{Float64}
    hessian::SparseMatrixCSC{Float64,Int64}
    solution::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
    line_search::BackTrackLineSearch
    use_line_search::Bool
end

mutable struct ExplicitSolver <: Explicit
    value::Float64
    gradient::Vector{Float64}
    lumped_hessian::Vector{Float64}
    solution::Vector{Float64}
    # The step of the free degrees of freedom, reused across steps.
    step_vector::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
end

mutable struct SteepestDescent <: MatrixFree
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    value::Float64
    gradient::Vector{Float64}
    solution::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
    line_search::BackTrackLineSearch
    # Energy stagnation exit for mesh smoothing (issue #220).  A window of 0
    # disables the criterion.  The history holds the energy at the end of each
    # iteration of the current solve; the streak counts consecutive iterations
    # whose windowed relative energy decrease fell below the tolerance.
    energy_stagnation_window::Int64
    energy_stagnation_tolerance::Float64
    energy_history::Vector{Float64}
    energy_stagnation_streak::Int64
    stagnated::Bool
    unconverged_action::String
end

struct NewtonStep <: Step
    step_length::Float64
end

struct ExplicitStep <: Step
    step_length::Float64
end

struct SteepestDescentStep <: Step
    step_length::Float64
end

# Limited-memory BFGS step for matrix-free minimization (e.g. EMS mesh
# smoothing).  Stores a short history of (s, y) = (Δsolution, Δgradient) pairs
# and reconstructs a quasi-Newton search direction via the two-loop recursion,
# giving curvature-aware steps without ever assembling a Hessian.  Mutable
# because the history evolves across iterations; `reset_step!` clears it at the
# start of each nonlinear solve.
mutable struct LBFGSStep <: Step
    step_length::Float64
    memory::Int64
    s_history::Vector{Vector{Float64}}
    y_history::Vector{Vector{Float64}}
    rho_history::Vector{Float64}
    prev_solution::Vector{Float64}
    prev_gradient::Vector{Float64}
    initialized::Bool
end

include("opinf/opinf_solver_types.jl")
