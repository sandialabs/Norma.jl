# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
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

mutable struct HessianMinimizer <: Minimizer
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    linear_solver_absolute_tolerance::Float64
    linear_solver_relative_tolerance::Float64
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

include("opinf/opinf_solver_types.jl")

