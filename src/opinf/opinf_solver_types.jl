# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

mutable struct RomHessianMinimizer <: Minimizer
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
    fom_solver::HessianMinimizer
end

mutable struct RomExplicitSolver <: Explicit
    value::Float64
    gradient::Vector{Float64}
    lumped_hessian::Vector{Float64}
    solution::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
    fom_solver::ExplicitSolver
end

