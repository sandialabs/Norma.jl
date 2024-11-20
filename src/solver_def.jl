abstract type Solver end
abstract type Minimizer <: Solver end
abstract type Step end
abstract type StepLength end

using SparseArrays

mutable struct HessianMinimizer <: Minimizer
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    value::Float64
    gradient::Vector{Float64}
    hessian::SparseMatrixCSC{Float64,Int64}
    solution::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
end

mutable struct ExplicitSolver <: Solver
    value::Float64
    gradient::Vector{Float64}
    lumped_hessian::Vector{Float64}
    solution::Vector{Float64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
end

mutable struct SteepestDescent <: Solver
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
end

struct FixedStepLength <: StepLength
    value::Float64
end

mutable struct AdaptiveStepLength <: StepLength
    value::Float64
end


struct NewtonStep <: Step
    step_length::FixedStepLength
end

struct ExplicitStep <: Step
    step_length::FixedStepLength
end

struct SteepestDescentStep <: Step
    step_length::StepLength
    max_step_length::Float64
    min_step_length::Float64
    reduction_factor::Float64
    increase_factor::Float64
end
