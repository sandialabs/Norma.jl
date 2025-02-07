# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type TimeIntegrator end
abstract type StaticTimeIntegrator <: TimeIntegrator end
abstract type DynamicTimeIntegrator <: TimeIntegrator end

mutable struct QuasiStatic <: StaticTimeIntegrator
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    runtime_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    prev_time::Float64
    time::Float64
    stop::Int64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    stored_energy::Float64
    initial_equilibrium::Bool
end

mutable struct Newmark <: DynamicTimeIntegrator
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    runtime_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    prev_time::Float64
    time::Float64
    stop::Int64
    β::Float64
    γ::Float64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    disp_pre::Vector{Float64}
    velo_pre::Vector{Float64}
    stored_energy::Float64
    kinetic_energy::Float64
end

mutable struct CentralDifference <: DynamicTimeIntegrator
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    runtime_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    user_time_step::Float64
    stable_time_step::Float64
    prev_time::Float64
    time::Float64
    stop::Int64
    CFL::Float64
    γ::Float64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    stored_energy::Float64
    kinetic_energy::Float64
end
