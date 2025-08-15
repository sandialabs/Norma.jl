# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

mutable struct RomNewmark <: DynamicTimeIntegrator
    prev_time::Float64
    time::Float64
    time_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    β::Float64
    γ::Float64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    disp_pre::Vector{Float64}
    velo_pre::Vector{Float64}
    prev_disp::Vector{Float64}
    prev_velo::Vector{Float64}
    prev_acce::Vector{Float64}
    prev_∂Ω_f::Vector{Float64}
    stored_energy::Float64
    kinetic_energy::Float64
    fom_integrator::Newmark
end


mutable struct RomCentralDifference <: ExplicitDynamicTimeIntegrator
    prev_time::Float64
    time::Float64
    time_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    CFL::Float64
    γ::Float64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    prev_disp::Vector{Float64}
    prev_velo::Vector{Float64}
    prev_acce::Vector{Float64}
    prev_∂Ω_f::Vector{Float64}
    stored_energy::Float64
    kinetic_energy::Float64
    fom_integrator::CentralDifference
end
