# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type TimeIntegrator end
abstract type StaticTimeIntegrator <: TimeIntegrator end
abstract type DynamicTimeIntegrator <: TimeIntegrator end
abstract type ExplicitDynamicTimeIntegrator <: DynamicTimeIntegrator end

mutable struct QuasiStatic <: StaticTimeIntegrator
    prev_time::Float64
    time::Float64
    time_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    displacement::Vector{Float64}
    velocity::Vector{Float64}
    acceleration::Vector{Float64}
    prev_disp::Vector{Float64}
    prev_velo::Vector{Float64}
    prev_acce::Vector{Float64}
    prev_∂Ω_f::Vector{Float64}
    stored_energy::Float64
    initial_equilibrium::Bool
end

mutable struct Newmark <: DynamicTimeIntegrator
    prev_time::Float64
    time::Float64
    time_step::Float64
    minimum_time_step::Float64
    maximum_time_step::Float64
    decrease_factor::Float64
    increase_factor::Float64
    β::Float64
    γ::Float64
    # HHT-α numerical dissipation: ᾱ ∈ [0, 1/3] (0 = plain Newmark). The
    # internal force is evaluated at the blended state
    # (1-ᾱ) u_{n+1} + ᾱ u_n, with γ = 1/2 + ᾱ and β = (1+ᾱ)²/4, giving
    # high-frequency dissipation ρ∞ = (1-ᾱ)/(1+ᾱ) with second-order
    # accuracy. External and interface forces are not blended.
    hht_alpha::Float64
    hht_disp_prev::Vector{Float64}
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
end

mutable struct CentralDifference <: ExplicitDynamicTimeIntegrator
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
end
include("opinf/opinf_time_integrator_types.jl")
