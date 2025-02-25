# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type SingleController end
abstract type SchwarzController end

mutable struct SolidSchwarzController <: SchwarzController
    num_domains::Int64
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    prev_time::Float64
    same_step::Bool
    num_stops::Int64
    stop::Int64
    converged::Bool
    iteration_number::Int64
    stop_disp::Vector{Vector{Float64}}
    stop_velo::Vector{Vector{Float64}}
    stop_acce::Vector{Vector{Float64}}
    stop_∂Ω_f::Vector{Vector{Float64}}
    schwarz_disp::Vector{Vector{Float64}}
    schwarz_velo::Vector{Vector{Float64}}
    schwarz_acce::Vector{Vector{Float64}}
    time_hist::Vector{Vector{Float64}}
    disp_hist::Vector{Vector{Vector{Float64}}}
    velo_hist::Vector{Vector{Vector{Float64}}}
    acce_hist::Vector{Vector{Vector{Float64}}}
    ∂Ω_f_hist::Vector{Vector{Vector{Float64}}}
    relaxation_parameter::Float64
    naive_stabilized::Bool
    lambda_disp::Vector{Vector{Float64}}
    lambda_velo::Vector{Vector{Float64}}
    lambda_acce::Vector{Vector{Float64}}
    schwarz_contact::Bool
    active_contact::Bool
    contact_hist::Vector{Bool}
    convergence_hist::Array{Float64}
end

mutable struct SolidSingleController <: SingleController
    time::Float64
    prev_time::Float64
    stop::Int64
    stop_disp::Vector{Float64}
    stop_velo::Vector{Float64}
    stop_acce::Vector{Float64}
    stop_∂Ω_f::Vector{Float64}
end
