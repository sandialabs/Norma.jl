# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type TimeController end
abstract type SingleTimeController <: TimeController end
abstract type MultiDomainTimeController <: TimeController end

# Identifier of one Schwarz interface: (own subdomain, partner subdomain, side
# set), which is what the relaxation state below is keyed by. The partner index
# alone does not identify an interface: a subdomain that is the partner of two
# interfaces (the gauge between the two holders of the laser weld example) made
# both of them blend and overwrite a single shared iterate.
const RelaxationKey = Tuple{Int64,Int64,Int64}

mutable struct SolidMultiDomainTimeController <: MultiDomainTimeController
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
    relaxation_method::Symbol
    aitken_N0::Int
    naive_stabilized::Bool
    # Relaxation state for the Schwarz fixed-point iteration, stored per
    # interface (RelaxationKey) AND per substep time slot within the current
    # stop: [key][slot], slots keyed by substep time in lambda_time and cleared
    # at every stop. In a windowed stop (controller step larger than the relaxed
    # side's time step) the relaxed side applies its coupling BC once per
    # substep; a single per-interface vector would blend iterates across TIME
    # instead of across Schwarz iterations — a causal low-pass on the exchanged data
    # that shifts the converged fixed point (measured on the cantilever
    # benchmark: monotone energy drain growing with the window under fixed θ,
    # and phase-lead energy injection under Aitken, whose residuals then compare
    # different times). With one substep per stop a single slot exists and this
    # reduces exactly to the classical per-interface iteration.
    lambda_time::Dict{RelaxationKey,Vector{Float64}}
    lambda_disp::Dict{RelaxationKey,Vector{Vector{Float64}}}
    lambda_velo::Dict{RelaxationKey,Vector{Vector{Float64}}}
    lambda_acce::Dict{RelaxationKey,Vector{Vector{Float64}}}
    aitken_prev_residual_disp::Dict{RelaxationKey,Vector{Vector{Float64}}}
    aitken_prev_residual_velo::Dict{RelaxationKey,Vector{Vector{Float64}}}
    aitken_prev_residual_acce::Dict{RelaxationKey,Vector{Vector{Float64}}}
    aitken_theta_disp::Dict{RelaxationKey,Vector{Float64}}
    # Previous interface displacement iterate g^(n-1), used by the Aitken-secant
    # (paper eq. 9) variant to form d^(n) = g^(n) - g^(n-1) directly.
    aitken_prev_lambda_disp::Dict{RelaxationKey,Vector{Vector{Float64}}}
    is_schwarz::Bool
    schwarz_contact::Bool
    active_contact::Bool
    contact_hist::Vector{Bool}
    schwarz_iters::Vector{Int64}
    use_interface_predictor::Bool
    predictor_disp::Vector{Vector{Float64}}
    predictor_velo::Vector{Vector{Float64}}
    predictor_acce::Vector{Vector{Float64}}
    predictor_∂Ω_f::Vector{Vector{Float64}}
    prev_stop_disp::Vector{Vector{Float64}}
    prev_stop_∂Ω_f::Vector{Vector{Float64}}
    # Set during a Schwarz iteration in which a relaxation factor of essentially
    # zero froze an interface iterate. Such an iteration re-solves every subdomain against
    # unchanged coupling data and reproduces the previous solution, so its
    # update is not evidence of convergence: the interface residual it leaves
    # behind is invisible to the displacement-based criterion. Cleared at the
    # start of each Schwarz iteration; see `frozen_relaxation_update!`.
    relaxation_frozen::Bool
end

mutable struct SolidSingleDomainTimeController <: SingleTimeController
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    prev_time::Float64
    num_stops::Int64
    stop::Int64
end

mutable struct SingleDomainSimulation <: Simulation
    name::String
    params::Parameters
    controller::SolidSingleDomainTimeController
    integrator::TimeIntegrator
    solver::Solver
    model::Model
    failed::Bool
    # Back-reference to the MultiDomainSimulation that owns this subsim, and a
    # stable slot id. Both are nothing when the simulation is run standalone.
    # Typed against the abstract Simulation because MultiDomainSimulation is
    # defined below.
    parent::Union{Nothing,Simulation}
    handle::Union{Nothing,DomainHandle}
    # Mid-run swap plans on this single-domain simulation.  Empty when the sim
    # is itself a subsim of a MultiDomainSimulation (which carries its own
    # swap list at the parent level).
    swaps::Vector{SwapPlan}
end

mutable struct MultiDomainSimulation <: Simulation
    name::String
    params::Parameters
    controller::MultiDomainTimeController
    num_domains::Int64
    subsims::Vector{SingleDomainSimulation}
    handle_by_name::Dict{String,DomainHandle}
    name_by_handle::Vector{String}
    swaps::Vector{SwapPlan}
    failed::Bool
end
