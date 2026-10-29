# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Schwarz coupling boundary conditions for multi-domain simulations.
#
# Includes: overlap, non-overlap (Dirichlet-Neumann), Robin-Robin,
# and contact Schwarz BCs. Projectors, time history interpolation,
# and interface force/displacement transfer.

using LinearAlgebra: dot

get_fom_model(sim::Simulation) = sim.model isa RomModel ? sim.model.fom_model : sim.model

# Resolve the partner (coupled) subsim of a Schwarz BC through its parent via
# the stable handle, so that swapping a subsim in its parent slot is transparent
# to every BC.
coupled_subsim_of(bc::SolidMechanicsSchwarzBoundaryCondition) = bc.parent.subsims[bc.coupled_handle.id]
self_subsim_of(bc::SolidMechanicsSchwarzBoundaryCondition) = bc.parent.subsims[bc.self_handle.id]

# The partner's coupling condition, read from the partner's own model: for a
# reduced order partner that is the RomModel, not its full order model.
function coupled_bc_of(bc::SolidMechanicsSchwarzBoundaryCondition)
    return coupled_subsim_of(bc).model.boundary_conditions[bc.coupled_bc_index]
end

# Relative floor of the differences that enter an Aitken factor. A difference
# δ of two successive residuals, or d of two successive iterates, is resolved
# when ‖δ‖ > AITKEN_RELATIVE_FLOOR ‖x‖, with x the current interface iterate;
# below that, ‖δ‖ is at most 45 unit roundoffs times ‖x‖, δ carries no secant
# information, and the factor falls back to the input relaxation parameter
# (secant form) or to the previous factor (recursive form). The
# Aitken factor itself is not clamped. The floor was absolute, 1e-20 on ‖δ‖²,
# until October 2026: on the matrix inclusion at 5% strain it was reached at
# interface differences of about 1e-10 m, and the fallback factor 0.5, at which
# that Dirichlet-Neumann iteration diverges, then prevented convergence below a
# relative tolerance of about 1e-10.
const AITKEN_RELATIVE_FLOOR = 1.0e-14

# True when a squared difference norm `difference_sq` is resolved relative to the
# squared norm `iterate_sq` of the iterate (see AITKEN_RELATIVE_FLOOR). A zero
# difference is never resolved, also when the iterate is zero.
function aitken_resolved(difference_sq::Float64, iterate_sq::Float64)
    return difference_sq > AITKEN_RELATIVE_FLOOR^2 * iterate_sq
end

# A relaxation factor at or below this magnitude leaves the interface iterate
# unchanged to working precision. Every subdomain then re-solves against the
# coupling data it already used and returns the solution it already had, so the
# displacement-based Schwarz criterion sees no update and reads it as
# convergence, however large the interface residual still is. Record such a
# Schwarz iteration on the controller so the criterion can refuse to convert it into
# convergence (see update_schwarz_convergence_criterion).
const FROZEN_THETA = 1.0e-12

function frozen_relaxation_update!(controller::MultiDomainTimeController, θ::Float64)
    abs(θ) <= FROZEN_THETA && (controller.relaxation_frozen = true)
    return θ
end

# Name a relaxation method the way the `relaxation` input value spells it, so a
# log line maps straight back to the input file. Defined once and used by both
# the start-up echo and the per-iteration factor lines, which previously spelled
# the same method two ways (issue #217).
# Default mixing parameter β of Anderson acceleration (see `relaxation: anderson`).
const ANDERSON_DEFAULT_MIXING = 0.5

function relaxation_method_name(method::Symbol)
    method === :anderson && return "Anderson"
    method === :aitken_recursive && return "Aitken recursive"
    method === :aitken_secant && return "Aitken secant"
    return "fixed"
end

# Identify the interface a Schwarz BC belongs to, for the relaxation state (see
# RelaxationKey): one key per (own subdomain, partner subdomain, side set), not
# one per partner subdomain, which aliases whenever a subdomain is the partner
# of more than one interface.
function relaxation_key(bc::SolidMechanicsSchwarzBoundaryCondition)
    return (bc.self_handle.id, bc.coupled_handle.id, bc.side_set_id)
end

# Fetch the per-slot relaxation state of one interface, creating it on first use.
function relaxation_slots!(state::Dict{RelaxationKey,Vector{Vector{Float64}}}, key::RelaxationKey)
    return get!(() -> Vector{Float64}[], state, key)
end

function relaxation_slots!(state::Dict{RelaxationKey,Vector{Float64}}, key::RelaxationKey)
    return get!(() -> Float64[], state, key)
end

# Resolve the relaxation time slot for interface `key` at substep time t,
# appending a new slot when t is not yet keyed. The relaxation state must be
# compared across Schwarz iterations AT THE SAME SUBSTEP TIME: with a windowed
# controller stop the relaxed side applies its BC once per substep, and a
# single per-interface state would blend iterates across time — a causal
# low-pass on the exchanged data that shifts the converged fixed point. Substep
# times repeat bitwise across Schwarz iterations (subcycle() restores the nominal step every
# pass), so the equality test hits; the isapprox fallback and the fresh-slot
# path only engage if the substep grid is perturbed mid-stop (e.g. adaptive
# stepping), where the new slot restarts that time's relaxation from scratch.
function relaxation_slot!(controller::MultiDomainTimeController, key::RelaxationKey, t::Float64)
    times = relaxation_slots!(controller.lambda_time, key)
    for (k, tk) in enumerate(times)
        if tk == t || isapprox(tk, t; rtol=1.0e-12)
            return k
        end
    end
    push!(times, t)
    if controller.iteration_number > 0
        norma_logf(1, :schwarz, "New relaxation slot at t = %.6e past iteration 0 (substep grid changed?)", t)
    end
    return length(times)
end

# Grow a per-slot state vector to hold slot k. New vector entries are empty —
# the "no previous iterate yet" sentinel used throughout the relaxation code;
# new scalar entries take the supplied default.
function ensure_slot!(state::Vector{Vector{Float64}}, k::Int)
    while length(state) < k
        push!(state, Float64[])
    end
    return nothing
end

function ensure_slot!(state::Vector{Float64}, k::Int, default::Float64)
    while length(state) < k
        push!(state, default)
    end
    return nothing
end

# Aitken acceleration applies only to single-slot (same-step) stops. In a
# windowed stop the Schwarz iteration map couples all time slots, and every Aitken policy
# measured on the 10 ms cantilever benchmark fails or loses there: per-slot
# θs take unclamped excursions that hand the solver a divergent interface
# force (deaths at 3.9/5.7 ms where fixed θ and θ = 1 ran clean), and pooling
# the residual inner products over the iteration's slots (waveform Aitken,
# Irons–Tuck base frozen per Schwarz iteration) survived on the implicit pair but cost
# 55 iterations/stop against 47.5 for fixed θ = 0.5 and 19.1 for θ = 1, while
# locking θ persistently negative on the explicit pair (dead at 0.33 ms).
# These measurements were made with the impedance transmission condition
# t + Z u̇ + α W u = g (Z the material impedance ρ c_p, u̇ the interface
# velocity), which has since been removed from Norma; the policy is retained for
# the remaining couplings. Windowed stops therefore use the configured
# relaxation parameter.
function aitken_applies(controller::MultiDomainTimeController, key::RelaxationKey)
    return length(relaxation_slots!(controller.lambda_time, key)) <= 1
end

# Returns the relaxation factor θ applied to interp_disp for this Schwarz
# iterate. Fixed mode returns the user-configured constant; Aitken-recursive
# mode uses Irons–Tuck with the previous residual stored on the controller,
# per interface and per substep time slot.
function relaxation_aitken_recursive_theta!(
    controller::MultiDomainTimeController,
    key::RelaxationKey,
    slot_k::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    if controller.relaxation_method !== :aitken_recursive || !aitken_applies(controller, key)
        return controller.relaxation_parameter
    end
    aitken_N0 = controller.aitken_N0
    residual_slots = relaxation_slots!(controller.aitken_prev_residual_disp, key)
    theta_slots = relaxation_slots!(controller.aitken_theta_disp, key)
    ensure_slot!(residual_slots, slot_k)
    ensure_slot!(theta_slots, slot_k, controller.relaxation_parameter)
    if iter < aitken_N0
        # Below N0 the input theta is the relaxation factor, as it is for the
        # secant form. This used to store the input theta as theta^(n-1) for the
        # recursion but apply and report 1.0, so the first N0 Schwarz iterations of every
        # run were unrelaxed whatever the input file asked for (issue #218).
        θ = controller.relaxation_parameter
        theta_slots[slot_k] = θ
        residual_slots[slot_k] = Float64[]
        norma_logf(1, :schwarz, "%s θ[iter=%d] = %.4e", relaxation_method_name(controller.relaxation_method), iter, θ)
        return θ
    end
    residual = interp_disp .- lambda_prev
    prev_residual = residual_slots[slot_k]
    θ_prev = theta_slots[slot_k]
    θ = controller.relaxation_parameter
    if !isempty(prev_residual) && length(prev_residual) == length(residual)
        δ = residual .- prev_residual
        δ_sq = dot(δ, δ)
        if aitken_resolved(δ_sq, dot(lambda_prev, lambda_prev))
            θ = -θ_prev * dot(prev_residual, δ) / δ_sq
        else
            θ = θ_prev
        end
    end
    residual_slots[slot_k] = residual
    theta_slots[slot_k] = θ
    norma_logf(1, :schwarz, "%s θ[iter=%d] = %.4e", relaxation_method_name(controller.relaxation_method), iter, θ)
    return θ
end

# Aitken-secant factor (non-recursive, the original-paper form):
# Sambataro-Tezaur eq. (9), equivalently Deparis-Discacciati-Quarteroni:
#
#   ρ^(n) = - (d^(n) · δ^(n)) / ‖δ^(n)‖²,
#
# with the interface jump (fixed-point residual) r^(n) = E^(n) = T(g^(n)) - g^(n),
# δ^(n) = r^(n) - r^(n-1), and d^(n) = g^(n) - g^(n-1) formed directly from the
# stored interface iterates. This minimizes ‖d^(n) + ρ δ^(n)‖² and, unlike the
# Aitken-recursive Irons-Tuck form in `relaxation_aitken_recursive_theta!`, does
# not carry θ^(n-1), so it is immune to the init/N0 bookkeeping required for the
# recursion to stay exact.
function relaxation_aitken_secant_theta!(
    controller::MultiDomainTimeController,
    key::RelaxationKey,
    slot_k::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    if !aitken_applies(controller, key)
        return controller.relaxation_parameter
    end
    residual_slots = relaxation_slots!(controller.aitken_prev_residual_disp, key)
    lambda_slots = relaxation_slots!(controller.aitken_prev_lambda_disp, key)
    ensure_slot!(residual_slots, slot_k)
    ensure_slot!(lambda_slots, slot_k)
    residual = interp_disp .- lambda_prev                 # r^(n) = T(g^(n)) - g^(n)
    prev_residual = residual_slots[slot_k]                # r^(n-1)
    prev_lambda = lambda_slots[slot_k]                    # g^(n-1)
    # For iter < N0 use the input ρ^(1) = relaxation_parameter (paper's N0 idea).
    θ = controller.relaxation_parameter
    if iter >= controller.aitken_N0 &&
        !isempty(prev_residual) &&
        length(prev_residual) == length(residual) &&
        !isempty(prev_lambda) &&
        length(prev_lambda) == length(lambda_prev)
        δ = residual .- prev_residual                     # δ^(n)
        d = lambda_prev .- prev_lambda                    # d^(n) = g^(n) - g^(n-1)
        δ_sq = dot(δ, δ)
        # A pair whose iterate did not move carries no secant information: d = 0
        # gives θ = 0 exactly, which freezes the interface. That is not a rare
        # numerical coincidence but the normal state of the first pair on the
        # DN and overlap path, where a fresh slot seeds g^(n-1) with the very
        # datum it is compared against (see the `isempty` fallback in apply_bc),
        # so r = 0, the written iterate equals the incoming one whatever θ is,
        # and the next Schwarz iteration finds g^(n) == g^(n-1). Fall back to the input
        # theta; the iteration after that has two genuine iterates to differentiate.
        # The test is relative to the iterate (AITKEN_RELATIVE_FLOOR), so that
        # differences far below the size of the iterate but far above its
        # rounding error still give a factor.
        iterate_sq = dot(lambda_prev, lambda_prev)
        if aitken_resolved(δ_sq, iterate_sq) && aitken_resolved(dot(d, d), iterate_sq)
            # Pure Aitken-secant factor, paper eq. (9): no value clamp, so the factor is
            # free to take the large/negative excursions that accelerate (or
            # damp) convergence. The floor above only guards the division when
            # successive residuals are equal to the rounding error of the iterate.
            θ = -dot(d, δ) / δ_sq
        end
    end
    residual_slots[slot_k] = residual
    lambda_slots[slot_k] = copy(lambda_prev)
    norma_logf(1, :schwarz, "%s θ[iter=%d] = %.4e", relaxation_method_name(controller.relaxation_method), iter, θ)
    return θ
end

# ---------------------------------------------------------------------------
# Global Aitken factor (`aitken scope: global`).
#
# The per-interface factors above give every relaxed interface its own θ from
# its own residuals. Under `aitken scope: global` one factor is formed for all
# relaxed interfaces from their residuals and iterates stacked into one vector,
# as the stacked Anderson step does (anderson_data!): the secant form
#
#   θ = -(Σ_i d_i · δ_i) / Σ_i ‖δ_i‖²,
#
# and the recursive (Irons-Tuck) form
#
#   θ = -θ_prev (Σ_i r_i^(n-1) · δ_i) / Σ_i ‖δ_i‖²,
#
# with the sums over the interfaces i, and each interface's r, δ, and d
# formed on its relaxed datum: the Robin datum, the projected trace of the
# constrained datum, or the projected trace of the partner displacement of
# an unconstrained pair. The first two are the vectors of the per-interface
# factors, so for a single Robin or constrained interface both scopes give
# the same factors; the per-interface factor of an unconstrained pair uses
# the whole partner field instead (see aitken_global_entry). The residual of
# an interface depends on the data of the others, which per-interface
# factors ignore.
#
# The residuals are gathered when the first Aitken-relaxed interface of a
# Schwarz iteration applies its datum, from the partner data at that moment,
# that is, at the end of the previous Gauss-Seidel pass (for an interface
# whose partner has already been solved in the current pass, of the current
# one). The factor is then applied by each interface at its own application,
# to the partner data it reads there, as the per-interface factor is. With
# all relaxed interfaces on one subdomain, as on the inclusion split, the two
# moments coincide. The interfaces are grouped by the relaxed quantity (the
# Robin datum, displacement, velocity), one factor per group, so that no sum
# mixes units. Couplings the global factor does not cover are rejected
# instead of falling back to other factors: the overlap matching coupling
# (both sides of a pair relaxed, each after its partner has solved) and pairs
# with a reduced-order partner at setup (reject_unsupported_global_aitken),
# and stops with several substeps of a relaxed side, where no Aitken factor
# applies (aitken_applies), when they occur (relaxation_theta!).
# ---------------------------------------------------------------------------

# Global secant factor of the stacked vectors: the iterate differences ds, the
# residual differences δs, and the current iterates. `fallback` when either
# difference is not resolved relative to the iterate (aitken_resolved). For
# one interface, the factor of relaxation_aitken_secant_theta!.
function aitken_secant_stacked(ds, δs, iterates, fallback::Float64)
    δ_sq = sum(δ -> dot(δ, δ), δs; init=0.0)
    d_sq = sum(d -> dot(d, d), ds; init=0.0)
    iterate_sq = sum(x -> dot(x, x), iterates; init=0.0)
    if aitken_resolved(δ_sq, iterate_sq) && aitken_resolved(d_sq, iterate_sq)
        return -sum(dot(d, δ) for (d, δ) in zip(ds, δs); init=0.0) / δ_sq
    end
    return fallback
end

# Global recursive (Irons-Tuck) factor of the stacked vectors: the previous
# residuals, the residual differences δs, and the current iterates. The
# previous factor θ_prev when the residual difference is not resolved. For
# one interface, the factor of relaxation_aitken_recursive_theta!.
function aitken_recursive_stacked(prev_residuals, δs, iterates, θ_prev::Float64)
    δ_sq = sum(δ -> dot(δ, δ), δs; init=0.0)
    iterate_sq = sum(x -> dot(x, x), iterates; init=0.0)
    if aitken_resolved(δ_sq, iterate_sq)
        return -θ_prev * sum(dot(r, δ) for (r, δ) in zip(prev_residuals, δs); init=0.0) / δ_sq
    end
    return θ_prev
end

# Per controller: the factor of every interface at the current Schwarz
# iteration, and from the previous one, per interface its residual and
# iterate, and per relaxed quantity the factor (recursive form). Cleared at
# every stop with the rest of the relaxation state (reset_relaxation_state!).
mutable struct AitkenGlobalState
    iteration::Int
    theta::Dict{RelaxationKey,Float64}
    prev_residual::Dict{RelaxationKey,Vector{Float64}}
    prev_iterate::Dict{RelaxationKey,Vector{Float64}}
    prev_theta::Dict{Symbol,Float64}
end

function AitkenGlobalState()
    return AitkenGlobalState(
        -1,
        Dict{RelaxationKey,Float64}(),
        Dict{RelaxationKey,Vector{Float64}}(),
        Dict{RelaxationKey,Vector{Float64}}(),
        Dict{Symbol,Float64}(),
    )
end

const AITKEN_GLOBAL_STATE = IdDict{Any,AitkenGlobalState}()

function reset_aitken_global_state!(controller)
    delete!(AITKEN_GLOBAL_STATE, controller)
    return nothing
end

# Whether the factor of the interface of `bc` is the global one.
function uses_global_aitken(bc::SolidMechanicsSchwarzBoundaryCondition)
    controller = bc.parent.controller
    controller.aitken_scope === :global || return false
    method = controller.relaxation_method
    method === :aitken_secant || method === :aitken_recursive || return false
    return !(bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
end

# The residual and iterate from which the global Aitken factor takes the
# contribution of the interface of `bc`, at the current partner data: the
# Robin datum, the projected trace of the constrained datum, or the projected
# trace of the partner displacement. An interface without a
# stored iterate takes the iterate that its application falls back to (zero
# for the Robin datum, the partner datum otherwise). Nothing when the
# interface is not Aitken-relaxed with the global factor at this stop: not
# its relaxed side, an inactive contact pair, a windowed stop, or not yet
# applied in this stop.
function aitken_global_entry(bc::SolidMechanicsSchwarzBoundaryCondition)
    uses_global_aitken(bc) || return nothing
    controller = bc.parent.controller
    key = relaxation_key(bc)
    if bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
        bc.self_handle.id < bc.coupled_handle.id || return nothing
    elseif is_swappable_dn_schwarz(bc)
        bc.is_dirichlet && !is_direct_dn(bc) || return nothing
        bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact && return nothing
    else
        return nothing
    end
    haskey(controller.lambda_time, key) && aitken_applies(controller, key) || return nothing
    times = controller.lambda_time[key]
    isempty(times) && return nothing
    time = times[1]
    disp_slots = relaxation_slots!(controller.lambda_disp, key)
    ensure_slot!(disp_slots, 1)
    if bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
        rhs = with_partner_state(() -> robin_datum(self_subsim_of(bc).model, bc), bc, time)
        g = isempty(disp_slots[1]) ? zeros(length(rhs)) : disp_slots[1]
        return (quantity=:robin, key=key, residual=rhs .- g, iterate=copy(g))
    end
    disp, velo, _, _ = partner_fields(bc, time)
    if is_constrained_dn(bc)
        relax_velocity = bc.constraint == :velocity
        slots = relax_velocity ? relaxation_slots!(controller.lambda_velo, key) : disp_slots
        ensure_slot!(slots, 1)
        g = relax_velocity ? velo : disp
        x = isempty(slots[1]) ? g : slots[1]
        trace_g = constrained_partner_trace(bc, g)
        trace_x = constrained_partner_trace(bc, x)
        quantity = relax_velocity ? :velocity : :displacement
        return (quantity=quantity, key=key, residual=trace_g .- trace_x, iterate=trace_x)
    end
    x = isempty(disp_slots[1]) ? disp : disp_slots[1]
    # The per-interface factor of an unconstrained pair is formed from the
    # whole partner field, so the side sets of one partner subdomain share it
    # and stacking them would only repeat that field. The global factor is
    # formed from the projected interface trace instead, as the Anderson
    # coefficients are (anderson_entry); a reduced-order partner, whose fields
    # are its reduced coordinates, keeps the whole field.
    coupled_subsim_of(bc).model isa RomModel &&
        return (quantity=:displacement, key=key, residual=disp .- x, iterate=copy(x))
    trace_g = constrained_partner_trace(bc, disp)
    trace_x = constrained_partner_trace(bc, x)
    return (quantity=:displacement, key=key, residual=trace_g .- trace_x, iterate=trace_x)
end

# The global Aitken factor of the interface of `bc`. The first call of a
# Schwarz iteration gathers the residuals of every Aitken-relaxed interface
# (aitken_global_entry) and forms one factor per relaxed quantity; the other
# calls of the iteration read it. As for the per-interface factors,
# `relaxation parameter` is the factor below `aitken N0 parameter` and
# wherever no previous residual exists; an interface not gathered (applied
# for the first time after the first call of its stop) takes it as well.
function aitken_global_theta!(bc::SolidMechanicsSchwarzBoundaryCondition)
    controller = bc.parent.controller
    state = get!(AitkenGlobalState, AITKEN_GLOBAL_STATE, controller)
    iter = controller.iteration_number
    θ_0 = controller.relaxation_parameter
    if state.iteration != iter
        empty!(state.theta)
        recursive = controller.relaxation_method === :aitken_recursive
        entries = []
        for subsim in bc.parent.subsims, other in subsim.model.boundary_conditions
            other isa SolidMechanicsSchwarzBoundaryCondition || continue
            entry = aitken_global_entry(other)
            entry === nothing || push!(entries, entry)
        end
        for quantity in unique(getfield.(entries, :quantity))
            group = [entry for entry in entries if entry.quantity === quantity]
            # Interfaces with a previous residual (and, for the secant form,
            # a previous iterate) of the same size.
            paired = [
                entry for entry in group if
                length(get(state.prev_residual, entry.key, Float64[])) == length(entry.residual) > 0 &&
                    (recursive || length(get(state.prev_iterate, entry.key, Float64[])) == length(entry.iterate))
            ]
            θ = if iter < controller.aitken_N0 || isempty(paired)
                θ_0
            elseif recursive
                aitken_recursive_stacked(
                    [state.prev_residual[entry.key] for entry in paired],
                    [entry.residual .- state.prev_residual[entry.key] for entry in paired],
                    [entry.iterate for entry in paired],
                    get(state.prev_theta, quantity, θ_0),
                )
            else
                aitken_secant_stacked(
                    [entry.iterate .- state.prev_iterate[entry.key] for entry in paired],
                    [entry.residual .- state.prev_residual[entry.key] for entry in paired],
                    [entry.iterate for entry in paired],
                    θ_0,
                )
            end
            for entry in group
                # Below N0 the recursive form starts afresh, as the
                # per-interface one does (relaxation_aitken_recursive_theta!).
                if recursive && iter < controller.aitken_N0
                    delete!(state.prev_residual, entry.key)
                else
                    state.prev_residual[entry.key] = entry.residual
                end
                state.prev_iterate[entry.key] = entry.iterate
                state.theta[entry.key] = θ
            end
            state.prev_theta[quantity] = θ
            norma_logf(
                1,
                :schwarz,
                "%s global θ[iter=%d] %s: %d interfaces = %.4e",
                relaxation_method_name(controller.relaxation_method),
                iter,
                string(quantity),
                length(group),
                θ,
            )
        end
        state.iteration = iter
    end
    return get(state.theta, relaxation_key(bc), θ_0)
end

# The relaxation factor of the interface of `bc` for this Schwarz iterate,
# from its unrelaxed candidate and its previous iterate: the global Aitken
# factor under `aitken scope: global` (aborting in a windowed stop), otherwise
# the per-interface factor of the selected Aitken form, or the fixed factor.
function relaxation_theta!(
    bc::SolidMechanicsSchwarzBoundaryCondition,
    key::RelaxationKey,
    slot_k::Int,
    iter::Int,
    candidate::AbstractVector{Float64},
    previous::AbstractVector{Float64},
)
    controller = bc.parent.controller
    if uses_global_aitken(bc)
        # A stop with several substeps of the relaxed side keys several time
        # slots, where no Aitken factor applies and the per-interface scope
        # falls back to `relaxation parameter`. Under the global scope that
        # fallback would silently replace the requested factor, so the run
        # stops instead. Whether a stop has several substeps is known only at
        # run time (stable step of an explicit integrator, adaptive stepping).
        aitken_applies(controller, key) || norma_abort(
            "`aitken scope: global` is not supported for stops with several substeps of a relaxed side " *
            "(interface $(key), controller time step $(controller.time_step)): the Aitken factor applies " *
            "only when the controller `time step` equals the time step of the relaxed subdomain. Use " *
            "`aitken scope: per interface`, which takes `relaxation parameter` in such stops.",
        )
        return aitken_global_theta!(bc)
    end
    if controller.relaxation_method === :aitken_secant
        return relaxation_aitken_secant_theta!(controller, key, slot_k, iter, candidate, previous)
    end
    return relaxation_aitken_recursive_theta!(controller, key, slot_k, iter, candidate, previous)
end

# ---------------------------------------------------------------------------
# Anderson acceleration (Walker and Ni 2011; the interface quasi-Newton form of
# Degroote, Bathe, and Vierendeels 2009) of the relaxed interface data.
#
# The fixed-point map takes the data x_k that the relaxed sides used to the
# data g_k = G(x_k) that the iteration returns. With the residuals
# f_j = g_j - x_j, their interface traces t_j, and the differences ΔT, ΔX, ΔF
# of the last m + 1 iterates (m the depth),
#   γ = argmin ‖t_k - ΔT γ‖,   x_{k+1} = x_k + β f_k - (ΔX + β ΔF) γ,
# with β the mixing parameter (`relaxation parameter`). One least-squares
# problem is solved for all relaxed interfaces, with their traces stacked
# into one vector, and its coefficients update every interface, as in
# Sambataro and Tezaur. Coefficients computed for each interface separately
# fit one secant model per interface although the residual of an interface
# depends on the data of the others: on an inclusion whose interface is
# split into four side sets they failed to converge, against 23 iterations
# for the stacked problem (issue #233). The interfaces are grouped by the
# relaxed quantity (interface displacement, interface velocity, or Robin
# datum), so that no least-squares problem combines quantities of different
# units. The step is computed when the first relaxed interface of a Schwarz
# iteration applies its datum, from the partner data at that moment; the
# other interfaces read their results. An interface whose partner is solved
# after that moment and before its own solve therefore receives the partner
# data of the previous iteration. The coefficients come from the traces, and
# the update is computed on the partner degrees of freedom that the transfer
# of the interface reads (anderson_footprint); the other entries of the
# relaxed partner field, which the exchange does not read, are those of the
# partner. The histories hold at most m + 1 columns, each with the iterates
# and residuals of the relaxed fields on the footprint and the trace once.
# The least-squares problem is solved by QR of ΔT, dropping the oldest column
# while the condition number of R exceeds ANDERSON_CONDITION_LIMIT. Without a
# difference (first iteration) the update is x + β f, the fixed factor β. For
# a linear map with β = 1 and a depth at least the iteration count,
# x_{k+1} = G of the k-th GMRES iterate for (I - A) x = b started from x_0
# (Walker and Ni, Theorem 2.2).
# ---------------------------------------------------------------------------

const ANDERSON_CONDITION_LIMIT = 1.0e10

# The history of one interface: per relaxed field, its iterates x_j and
# residuals f_j = g_j - x_j, restricted to the partner degrees of freedom that
# the transfer of the interface reads (anderson_footprint); and the traces t_j
# of the residual, stored once for all fields. A field can enter the history
# later than the first (the window-end reaction), so each field has its own
# number of columns; the traces are those of the columns of the first field.
mutable struct AndersonHistory
    X::Vector{Vector{Vector{Float64}}}   # per field, the iterates x_j
    F::Vector{Vector{Vector{Float64}}}   # per field, the residuals f_j
    T::Vector{Vector{Float64}}           # the residual traces t_j
end

AndersonHistory() = AndersonHistory(Vector{Vector{Float64}}[], Vector{Vector{Float64}}[], Vector{Float64}[])

# Per controller: one history per interface, and the data computed at the
# current Schwarz iteration, per interface. Cleared at every stop with the
# rest of the relaxation state.
mutable struct AndersonState
    histories::Dict{RelaxationKey,AndersonHistory}
    iteration::Int
    data::Dict{RelaxationKey,Vector{Vector{Float64}}}
end

function AndersonState()
    return AndersonState(Dict{RelaxationKey,AndersonHistory}(), -1, Dict{RelaxationKey,Vector{Vector{Float64}}}())
end

const ANDERSON_STATE = IdDict{Any,AndersonState}()

# Memory held by the Anderson histories of a controller: the bytes of the
# Float64 data of the iterates X, residuals F, and traces T, over every
# interface and relaxed field, without the headers of the Julia arrays. X and
# F hold the partner data on the footprint of the interface and T its trace,
# once per interface. The peak over a stop is logged at its end: columns
# dropped by the condition limit make the history shrink within a stop.
const ANDERSON_PEAK_BYTES = IdDict{Any,Int}()

function anderson_history_bytes(controller)
    bytes = 0
    groups = Any[]
    haskey(ANDERSON_STATE, controller) && push!(groups, values(ANDERSON_STATE[controller].histories))
    haskey(OVERLAP_ANDERSON_HISTORIES, controller) && push!(groups, values(OVERLAP_ANDERSON_HISTORIES[controller]))
    for histories in groups, history in histories
        for columns in Iterators.flatten((history.X, history.F, (history.T,))), column in columns
            bytes += sizeof(column)
        end
    end
    return bytes
end

function record_anderson_memory!(controller)
    bytes = anderson_history_bytes(controller)
    ANDERSON_PEAK_BYTES[controller] = max(get(ANDERSON_PEAK_BYTES, controller, 0), bytes)
    return bytes
end

# Logs the peak of the stop that ended, records it in anderson_peak_bytes for
# that stop, and clears it. Called after the Schwarz loop, since the state
# itself can be cleared within the stop.
function log_anderson_memory_peak!(controller)
    peak = get(ANDERSON_PEAK_BYTES, controller, 0)
    if peak > 0
        norma_logf(0, :schwarz, "Anderson history memory peak: %.1f KB", peak / 1024)
        controller.anderson_peak_bytes[controller.stop] = peak
    end
    delete!(ANDERSON_PEAK_BYTES, controller)
    return nothing
end

function reset_anderson_state!(controller)
    delete!(ANDERSON_STATE, controller)
    delete!(OVERLAP_ANDERSON_HISTORIES, controller)
    delete!(OVERLAP_WINDOW_STEPS, controller)
    delete!(ANDERSON_FOOTPRINTS, controller)
    return nothing
end

# Keep at most the last n columns of every field and of the traces.
function keep_last!(history::AndersonHistory, n::Int)
    for columns in Iterators.flatten((history.X, history.F, (history.T,)))
        length(columns) > n && deleteat!(columns, 1:(length(columns) - n))
    end
    return nothing
end

# Append the column of the current iterates xs, data gs, and trace t; a field
# seen for the first time gets its own lists.
function push_column!(history::AndersonHistory, xs, gs, t::AbstractVector{Float64})
    while length(history.X) < length(xs)
        push!(history.X, Vector{Float64}[])
        push!(history.F, Vector{Float64}[])
    end
    for (field, (x, g)) in enumerate(zip(xs, gs))
        push!(history.X[field], copy(x))
        push!(history.F[field], g .- x)
    end
    push!(history.T, copy(t))
    return history
end

# The degrees of freedom of the relaxed field of the interface of `bc` that
# the exchange reads: the three components of the partner nodes of
# partner_node_map; for an overlap matching interface, of the partner nodes
# with a nonzero column of its projector; and for a Robin interface, of the
# nodes of its own side set, on which alone its datum is nonzero. The relaxed
# fields are full fields, but the exchange depends only on these entries, and
# the Anderson coefficients come from the traces, so the histories hold the
# iterates and residuals on these entries alone. Nothing (all entries) for a
# reduced-order partner, whose relaxed fields are its reduced coordinates.
# Cached per controller and cleared at every stop.
const ANDERSON_FOOTPRINTS = IdDict{Any,Dict{RelaxationKey,Union{Vector{Int},Nothing}}}()

function anderson_footprint(bc::SolidMechanicsSchwarzBoundaryCondition)
    controller = bc.parent.controller
    footprints = get!(Dict{RelaxationKey,Union{Vector{Int},Nothing}}, ANDERSON_FOOTPRINTS, controller)
    return get!(footprints, relaxation_key(bc)) do
        nodes = if bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
            bc.global_from_local_map
        elseif coupled_subsim_of(bc).model isa RomModel
            return nothing
        elseif bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
            P = bc.dirichlet_projector
            [node for node in partner_node_map(bc) if any(!iszero, view(P, :, node))]
        else
            partner_node_map(bc)
        end
        return vec([3 * node - c for c in (2, 1, 0), node in nodes])
    end
end

restrict(v::AbstractVector{Float64}, footprint::Nothing) = v
restrict(v::AbstractVector{Float64}, footprint::Vector{Int}) = v[footprint]

# The full field whose footprint entries are the relaxed values and whose other
# entries, which the exchange does not read, are those of `base`.
expand(relaxed::Vector{Float64}, base::AbstractVector{Float64}, footprint::Nothing) = copy(relaxed)
function expand(relaxed::Vector{Float64}, base::AbstractVector{Float64}, footprint::Vector{Int})
    full = copy(base)
    full[footprint] = relaxed
    return full
end

# One Anderson step for a group of interfaces. Each entry holds the history of
# one interface, the current iterates xs and the data gs of its relaxed
# fields, restricted to its footprint, and the trace t of the residual. All
# fields of an interface share its trace, and all interfaces share the
# coefficients. Returns the next iterates, per entry and field, on the
# footprint. The histories are trimmed to `depth` columns before the new
# column is appended, so that none holds more than depth + 1.
function anderson_stacked_step!(entries, β::Float64, depth::Int)
    for entry in entries
        keep_last!(entry.history, depth)
        push_column!(entry.history, entry.xs, entry.gs, entry.t)
    end
    # An interface that entered the group later has fewer columns; the
    # stacked differences use the columns that all interfaces have.
    n = minimum(length(entry.history.T) for entry in entries)
    for entry in entries
        keep_last!(entry.history, n)
    end
    γ = Float64[]
    while n >= 2
        ΔT = reduce(
            vcat, [reduce(hcat, [h.T[j + 1] .- h.T[j] for j in 1:(n - 1)]) for h in getfield.(entries, :history)]
        )
        R = qr(ΔT).R
        if any(iszero, diag(R)) || cond(UpperTriangular(R)) > ANDERSON_CONDITION_LIMIT
            n -= 1
            for entry in entries
                keep_last!(entry.history, n)
            end
            continue
        end
        γ = ΔT \ reduce(vcat, [entry.t for entry in entries])
        break
    end
    return map(entries) do entry
        h = entry.history
        map(eachindex(entry.xs)) do field
            X, F = h.X[field], h.F[field]
            x_new = X[end] .+ β .* F[end]
            for j in eachindex(γ)
                x_new .-= γ[j] .* ((X[j + 1] .- X[j]) .+ β .* (F[j + 1] .- F[j]))
            end
            return x_new
        end
    end
end

# Anderson step of a single datum with its own history (coupled initial
# acceleration of matching pairs).
function anderson_step!(
    history::AndersonHistory,
    x::AbstractVector{Float64},
    g::AbstractVector{Float64},
    t::AbstractVector{Float64},
    β::Float64,
    depth::Int,
)
    entry = (history=history, xs=[x], gs=[g], t=t)
    return anderson_stacked_step!([entry], β, depth)[1][1]
end

# The relaxed data of one interface for the stacked step: the relaxed quantity,
# the current iterates and returned data of its relaxed fields, and the trace
# of the residual. Nothing when the interface is not relaxed by Anderson
# acceleration at this iteration: not its relaxed side, an inactive contact
# pair, a stop of several time slots, or no previous iterate yet.
# Scale of the window-end reaction residual in the joint Anderson step of a
# window-end-constrained side, set at its first step and kept for the run.
const WINDOW_REACTION_SCALE = Dict{RelaxationKey,Float64}()

# Anderson histories of the overlap matching interfaces of a controller, keyed
# as the relaxation state is, with the datum as field 1 and the window-end
# reaction as field 2, and the iteration at which a window-end-constrained
# interface took its joint step; both cleared with the Anderson state at every
# stop (reset_anderson_state!).
const OVERLAP_ANDERSON_HISTORIES = IdDict{Any,Dict{RelaxationKey,AndersonHistory}}()
const OVERLAP_WINDOW_STEPS = IdDict{Any,Dict{RelaxationKey,Int}}()
function overlap_anderson_histories(controller)
    return get!(Dict{RelaxationKey,AndersonHistory}, OVERLAP_ANDERSON_HISTORIES, controller)
end
overlap_window_steps(controller) = get!(Dict{RelaxationKey,Int}, OVERLAP_WINDOW_STEPS, controller)

# Joint Anderson step of a window-end-constrained overlap matching side at its
# first interior substep of an iteration: on the datum of the window end and
# on the window-end reaction that forms the ramp, with the reaction residual
# scaled once per interface to the datum residual. The datum goes into the
# slot of the window end, which the application at the end substep reads.
function overlap_matching_window_step!(bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    controller = bc.parent.controller
    key = relaxation_key(bc)
    iter = controller.iteration_number
    get(overlap_window_steps(controller), key, -1) == iter && return nothing
    isempty(bc.window_end_reaction) && return nothing
    length(bc.window_end_reaction_raw) == length(bc.window_end_reaction) || return nothing
    t_end = controller.time
    slot_k = relaxation_slot!(controller, key, t_end)
    relax_velocity = bc.constraint == :velocity
    c_slots = relaxation_slots!(relax_velocity ? controller.lambda_velo : controller.lambda_disp, key)
    ensure_slot!(c_slots, slot_k)
    isempty(c_slots[slot_k]) && return nothing
    disp, velo, _, _ = partner_fields(bc, t_end)
    interp_c = relax_velocity ? velo : disp
    λ_c_prev = c_slots[slot_k]
    t = constrained_partner_trace(bc, interp_c) .- constrained_partner_trace(bc, λ_c_prev)
    x_λ = bc.window_end_reaction
    g_λ = bc.window_end_reaction_raw
    t_λ = g_λ .- x_λ
    scale = get!(WINDOW_REACTION_SCALE, key) do
        n_t, n_λ = norm(t), norm(t_λ)
        return n_t > 0.0 && n_λ > 0.0 ? n_t / n_λ : 1.0
    end
    footprint = anderson_footprint(bc)
    entry = (
        history=get!(AndersonHistory, overlap_anderson_histories(controller), key),
        xs=[restrict(λ_c_prev, footprint), x_λ],
        gs=[restrict(interp_c, footprint), g_λ],
        t=vcat(t, scale .* t_λ),
    )
    data = anderson_stacked_step!([entry], controller.relaxation_parameter, controller.anderson_depth)[1]
    c_slots[slot_k] = expand(data[1], interp_c, footprint)
    bc.window_end_reaction = copy(data[2])
    overlap_window_steps(controller)[key] = iter
    record_anderson_memory!(controller)
    return nothing
end

function anderson_entry(bc::SolidMechanicsSchwarzBoundaryCondition, time::Float64)
    controller = bc.parent.controller
    key = relaxation_key(bc)
    bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition && return nothing
    if bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
        bc.self_handle.id < bc.coupled_handle.id || return nothing
    elseif is_swappable_dn_schwarz(bc)
        bc.is_dirichlet || return nothing
        bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact && return nothing
    else
        return nothing
    end
    haskey(controller.lambda_time, key) && aitken_applies(controller, key) || return nothing
    disp_slots = relaxation_slots!(controller.lambda_disp, key)
    velo_slots = relaxation_slots!(controller.lambda_velo, key)
    acce_slots = relaxation_slots!(controller.lambda_acce, key)
    if bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
        (isempty(disp_slots) || isempty(disp_slots[1])) && return nothing
        x = disp_slots[1]
        g = with_partner_state(() -> robin_datum(self_subsim_of(bc).model, bc), bc, time)
        # The datum is a field of this side, nonzero on its interface degrees
        # of freedom only, and its residual is its own trace.
        footprint = anderson_footprint(bc)
        x_Γ, g_Γ = restrict(x, footprint), restrict(g, footprint)
        return (quantity=:robin, key=key, xs=[x_Γ], gs=[g_Γ], t=g_Γ .- x_Γ)
    end
    # The datum of a single-slot interface is that of the window end, whatever
    # the time of the application that triggered the step.
    disp, velo, acce, _ = partner_fields(bc, controller.time)
    if is_constrained_dn(bc)
        relax_velocity = bc.constraint == :velocity
        slots = relax_velocity ? velo_slots : disp_slots
        (isempty(slots) || isempty(slots[1])) && return nothing
        x = slots[1]
        g = relax_velocity ? velo : disp
        t = constrained_partner_trace(bc, g) .- constrained_partner_trace(bc, x)
        quantity = relax_velocity ? :velocity : :displacement
        if bc.window_end_constraint &&
            !isempty(bc.window_end_reaction) &&
            length(bc.window_end_reaction_raw) == length(bc.window_end_reaction)
            # Joint step on the datum and the window-end reaction: the residual
            # of the reaction is scaled once per interface to the size of the
            # datum's residual, so that both enter the least squares.
            x_λ = bc.window_end_reaction
            g_λ = bc.window_end_reaction_raw
            t_λ = g_λ .- x_λ
            scale = get!(WINDOW_REACTION_SCALE, key) do
                n_t, n_λ = norm(t), norm(t_λ)
                return n_t > 0.0 && n_λ > 0.0 ? n_t / n_λ : 1.0
            end
            footprint = anderson_footprint(bc)
            return (
                quantity=quantity,
                key=key,
                xs=[restrict(x, footprint), x_λ],
                gs=[restrict(g, footprint), g_λ],
                t=vcat(t, scale .* t_λ),
            )
        end
        footprint = anderson_footprint(bc)
        return (quantity=quantity, key=key, xs=[restrict(x, footprint)], gs=[restrict(g, footprint)], t=t)
    end
    (isempty(disp_slots) || isempty(disp_slots[1])) && return nothing
    t = constrained_partner_trace(bc, disp) .- constrained_partner_trace(bc, disp_slots[1])
    footprint = anderson_footprint(bc)
    # A static partner does not change its velocity and acceleration during the
    # Schwarz iteration, so they are passed through unrelaxed (apply_bc), and
    # the displacement is the only relaxed field.
    fields = relaxes_displacement_only(bc) ? (1:1) : (1:3)
    xs = [restrict(slots[1], footprint) for slots in (disp_slots, velo_slots, acce_slots)[fields]]
    gs = [restrict(field, footprint) for field in (disp, velo, acce)[fields]]
    return (quantity=:displacement, key=key, xs=xs, gs=gs, t=t)
end

# Whether the Anderson step of an unconstrained Dirichlet-Neumann or contact
# interface relaxes the partner displacement alone: when the partner's time
# integrator is static, whose velocity and acceleration stay constant.
function relaxes_displacement_only(bc::SolidMechanicsSchwarzBoundaryCondition)
    return coupled_subsim_of(bc).integrator isa StaticTimeIntegrator
end

# Next data of the interface of `bc` under Anderson acceleration, one vector
# per relaxed field. The first call of a Schwarz iteration computes the
# stacked step for every relaxed interface of the simulation.
function anderson_data!(bc::SolidMechanicsSchwarzBoundaryCondition, time::Float64)
    controller = bc.parent.controller
    state = get!(AndersonState, ANDERSON_STATE, controller)
    if state.iteration != controller.iteration_number
        empty!(state.data)
        entries = []
        for subsim in bc.parent.subsims, other in subsim.model.boundary_conditions
            other isa SolidMechanicsSchwarzBoundaryCondition || continue
            entry = anderson_entry(other, time)
            entry === nothing || push!(entries, entry)
        end
        for quantity in unique(getfield.(entries, :quantity))
            group = [
                merge(e, (history=get!(AndersonHistory, state.histories, e.key),)) for
                e in entries if e.quantity === quantity
            ]
            for (entry, data) in
                zip(group, anderson_stacked_step!(group, controller.relaxation_parameter, controller.anderson_depth))
                state.data[entry.key] = data
            end
            norma_logf(
                1,
                :schwarz,
                "Anderson [iter=%d] %s: %d interfaces, %d columns",
                controller.iteration_number,
                quantity,
                length(group),
                length(group[1].history.T) - 1,
            )
        end
        norma_logf(
            1,
            :schwarz,
            "Anderson [iter=%d] history memory: %.1f KB",
            controller.iteration_number,
            record_anderson_memory!(controller) / 1024,
        )
        state.iteration = controller.iteration_number
    end
    key = relaxation_key(bc)
    haskey(state.data, key) || norma_abort("Anderson acceleration: no stacked data for interface $(key).")
    return state.data[key]
end

# Recover interface velocity and acceleration consistent with a relaxed interface
# displacement, instead of relaxing them independently. Norma's implicit Newmark
# is displacement-form (the solver unknown is u; see `correct`, time_integrator.jl),
# so the single interface unknown to relax is the displacement, and v, a follow
# from the integrator's own recovery relations using the stored predictors:
#   a = (u - u_pre) / (β Δt²),   v = v_pre + γ Δt a.
# For integrators without these predictors (e.g. quasi-static, explicit), the
# interpolated kinematics are passed through unrelaxed.
function recover_interface_kinematics!(controller, key, slot_k, integrator, interp_velo, interp_acce)
    velo_slots = relaxation_slots!(controller.lambda_velo, key)
    acce_slots = relaxation_slots!(controller.lambda_acce, key)
    if integrator isa Newmark
        Δt = integrator.time_step
        β = integrator.β
        γ = integrator.γ
        u = relaxation_slots!(controller.lambda_disp, key)[slot_k]
        a = (u .- integrator.disp_pre) ./ (β * Δt * Δt)
        acce_slots[slot_k] = a
        velo_slots[slot_k] = integrator.velo_pre .+ (γ * Δt) .* a
    else
        velo_slots[slot_k] = interp_velo
        acce_slots[slot_k] = interp_acce
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Robin-Robin nonoverlap Schwarz: t + α W u = g
# ---------------------------------------------------------------------------
#
# Each side receives the datum g = -t_src + α W u_src, where t_src is the
# interface reaction of the partner transferred by the Neumann projector, u_src
# the interface displacement of the partner transferred by the Dirichlet
# projector, W the boundary mass matrix of this side, and α its Robin
# parameter. The self term α W u is added to the internal force
# (build_robin_schwarz_force) and its tangent α W to the stiffness
# (build_robin_schwarz_stiffness). The side listed first in `domains` relaxes
# its datum: in the Gauss-Seidel pass over the subdomains it is solved before
# its partner and so receives the partner's datum of the previous iteration,
# the datum that carries the state from one Schwarz iteration to the next.
# The side listed later receives a datum computed in the same pass. Anderson
# acceleration needs the iteration to be a function of the accelerated data
# alone, which holds only when the lagged data are the accelerated ones; with
# the side listed later relaxed it stalled on the five-subdomain inclusion of
# issue #233 (issue #245). This is the arrangement of a Dirichlet-Neumann pair
# whose Dirichlet side is listed first.
# Unrelaxed Robin datum -t_src + α W u_src of this side, on its degrees of
# freedom, from the partner state currently set on the partner integrator.
function robin_datum(model::SolidMechanics, bc::SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition)
    α = bc.robin_parameter
    W = bc.square_projector
    # Neumann part: -t_src projected
    neumann_force = get_dst_force(bc)

    # Source displacement, projected to destination
    src_sim = coupled_subsim_of(bc)
    src_model = get_fom_model(src_sim)
    src_bc = src_model.boundary_conditions[bc.coupled_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
    num_src_nodes = length(src_global_from_local_map)
    src_disp = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_disp[:, i_local] = src_model.displacement[:, i_global]
    end
    dst_disp = project_components(bc.dirichlet_projector, src_disp)
    global_from_local_map = bc.global_from_local_map

    # Datum (partner contribution only), on the interface degrees of freedom:
    # rhs = -t_src + α W u_src. The self term α W u_self is added by
    # build_robin_schwarz_force at each evaluation. Other loads already in
    # boundary_force (Neumann and pressure conditions on interface nodes) are
    # not part of the datum; relaxing them would scale them by 1/θ at the
    # fixed point.
    rhs = zeros(length(model.boundary_force))
    for comp in 1:3
        α_W_u = α * (W * dst_disp[comp, :])
        for (i_local, i_global) in enumerate(global_from_local_map)
            dof_i = 3 * (i_global - 1) + comp
            rhs[dof_i] = neumann_force[3 * (i_local - 1) + comp] + α_W_u[i_local]
        end
    end
    return rhs
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition)
    parent_sim = bc.parent
    controller = parent_sim.controller
    iter = controller.iteration_number
    coupled_index = bc.coupled_handle.id
    this_index = bc.self_handle.id
    global_from_local_map = bc.global_from_local_map
    rhs = robin_datum(model, bc)
    if (this_index > coupled_index)  # side listed later: not relaxed
        datum = rhs
    else  # side listed first: relaxed
        # The relaxation state is per substep time slot (see relaxation_slot!):
        # relaxing against the datum of the previous Schwarz iteration at the
        # same time keeps a windowed exchange consistent in time. A fresh slot
        # (every slot on iteration 0, since the state is reset at each stop) has
        # no previous iterate and starts from zero.
        key = relaxation_key(bc)
        slot_k = relaxation_slot!(controller, key, model.time)
        g_slots = relaxation_slots!(controller.lambda_disp, key)
        ensure_slot!(g_slots, slot_k)
        g_stored = g_slots[slot_k]
        g = isempty(g_stored) ? zeros(length(model.boundary_force)) : g_stored
        # Optional Aitken relaxation factor. The fixed-point iterate is the datum
        # stored in lambda_disp; its unrelaxed (theta = 1) candidate is rhs, from
        # which the residual r = rhs - g is formed.
        if controller.relaxation_method === :anderson && aitken_applies(controller, key) && !isempty(g_stored)
            # Anderson acceleration of the Robin datum, stacked with the other
            # relaxed Robin interfaces (anderson_data!); the datum lives on the
            # interface degrees of freedom only, so it is its own trace, and the
            # step is computed on those degrees of freedom (anderson_footprint).
            # A fresh slot (g = 0) is not entered into the history.
            datum = expand(anderson_data!(bc, model.time)[1], rhs, anderson_footprint(bc))
        else
            theta = controller.relaxation_parameter
            if controller.relaxation_method === :aitken_recursive || controller.relaxation_method === :aitken_secant
                theta = relaxation_theta!(bc, key, slot_k, iter, rhs, g)
            end
            frozen_relaxation_update!(controller, theta)
            datum = (1 - theta) .* g .+ theta .* rhs
        end
        g_slots[slot_k] = datum
    end
    for i_global in global_from_local_map, comp in 1:3
        dof_i = 3 * (i_global - 1) + comp
        model.boundary_force[dof_i] += datum[dof_i]
    end
end

# Tangent contribution α W of the Robin self term, for the implicit
# integrators. Explicit integrators assemble no tangent.
function build_robin_schwarz_stiffness(model::SolidMechanics)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_rs = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition || continue
        α = bc.robin_parameter
        W = bc.square_projector
        global_from_local_map = bc.global_from_local_map
        for (i_local, i_global) in enumerate(global_from_local_map)
            for (j_local, j_global) in enumerate(global_from_local_map)
                w_ij = α * W[i_local, j_local]
                for comp in 1:3
                    dof_i = 3 * (i_global - 1) + comp
                    dof_j = 3 * (j_global - 1) + comp
                    K_rs[dof_i, dof_j] += w_ij
                end
            end
        end
    end
    return K_rs
end

# Robin self term W α u_self, returned as a separate force vector rather than
# added to model.internal_force, because get_dst_force reads
# model.internal_force for the traction transfer and must see only the elastic
# internal force.
function build_robin_schwarz_force(model::SolidMechanics)
    num_dofs = 3 * size(model.reference, 2)
    f = zeros(num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition || continue
        α = bc.robin_parameter
        W = bc.square_projector
        global_from_local_map = bc.global_from_local_map
        num_nodes = length(global_from_local_map)
        self_disp = zeros(num_nodes)
        for comp in 1:3
            for (i_local, i_global) in enumerate(global_from_local_map)
                self_disp[i_local] = model.displacement[comp, i_global]
            end
            f_rs = W * (α * self_disp)
            for (i_local, i_global) in enumerate(global_from_local_map)
                dof_i = 3 * (i_global - 1) + comp
                f[dof_i] += f_rs[i_local]
            end
        end
    end
    return f
end

function pair_bc(bc::SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition, bc_index::Int64)
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = coupled_subsim_of(bc).model
    coupled_bcs = coupled_model.boundary_conditions
    # The partner is the Robin condition on the named side set; other
    # conditions on the same side set (a Neumann load, for example) share its
    # name and are not partners.
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name && coupled_bc isa SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
        end
    end
    return nothing
end

# Per-side transfer operators of a Robin-Robin interface: each side builds its
# own Dirichlet and Neumann projectors and its boundary mass matrix W.
function compute_robin_schwarz_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition
)
    cache = RectangularProjectionCache()
    compute_dirichlet_projector(dst_model, dst_bc; cache=cache)
    compute_neumann_projector(dst_model, dst_bc; cache=cache)
    dst_bc.square_projector = get_square_projection_matrix(dst_model, dst_bc)
    return nothing
end

# Store the boundary mass matrix W of one side of a Dirichlet-Neumann pair
# together with its Cholesky factorization. W is symmetric positive definite and
# does not change during a run. The interface force residual evaluates the W⁻¹
# norm at every Schwarz iteration and at every stop of the energy output; before
# the factorization was stored, each evaluation factored the dense matrix anew,
# six times per residual.
function set_square_projector!(bc::SolidMechanicsConstrainedSchwarzBoundaryCondition, W::Matrix{Float64})
    bc.square_projector = W
    bc.square_factor = cholesky(Symmetric(W))
    return nothing
end

# Transfer operators of a constrained Dirichlet-Neumann pair, built from one
# cross mass matrix B_mn = ∫_Γ φ¹_m φ²_n dS, with φ¹ and φ² the trace shape
# functions of the two sides: Π₁ = W₁⁻¹ B, Π₂ = W₂⁻¹ Bᵀ, and each side's force transfer is the
# transpose of the partner's kinematic transfer, N₁ = Π₂ᵀ and N₂ = Π₁ᵀ. When
# side 1 is the Dirichlet side, the Neumann side receives N₂ = Π_Dᵀ, so the
# work of the transferred reaction on the Neumann velocity equals the work of
# the reaction on the projected velocity imposed on the Dirichlet side. W_k is
# stored as the square projector of side k and is the same matrix in the
# constraint, in the projector, and in the force transfer. B is assembled
# with the facet quadrature of either side. Integrated on the facets of the
# Neumann side, it is exact when those facets each lie within a facet of the
# Dirichlet side (nested facets, the Neumann side the finer), and Π_D then
# reproduces constants to roundoff; that integration is used whenever Π_D
# reproduces constants to NESTED_PARTITION_TOLERANCE with it. Otherwise B is
# integrated over the facets of the Dirichlet side with the same quadrature as
# W_D, so that its row sums equal those of W_D exactly: Π_D then reproduces
# constants to roundoff and the force Π_Dᵀλ transferred to the Neumann side
# has the total of λ. On facets that are not nested, B integrated on the
# Neumann side carries the quadrature error of the product of the two trace
# bases into Π_D (1.2e-2 on the beam at mesh ratio 1:0.75 and 0.72 on the
# cylinders). On nested facets with the coarser mesh as the Dirichlet side,
# B integrated on the Dirichlet side is not exact: Π_D then differs from the
# L2 projection, and on the beam at 1:0.5 with 4:1 steps (clamped part, coarse
# in mesh and step, the Dirichlet side) the trapezoidal pair needed 61.9
# Schwarz iterations per window at a fixed factor of 0.5 against 44.5 with the
# exact B. The projection onto the Neumann side, W_N⁻¹Bᵀ, is never applied;
# its partition-of-unity error is logged for information only.
const NESTED_PARTITION_TOLERANCE = 1.0e-12

# Transfer operators of an overlap matching pair: Π_k = W_k⁻¹ L_k on each
# Schwarz boundary, with L_k the cross-mass matrix of Γ_k against the partner's
# volume shape functions and W_k the facet mass matrix of Γ_k by the same
# segmented quadrature (get_overlap_rectangular_projection_matrix), and each
# side's force transfer the transpose of the partner's projector. Both sides
# are built at the first call.
function compute_constrained_dn_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
)
    size(dst_bc.dirichlet_projector, 1) > 0 && return nothing
    dst_sim = self_subsim_of(dst_bc)
    src_sim = coupled_subsim_of(dst_bc)
    if !(dst_sim.model isa SolidMechanics) || !(src_sim.model isa SolidMechanics)
        norma_abort("`Schwarz overlap matching` requires full order (solid mechanics) models on both sides.")
    end
    src_model = src_sim.model
    src_bc = coupled_bc_of(dst_bc)
    constraint = resolve_constraint(dst_bc, src_bc)
    check_constrained_integrators(dst_sim, src_sim, constraint)
    check_direct_interface_solve(dst_bc, src_bc, dst_sim, src_sim)
    if dst_bc.direct_solve && (dst_bc.window_end_constraint || src_bc.window_end_constraint)
        norma_abort("`interface solve: direct` and `matching at: window end` exclude each other.")
    end
    Δt_1 = dst_sim.integrator.time_step
    Δt_2 = src_sim.integrator.time_step
    if !isapprox(Δt_1, Δt_2; rtol=1.0e-12)
        fine_bc = Δt_1 < Δt_2 ? dst_bc : src_bc
        fine_sim = Δt_1 < Δt_2 ? dst_sim : src_sim
        fine_bc.window_end_constraint || norma_abort(
            "`Schwarz overlap matching` with different time steps requires `matching at: window end` on the " *
            "fine-step side ('$(fine_sim.name)', side set '$(fine_bc.name)'): both sides are constrained, " *
            "and a side that takes several steps per window matches at the window end and carries the ramp " *
            "of its own reaction at its interior substeps.",
        )
    end
    for (model, bc) in ((dst_model, dst_bc), (src_model, src_bc))
        partner = get_fom_model(coupled_subsim_of(bc))
        L, W = get_overlap_rectangular_projection_matrix(model, bc, partner, bc.coupled_block_name, bc.search_tolerance)
        set_square_projector!(bc, W)
        bc.dirichlet_projector = bc.square_factor \ L
    end
    dst_bc.neumann_projector = Matrix(transpose(src_bc.dirichlet_projector))
    src_bc.neumann_projector = Matrix(transpose(dst_bc.dirichlet_projector))
    for bc in (dst_bc, src_bc)
        P = bc.dirichlet_projector
        norma_logf(
            0,
            :info,
            "Overlap matching interface '%s' of %s: %s matching, %d constrained nodes against %d partner nodes, " *
            "partition-of-unity error of Π %.2e.",
            bc.name,
            self_subsim_of(bc).name,
            String(constraint),
            size(P, 1),
            size(P, 2),
            maximum(abs.(P * ones(size(P, 2)) .- 1.0)),
        )
    end
    return nothing
end

function compute_constrained_dn_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition
)
    # The pair is processed once; the partner's pass finds its operators set.
    if size(dst_bc.dirichlet_projector, 1) > 0
        return nothing
    end
    dst_sim = self_subsim_of(dst_bc)
    src_sim = coupled_subsim_of(dst_bc)
    if !(dst_sim.model isa SolidMechanics) || !(src_sim.model isa SolidMechanics)
        norma_abort("`matching` requires full order (solid mechanics) models on both sides of the pair.")
    end
    src_model = src_sim.model
    src_bc = coupled_bc_of(dst_bc)
    if !(src_bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition) || !src_bc.constrained
        norma_abort(
            "`matching` must be set on BOTH sides of a Schwarz DN nonoverlap pair " *
            "(missing on the side coupled to '$(dst_bc.name)').",
        )
    end
    constraint = resolve_constraint(dst_bc, src_bc)
    dst_bc.constraint = src_bc.constraint = constraint
    check_constrained_integrators(dst_sim, src_sim, constraint)
    check_direct_interface_solve(dst_bc, src_bc, dst_sim, src_sim)
    # B in the orientation (rows of side 1, columns of side 2), integrated over
    # the facets of side 1 or of side 2 (the first argument of
    # get_rectangular_projection_matrix is the side whose facets are integrated),
    # each with the boundary mass matrix of the integrated side by the same
    # quadrature, so that the row sums of B_on_1 equal those of W1 and the
    # column sums of B_on_2 those of W2 exactly.
    B_on_1, W1 = get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    B_on_2_transposed, W2 = get_rectangular_projection_matrix(src_model, src_bc, dst_model, dst_bc)
    B_on_2 = Matrix(transpose(B_on_2_transposed))
    n1 = size(W1, 1)
    n2 = size(W2, 1)
    B_on_D, B_on_N = dst_bc.is_dirichlet ? (B_on_1, B_on_2) : (B_on_2, B_on_1)
    set_square_projector!(dst_bc, W1)
    set_square_projector!(src_bc, W2)
    # Partition-of-unity error of Π_D with B integrated on the Neumann side.
    pu_error_D_on_N = if dst_bc.is_dirichlet
        maximum(abs.(dst_bc.square_factor \ (B_on_N * ones(n2)) .- 1.0))
    else
        maximum(abs.(src_bc.square_factor \ (transpose(B_on_N) * ones(n1)) .- 1.0))
    end
    nested = pu_error_D_on_N <= NESTED_PARTITION_TOLERANCE
    B = nested ? B_on_N : B_on_D
    if dst_bc.is_dirichlet
        warn_dirichlet_role(dst_sim, src_sim, dst_bc.name, src_bc.name, n1, n2; window_end=dst_bc.window_end_constraint)
    else
        warn_dirichlet_role(src_sim, dst_sim, src_bc.name, dst_bc.name, n2, n1; window_end=src_bc.window_end_constraint)
    end
    P1 = dst_bc.square_factor \ B
    P2 = src_bc.square_factor \ Matrix(transpose(B))
    dst_bc.dirichlet_projector = P1
    dst_bc.neumann_projector = Matrix(transpose(P2))
    src_bc.dirichlet_projector = P2
    src_bc.neumann_projector = Matrix(transpose(P1))
    pu_error_1 = maximum(abs.(P1 * ones(n2) .- 1.0))
    pu_error_2 = maximum(abs.(P2 * ones(n1) .- 1.0))
    pu_error_D, pu_error_N = dst_bc.is_dirichlet ? (pu_error_1, pu_error_2) : (pu_error_2, pu_error_1)
    norma_logf(
        0,
        :info,
        "Matching DN interface '%s'/'%s': %s matching, B integrated on the facets of the %s side, " *
        "partition-of-unity error of Π_D %.2e (unused projection onto N: %.2e).",
        dst_bc.name,
        src_bc.name,
        String(constraint),
        nested ? "Neumann" : "Dirichlet",
        pu_error_D,
        pu_error_N,
    )
    return nothing
end

# Warns when the Dirichlet role of a matching pair is assigned against the
# role rules (coupling note, assignment of the Dirichlet role), in their order
# of precedence: with different time steps the coarse-step side is the
# Dirichlet side under the constraint at every substep (the impulse residual
# vanishes only then; the window-end constraint balances it with either side as
# the Dirichlet side, so no warning is given there); with different
# integrators the implicit side is the Dirichlet side (with the explicit side
# as the Dirichlet side the gain of the iteration reaches -13 on the beam and
# no fixed factor converges); with the same integrator and step the finer
# mesh is the Dirichlet side (the gain then lies in [-1, 0], and grows with
# the cube of the mesh ratio otherwise). The node count stands in for the
# mesh size on the interface. The roles are not changed: they are the input's.
function warn_dirichlet_role(
    sim_D::SingleDomainSimulation,
    sim_N::SingleDomainSimulation,
    name_D::String,
    name_N::String,
    n_D::Int,
    n_N::Int;
    window_end::Bool=false,
)
    Δt_D = sim_D.integrator.time_step
    Δt_N = sim_N.integrator.time_step
    pair = "Matching DN interface '$name_D'/'$name_N' ($(sim_D.name) is the Dirichlet side)"
    if !isapprox(Δt_D, Δt_N; rtol=1.0e-12)
        if Δt_D < Δt_N && !window_end
            norma_logf(
                0,
                :warn,
                "%s: the Dirichlet side has the finer time step (%.3e against %.3e). The coarse-step side " *
                "should be the Dirichlet side: with the fine-step side as the Dirichlet side the rate invariant " *
                "is conserved but the interface impulse is not balanced, and the energy varies by the order of " *
                "the impulse residual.",
                pair,
                Δt_D,
                Δt_N,
            )
        end
        return nothing
    end
    explicit_D = sim_D.integrator isa CentralDifference
    explicit_N = sim_N.integrator isa CentralDifference
    if explicit_D != explicit_N
        if explicit_D
            norma_logf(
                0,
                :warn,
                "%s: the Dirichlet side is the explicit member of a mixed pair. The implicit side should be " *
                "the Dirichlet side: with the explicit side as the Dirichlet side the gain of the Schwarz " *
                "iteration reaches -13 on the cantilever beam and no fixed relaxation factor converges.",
                pair,
            )
        end
        return nothing
    end
    if n_D < n_N
        norma_logf(
            0,
            :warn,
            "%s: the Dirichlet side has fewer interface nodes (%d against %d). With the same integrator and " *
            "step on both sides the finer mesh should be the Dirichlet side: the gain of the Schwarz iteration " *
            "then lies in [-1, 0], and with the coarser mesh as the Dirichlet side it grows with the cube of " *
            "the mesh ratio.",
            pair,
            n_D,
            n_N,
        )
    end
    return nothing
end

# The matched quantity of a pair, the value of `matching`, which both sides
# give and which must agree.
function resolve_constraint(
    bc::SolidMechanicsConstrainedSchwarzBoundaryCondition, partner::SolidMechanicsConstrainedSchwarzBoundaryCondition
)
    if bc.constraint != partner.constraint
        norma_abort(
            "The sides '$(bc.name)' and '$(partner.name)' of a matching DN pair give different " *
            "values of `matching` ($(bc.constraint) and $(partner.constraint)); a pair has one matched quantity.",
        )
    end
    return bc.constraint
end

# The displacement constraint derives the velocity and the acceleration of the
# Dirichlet side from the Newmark relations a = (u - u_pre)/(βΔt²), which
# requires β > 0, and it reproduces the coupled problem only when both sides
# advance with the same relations and step, so it is restricted to equal steps.
# The velocity constraint also admits different steps: the side with the finer
# step receives the partner's velocity or reaction interpolated linearly in
# time from the partner's substep history (apply_bc). With the coarse-step side
# as the Dirichlet side this is the r = 1 multirate scheme of Connors, Owen,
# Kuberry, and Bochev (2024) and the scheme of Prakash and Hjelmstad (2004),
# whose interface terms cancel and which conserves the rate invariant Ẽ (coupling note).
function check_constrained_integrators(sim_1::SingleDomainSimulation, sim_2::SingleDomainSimulation, constraint::Symbol)
    for sim_k in (sim_1, sim_2)
        integrator = sim_k.integrator
        if !(integrator isa Newmark || integrator isa CentralDifference)
            norma_abort(
                "`matching` requires Newmark or central difference integrators; " *
                "subdomain '$(sim_k.name)' uses $(typeof(integrator)).",
            )
        end
        if integrator isa Newmark && integrator.hht_alpha > 0.0
            norma_abort("`matching` does not support HHT-α (subdomain '$(sim_k.name)').")
        end
    end
    if constraint == :displacement
        Δt_1 = sim_1.integrator.time_step
        Δt_2 = sim_2.integrator.time_step
        if !isapprox(Δt_1, Δt_2; rtol=1.0e-12)
            norma_abort(
                "`matching: displacement` is implemented for equal time steps; subdomains " *
                "'$(sim_1.name)' and '$(sim_2.name)' use $(Δt_1) and $(Δt_2). Use `matching: velocity`.",
            )
        end
        i_1 = sim_1.integrator
        i_2 = sim_2.integrator
        same_newmark = i_1 isa Newmark && i_2 isa Newmark && i_1.β == i_2.β && i_1.γ == i_2.γ && i_1.β > 0.0
        if !same_newmark
            norma_abort(
                "`matching: displacement` requires both members of the pair to be Newmark " *
                "integrators with the same β > 0, γ, and time step; subdomains '$(sim_1.name)' " *
                "and '$(sim_2.name)' do not satisfy this. Use `matching: velocity`.",
            )
        end
    end
    return nothing
end

# The direct interface solve applies to a constrained pair of two central
# difference subdomains at equal steps under the velocity constraint; both
# sides must request it.
function check_direct_interface_solve(
    bc::SolidMechanicsConstrainedSchwarzBoundaryCondition,
    partner::SolidMechanicsConstrainedSchwarzBoundaryCondition,
    sim_1::SingleDomainSimulation,
    sim_2::SingleDomainSimulation,
)
    bc.direct_solve || partner.direct_solve || return nothing
    if bc.direct_solve != partner.direct_solve
        norma_abort("`interface solve: direct` must be set on both sides of the pair '$(bc.name)'/'$(partner.name)'.")
    end
    explicit = sim_1.integrator isa CentralDifference && sim_2.integrator isa CentralDifference
    equal = isapprox(sim_1.integrator.time_step, sim_2.integrator.time_step; rtol=1.0e-12)
    if !explicit || !equal || bc.constraint != :velocity
        norma_abort(
            "`interface solve: direct` is implemented for a matching pair of two central difference " *
            "subdomains at equal time steps under `matching: velocity`; subdomains '$(sim_1.name)' and " *
            "'$(sim_2.name)' do not satisfy this. Use `interface solve: iterative`.",
        )
    end
    return nothing
end

# True for the Dirichlet-Neumann condition with the constrained exchange.
function is_constrained_dn(bc::SolidMechanicsBoundaryCondition)
    return (bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition && bc.constrained) ||
           bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
end

# True for a constrained pair solved by the direct interface solve.
is_direct_dn(bc::SolidMechanicsBoundaryCondition) = is_constrained_dn(bc) && bc.direct_solve

# ---------------------------------------------------------------------------
# Direct interface solve for constrained pairs of two central difference
# subdomains at equal steps.
#
# Let D be the Dirichlet side and N the Neumann side, λ the interface force
# applied on the interface rows of D and -Π_Dᵀ λ the force on the interface
# rows of N (the force transfer of the constrained exchange), and m_D, m_N the
# lumped masses of the interface rows. Central difference computes the
# displacement u_{n+1} before the acceleration, and the internal force at
# u_{n+1} does not depend on λ. Advancing both sides with no interface force
# gives the free velocities v_free; the interface force adds
#   v_D = v_D,free + γΔt m_D⁻¹ λ,   v_N = v_N,free - γΔt m_N⁻¹ Π_Dᵀ λ.
# The velocity constraint v_D = Π_D v_N at the end of the step then gives
#   γΔt (m_D⁻¹ + Π_D m_N⁻¹ Π_Dᵀ) λ = Π_D v_N,free - v_D,free,
# per Cartesian component. The matrix H₀ = m_D⁻¹ + Π_D m_N⁻¹ Π_Dᵀ is symmetric
# positive definite, of the size of the Dirichlet interface, and constant
# while the lumped masses and Π_D are, so it is factored once. Interface rows
# held by other Dirichlet conditions do not move and take a zero inverse mass.
# The solve is exact for nonlinear materials too, because the internal force
# is evaluated at the predictor displacement. The same matrix solves the
# acceleration constraint a_D = Π_D a_N at t = 0:
#   H₀ λ = Π_D a_N,free - a_D,free.
# ---------------------------------------------------------------------------

struct DirectInterfaceFactor
    factors::Vector{Any}              # Cholesky factor of H₀ per component
    d_inverse_mass::Vector{Vector{Float64}}  # interface_inverse_mass of D per component
    n_inverse_mass::Vector{Vector{Float64}}  # interface_inverse_mass of N per component
end

const DIRECT_INTERFACE_CACHE = IdDict{Any,DirectInterfaceFactor}()

# Inverse lumped masses of the interface rows of one side for one component,
# zero where another Dirichlet condition holds the degree of freedom.
function interface_inverse_mass(model::SolidMechanics, map::Vector{Int64}, comp::Int, fixed::BitVector)
    m = model.lumped_mass
    return [fixed[3 * (n - 1) + comp] ? 0.0 : 1.0 / m[3 * (n - 1) + comp] for n in map]
end

# Degrees of freedom held by Dirichlet conditions other than the Schwarz ones.
function prescribed_dofs(model::SolidMechanics)
    fixed = falses(length(model.free_dofs))
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsDirichletBoundaryCondition
            for n in bc.node_set_node_indices
                fixed[3 * (n - 1) + bc.offset] = true
            end
        elseif bc isa SolidMechanicsSideSetDirichletBoundaryCondition
            for n in bc.side_set_node_indices
                fixed[3 * (n - 1) + bc.offset] = true
            end
        end
    end
    return fixed
end

function direct_interface_factor(bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    return get!(DIRECT_INTERFACE_CACHE, bc) do
        d_model = self_subsim_of(bc).model
        n_model = coupled_subsim_of(bc).model
        n_bc = coupled_bc_of(bc)
        P = bc.dirichlet_projector
        d_fixed = prescribed_dofs(d_model)
        n_fixed = prescribed_dofs(n_model)
        t0 = time()
        factors = Any[]
        d_inverse_mass = Vector{Float64}[]
        n_inverse_mass = Vector{Float64}[]
        for comp in 1:3
            d_inv = interface_inverse_mass(d_model, bc.global_from_local_map, comp, d_fixed)
            n_inv = interface_inverse_mass(n_model, n_bc.global_from_local_map, comp, n_fixed)
            push!(d_inverse_mass, d_inv)
            push!(n_inverse_mass, n_inv)
            H = P * (n_inv .* transpose(P))
            for i in eachindex(d_inv)
                H[i, i] += d_inv[i]
            end
            F = cholesky(Symmetric(H); check=false)
            if !issuccess(F)
                norma_abort(
                    "The direct interface matrix of '$(bc.name)' is not positive definite (component $comp): " *
                    "an interface degree of freedom is held by Dirichlet conditions on both sides.",
                )
            end
            push!(factors, F)
        end
        norma_logf(
            0,
            :setup,
            "Direct interface solve '%s': %d interface nodes, three factorizations in %.2f s.",
            bc.name,
            size(P, 1),
            time() - t0,
        )
        return DirectInterfaceFactor(factors, d_inverse_mass, n_inverse_mass)
    end
end

# Apply the interface force λ (3 × n_D) to the two sides: change of the
# acceleration on the interface rows by m⁻¹ times the force, and of the
# velocity by `velocity_factor` times that (γΔt within a step, 0 at t = 0).
# The inverse masses are those of the factorization, which direct_interface_force
# computed before.
function apply_direct_interface_force!(
    bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition, λ::Matrix{Float64}, velocity_factor::Float64
)
    factor = direct_interface_factor(bc)
    d_model = self_subsim_of(bc).model
    n_model = coupled_subsim_of(bc).model
    n_bc = coupled_bc_of(bc)
    P = bc.dirichlet_projector
    f_N = zeros(3, length(n_bc.global_from_local_map))
    for comp in 1:3
        d_inv = factor.d_inverse_mass[comp]
        n_inv = factor.n_inverse_mass[comp]
        f_N[comp, :] = -(transpose(P) * λ[comp, :])
        for (i, node) in enumerate(bc.global_from_local_map)
            Δa = d_inv[i] * λ[comp, i]
            d_model.acceleration[comp, node] += Δa
            d_model.velocity[comp, node] += velocity_factor * Δa
        end
        for (i, node) in enumerate(n_bc.global_from_local_map)
            Δa = n_inv[i] * f_N[comp, i]
            n_model.acceleration[comp, node] += Δa
            n_model.velocity[comp, node] += velocity_factor * Δa
        end
    end
    # Record the force on the Neumann side, as the iterated exchange does, so
    # that its boundary force and the force residual describe the same state.
    n_bc.transferred_force = vec(f_N)
    for (i, node) in enumerate(n_bc.global_from_local_map)
        n_model.boundary_force[(3 * node - 2):(3 * node)] .+= f_N[:, i]
    end
    return nothing
end

# Solve H₀ λ = Π_D q_N - q_D for the interface force, with q the velocity or
# the acceleration of the free step, divided by `scale` (γΔt within a step,
# 1 at t = 0).
function direct_interface_force(bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition, field::Symbol, scale::Float64)
    factor = direct_interface_factor(bc)
    d_model = self_subsim_of(bc).model
    n_model = coupled_subsim_of(bc).model
    n_bc = coupled_bc_of(bc)
    P = bc.dirichlet_projector
    q_D = getfield(d_model, field)[:, bc.global_from_local_map]
    q_N = getfield(n_model, field)[:, n_bc.global_from_local_map]
    λ = zeros(size(q_D))
    for comp in 1:3
        rhs = (P * q_N[comp, :] .- q_D[comp, :]) ./ scale
        λ[comp, :] = factor.factors[comp] \ rhs
    end
    return λ
end

# True when the simulation's Schwarz couplings are direct pairs. A mixture of
# direct pairs and other couplings aborts.
function uses_direct_interface_solve(sim::MultiDomainSimulation)
    direct = 0
    other = 0
    for subsim in sim.subsims, bc in subsim.model.boundary_conditions
        bc isa SolidMechanicsSchwarzBoundaryCondition || continue
        is_direct_dn(bc) ? (direct += 1) : (other += 1)
    end
    direct == 0 && return false
    other == 0 || norma_abort("`interface solve: direct` requires that every Schwarz coupling be a direct pair.")
    return true
end

# One controller stop with the direct interface solve: each subdomain takes its
# step with the interface rows free and no interface force, the interface force
# is solved from the velocity constraint, and the acceleration and velocity of
# the interface rows are corrected. No Schwarz iteration.
function direct_interface_stop!(sim::MultiDomainSimulation)
    controller = sim.controller
    controller.is_schwarz = false
    save_stop_state(sim)
    set_initial_subcycle_time(sim)
    for subsim in sim.subsims
        integrator = subsim.integrator
        if integrator.minimum_time_step == integrator.maximum_time_step
            integrator.time_step = integrator.maximum_time_step
        end
        advance_time(subsim)
        # The step is clipped to land on the stop, so it differs from the nominal
        # step by the rounding of the accumulated time.
        if !isapprox(integrator.time_step, controller.time_step; rtol=1.0e-6) || !stop_subcyle(subsim)
            norma_abort(
                "`interface solve: direct` needs one step per stop, but subdomain '$(subsim.name)' took a step of " *
                "$(integrator.time_step) for the stop of $(controller.time_step) (the explicit stable step can " *
                "shorten it; raise CFL or lower the step).",
            )
        end
        advance_one_step(subsim)
    end
    for subsim in sim.subsims, bc in subsim.model.boundary_conditions
        is_direct_dn(bc) && bc.is_dirichlet || continue
        γΔt = subsim.integrator.γ * subsim.integrator.time_step
        if bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
            # Both sides are constrained: the pair is solved once, from the
            # side listed first.
            bc.self_handle.id < bc.coupled_handle.id || continue
            direct_overlap_pair_step!(bc, :velocity, γΔt)
            continue
        end
        λ = direct_interface_force(bc, :velocity, γΔt)
        apply_direct_interface_force!(bc, λ, γΔt)
    end
    controller.schwarz_iters[controller.stop] = 0
    controller.is_schwarz = true
    return nothing
end

# Coupled initial acceleration of direct pairs: one solve of the acceleration
# constraint a_D = Π_D a_N with the same matrix, after the interface
# displacement and velocity of the Dirichlet side are set to the projected
# initial values of the partner.
function direct_initial_acceleration!(bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    bc.self_handle.id < bc.coupled_handle.id || return nothing
    # The displacement and velocity of both boundaries are set to the projected
    # initial values of the partner before the acceleration constraint is solved.
    for side in (bc, coupled_bc_of(bc))
        model = self_subsim_of(side).model
        partner = coupled_subsim_of(side).model
        P = side.dirichlet_projector
        for field in (:displacement, :velocity)
            q = getfield(model, field)
            q_partner = getfield(partner, field)
            for comp in 1:3
                q[comp, side.global_from_local_map] = P * q_partner[comp, :]
            end
        end
    end
    direct_overlap_pair_step!(bc, :acceleration, 1.0)
    for side in (bc, coupled_bc_of(bc))
        r = dn_interface_residuals(self_subsim_of(side).model, side)
        norma_logf(
            0,
            :acceleration,
            "Direct initial acceleration %s/%s: acceleration jump %.2e, force residual %.2e",
            r.dirichlet_name,
            r.neumann_name,
            r.acceleration_jump,
            r.force_residual,
        )
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Direct interface solve of an overlap matching pair of two central difference
# subdomains at equal steps. Both sides are constrained: λ₁ on the rows of Γ₁
# of side 1 and λ₂ on the rows of Γ₂ of side 2 are the unknowns, side k also
# receives -Π_jᵀ λ_j on all of its nodes, and the velocity constraints at the
# end of the step on both boundaries give the symmetric block system
#   [ A₁  -B ] [λ₁]   1   [Π₁ v₂ - C₁ v₁]
#   [-Bᵀ  A₂ ] [λ₂] = --- [Π₂ v₁ - C₂ v₂],   A_k = m_Γk⁻¹ + Π_k m_j⁻¹ Π_kᵀ,
#                     γΔt
#   B = C₁ m₁⁻¹ Π₂ᵀ + Π₁ m₂⁻¹ C₂ᵀ,
# with v the free velocities of the step and m the lumped masses (zero
# inverse on degrees of freedom that other Dirichlet conditions hold). B
# couples the rows of Γ₁ that Π₂ reaches and those of Γ₂ that Π₁ reaches, and
# vanishes when the two boundaries are farther apart than one element. The
# matrix is constant and factored once per component.
# ---------------------------------------------------------------------------

struct DirectOverlapFactor
    factors::Vector{Any}                      # Cholesky factor of the block matrix per component
    inverse_mass_1::Vector{Vector{Float64}}   # of all nodes of side 1, per component
    inverse_mass_2::Vector{Vector{Float64}}   # of all nodes of side 2, per component
end

const DIRECT_OVERLAP_CACHE = IdDict{Any,DirectOverlapFactor}()

# Inverse lumped mass of every node of a model for one component, zero where
# a Dirichlet condition holds the degree of freedom.
function full_inverse_mass(model::SolidMechanics, comp::Int, fixed::BitVector)
    m = model.lumped_mass
    n = size(model.reference, 2)
    return [fixed[3 * (k - 1) + comp] ? 0.0 : 1.0 / m[3 * (k - 1) + comp] for k in 1:n]
end

function direct_overlap_factor(bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    return get!(DIRECT_OVERLAP_CACHE, bc) do
        bc_1, bc_2 = bc, coupled_bc_of(bc)
        model_1 = self_subsim_of(bc_1).model
        model_2 = self_subsim_of(bc_2).model
        P_1, P_2 = bc_1.dirichlet_projector, bc_2.dirichlet_projector
        Γ_1, Γ_2 = bc_1.global_from_local_map, bc_2.global_from_local_map
        fixed_1, fixed_2 = prescribed_dofs(model_1), prescribed_dofs(model_2)
        t0 = time()
        factors = Any[]
        inverse_mass_1 = Vector{Float64}[]
        inverse_mass_2 = Vector{Float64}[]
        for comp in 1:3
            inv_1 = full_inverse_mass(model_1, comp, fixed_1)
            inv_2 = full_inverse_mass(model_2, comp, fixed_2)
            push!(inverse_mass_1, inv_1)
            push!(inverse_mass_2, inv_2)
            n_1, n_2 = length(Γ_1), length(Γ_2)
            A_1 = P_1 * (inv_2 .* transpose(P_1))
            A_2 = P_2 * (inv_1 .* transpose(P_2))
            for (i, node) in enumerate(Γ_1)
                A_1[i, i] += inv_1[node]
            end
            for (i, node) in enumerate(Γ_2)
                A_2[i, i] += inv_2[node]
            end
            # B = C₁ m₁⁻¹ Π₂ᵀ + Π₁ m₂⁻¹ C₂ᵀ: the rows Γ₁ of m₁⁻¹ Π₂ᵀ and the
            # columns Γ₂ of Π₁ m₂⁻¹.
            B = (inv_1[Γ_1] .* transpose(P_2)[Γ_1, :]) + (P_1[:, Γ_2] .* transpose(inv_2[Γ_2]))
            H = [A_1 -B; -transpose(B) A_2]
            F = cholesky(Symmetric(H); check=false)
            if !issuccess(F)
                norma_abort(
                    "The direct interface matrix of the overlap matching pair '$(bc_1.name)'/'$(bc_2.name)' is not " *
                    "positive definite (component $comp): a constrained degree of freedom is held by Dirichlet " *
                    "conditions on both sides.",
                )
            end
            push!(factors, F)
        end
        norma_logf(
            0,
            :setup,
            "Direct interface solve of the overlap matching pair '%s'/'%s': %d + %d constrained nodes, three " *
            "factorizations in %.2f s.",
            bc_1.name,
            bc_2.name,
            length(Γ_1),
            length(Γ_2),
            time() - t0,
        )
        return DirectOverlapFactor(factors, inverse_mass_1, inverse_mass_2)
    end
end

# One step of the direct solve of the pair: the reactions from the velocity
# (or acceleration) constraints on both boundaries after the free step, and
# the correction of both sides: Δa_k = m_k⁻¹ (C_kᵀ λ_k - Π_jᵀ λ_j) on all of
# its nodes and Δv_k = velocity_factor Δa_k. The forces received from the
# partner are recorded as the iterated exchange does.
function direct_overlap_pair_step!(
    bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition, field::Symbol, scale::Float64
)
    factor = direct_overlap_factor(bc)
    bc_1, bc_2 = bc, coupled_bc_of(bc)
    model_1 = self_subsim_of(bc_1).model
    model_2 = self_subsim_of(bc_2).model
    P_1, P_2 = bc_1.dirichlet_projector, bc_2.dirichlet_projector
    Γ_1, Γ_2 = bc_1.global_from_local_map, bc_2.global_from_local_map
    n_1, n_2 = length(Γ_1), length(Γ_2)
    q_1 = getfield(model_1, field)
    q_2 = getfield(model_2, field)
    f_1 = zeros(3, size(model_1.reference, 2))   # force received by side 1 from λ₂
    f_2 = zeros(3, size(model_2.reference, 2))   # force received by side 2 from λ₁
    velocity_factor = field == :velocity ? scale : 0.0
    for comp in 1:3
        rhs = vcat(P_1 * q_2[comp, :] .- q_1[comp, Γ_1], P_2 * q_1[comp, :] .- q_2[comp, Γ_2]) ./ scale
        λ = factor.factors[comp] \ rhs
        λ_1 = λ[1:n_1]
        λ_2 = λ[(n_1 + 1):(n_1 + n_2)]
        f_1[comp, :] = -(transpose(P_2) * λ_2)
        f_2[comp, :] = -(transpose(P_1) * λ_1)
        force_1 = copy(f_1[comp, :])
        force_2 = copy(f_2[comp, :])
        force_1[Γ_1] .+= λ_1
        force_2[Γ_2] .+= λ_2
        Δa_1 = factor.inverse_mass_1[comp] .* force_1
        Δa_2 = factor.inverse_mass_2[comp] .* force_2
        model_1.acceleration[comp, :] .+= Δa_1
        model_2.acceleration[comp, :] .+= Δa_2
        model_1.velocity[comp, :] .+= velocity_factor .* Δa_1
        model_2.velocity[comp, :] .+= velocity_factor .* Δa_2
    end
    bc_1.transferred_force = vec(f_1)
    bc_2.transferred_force = vec(f_2)
    model_1.boundary_force .+= vec(f_1)
    model_2.boundary_force .+= vec(f_2)
    return nothing
end

function direct_initial_acceleration!(bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    d_model = self_subsim_of(bc).model
    n_model = coupled_subsim_of(bc).model
    n_bc = coupled_bc_of(bc)
    P = bc.dirichlet_projector
    for field in (:displacement, :velocity)
        q_N = getfield(n_model, field)[:, n_bc.global_from_local_map]
        q_D = getfield(d_model, field)
        for comp in 1:3
            q_D[comp, bc.global_from_local_map] = P * q_N[comp, :]
        end
    end
    λ = direct_interface_force(bc, :acceleration, 1.0)
    apply_direct_interface_force!(bc, λ, 0.0)
    r = dn_interface_residuals(d_model, bc)
    norma_logf(
        0,
        :acceleration,
        "Direct initial acceleration %s/%s: acceleration jump %.2e, force residual %.2e",
        r.dirichlet_name,
        r.neumann_name,
        r.acceleration_jump,
        r.force_residual,
    )
    return nothing
end

# --- Dirichlet-Neumann interface residuals ---------------------------------
#
# For a DN pair with Dirichlet side D and Neumann side N, Π_D the Dirichlet
# projector (interface values of N to interface values of D) and W_D, W_N the
# square projection (boundary mass) matrices of the two interfaces:
#
#   jump of q:       ‖q_D - Π_D q_N‖_{W_D} / ‖q_D‖_{W_D},  q = velocity, displacement,
#   force residual:  ‖Π_Dᵀ r_D - f_N‖_{W_N⁻¹} / ‖Π_Dᵀ r_D‖_{W_N⁻¹},
#
# with ‖x‖²_W = Σ_c x_cᵀ W x_c over the three Cartesian components, r_D the
# interface reaction of D, the restriction to its interface rows of
# -(M a + f_int - f_body - f_boundary) evaluated on its current state, and f_N
# the interface force applied on N in the same Schwarz iteration. The jump is
# one sided: it is measured on the Dirichlet side. Its root mean square value
# over the interface, ‖q_D - Π_D q_N‖_{W_D} / √|Γ_D|, is also returned, so that
# a jump of a quantity that is itself near zero can be compared with the
# absolute tolerance.
struct DNInterfaceResiduals
    dirichlet_name::String
    neumann_name::String
    velocity_jump::Float64
    displacement_jump::Float64
    velocity_jump_rms::Float64
    displacement_jump_rms::Float64
    force_residual::Float64
    acceleration_jump::Float64
    acceleration_jump_rms::Float64
end

function weighted_norm(W::AbstractMatrix{Float64}, x::Matrix{Float64})
    total = 0.0
    for comp in axes(x, 1)
        xc = x[comp, :]
        total += dot(xc, W * xc)
    end
    return sqrt(max(total, 0.0))
end

# ‖x‖_{W⁻¹}, with W given by its Cholesky factorization.
function inverse_weighted_norm(W_factor::Cholesky{Float64,Matrix{Float64}}, x::Matrix{Float64})
    total = 0.0
    for comp in axes(x, 1)
        xc = x[comp, :]
        total += dot(xc, W_factor \ xc)
    end
    return sqrt(max(total, 0.0))
end

# The operator A applied to each of the three components of a nodal field q
# (3 × n), row by row.
function project_components(A::AbstractMatrix{Float64}, q::AbstractMatrix{Float64})
    result = zeros(3, size(A, 1))
    for comp in 1:3
        result[comp, :] = A * q[comp, :]
    end
    return result
end

# Force residuals of the matching exchanges. Each is the difference between
# the reaction transferred to a side and the force applied there, relative to
# the larger of the transferred force and the momentum rate ‖M v‖/Δt of that
# side on the rows where the force acts, both in the dual norm of the mass
# that carries the force: W⁻¹ of the boundary mass on an interface, the
# inverse lumped mass on the volume rows of an overlap. A residual of a side
# that receives the partner's reaction of the previous iteration is the
# change of that reaction between iterations, and the window-end ramp
# residual compares two reactions; when the interface forces are at roundoff,
# as in a body that still moves rigidly, their relative change is not a
# measure of convergence, and the momentum rate, the force that would stop the
# rows in one step, is the scale of the forces that act on them.

# Momentum M v of a subdomain with the mass its integrator uses, and the
# row-sum lumped mass, per degree of freedom; nothing before the first
# evaluation assembles a mass.
function momentum_and_lumped_mass(model::SolidMechanics)
    v = vec(model.velocity)
    if length(model.lumped_mass) == length(v)
        return model.lumped_mass .* v, model.lumped_mass
    elseif size(model.mass, 1) == length(v)
        return model.mass * v, vec(sum(model.mass; dims=2))
    end
    return nothing
end

# ‖C (M v)‖_{W⁻¹}/Δt of a side on the interface rows of `bc`, with W given by
# its Cholesky factorization.
function interface_momentum_rate(
    sim::SingleDomainSimulation, bc::SolidMechanicsSchwarzBoundaryCondition, W_factor::Cholesky{Float64,Matrix{Float64}}
)
    model = get_fom_model(sim)
    p = mass_rows_product(model, interface_dofs(bc), vec(model.velocity), true)
    p === nothing && return 0.0
    Δt = sim.integrator.time_step
    Δt > 0.0 || return 0.0
    return inverse_weighted_norm(W_factor, reshape(p, 3, :)) / Δt
end

# Relative difference, in the norm of W_N⁻¹, between the reaction r_D of the
# Dirichlet side transferred by Pᵀ and the force `applied` to the Neumann side,
# with W_N given by its Cholesky factorization, relative to the larger of the
# transferred force and the momentum rate of the Neumann side's interface.
function dn_force_residual(
    W_N_factor::Cholesky{Float64,Matrix{Float64}},
    P::AbstractMatrix{Float64},
    r_D::AbstractMatrix{Float64},
    applied::Matrix{Float64},
    momentum_rate::Float64,
)
    transferred = project_components(transpose(P), r_D)
    scale = max(inverse_weighted_norm(W_N_factor, transferred), momentum_rate)
    difference = inverse_weighted_norm(W_N_factor, transferred - applied)
    return scale > 0.0 ? difference / scale : difference
end

# Nodes of the receiving side on which the transferred force of an overlap
# matching side acts: the columns of P with a nonzero entry, cached per
# projector.
const OVERLAP_FORCE_ROWS = IdDict{Any,BitVector}()
overlap_force_rows(P::AbstractMatrix{Float64}) = get!(() -> vec(any(!iszero, P; dims=1)), OVERLAP_FORCE_ROWS, P)

# Relative difference between the reaction r_D of an overlap matching side
# transferred by Pᵀ to the partner's nodes and the force `applied` there, in
# the dual norm of the receiving side's lumped mass on the rows where the
# force acts, relative to the larger of the transferred force and the
# momentum rate of those rows.
function overlap_force_residual(
    P::AbstractMatrix{Float64}, r_D::AbstractMatrix{Float64}, applied::Matrix{Float64}, receiver::SingleDomainSimulation
)
    transferred = project_components(transpose(P), r_D)
    rows = overlap_force_rows(P)
    pm = momentum_and_lumped_mass(get_fom_model(receiver))
    if pm === nothing
        scale = norm(transferred[:, rows])
        difference = norm((transferred - applied)[:, rows])
        return scale > 0.0 ? difference / scale : difference
    end
    m = reshape(pm[2], 3, :)[:, rows]
    v = reshape(vec(get_fom_model(receiver).velocity), 3, :)[:, rows]
    Δt = receiver.integrator.time_step
    dual(f) = sqrt(sum(f[:, rows] .^ 2 ./ m))
    momentum_rate = Δt > 0.0 ? sqrt(sum(m .* v .^ 2)) / Δt : 0.0
    scale = max(dual(transferred), momentum_rate)
    difference = dual(transferred - applied)
    return scale > 0.0 ? difference / scale : difference
end

function dn_one_sided_jump(
    W::AbstractMatrix{Float64}, P::AbstractMatrix{Float64}, q_D::Matrix{Float64}, q_N::Matrix{Float64}
)
    transferred = project_components(P, q_N)
    jump = weighted_norm(W, q_D - transferred)
    scale = weighted_norm(W, q_D)
    area = sum(W)
    rms = area > 0.0 ? jump / sqrt(area) : jump
    relative = scale > 0.0 ? jump / scale : (jump > 0.0 ? Inf : 0.0)
    return relative, rms
end

# Residuals of the DN pair whose Dirichlet side is `bc`, with `model` the full
# order model of the Dirichlet side. The partner's boundary condition belongs to
# the partner's own model (n_sim.model), as in get_dst_force: for a reduced
# order partner that is the RomModel, whose full order model carries the
# kinematic fields but an empty list of boundary conditions. The fields are
# therefore read from get_fom_model(n_sim) and the condition from n_sim.model.
function dn_interface_residuals(model::SolidMechanics, bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    n_sim = coupled_subsim_of(bc)
    n_model = get_fom_model(n_sim)
    n_bc = coupled_bc_of(bc)
    W_D = bc.square_projector
    P = bc.dirichlet_projector
    d_map = bc.global_from_local_map
    n_map = partner_node_map(bc)
    jv, jv_rms = dn_one_sided_jump(W_D, P, model.velocity[:, d_map], n_model.velocity[:, n_map])
    ju, ju_rms = dn_one_sided_jump(W_D, P, model.displacement[:, d_map], n_model.displacement[:, n_map])
    ja, ja_rms = dn_one_sided_jump(W_D, P, model.acceleration[:, d_map], n_model.acceleration[:, n_map])
    force_residual = NaN
    # A reduced order partner's coupling condition may be of another type,
    # without the transferred force and the square projector the residual needs.
    if n_bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
        f_N = n_bc.transferred_force
        W_N_factor = n_bc.square_factor
    elseif n_bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
        f_N = n_bc.transferred_force
        W_N_factor = nothing
    else
        f_N = Float64[]
        W_N_factor = nothing
    end
    if !isempty(f_N) && length(model.internal_force) == length(model.displacement)
        dofs = interface_dofs(bc)
        r_D = reshape(-(model.internal_force[dofs] + dalembert_inertia_minus_loads(model, dofs)), 3, :)
        if bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
            force_residual = overlap_force_residual(P, r_D, reshape(f_N, 3, :), n_sim)
        elseif W_N_factor !== nothing && size(W_N_factor, 1) == length(n_map)
            force_residual = dn_force_residual(
                W_N_factor, P, r_D, reshape(f_N, 3, :), interface_momentum_rate(n_sim, n_bc, W_N_factor)
            )
        end
    end
    return DNInterfaceResiduals(self_subsim_of(bc).name, n_sim.name, jv, ju, jv_rms, ju_rms, force_residual, ja, ja_rms)
end

# The Dirichlet side of every DN pair of a multidomain simulation whose
# projectors are built, with the full order model of that side, in the order of
# the subdomains that hold it.
function dn_dirichlet_sides(sim::MultiDomainSimulation)
    sides = Tuple{SolidMechanics,SolidMechanicsConstrainedSchwarzBoundaryCondition}[]
    for subsim in sim.subsims
        # The conditions belong to subsim.model (a RomModel for a reduced order
        # subdomain); the fields to its full order model.
        model = get_fom_model(subsim)
        model isa SolidMechanics || continue
        for bc in subsim.model.boundary_conditions
            bc isa SolidMechanicsConstrainedSchwarzBoundaryCondition || continue
            bc.is_dirichlet || continue
            size(bc.dirichlet_projector, 1) > 0 || continue
            size(bc.square_projector, 1) > 0 || continue
            push!(sides, (model, bc))
        end
    end
    return sides
end

# Residuals of every DN pair of a multidomain simulation, one entry per pair,
# in the order of the subdomains that hold the Dirichlet side.
function dn_interface_residuals(sim::MultiDomainSimulation)
    return [(bc, dn_interface_residuals(model, bc)) for (model, bc) in dn_dirichlet_sides(sim)]
end

# Residuals of a constrained pair at every substep of the stop, from the
# substep histories of the Schwarz iteration just completed: the one-sided
# jumps at each substep time of the Dirichlet side, against the Neumann side's
# history interpolated linearly to that time, and the force residual at each
# substep time of the Neumann side, against the Dirichlet side's d'Alembert
# reaction interpolated linearly to that time. Each reported value is the
# largest over the substeps. At equal steps there is one substep and the values
# are those of dn_interface_residuals.
function dn_substep_residuals(sim::MultiDomainSimulation, bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    controller = sim.controller
    d_id = bc.self_handle.id
    n_id = bc.coupled_handle.id
    d_model = get_fom_model(self_subsim_of(bc))
    n_sim = coupled_subsim_of(bc)
    n_model = get_fom_model(n_sim)
    n_bc = coupled_bc_of(bc)
    W_D = bc.square_projector
    P = bc.dirichlet_projector
    d_map = bc.global_from_local_map
    n_map = partner_node_map(bc)
    t_start = controller.prev_time
    tol_t = 1.0e-12 * max(1.0, abs(controller.time))
    d_times = controller.time_hist[d_id]
    n_times = controller.time_hist[n_id]
    if isempty(d_times) || isempty(n_times)
        return dn_interface_residuals(d_model, bc)
    end
    nodal(v) = reshape(v, 3, :)
    jv = ju = jv_rms = ju_rms = 0.0
    # The partner's fields are interpolated on its interface rows only.
    n_dofs = node_dofs(n_map)
    for (k, t) in enumerate(d_times)
        t > t_start + tol_t || continue
        # A window-end-constrained side matches at the end of the window only.
        bc.window_end_constraint && k < length(d_times) && continue
        v_N = nodal(interpolate(n_times, controller.velo_hist[n_id], t, n_dofs))
        u_N = nodal(interpolate(n_times, controller.disp_hist[n_id], t, n_dofs))
        r_v, rms_v = dn_one_sided_jump(W_D, P, nodal(controller.velo_hist[d_id][k])[:, d_map], v_N)
        r_u, rms_u = dn_one_sided_jump(W_D, P, nodal(controller.disp_hist[d_id][k])[:, d_map], u_N)
        jv = max(jv, r_v)
        jv_rms = max(jv_rms, rms_v)
        ju = max(ju, r_u)
        ju_rms = max(ju_rms, rms_u)
    end
    force_residual = 0.0
    found = false
    W_N_factor = n_bc.square_factor
    d_dofs = interface_dofs(bc)
    for (t, f_N) in n_bc.transferred_force_history
        t > t_start + tol_t || continue
        found = true
        # The inertial rows need the whole acceleration (the rows of the mass
        # reach the neighbors of the interface nodes); the force only its rows.
        a_D = interpolate(d_times, controller.acce_hist[d_id], t)
        f_int = interpolate(d_times, controller.∂Ω_f_hist[d_id], t, d_dofs)
        r_D = nodal(-(f_int + dalembert_inertia_minus_loads(d_model, d_dofs, a_D)))
        if bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
            force_residual = max(force_residual, overlap_force_residual(P, r_D, nodal(f_N), n_sim))
        else
            force_residual = max(
                force_residual,
                dn_force_residual(W_N_factor, P, r_D, nodal(f_N), interface_momentum_rate(n_sim, n_bc, W_N_factor)),
            )
        end
    end
    found || (force_residual = NaN)
    # A window-end-constrained side: the end value of the ramp that its
    # interior substeps applied must equal the reaction at the end of the
    # window, or the window balance of the rate invariant acquires the term
    # (m C_D Δv_D^(m) - C_D Δv_D)ᵀ(r_D - ℓ)/ΔT. Measured in the W_D⁻¹ norm
    # relative to the larger of the reaction and the momentum rate of the
    # interface, and folded into the force residual.
    if bc.window_end_constraint &&
        !isempty(bc.window_ramp_end) &&
        length(bc.window_ramp_end) == length(bc.window_end_reaction_raw) &&
        bc.square_factor !== nothing
        r_end = nodal(bc.window_end_reaction_raw)
        scale = max(
            inverse_weighted_norm(bc.square_factor, r_end),
            interface_momentum_rate(self_subsim_of(bc), bc, bc.square_factor),
        )
        difference = inverse_weighted_norm(bc.square_factor, r_end - nodal(bc.window_ramp_end))
        ramp_residual = scale > 0.0 ? difference / scale : difference
        force_residual = isnan(force_residual) ? ramp_residual : max(force_residual, ramp_residual)
    end
    return DNInterfaceResiduals(self_subsim_of(bc).name, n_sim.name, jv, ju, jv_rms, ju_rms, force_residual, NaN, NaN)
end

# Interface impulse residual of a constrained pair over the stop just
# completed: the trapezoid-in-time sum of the force applied on the Neumann side
# over its substeps, from the start of the stop, minus Π_Dᵀ times the
# trapezoid-in-time sum of the Dirichlet side's reaction over its own
# substeps, from its stop-start anchor. Returns the sum of the residual over
# the interface nodes per component (the net impulse error, in units of force times time) and the
# residual in the W_N⁻¹ norm relative to the transferred impulse. At the fixed
# point of the exchange the residual vanishes when the force applied on the
# Neumann side is linear in time between the shared end reactions (Connors et
# al. 2024, Eq. (80)).
function dn_impulse_residual(sim::MultiDomainSimulation, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    controller = sim.controller
    nan = (fill(NaN, 3), NaN)
    d_id = bc.self_handle.id
    d_times = controller.time_hist[d_id]
    n_bc = coupled_bc_of(bc)
    n_bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition || return nan
    t_start = controller.prev_time
    tol_t = 1.0e-12 * max(1.0, abs(controller.time))
    d_model = get_fom_model(self_subsim_of(bc))
    P = bc.dirichlet_projector
    W_N_factor = n_bc.square_factor
    W_N_factor === nothing && return nan
    length(d_times) >= 2 || return nan
    # Dirichlet side: reaction at each of its snapshots, trapezoid in time. A
    # window-end-constrained side carries the ramp at its interior substeps.
    function reaction(k)
        if bc.window_end_constraint && k < length(d_times)
            # The ramp is the force on the rows; the reaction here is its negative.
            return reshape(-window_ramp(bc, d_times[k]), 3, :)
        end
        a_k = controller.acce_hist[d_id][k]
        dofs = interface_dofs(bc)
        return reshape(-(controller.∂Ω_f_hist[d_id][k][dofs] + dalembert_inertia_minus_loads(d_model, dofs, a_k)), 3, :)
    end
    impulse_D = zeros(3, size(P, 1))
    for k in 2:length(d_times)
        impulse_D .+= 0.5 * (d_times[k] - d_times[k - 1]) .* (reaction(k - 1) .+ reaction(k))
    end
    # Neumann side: applied force from the window start, trapezoid in time.
    entries = filter(e -> e[1] >= t_start - tol_t, n_bc.transferred_force_history)
    (length(entries) >= 2 && abs(entries[1][1] - t_start) <= tol_t) || return nan
    impulse_N = zeros(3, size(P, 2))
    for k in 2:length(entries)
        impulse_N .+=
            0.5 * (entries[k][1] - entries[k - 1][1]) .*
            (reshape(entries[k - 1][2], 3, :) .+ reshape(entries[k][2], 3, :))
    end
    transferred = project_components(transpose(P), impulse_D)
    residual = impulse_N .- transferred
    net = vec(sum(residual; dims=2))
    scale = inverse_weighted_norm(W_N_factor, transferred)
    relative = scale > 0.0 ? inverse_weighted_norm(W_N_factor, residual) / scale : NaN
    return net, relative
end

# A pair under the displacement constraint requires both members to have taken
# the same step; the explicit stable-step cap can shorten one of them.
function check_constrained_time_steps(bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    bc.constraint == :displacement || return nothing
    d_sim = self_subsim_of(bc)
    n_sim = coupled_subsim_of(bc)
    Δt_D = d_sim.integrator.time_step
    Δt_N = n_sim.integrator.time_step
    if !isapprox(Δt_D, Δt_N; rtol=1.0e-12)
        norma_abort(
            "`matching: displacement` is implemented for equal time steps, but subdomains " *
            "'$(d_sim.name)' and '$(n_sim.name)' took steps $(Δt_D) and $(Δt_N) " *
            "(the explicit stable step can shorten the requested step; raise CFL or lower the step).",
        )
    end
    return nothing
end

# ---------------------------------------------------------------------------

# Locate `point` in block `block_id` of `model` on the reference configuration,
# returning the element's node indices and the parametric coordinates. Only
# the elements binned in the point's grid cell are tested (see
# spatial_search.jl); the bins are padded by the same margin as is_inside's
# bounding-box prefilter, so the result matches a scan of every element.
function find_point_in_mesh(point::Vector{Float64}, model::SolidMechanics, block_id::Int, tol::Float64)
    block_index = findfirst(block -> block.id == block_id, model.blocks)
    block_index === nothing && return Int64[], zeros(length(point)), false
    block = model.blocks[block_index]
    elements = get_element_grid(model)
    for item in box_grid_items(elements.grid, SVector{3,Float64}(point))
        elements.block_index[item] == block_index || continue
        node_indices = block.connectivity[:, elements.element_index[item]]
        ξ, found = is_inside(block.element_type, model.reference[:, node_indices], point, tol)
        found && return node_indices, ξ, true
    end
    return Int64[], zeros(length(point)), false
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    if bc.is_dirichlet == true
        contact_weak_dbc(model, bc)
    else
        contact_weak_nbc(model, bc)
    end
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    if bc.use_weak
        coupling_weak_overlap_dbc(model, bc)
    else
        coupling_strong_dbc(model, bc)
    end
    return nothing
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    if bc.is_dirichlet == true
        coupling_weak_dbc(model, bc)
    else
        coupling_weak_nbc(model, bc)
    end
end

# Overlap matching coupling: the side is constrained on its own Schwarz
# boundary as the Dirichlet side of a nonoverlap matching pair, and receives on
# all of its nodes the partner's reaction on the partner's Schwarz boundary
# through the transpose of the partner's projector.
function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    coupling_weak_dbc(model, bc)
    overlap_matching_nbc(model, bc)
    return nothing
end

function overlap_matching_nbc(model::SolidMechanics, bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    size(bc.neumann_projector, 1) > 0 || return nothing
    nodal_force = get_dst_force(bc)
    bc.transferred_force = nodal_force
    history = bc.transferred_force_history
    t_start = bc.parent.controller.prev_time
    tol_t = 1.0e-12 * max(1.0, abs(model.time))
    if !isempty(history) && model.time <= history[end][1]
        filter!(entry -> entry[1] <= t_start + tol_t, history)
    end
    while length(history) >= 2 && history[2][1] <= t_start + tol_t
        popfirst!(history)
    end
    push!(history, (model.time, nodal_force))
    model.boundary_force .+= nodal_force
    return nothing
end

function coupling_strong_dbc(model::SolidMechanics, bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    coupled_model_obj = coupled_subsim_of(bc).model
    get_coupled_field = if coupled_model_obj isa SolidMechanics
        (field -> getfield(coupled_model_obj, field))
    else
        (field -> getfield(coupled_model_obj.fom_model, field))
    end

    coupled_reference = get_coupled_field(:reference)
    coupled_displacement = get_coupled_field(:displacement)
    velocity = get_coupled_field(:velocity)
    acceleration = get_coupled_field(:acceleration)

    unique_node_indices = unique(bc.side_set_node_indices)

    for i in eachindex(unique_node_indices)
        node_index = unique_node_indices[i]
        coupled_node_indices = bc.coupled_nodes_indices[i]
        N = bc.interpolation_function_values[i]

        model.displacement[:, node_index] =
            (coupled_reference[:, coupled_node_indices] + coupled_displacement[:, coupled_node_indices]) * N -
            model.reference[:, node_index]
        model.velocity[:, node_index] = velocity[:, coupled_node_indices] * N
        model.acceleration[:, node_index] = acceleration[:, coupled_node_indices] * N

        dof_index = (3 * node_index - 2):(3 * node_index)
        model.free_dofs[dof_index] .= false
    end
end

function coupling_weak_overlap_dbc(model::SolidMechanics, bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    coupled_model_obj = coupled_subsim_of(bc).model
    src_model = coupled_model_obj isa SolidMechanics ? coupled_model_obj : coupled_model_obj.fom_model

    P = bc.dirichlet_projector
    global_from_local_map = bc.global_from_local_map

    for comp in 1:3
        src_current = src_model.reference[comp, :] + src_model.displacement[comp, :]
        src_velocity = src_model.velocity[comp, :]
        src_acceleration = src_model.acceleration[comp, :]
        proj_current = P * src_current
        proj_velocity = P * src_velocity
        proj_acceleration = P * src_acceleration
        for (i_local, i_global) in enumerate(global_from_local_map)
            model.displacement[comp, i_global] = proj_current[i_local] - model.reference[comp, i_global]
            model.velocity[comp, i_global] = proj_velocity[i_local]
            model.acceleration[comp, i_global] = proj_acceleration[i_local]
        end
    end

    for i_global in global_from_local_map
        dof_index = (3 * i_global - 2):(3 * i_global)
        model.free_dofs[dof_index] .= false
    end
end

function coupling_weak_dbc(model::SolidMechanics, bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    # Transfer the projected DISPLACEMENT, not the projected current position.
    # On a curved interface the two sides discretize the geometry as different
    # facet polyhedra, so the projector does not map the source reference onto
    # the destination reference (P x_src ≠ x_dst; the gap is the coarser
    # side's facet sagitta). Position transfer injects that geometric mismatch
    # as a spurious scalloped Dirichlet displacement at the coarse-facet
    # frequency; displacement transfer leaves each side its own reference
    # geometry and its error is O(h²) in the field. On flat interfaces the two
    # forms coincide because the L2 projection reproduces linear functions.
    # (Contact keeps position transfer: closure there is coincidence of the
    # CURRENT surfaces, see contact_weak_dbc.)
    if interior_substep(bc)
        apply_window_ramp!(model, bc)
        return nothing
    end
    _, nodal_disp, nodal_velo, nodal_acce = get_dst_curr_disp_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    if bc.constrained && constrained_step_in_progress(bc)
        impose_constrained_kinematics!(model, bc, nodal_disp, nodal_velo)
        return nothing
    end
    for (i_local, i_global) in enumerate(global_from_local_map)
        @inbounds model.displacement[:, i_global] = nodal_disp[:, i_local]
        @inbounds model.velocity[:, i_global] = nodal_velo[:, i_local]
        @inbounds model.acceleration[:, i_global] = nodal_acce[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        model.free_dofs[global_range] .= false
    end
end

# The constrained imposition applies within a time step, where the receiver's
# fields at the interface hold its state at the start of the step. At
# initialization, before any step, the partner's projected fields are imposed
# as in the unconstrained exchange.
function constrained_step_in_progress(bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    integrator = self_subsim_of(bc).integrator
    return isfinite(integrator.prev_time) && integrator.time > integrator.prev_time
end

# True at an interior substep of the window of a side with `matching at:
# window end`: the step in progress does not end at the controller's stop time.
function interior_substep(bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    bc.window_end_constraint || return false
    constrained_step_in_progress(bc) || return false
    integrator = self_subsim_of(bc).integrator
    stop_time = bc.parent.controller.time
    return !isapprox(integrator.time, stop_time; rtol=1.0e-6, atol=1.0e-12)
end

# Reaction of a constrained side on its interface rows at its current state:
# the rows of M a + f_int - f_body - f_boundary, in interleaved components of
# the interface nodes.
function window_reaction(model::SolidMechanics, bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    length(model.internal_force) == length(model.displacement) || return Float64[]
    dofs = interface_dofs(bc)
    return model.internal_force[dofs] + dalembert_inertia_minus_loads(model, dofs)
end

# Interface force of a window-end-constrained side at time t of the window:
# the linear ramp between its window-start and window-end reactions.
function window_ramp(bc::SolidMechanicsConstrainedSchwarzBoundaryCondition, t::Float64)
    controller = bc.parent.controller
    λ0 = bc.window_start_reaction
    λ1 = isempty(bc.window_end_reaction) ? λ0 : bc.window_end_reaction
    t0, t1 = controller.prev_time, controller.time
    s = t1 > t0 ? clamp((t - t0) / (t1 - t0), 0.0, 1.0) : 1.0
    return (1.0 - s) .* λ0 .+ s .* λ1
end

# Interior substep of a window-end-constrained side: its interface rows are
# free and carry the ramp of its own window reaction as an applied force.
function apply_window_ramp!(model::SolidMechanics, bc::SolidMechanicsConstrainedSchwarzBoundaryCondition)
    isempty(bc.window_start_reaction) && return nothing
    force = window_ramp(bc, self_subsim_of(bc).integrator.time)
    bc.window_ramp_end = copy(isempty(bc.window_end_reaction) ? bc.window_start_reaction : bc.window_end_reaction)
    for (i_local, i_global) in enumerate(bc.global_from_local_map)
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        local_range = (3 * (i_local - 1) + 1):(3 * i_local)
        @inbounds model.boundary_force[global_range] += force[local_range]
    end
    return nothing
end

# Record the reaction of every window-end-constrained side of a subdomain at
# the end of its pass over the window (subcycle), for the ramp of the next
# Schwarz iteration and the window start of the next stop.
function record_window_end_reactions!(subsim::SingleDomainSimulation)
    model = subsim.model
    model isa SolidMechanics || return nothing
    for bc in model.boundary_conditions
        bc isa SolidMechanicsConstrainedSchwarzBoundaryCondition && bc.window_end_constraint || continue
        bc.is_dirichlet || continue
        reaction = window_reaction(model, bc)
        bc.window_end_reaction_raw = copy(reaction)
        # The ramp of the next iteration is formed from the relaxed reaction:
        # the interior substeps respond to the ramp and the end constraint
        # removes that response through the reaction, a dependence with a
        # gain of about -(m - 1) for m substeps that the relaxation of the
        # datum does not act on. Under Anderson acceleration the end value of
        # the ramp is a second field of the interface's entry: it keeps the
        # value applied in this pass, so that the joint step of the next
        # iteration forms its residual from the measured reaction against
        # that value, and the step sets the new end value with the datum.
        # Overwriting it here with the measured reaction made that residual
        # identically zero.
        controller = bc.parent.controller
        θ = controller.relaxation_parameter
        anderson = controller.relaxation_method === :anderson
        same_size = length(bc.window_end_reaction) == length(reaction)
        if anderson
            same_size || (bc.window_end_reaction = reaction)
            continue
        end
        if controller.iteration_number > 0 && same_size
            reaction = θ .* reaction .+ (1.0 - θ) .* bc.window_end_reaction
        end
        bc.window_end_reaction = reaction
    end
    return nothing
end

# At the start of a stop the converged reaction of the previous stop becomes
# the window-start reaction; at the first stop it is the reaction of the
# initial state after the coupled initial acceleration.
function start_window_reactions!(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        model = subsim.model
        model isa SolidMechanics || continue
        for bc in model.boundary_conditions
            bc isa SolidMechanicsConstrainedSchwarzBoundaryCondition && bc.window_end_constraint || continue
            bc.is_dirichlet || continue
            # The converged reaction of the previous stop, as measured at the
            # end of its last pass, starts the new window and is the first
            # end value of its ramp.
            start = if !isempty(bc.window_end_reaction_raw)
                bc.window_end_reaction_raw
            elseif !isempty(bc.window_end_reaction)
                bc.window_end_reaction
            else
                window_reaction(model, bc)
            end
            bc.window_start_reaction = copy(start)
            bc.window_end_reaction = copy(start)
        end
    end
    return nothing
end

# Newmark parameters (β, γ) of a receiver; central difference is β = 0.
newmark_parameters(integrator::Newmark) = (integrator.β, integrator.γ)
newmark_parameters(integrator::CentralDifference) = (0.0, integrator.γ)

# Constrained Dirichlet side: impose the projected constrained quantity of the
# partner and derive the other two fields from the receiver's own Newmark
# relations over the current step. apply_bcs runs before the predictor, so the
# interface fields of the receiver still hold its state (u_n, v_n, a_n) at the
# start of the step, from which
#   u_pre = u_n + Δt v_n + (1/2 - β) Δt² a_n,   v_pre = v_n + (1 - γ) Δt a_n.
# Velocity constraint: v = Π v_N, a = (v - v_pre)/(γ Δt), u = u_pre + β Δt² a
# (central difference: β = 0, so u is the predictor value).
# Displacement constraint: u = Π u_N, a = (u - u_pre)/(β Δt²), v = v_pre + γ Δt a.
# The imposed values are those of the end of the step; the Newmark predictor
# leaves fixed degrees of freedom unchanged and the corrector updates only the
# free ones, so the receiver ends the step with these values.
function impose_constrained_kinematics!(
    model::SolidMechanics,
    bc::SolidMechanicsConstrainedSchwarzBoundaryCondition,
    nodal_disp::Matrix{Float64},
    nodal_velo::Matrix{Float64},
)
    integrator = self_subsim_of(bc).integrator
    Δt = integrator.time_step
    β, γ = newmark_parameters(integrator)
    for (i_local, i_global) in enumerate(bc.global_from_local_map)
        for comp in 1:3
            u_n = model.displacement[comp, i_global]
            v_n = model.velocity[comp, i_global]
            a_n = model.acceleration[comp, i_global]
            u_pre = u_n + Δt * v_n + (0.5 - β) * Δt * Δt * a_n
            v_pre = v_n + (1.0 - γ) * Δt * a_n
            if bc.constraint == :displacement
                u = nodal_disp[comp, i_local]
                a = (u - u_pre) / (β * Δt * Δt)
                v = v_pre + γ * Δt * a
            else
                v = nodal_velo[comp, i_local]
                a = (v - v_pre) / (γ * Δt)
                u = u_pre + β * Δt * Δt * a
            end
            model.displacement[comp, i_global] = u
            model.velocity[comp, i_global] = v
            model.acceleration[comp, i_global] = a
        end
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        model.free_dofs[global_range] .= false
    end
    return nothing
end

function coupling_weak_nbc(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    nodal_force = get_dst_force(bc)
    bc.transferred_force = nodal_force
    # Keep one entry per substep of the current Schwarz iteration of the current
    # stop, preceded by the last entry at or before the start of the stop (the
    # force at the window start, for the impulse residual): a new iteration
    # restarts at a time not after the last entry and drops the entries of the
    # previous iteration.
    history = bc.transferred_force_history
    t_start = bc.parent.controller.prev_time
    tol_t = 1.0e-12 * max(1.0, abs(model.time))
    if !isempty(history) && model.time <= history[end][1]
        filter!(entry -> entry[1] <= t_start + tol_t, history)
    end
    while length(history) >= 2 && history[2][1] <= t_start + tol_t
        popfirst!(history)
    end
    push!(history, (model.time, nodal_force))
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        local_range = (3 * (i_local - 1) + 1):(3 * i_local)
        @inbounds model.boundary_force[global_range] += nodal_force[local_range]
    end
end

function get_internal_force(model::SolidMechanics)
    return model.internal_force
end

# The model keeps one internal force vector, which the assembly writes in
# place (reduce_element_vector!); setting the force copies into it, so that a
# vector handed in is never written by a later assembly. Callers that keep the
# force across a set copy it first (get_internal_force returns the vector).
function set_internal_force!(model::SolidMechanics, force)
    current = model.internal_force
    if current === force
        return current
    elseif length(current) == length(force)
        copyto!(current, force)
    else
        model.internal_force = copy(force)
    end
    return model.internal_force
end

# Projected interface trace Π_D q_N of a partner field q (all degrees of
# freedom of the Neumann side, interleaved components), as a vector.
function constrained_partner_trace(bc::SolidMechanicsSchwarzBoundaryCondition, field::AbstractVector{Float64})
    q_N = reshape(field, 3, :)[:, partner_node_map(bc)]
    return vec(transpose(bc.dirichlet_projector * transpose(q_N)))
end

# Partner kinematic fields and boundary force at time t, as apply_bc reads them:
# interpolated from the partner history, the interface predictor at iteration 0,
# the stop state, or zero.
function partner_fields(bc::SolidMechanicsSchwarzBoundaryCondition, time::Float64)
    controller = bc.parent.controller
    coupled_model = coupled_subsim_of(bc).model
    coupled_index = bc.coupled_handle.id
    coupled_num_dofs = length(coupled_model.free_dofs)
    saved_∂Ω_f = copy(get_internal_force(coupled_model))
    time_hist = controller.time_hist[coupled_index]
    disp_hist = controller.disp_hist[coupled_index]
    velo_hist = controller.velo_hist[coupled_index]
    acce_hist = controller.acce_hist[coupled_index]
    ∂Ω_f_hist = controller.∂Ω_f_hist[coupled_index]

    # Interpolate or use fallback
    use_predictor =
        bc isa SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition &&
        controller.use_interface_predictor &&
        controller.iteration_number == 0 &&
        !isempty(controller.predictor_disp[coupled_index])
    if !isempty(time_hist)
        # Piecewise-linear interpolation of the partner trajectory between
        # the stop-start anchor and the partner's substep snapshots. For a
        # constrained pair whose coarse-step side is the Dirichlet side, the
        # fine-step Neumann side receives the coarse side's d'Alembert reaction
        # interpolated linearly between the two end values of the window, and
        # the coarse side imposes the fine side's velocity at the window end:
        # the r = 1 multirate scheme of Connors, Owen, Kuberry, and Bochev
        # (2024, Eqs. (53)-(55), (89)-(90), (104)-(105)), which is also that of
        # Prakash and Hjelmstad (2004, Eqs. (25), (40)). Its interface terms
        # cancel at the fixed point (Connors et al. Eq. (103)), so the
        # subcycled exchange conserves the rate invariant Ẽ without a direct interface solve.
        # The cancellation needs the window-start reaction to be the previous
        # window's converged one, which holds because restore_stop_state runs
        # before subcycle pushes the anchor snapshot, and loads on the
        # Dirichlet interface rows that are constant within the window
        # (dalembert_inertia_minus_loads reads their current value).
        interp_disp = interpolate(time_hist, disp_hist, time)
        interp_velo = interpolate(time_hist, velo_hist, time)
        interp_acce = interpolate(time_hist, acce_hist, time)
        interp_∂Ω_f = interpolate(time_hist, ∂Ω_f_hist, time)
    elseif use_predictor
        interp_disp = controller.predictor_disp[coupled_index]
        interp_velo = controller.predictor_velo[coupled_index]
        interp_acce = controller.predictor_acce[coupled_index]
        interp_∂Ω_f = if !isempty(controller.predictor_∂Ω_f[coupled_index])
            controller.predictor_∂Ω_f[coupled_index]
        else
            controller.stop_∂Ω_f[coupled_index]
        end
    elseif !isempty(controller.stop_disp[coupled_index])
        interp_disp = controller.stop_disp[coupled_index]
        interp_velo = controller.stop_velo[coupled_index]
        interp_acce = controller.stop_acce[coupled_index]
        interp_∂Ω_f = controller.stop_∂Ω_f[coupled_index]
    else
        interp_disp = zeros(coupled_num_dofs)
        interp_velo = zeros(coupled_num_dofs)
        interp_acce = zeros(coupled_num_dofs)
        interp_∂Ω_f = zeros(size(saved_∂Ω_f))
    end

    return interp_disp, interp_velo, interp_acce, interp_∂Ω_f
end

# Evaluate f() with the partner of `bc` set to its state at time t, as
# apply_bc sets it for an unrelaxed exchange, and restore the partner state.
function with_partner_state(f, bc::SolidMechanicsSchwarzBoundaryCondition, time::Float64)
    coupled_subsim = coupled_subsim_of(bc)
    coupled_integrator = coupled_subsim.integrator
    coupled_model = coupled_subsim.model
    saved_disp = copy(coupled_integrator.displacement)
    saved_velo = copy(coupled_integrator.velocity)
    saved_acce = copy(coupled_integrator.acceleration)
    saved_∂Ω_f = copy(get_internal_force(coupled_model))
    interp_disp, interp_velo, interp_acce, interp_∂Ω_f = partner_fields(bc, time)
    set_internal_force!(coupled_model, interp_∂Ω_f)
    coupled_integrator.displacement .= interp_disp
    coupled_integrator.velocity .= interp_velo
    coupled_integrator.acceleration .= interp_acce
    if coupled_model isa RomModel
        reconstruct_fom_fields!(coupled_integrator, coupled_subsim.solver, coupled_model)
    end
    result = f()
    coupled_integrator.displacement .= saved_disp
    coupled_integrator.velocity .= saved_velo
    coupled_integrator.acceleration .= saved_acce
    set_internal_force!(coupled_model, saved_∂Ω_f)
    if coupled_model isa RomModel
        reconstruct_fom_fields!(coupled_integrator, coupled_subsim.solver, coupled_model)
    end
    return result
end

function apply_bc(model::Model, bc::SolidMechanicsSchwarzBoundaryCondition)
    # A direct pair leaves its interface rows free and unloaded during the step;
    # direct_interface_stop! applies the interface force afterwards.
    is_direct_dn(bc) && return nothing
    parent_sim = bc.parent
    controller = parent_sim.controller

    # Skip application if contact is inactive
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact
        return nothing
    end

    coupled_subsim = coupled_subsim_of(bc)
    coupled_integrator = coupled_subsim.integrator
    coupled_model = coupled_subsim.model

    # Save current state (copy data, not reference — aliased integrators share memory with model)
    # For coupled RomModel, these are reduced states
    saved_disp = copy(coupled_integrator.displacement)
    saved_velo = copy(coupled_integrator.velocity)
    saved_acce = copy(coupled_integrator.acceleration)

    # Even for RomModel, this is full-dimensional
    saved_∂Ω_f = copy(get_internal_force(coupled_model))

    # Fetch interpolation inputs
    time = model.time
    coupled_index = bc.coupled_handle.id
    coupled_num_dofs = length(coupled_model.free_dofs)

    interp_disp, interp_velo, interp_acce, interp_∂Ω_f = partner_fields(bc, time)

    # Assign interpolated force
    set_internal_force!(coupled_model, interp_∂Ω_f)

    # Apply relaxed update if needed. Only the Dirichlet side of the interface
    # is relaxed: the Neumann side's traction datum is the interpolated partner
    # reaction (set_internal_force! above), which the kinematic blend below
    # never touches, so relaxation state kept there is dead — it would only
    # compute and log a spurious second θ per iteration.
    relaxed_application = is_swappable_dn_schwarz(bc) && bc.is_dirichlet
    if relaxed_application && is_constrained_dn(bc) && interior_substep(bc)
        # An interior substep of a window-end-constrained side imposes no
        # datum. Under Anderson acceleration the window-end reaction that
        # forms its ramp is the second field of the interface's entry, set
        # with the datum from the joint step of this iteration.
        relaxed_application = false
        key = relaxation_key(bc)
        if controller.relaxation_method === :anderson &&
            controller.iteration_number > 0 &&
            aitken_applies(controller, key) &&
            !isempty(relaxation_slots!(controller.lambda_velo, key)) &&
            !isempty(relaxation_slots!(controller.lambda_velo, key)[1])
            if bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
                overlap_matching_window_step!(bc)
            else
                fields = anderson_data!(bc, time)
                length(fields) >= 2 && (bc.window_end_reaction = copy(fields[2]))
            end
        end
    end
    if relaxed_application
        iter = controller.iteration_number
        # Per-substep-time relaxation slot (see relaxation_slot!): the previous
        # iterate must come from the previous Schwarz iteration at the SAME time. A
        # fresh slot (every slot on iteration 0) has no previous iterate and
        # falls back to the interpolated partner state, as before.
        key = relaxation_key(bc)
        slot_k = relaxation_slot!(controller, key, time)
        use_anderson = controller.relaxation_method === :anderson && aitken_applies(controller, key)
        disp_slots = relaxation_slots!(controller.lambda_disp, key)
        velo_slots = relaxation_slots!(controller.lambda_velo, key)
        acce_slots = relaxation_slots!(controller.lambda_acce, key)
        ensure_slot!(disp_slots, slot_k)
        ensure_slot!(velo_slots, slot_k)
        ensure_slot!(acce_slots, slot_k)
        λ_u_stored = disp_slots[slot_k]
        λ_u_prev = isempty(λ_u_stored) ? interp_disp : λ_u_stored
        λ_v_stored = velo_slots[slot_k]
        λ_v_prev = isempty(λ_v_stored) ? interp_velo : λ_v_stored
        λ_a_stored = acce_slots[slot_k]
        λ_a_prev = isempty(λ_a_stored) ? interp_acce : λ_a_stored

        if is_constrained_dn(bc)
            # Constrained exchange: relax the constrained datum only. The other
            # two fields are derived from it by the receiver's Newmark relations
            # in coupling_weak_dbc, after relaxation, so they are passed through
            # unrelaxed here and never read.
            relax_velocity = bc.constraint == :velocity
            interp_c = relax_velocity ? interp_velo : interp_disp
            c_slots = relax_velocity ? velo_slots : disp_slots
            λ_c_prev = relax_velocity ? λ_v_prev : λ_u_prev
            # The Aitken factor is formed from the datum the Dirichlet side
            # receives, the projected interface trace Π_D q_N, not from the
            # partner's whole field. The interior values of the partner field do
            # not enter the exchange and follow a map of zero gain in the relaxed
            # field, so on the whole field the factor tends to 1 as soon as they
            # dominate the residual, which leaves the interface iteration, of gain
            # near -1, without contraction (measured on the conforming beam).
            trace_c = constrained_partner_trace(bc, interp_c)
            trace_prev = constrained_partner_trace(bc, λ_c_prev)
            fresh_slot = isempty(relax_velocity ? λ_v_stored : λ_u_stored)
            if use_anderson && fresh_slot
                # A fresh slot seeds the previous iterate with the incoming datum,
                # so its residual is exactly zero; entered into the history, that
                # pair makes the least-squares coefficient cancel the next update.
                c_slots[slot_k] = copy(interp_c)
            elseif use_anderson && bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
                # Both sides of the overlap matching coupling are relaxed, and
                # each is applied after its partner has solved in the iteration
                # (Gauss-Seidel), so each side takes its own Anderson step at its
                # application, on the partner data of the current iteration; the
                # stacked step would read the second side's partner data before
                # that partner solved. A window-end-constrained side took its
                # joint step at its first interior substep, and its slot holds
                # the datum already.
                if get(overlap_window_steps(controller), key, -1) != iter
                    history = get!(AndersonHistory, overlap_anderson_histories(controller), key)
                    footprint = anderson_footprint(bc)
                    relaxed = anderson_step!(
                        history,
                        restrict(λ_c_prev, footprint),
                        restrict(interp_c, footprint),
                        trace_c .- trace_prev,
                        controller.relaxation_parameter,
                        controller.anderson_depth,
                    )
                    c_slots[slot_k] = expand(relaxed, interp_c, footprint)
                    record_anderson_memory!(controller)
                end
            elseif use_anderson
                c_slots[slot_k] = expand(anderson_data!(bc, time)[1], interp_c, anderson_footprint(bc))
            else
                θ = relaxation_theta!(bc, key, slot_k, iter, trace_c, trace_prev)
                frozen_relaxation_update!(controller, θ)
                c_slots[slot_k] = θ * interp_c + (1 - θ) * λ_c_prev
            end
            if relax_velocity
                disp_slots[slot_k] = interp_disp
            else
                velo_slots[slot_k] = interp_velo
            end
            acce_slots[slot_k] = interp_acce
        elseif use_anderson
            # Anderson acceleration of the three imposed kinematic fields. The
            # least-squares coefficients come from the projected interface trace
            # of the displacement residual, stacked with the traces of the other
            # relaxed interfaces (anderson_data!), and the same coefficients update the
            # displacement, the velocity, and the acceleration, each with its own
            # iterates on the footprint of the interface (anderson_footprint).
            # The velocity and acceleration of a static partner are constant and
            # pass through unrelaxed (relaxes_displacement_only). A fresh slot
            # has a residual of exactly zero and is not entered into the history.
            if isempty(λ_u_stored)
                disp_slots[slot_k] = copy(interp_disp)
                velo_slots[slot_k] = copy(interp_velo)
                acce_slots[slot_k] = copy(interp_acce)
            else
                fields = anderson_data!(bc, time)
                footprint = anderson_footprint(bc)
                disp_slots[slot_k] = expand(fields[1], interp_disp, footprint)
                velo_slots[slot_k] = length(fields) ≥ 2 ? expand(fields[2], interp_velo, footprint) : copy(interp_velo)
                acce_slots[slot_k] = length(fields) ≥ 3 ? expand(fields[3], interp_acce, footprint) : copy(interp_acce)
            end
        elseif controller.relaxation_method === :aitken_secant
            # Relax the single d-form interface unknown (displacement); recover
            # velocity and acceleration consistently from it (see functions above).
            θ = relaxation_theta!(bc, key, slot_k, iter, interp_disp, λ_u_prev)
            frozen_relaxation_update!(controller, θ)
            disp_slots[slot_k] = θ * interp_disp + (1 - θ) * λ_u_prev
            recover_interface_kinematics!(controller, key, slot_k, coupled_integrator, interp_velo, interp_acce)
        else
            θ = relaxation_theta!(bc, key, slot_k, iter, interp_disp, λ_u_prev)
            frozen_relaxation_update!(controller, θ)

            disp_slots[slot_k] = θ * interp_disp + (1 - θ) * λ_u_prev
            velo_slots[slot_k] = θ * interp_velo + (1 - θ) * λ_v_prev
            acce_slots[slot_k] = θ * interp_acce + (1 - θ) * λ_a_prev
        end

        coupled_integrator.displacement .= disp_slots[slot_k]
        coupled_integrator.velocity .= velo_slots[slot_k]
        coupled_integrator.acceleration .= acce_slots[slot_k]
    else
        coupled_integrator.displacement .= interp_disp
        coupled_integrator.velocity .= interp_velo
        coupled_integrator.acceleration .= interp_acce
    end

    # For ROM coupled subdomains, reconstruct FOM displacement from reduced state before reading it
    if coupled_subsim.model isa RomModel
        reconstruct_fom_fields!(coupled_integrator, coupled_subsim.solver, coupled_subsim.model)
    end
    apply_bc_detail(model, bc)

    # Restore previous state (in-place to keep alias intact)
    coupled_integrator.displacement .= saved_disp
    coupled_integrator.velocity .= saved_velo
    coupled_integrator.acceleration .= saved_acce
    set_internal_force!(coupled_model, saved_∂Ω_f)
    if coupled_subsim.model isa RomModel
        reconstruct_fom_fields!(coupled_integrator, coupled_subsim.solver, coupled_subsim.model)
    end
    return nothing
end

function transfer_normal_component(source::Vector{Float64}, target::Vector{Float64}, normal::Vector{Float64})
    normal_projection = normal * normal'
    tangent_projection = I(length(normal)) - normal_projection
    return tangent_projection * target + normal_projection * source
end

# Largest departure of |n_x| from 1 accepted for a frictionless contact normal.
const FRICTIONLESS_NORMAL_TOLERANCE = 1.0e-6

function contact_weak_dbc(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    nodal_curr, _, nodal_velo, nodal_acce = get_dst_curr_disp_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    for (i_local, i_global) in enumerate(global_from_local_map)
        normal = normals[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        if bc.friction_type == 0
            # The constrained degree of freedom is the Cartesian x component, so
            # the normal must lie along x; a contact surface of another
            # orientation would need the constraint in a rotated nodal frame.
            if abs(normal[1]) < 1.0 - FRICTIONLESS_NORMAL_TOLERANCE
                norma_abort(
                    "Frictionless Schwarz contact on side set $(bc.name) has a normal " *
                    "($(normal[1]), $(normal[2]), $(normal[3])) at node $(i_global) that is not along x. " *
                    "Frictionless contact constrains the x component only and supports contact " *
                    "surfaces normal to x; use friction type: tied, or orient the contact along x.",
                )
            end
            @inbounds model.displacement[:, i_global] =
                transfer_normal_component(
                    nodal_curr[:, i_local], model.reference[:, i_global] + model.displacement[:, i_global], normal
                ) - model.reference[:, i_global]
            @inbounds model.velocity[:, i_global] = transfer_normal_component(
                nodal_velo[:, i_local], model.velocity[:, i_global], normal
            )
            @inbounds model.acceleration[:, i_global] = transfer_normal_component(
                nodal_acce[:, i_local], model.acceleration[:, i_global], normal
            )
            model.free_dofs[[3 * i_global - 2]] .= false
            model.free_dofs[[3 * i_global - 1]] .= true
            model.free_dofs[[3 * i_global]] .= true
        elseif bc.friction_type == 1
            @inbounds model.displacement[:, i_global] = nodal_curr[:, i_local] - model.reference[:, i_global]
            @inbounds model.velocity[:, i_global] = nodal_velo[:, i_local]
            @inbounds model.acceleration[:, i_global] = nodal_acce[:, i_local]
            model.free_dofs[global_range] .= false
        else
            norma_abort("Unknown or not implemented friction type.")
        end
    end
end

function apply_naive_stabilized_bcs(subsim::SingleDomainSimulation)
    bcs = subsim.model.boundary_conditions
    for bc in bcs
        if bc isa SolidMechanicsContactSchwarzBoundaryCondition
            unique_node_indices = unique(bc.side_set_node_indices)
            for node_index in unique_node_indices
                subsim.model.acceleration[:, node_index] .= 0.0
            end
        end
    end
    return nothing
end

function contact_weak_nbc(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    friction_type = bc.friction_type
    nodal_force = get_dst_force(bc)
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        local_range = (3 * (i_local - 1) + 1):(3 * i_local)
        normal = normals[:, i_local]
        node_force = nodal_force[local_range]
        if friction_type == 0
            target = model.boundary_force[global_range]
            eff_node_force = transfer_normal_component(node_force, target, normal)
        else
            eff_node_force = node_force
        end
        @inbounds model.boundary_force[global_range] += eff_node_force
    end
end

function extract_local_vector(global_vector::Vector{Float64}, global_from_local_map::Vector{Int64}, dim::Int64)
    num_local_nodes = length(global_from_local_map)
    local_vector = Vector{Float64}(undef, dim * num_local_nodes)
    for (i_local, i_global) in enumerate(global_from_local_map)
        global_range = (dim * (i_global - 1) + 1):(dim * i_global)
        local_range = (dim * (i_local - 1) + 1):(dim * i_local)
        @inbounds local_vector[local_range] = global_vector[global_range]
    end
    return local_vector
end

function extract_local_vector(bc::SolidMechanicsSchwarzBoundaryCondition, global_vector::Vector{Float64}, dim::Int64)
    global_from_local_map = bc.global_from_local_map
    return extract_local_vector(global_vector, global_from_local_map, dim)
end

# The Dirichlet and Neumann projectors of one interface use the same two
# rectangular projection matrices; a cache shared between the two calls keeps
# each from being computed twice.
const RectangularProjectionCache = Dict{Symbol,Tuple{Matrix{Float64},Matrix{Float64}}}

function cached_rectangular_projection(
    cache::RectangularProjectionCache,
    which::Symbol,
    dst_fom::SolidMechanics,
    dst_bc::SolidMechanicsSchwarzBoundaryCondition,
    src_fom::SolidMechanics,
    src_bc::SolidMechanicsSchwarzBoundaryCondition,
)
    return get!(cache, which) do
        if which == :destination_integrated
            get_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
        else
            get_src_integrated_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
        end
    end
end

function compute_neumann_projector(
    dst_model::Model,
    dst_bc::SolidMechanicsSchwarzBoundaryCondition;
    cache::RectangularProjectionCache=RectangularProjectionCache(),
)
    src_model = coupled_subsim_of(dst_bc).model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_fom = src_model isa RomModel ? src_model.fom_model : src_model
    dst_fom = dst_model isa RomModel ? dst_model.fom_model : dst_model
    H = get_square_projection_matrix(src_fom, src_bc)
    # Prefer the destination-integrated L: where the destination facets
    # resolve the source trace mesh (nested refinement), it is exactly
    # integrated, hence also conservative. Conservation of the transferred
    # force totals is the certificate of that exactness — its column sums
    # must match the source mass row sums — and when it fails (coarser
    # destination facets, non-nested interfaces), fall back to the
    # source-integrated L, which is conservative by construction.
    L, _ = cached_rectangular_projection(cache, :destination_integrated, dst_fom, dst_bc, src_fom, src_bc)
    src_lumped = H * ones(size(H, 2))
    conservation_error = maximum(abs.(vec(sum(L; dims=1)) - src_lumped)) / maximum(abs.(src_lumped))
    if conservation_error > 1.0e-10
        # The source-integrated L comes with the source mass matrix of its
        # own quadrature, whose row sums its column sums match exactly.
        L, H = cached_rectangular_projection(cache, :source_integrated, dst_fom, dst_bc, src_fom, src_bc)
    end
    dst_bc.neumann_projector = L * (H \ I)
    return nothing
end

function compute_dirichlet_projector(
    dst_model::Model,
    dst_bc::SolidMechanicsSchwarzBoundaryCondition;
    cache::RectangularProjectionCache=RectangularProjectionCache(),
)
    src_model = coupled_subsim_of(dst_bc).model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_fom = src_model isa RomModel ? src_model.fom_model : src_model
    dst_fom = dst_model isa RomModel ? dst_model.fom_model : dst_model
    W = get_square_projection_matrix(dst_fom, dst_bc)
    # Prefer the source-integrated L: on interfaces where the source trace
    # mesh resolves the destination facets (nested refinement), it is exactly
    # integrated, so the transferred kinematics are exact. Its partition of
    # unity is the certificate of that exactness — verify it numerically, and
    # when it fails (non-nested facets, partial coverage), fall back to the
    # destination-integrated L, whose Dirichlet projector reproduces constants
    # for any quadrature by construction.
    L, _ = cached_rectangular_projection(cache, :source_integrated, dst_fom, dst_bc, src_fom, src_bc)
    P = (W \ I) * L
    pu_error = maximum(abs.(P * ones(size(P, 2)) .- 1.0))
    if pu_error > 1.0e-10
        # The destination-integrated L comes with the destination mass matrix
        # of its own quadrature, whose row sums equal those of L exactly.
        L, W = cached_rectangular_projection(cache, :destination_integrated, dst_fom, dst_bc, src_fom, src_bc)
        P = (W \ I) * L
    end
    dst_bc.dirichlet_projector = P
    return nothing
end

function get_dst_force(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = coupled_subsim_of(dst_bc)
    src_model = src_sim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_dofs = interface_dofs(src_bc)
    src_force = -get_internal_force(src_model)[src_dofs]
    if is_constrained_dn(dst_bc)
        # Constrained DN exchange: the Neumann side receives the d'Alembert
        # reaction of the Dirichlet side, -(M a + f_int - f_body - f_boundary)
        # on its interface rows, so that the two interface rows add to the row
        # of the undecomposed problem. f_boundary is the Dirichlet side's own
        # applied surface load on the interface nodes, if any.
        src_force .-= dalembert_inertia_minus_loads(get_fom_model(src_sim), src_dofs)
    end
    neumann_projector = dst_bc.neumann_projector
    num_dst_nodes = size(neumann_projector, 1)
    dst_force = zeros(3 * num_dst_nodes)
    dst_force[1:3:end] = neumann_projector * src_force[1:3:end]
    dst_force[2:3:end] = neumann_projector * src_force[2:3:end]
    dst_force[3:3:end] = neumann_projector * src_force[3:3:end]
    return dst_force
end

# M a - f_body - f_boundary of a subdomain, with the mass matrix its integrator
# uses: the consistent mass for Newmark, the lumped mass for central difference.
# Added to the internal force it gives the d'Alembert residual whose interface
# rows are the reaction of the coupling constraint.
function dalembert_inertia_minus_loads(model::SolidMechanics, a::AbstractVector{Float64}=vec(model.acceleration))
    inertial_force = if size(model.mass, 1) == length(a)
        model.mass * a
    elseif length(model.lumped_mass) == length(a)
        model.lumped_mass .* a
    else
        # Before the first evaluation no mass is assembled and the acceleration
        # is still zero, so the inertial term vanishes.
        zeros(length(a))
    end
    if length(model.body_force) == length(a)
        inertial_force = inertial_force - model.body_force
    end
    if length(model.boundary_force) == length(a)
        inertial_force = inertial_force - model.boundary_force
    end
    return inertial_force
end

# Degrees of freedom of `nodes`, interleaved by component, in the order of
# extract_local_vector; `interface_dofs` gives those of the interface nodes
# of `bc`.
function node_dofs(nodes::Vector{Int64})
    dofs = Vector{Int64}(undef, 3 * length(nodes))
    for (i_local, i_global) in enumerate(nodes)
        @inbounds for c in 1:3
            dofs[3 * (i_local - 1) + c] = 3 * (i_global - 1) + c
        end
    end
    return dofs
end

interface_dofs(bc::SolidMechanicsSchwarzBoundaryCondition) = node_dofs(bc.global_from_local_map)

# Rows `dofs` of M x with the mass of the model, the consistent mass or the
# lumped one in the order of preference given, or nothing when none is
# assembled. The consistent mass is symmetric, so its row is its column, and
# the rows of the interface cost a few hundred dot products where the full
# product cost 1.1 ms on 46 875 degrees of freedom.
function mass_rows_product(model::SolidMechanics, dofs::Vector{Int64}, x::AbstractVector{Float64}, lumped_first::Bool)
    M = model.mass
    lumped = model.lumped_mass
    has_consistent = size(M, 1) == length(x)
    has_lumped = length(lumped) == length(x)
    if has_lumped && (lumped_first || !has_consistent)
        return [lumped[i] * x[i] for i in dofs]
    end
    has_consistent || return nothing
    colptr, rowval, nzval = M.colptr, M.rowval, M.nzval
    p = Vector{Float64}(undef, length(dofs))
    @inbounds for (k, i) in enumerate(dofs)
        s = 0.0
        for q in colptr[i]:(colptr[i + 1] - 1)
            s += nzval[q] * x[rowval[q]]
        end
        p[k] = s
    end
    return p
end

# dalembert_inertia_minus_loads on the rows `dofs`, without the full product.
function dalembert_inertia_minus_loads(
    model::SolidMechanics, dofs::Vector{Int64}, a::AbstractVector{Float64}=vec(model.acceleration)
)
    inertial = mass_rows_product(model, dofs, a, false)
    r = inertial === nothing ? zeros(length(dofs)) : inertial
    if length(model.body_force) == length(a)
        @inbounds for (k, i) in enumerate(dofs)
            r[k] -= model.body_force[i]
        end
    end
    if length(model.boundary_force) == length(a)
        @inbounds for (k, i) in enumerate(dofs)
            r[k] -= model.boundary_force[i]
        end
    end
    return r
end

# Partner nodes that the columns of the Dirichlet projector refer to: the
# partner's interface nodes of a nonoverlap or contact pair, every node of the
# partner for the overlap matching coupling, whose projector acts on the
# partner's volume shape functions.
partner_node_map(bc::SolidMechanicsSchwarzBoundaryCondition) = coupled_bc_of(bc).global_from_local_map
function partner_node_map(bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition)
    return collect(1:size(get_fom_model(coupled_subsim_of(bc)).reference, 2))
end

function get_dst_curr_disp_velo_acce(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = coupled_subsim_of(dst_bc)
    src_model = src_sim.model
    src_fom = src_model isa RomModel ? src_model.fom_model : src_model
    src_global_from_local_map = partner_node_map(dst_bc)
    num_src_nodes = length(src_global_from_local_map)
    src_curr = zeros(3, num_src_nodes)
    src_refe = zeros(3, num_src_nodes)
    src_velo = zeros(3, num_src_nodes)
    src_acce = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_curr[:, i_local] = src_fom.reference[:, i_global] + src_fom.displacement[:, i_global]
        src_refe[:, i_local] = src_fom.reference[:, i_global]
        src_velo[:, i_local] = src_fom.velocity[:, i_global]
        src_acce[:, i_local] = src_fom.acceleration[:, i_global]
    end
    dirichlet_projector = dst_bc.dirichlet_projector
    num_dst_nodes = size(dirichlet_projector, 1)
    dst_curr = zeros(3, num_dst_nodes)
    dst_disp = zeros(3, num_dst_nodes)
    dst_velo = zeros(3, num_dst_nodes)
    dst_acce = zeros(3, num_dst_nodes)
    for i in 1:3
        dst_curr[i, :] = dirichlet_projector * src_curr[i, :]
        dst_disp[i, :] = dirichlet_projector * (src_curr[i, :] - src_refe[i, :])
        dst_velo[i, :] = dirichlet_projector * src_velo[i, :]
        dst_acce[i, :] = dirichlet_projector * src_acce[i, :]
    end
    return dst_curr, dst_disp, dst_velo, dst_acce
end

function set_id_from_name(name::String, mesh::ExodusDatabase, ::Type{T}) where {T}
    names = Exodus.read_names(mesh, T)
    idx = findfirst(==(name), names)
    if idx === nothing
        type_str = if T === NodeSet
            "node set"
        elseif T === SideSet
            "side set"
        else
            "block"
        end
        norma_abort("$type_str $name cannot be found in mesh")
    end
    return Int64(Exodus.read_ids(mesh, T)[idx])
end

node_set_id_from_name(name::String, mesh::ExodusDatabase) = set_id_from_name(name, mesh, NodeSet)
side_set_id_from_name(name::String, mesh::ExodusDatabase) = set_id_from_name(name, mesh, SideSet)
block_id_from_name(name::String, mesh::ExodusDatabase) = set_id_from_name(name, mesh, Block)

function component_offset_from_string(name::String)
    offset = 0
    if name == "x"
        offset = 1
    elseif name == "y"
        offset = 2
    elseif name == "z"
        offset = 3
    else
        norma_abort("invalid component name $name")
    end
    return offset
end

# Boundary conditions are applied in creation order, and a node can belong to
# more than one of them: a Schwarz interface meets the clamped or symmetry
# faces along its edges. The last condition applied to such a node wins. The
# input dictionary iterates in hash order, which changed between Julia 1.12
# and 1.13 and flipped which condition won on those nodes. The types are
# therefore created in a fixed order: prescribed (Dirichlet) conditions last,
# so that they hold exactly on shared nodes, with the others before them in
# alphabetical order. Entries of one type keep their order in the input file.
bc_type_rank(bc_type::AbstractString) = (endswith(bc_type, "Dirichlet"), bc_type)

function _create_bcs(subsim::SingleDomainSimulation)
    boundary_conditions = Vector{BoundaryCondition}()
    params = subsim.params
    if haskey(params, "boundary conditions") == false
        return boundary_conditions
    end
    input_mesh = params["input_mesh"]
    bc_params = params["boundary conditions"]
    for bc_type in sort(collect(keys(bc_params)); by=bc_type_rank)
        bc_type_params = bc_params[bc_type]
        for bc_setting_params in bc_type_params
            if bc_type == "Dirichlet"
                # Same "Dirichlet" syntax for both: a "node set" entry applies on
                # a node set, a "side set" entry applies on the nodes of a side set.
                if haskey(bc_setting_params, "side set")
                    boundary_condition = SolidMechanicsSideSetDirichletBoundaryCondition(input_mesh, bc_setting_params)
                elseif haskey(bc_setting_params, "node set")
                    boundary_condition = SolidMechanicsDirichletBoundaryCondition(input_mesh, bc_setting_params)
                else
                    norma_abort(
                        "A Dirichlet boundary condition requires either a \"node set\" or a \"side set\" entry."
                    )
                end
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "OpInf Dirichlet"
                boundary_condition = SolidMechanicsOpInfDirichletBC(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann"
                boundary_condition = SolidMechanicsNeumannBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann pressure"
                boundary_condition = SolidMechanicsNeumannPressureBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Surface"
                boundary_condition = SolidMechanicsSurfaceBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Robin"
                boundary_condition = SolidMechanicsRobinBoundaryCondition(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Periodic"
                if length(bc_type_params) > 1
                    norma_abort("A model has at most one Periodic boundary condition, with up to three node set pairs.")
                end
                boundary_condition = SolidMechanicsPeriodicBoundaryCondition(subsim, input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz contact"
                sim = subsim.parent
                sim.controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SolidMechanicsContactSchwarzBoundaryCondition(
                    subsim, coupled_subsim, input_mesh, bc_setting_params
                )
                push!(boundary_conditions, boundary_condition)
            elseif bc_type in
                ("Schwarz overlap", "Schwarz overlap matching", "Schwarz DN nonoverlap", "Schwarz RR nonoverlap")
                sim = subsim.parent
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SMCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "OpInf Schwarz overlap"
                sim = subsim.parent
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SMOpInfCouplingSchwarzBC(
                    subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params
                )
                push!(boundary_conditions, boundary_condition)
            else
                norma_abort("Unknown boundary condition type : $bc_type")
            end
        end
    end
    return boundary_conditions
end

function apply_bcs(model::SolidMechanics)
    model.boundary_force .= 0.0
    model.free_dofs .= true
    for boundary_condition in model.boundary_conditions
        apply_bc(model, boundary_condition)
    end
    # A periodic cell sets its corners and dependent nodes from the final
    # displacement of its origin. A no-op without a Periodic BC.
    finalize_periodic_bcs!(model)
    return nothing
end

function assign_velocity!(
    velocity::Matrix{Float64}, offset::Int64, node_index::Int32, velo_val::Float64, context::String
)
    current_val = velocity[offset, node_index]
    velocity_already_defined = !(current_val ≈ 0.0)
    dissimilar_velocities = !(current_val ≈ velo_val)
    if velocity_already_defined && dissimilar_velocities
        norma_abortf(
            "Inconsistent velocity initial conditions for node %d: " *
            "attempted to assign velocity %s (v = (%.4e, %.4e, %.4e)), " *
            "which conflicts with an already assigned value (v = (%.4e, %.4e, %.4e)).",
            node_index,
            context,
            velo_val[1],
            velo_val[2],
            velo_val[3],
            current_val[1],
            current_val[2],
            current_val[3],
        )
    else
        velocity[offset, node_index] = velo_val
    end
    return nothing
end

function apply_ics(params::Parameters, model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    if haskey(params, "initial conditions") == false
        return nothing
    end
    input_mesh = params["input_mesh"]
    ic_params = params["initial conditions"]
    for (ic_type, ic_type_params) in ic_params
        for ic in ic_type_params
            node_set_name = ic["node set"]
            expression = ic["function"]
            component = ic["component"]
            offset = component_offset_from_string(component)
            node_set_id = node_set_id_from_name(node_set_name, input_mesh)
            node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
            # expression is an arbitrary function of t, x, y, z in the input file
            # Compile to a Float64-valued callable once, then evaluate per node.
            if ic_type == "displacement"
                disp_num = eval(Meta.parse(expression))
                velo_num = expand_derivatives(D(disp_num))
                disp_fn = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
                velo_fn = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
            elseif ic_type == "velocity"
                disp_fn = nothing
                velo_num = eval(Meta.parse(expression))
                velo_fn = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
            else
                norma_abort(
                    "Invalid initial condition type: '$ic_type'. Supported types are: displacement or velocity."
                )
            end
            for node_index in node_set_node_indices
                txzy = (
                    model.time,
                    model.reference[1, node_index],
                    model.reference[2, node_index],
                    model.reference[3, node_index],
                )
                disp_val = disp_fn === nothing ? 0.0 : Float64(disp_fn(txzy))
                velo_val = Float64(velo_fn(txzy))
                if ic_type == "displacement"
                    model.displacement[offset, node_index] = disp_val
                    non_zero_velocity = !(velo_val ≈ 0.0)
                    if non_zero_velocity
                        assign_velocity!(model.velocity, offset, node_index, velo_val, "derived from displacement")
                    end
                end
                if ic_type == "velocity"
                    assign_velocity!(model.velocity, offset, node_index, velo_val, "directly from velocity IC")
                end
            end
        end
    end
    return nothing
end

function pair_schwarz_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        model = subsim.model
        bcs = model.boundary_conditions
        for (bc_index, bc) in enumerate(bcs)
            pair_bc(bc, bc_index)
        end
    end
end

function pair_bc(_::SolidMechanicsRegularBoundaryCondition, _::Int64) end

function pair_bc(bc::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition, bc_index::Int64)
    # The partner is the overlap matching condition of the coupled subdomain
    # that couples back to this one.
    coupled_bcs = coupled_subsim_of(bc).model.boundary_conditions
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        coupled_bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition || continue
        coupled_bc.coupled_handle.id == bc.self_handle.id || continue
        bc.coupled_bc_index = coupled_bc_index
        coupled_bc.coupled_bc_index = bc_index
        return nothing
    end
    return norma_abort(
        "`Schwarz overlap matching` on side set \"$(bc.name)\" has no partner condition of the same type " *
        "on the coupled subdomain; both sides of the overlap must set it.",
    )
end

function pair_bc(bc::SolidMechanicsSchwarzBoundaryCondition, bc_index::Int64)
    if bc isa SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
        return nothing
    end
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = coupled_subsim_of(bc).model
    coupled_bcs = coupled_model.boundary_conditions
    # The partner is the Schwarz condition on the named side set; other
    # conditions on the same side set (a Neumann load, for example) share its
    # name and are not partners.
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name && coupled_bc isa SolidMechanicsSchwarzBoundaryCondition
            if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition &&
                coupled_bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                if bc.is_dirichlet == coupled_bc.is_dirichlet
                    norma_abort(
                        "Nonoverlap Schwarz BCs must specify different default BC types (Dirichlet vs Neumann)."
                    )
                end
            end
            coupled_bc.is_dirichlet = !bc.is_dirichlet
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
        end
    end
    return nothing
end
