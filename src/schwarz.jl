# Schwarz coupling boundary conditions for multi-domain simulations.
#
# Includes: overlap, non-overlap (Dirichlet-Neumann), Robin-Robin,
# and contact Schwarz BCs. Projectors, time history interpolation,
# and interface force/displacement transfer.

using LinearAlgebra: dot

# ---------------------------------------------------------------------------
# Impedance-matching Robin Schwarz: t + Z u̇ + α W u = g
# ---------------------------------------------------------------------------
#
# g = -t_src_projected + Z * u̇_src_projected + α * W * u_src_projected
#
# The impedance term Z * u̇ absorbs outgoing waves at the interface,
# preventing reflections that cause energy growth with mixed integrators.

get_fom_model(sim::Simulation) = sim.model isa RomModel ? sim.model.fom_model : sim.model

# Resolve the partner (coupled) subsim of a Schwarz BC through its parent via
# the stable handle, so that swapping a subsim in its parent slot is transparent
# to every BC.
coupled_subsim_of(bc::SolidMechanicsSchwarzBoundaryCondition) = bc.parent.subsims[bc.coupled_handle.id]
self_subsim_of(bc::SolidMechanicsSchwarzBoundaryCondition)    = bc.parent.subsims[bc.self_handle.id]

# Floor on ‖δ‖² guarding the division when successive residuals are essentially
# identical (a stale residual direction). This is numerical safety only; the
# Aitken factor itself is left unclamped so it is free to take the large or
# negative excursions that accelerate convergence.
const AITKEN_DELTA_SQ_FLOOR = 1.0e-20

# Resolve the relaxation time slot for coupling pair `pair` at substep time t,
# appending a new slot when t is not yet keyed. The relaxation state must be
# compared across Schwarz sweeps AT THE SAME SUBSTEP TIME: with a windowed
# controller stop the relaxed side applies its BC once per substep, and a
# single per-pair state would blend iterates across time — a causal low-pass
# on the exchanged data that shifts the converged fixed point. Substep times
# repeat bitwise across sweeps (subcycle() restores the nominal step every
# pass), so the equality test hits; the isapprox fallback and the fresh-slot
# path only engage if the substep grid is perturbed mid-stop (e.g. adaptive
# stepping), where the new slot restarts that time's relaxation from scratch.
function relaxation_slot!(controller::MultiDomainTimeController, pair::Int, t::Float64)
    times = controller.lambda_time[pair]
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
# windowed stop the sweep map couples all time slots, and every Aitken policy
# measured on the 10 ms cantilever benchmark fails or loses there: per-slot
# θs take unclamped excursions that hand the solver a divergent interface
# force (deaths at 3.9/5.7 ms where fixed θ and θ = 1 ran clean), and pooling
# the residual inner products over the sweep's slots (waveform Aitken,
# Irons–Tuck base frozen per sweep) survived on the implicit pair but cost
# 55 sweeps/stop against 47.5 for fixed θ = 0.5 and 19.1 for θ = 1, while
# locking θ persistently negative on the explicit pair (dead at 0.33 ms).
# Windowed stops therefore use the configured relaxation parameter; for the
# dashpot-stabilized impedance exchange θ = 1 is the measured optimum there.
function aitken_applies(controller::MultiDomainTimeController, pair::Int)
    return length(controller.lambda_time[pair]) <= 1
end

# Returns the relaxation factor θ applied to interp_disp for this Schwarz
# iterate. Fixed mode returns the user-configured constant; Aitken-recursive
# mode uses Irons–Tuck with the previous residual stored on the controller,
# per pair and per substep time slot.
function relaxation_aitken_recursive_theta!(
    controller::MultiDomainTimeController,
    pair::Int,
    slot_k::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    if controller.relaxation_method !== :aitken_recursive || !aitken_applies(controller, pair)
        return controller.relaxation_parameter
    end
    aitken_N0 = controller.aitken_N0
    ensure_slot!(controller.aitken_prev_residual_disp[pair], slot_k)
    ensure_slot!(controller.aitken_theta_disp[pair], slot_k, controller.relaxation_parameter)
    if iter < aitken_N0
        controller.aitken_theta_disp[pair][slot_k] = controller.relaxation_parameter
        controller.aitken_prev_residual_disp[pair][slot_k] = Float64[]
        norma_logf(1, :schwarz, "Aitken-recursive θ[pair=%d, slot=%d, iter=%d] = %.4f", pair, slot_k, iter, 1.0)
        return 1.0
    end
    residual = interp_disp .- lambda_prev
    prev_residual = controller.aitken_prev_residual_disp[pair][slot_k]
    θ_prev = controller.aitken_theta_disp[pair][slot_k]
    θ = controller.relaxation_parameter
    if !isempty(prev_residual) && length(prev_residual) == length(residual)
        δ = residual .- prev_residual
        δ_sq = dot(δ, δ)
        if δ_sq > AITKEN_DELTA_SQ_FLOOR
            θ = -θ_prev * dot(prev_residual, δ) / δ_sq
        else
            θ = θ_prev
        end
    end
    controller.aitken_prev_residual_disp[pair][slot_k] = residual
    controller.aitken_theta_disp[pair][slot_k] = θ
    norma_logf(1, :schwarz, "Aitken-recursive θ[pair=%d, slot=%d, iter=%d] = %.4f", pair, slot_k, iter, θ)
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
    pair::Int,
    slot_k::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    if !aitken_applies(controller, pair)
        return controller.relaxation_parameter
    end
    ensure_slot!(controller.aitken_prev_residual_disp[pair], slot_k)
    ensure_slot!(controller.aitken_prev_lambda_disp[pair], slot_k)
    residual = interp_disp .- lambda_prev                 # r^(n) = T(g^(n)) - g^(n)
    prev_residual = controller.aitken_prev_residual_disp[pair][slot_k]   # r^(n-1)
    prev_lambda = controller.aitken_prev_lambda_disp[pair][slot_k]       # g^(n-1)
    # For iter < N0 use the input ρ^(1) = relaxation_parameter (paper's N0 idea).
    θ = controller.relaxation_parameter
    if iter >= controller.aitken_N0 &&
       !isempty(prev_residual) && length(prev_residual) == length(residual) &&
       !isempty(prev_lambda) && length(prev_lambda) == length(lambda_prev)
        δ = residual .- prev_residual                     # δ^(n)
        d = lambda_prev .- prev_lambda                    # d^(n) = g^(n) - g^(n-1)
        δ_sq = dot(δ, δ)
        if δ_sq > AITKEN_DELTA_SQ_FLOOR
            # Pure Aitken-secant factor, paper eq. (9): no value clamp, so the factor is
            # free to take the large/negative excursions that accelerate (or
            # damp) convergence. The δ_sq floor above is kept only to guard the
            # division when successive residuals are essentially identical.
            θ = -dot(d, δ) / δ_sq
        end
    end
    controller.aitken_prev_residual_disp[pair][slot_k] = residual
    controller.aitken_prev_lambda_disp[pair][slot_k] = copy(lambda_prev)
    norma_logf(1, :schwarz, "Aitken-secant θ[pair=%d, slot=%d, iter=%d] = %.4f", pair, slot_k, iter, θ)
    return θ
end

# Recover interface velocity and acceleration consistent with a relaxed interface
# displacement, instead of relaxing them independently. Norma's implicit Newmark
# is displacement-form (the solver unknown is u; see `correct`, time_integrator.jl),
# so the single interface unknown to relax is the displacement, and v, a follow
# from the integrator's own recovery relations using the stored predictors:
#   a = (u - u_pre) / (β Δt²),   v = v_pre + γ Δt a.
# For integrators without these predictors (e.g. quasi-static, explicit), the
# interpolated kinematics are passed through unrelaxed.
function recover_interface_kinematics!(controller, pair, slot_k, integrator, interp_velo, interp_acce)
    if integrator isa Newmark
        Δt = integrator.time_step
        β = integrator.β
        γ = integrator.γ
        u = controller.lambda_disp[pair][slot_k]
        a = (u .- integrator.disp_pre) ./ (β * Δt * Δt)
        controller.lambda_acce[pair][slot_k] = a
        controller.lambda_velo[pair][slot_k] = integrator.velo_pre .+ (γ * Δt) .* a
    else
        controller.lambda_velo[pair][slot_k] = interp_velo
        controller.lambda_acce[pair][slot_k] = interp_acce
    end
    return nothing
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition)
    Z = bc.impedance
    α = bc.robin_parameter
    W = bc.square_projector
    parent_sim = bc.parent
    controller = parent_sim.controller
    iter = controller.iteration_number
    coupled_index = bc.coupled_handle.id
    this_index = bc.self_handle.id

    # Neumann part: -t_src projected
    neumann_force = get_dst_force(bc)

    # Source displacement and velocity, projected to destination
    src_sim = coupled_subsim_of(bc)
    src_model = get_fom_model(src_sim)
    src_bc = src_model.boundary_conditions[bc.coupled_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
    num_src_nodes = length(src_global_from_local_map)
    src_disp = zeros(3, num_src_nodes)
    src_velo = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_disp[:, i_local] = src_model.displacement[:, i_global]
        src_velo[:, i_local] = src_model.velocity[:, i_global]
    end
    dirichlet_projector = bc.dirichlet_projector
    num_dst_nodes = size(dirichlet_projector, 1)
    dst_disp = zeros(3, num_dst_nodes)
    dst_velo = zeros(3, num_dst_nodes)
    for i in 1:3
        dst_disp[i, :] = dirichlet_projector * src_disp[i, :]
        dst_velo[i, :] = dirichlet_projector * src_velo[i, :]
    end
    global_from_local_map = bc.global_from_local_map
    theta = controller.relaxation_parameter

    # RHS (partner contribution only): g = -t_src + Z*W*u̇_src + α*W*u_src
    # The self terms (Z*W*u̇_self + α*W*u_self) are handled by
    # apply_impedance_bcs_internal_force! called at each Newton iteration.
    if (this_index < coupled_index)  # right subdomain
        for comp in 1:3
            Z_W_vdot = Z * (W * dst_velo[comp, :])
            α_W_u = α * (W * dst_disp[comp, :])
            for (i_local, i_global) in enumerate(global_from_local_map)
                dof_i = 3 * (i_global - 1) + comp
                model.boundary_force[dof_i] += neumann_force[3 * (i_local - 1) + comp] +
                                                Z_W_vdot[i_local] + α_W_u[i_local]
            end
        end
    else  # left subdomain (with relaxation)
        # The relaxation state is per substep time slot (see relaxation_slot!):
        # relaxing against the previous sweep's RHS at the SAME time keeps the
        # windowed exchange waveform-consistent. A fresh slot (every slot on
        # iteration 0, since the state resets at each stop) has no previous
        # iterate and starts from zero, as before.
        slot_k = relaxation_slot!(controller, coupled_index, model.time)
        ensure_slot!(controller.lambda_disp[coupled_index], slot_k)
        g_stored = controller.lambda_disp[coupled_index][slot_k]
        g = isempty(g_stored) ? zeros(length(model.boundary_force)) : g_stored
        # Optional dynamic Aitken relaxation factor. The fixed-point iterate is
        # the impedance RHS stored in lambda_disp; its unrelaxed (theta = 1)
        # candidate is T(g) = boundary_force + impedance coupling term, from which
        # the residual r = T(g) - g is formed. The fixed-theta blend is unchanged.
        if controller.relaxation_method === :aitken_recursive || controller.relaxation_method === :aitken_secant
            candidate = copy(model.boundary_force)
            for comp in 1:3
                Z_W_vdot = Z * (W * dst_velo[comp, :])
                α_W_u = α * (W * dst_disp[comp, :])
                for (i_local, i_global) in enumerate(global_from_local_map)
                    dof_i = 3 * (i_global - 1) + comp
                    candidate[dof_i] += neumann_force[3 * (i_local - 1) + comp] + Z_W_vdot[i_local] + α_W_u[i_local]
                end
            end
            theta = controller.relaxation_method === :aitken_secant ?
                relaxation_aitken_secant_theta!(controller, coupled_index, slot_k, iter, candidate, g) :
                relaxation_aitken_recursive_theta!(controller, coupled_index, slot_k, iter, candidate, g)
        end
        controller.lambda_disp[coupled_index][slot_k] = copy(model.boundary_force)
        for comp in 1:3
            Z_W_vdot = Z * (W * dst_velo[comp, :])
            α_W_u = α * (W * dst_disp[comp, :])
            for (i_local, i_global) in enumerate(global_from_local_map)
                dof_i = 3 * (i_global - 1) + comp
                rhs_i = neumann_force[3 * (i_local - 1) + comp] + Z_W_vdot[i_local] + α_W_u[i_local]
                controller.lambda_disp[coupled_index][slot_k][dof_i] += (1 - theta) * g[dof_i] + theta * rhs_i
                model.boundary_force[dof_i] = controller.lambda_disp[coupled_index][slot_k][dof_i]
            end
        end
    end
end

# Tangent contribution: K_impedance = (Z * c_v + α) * W
# where c_v = γ/(β·Δt) is the Newmark velocity coefficient.
# For quasi-static (no velocity), only the α * W term contributes.
# For explicit integrators, this is not called (no tangent assembly).
function build_impedance_schwarz_stiffness(model::SolidMechanics, integrator::TimeIntegrator)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_is = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition || continue
        Z = bc.impedance
        α = bc.robin_parameter
        W = bc.square_projector
        # Velocity coefficient from Newmark: c_v = γ/(β·Δt)
        c_v = if integrator isa Newmark
            integrator.γ / (integrator.β * integrator.time_step)
        else
            0.0  # quasi-static: no velocity contribution
        end
        coeff = Z * c_v + α
        global_from_local_map = bc.global_from_local_map
        for (i_local, i_global) in enumerate(global_from_local_map)
            for (j_local, j_global) in enumerate(global_from_local_map)
                w_ij = coeff * W[i_local, j_local]
                for comp in 1:3
                    dof_i = 3 * (i_global - 1) + comp
                    dof_j = 3 * (j_global - 1) + comp
                    K_is[dof_i, dof_j] += w_ij
                end
            end
        end
    end
    return K_is
end

has_paired_impedance_bcs(model::SolidMechanics) = any(
    bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition && bc.adjoint_pairing for
    bc in model.boundary_conditions
)

# IMEX treatment of the paired impedance interface for explicit (central
# difference) subdomains. The interface dashpot must act on the END-of-step
# velocity v_{n+1} = ṽ + γΔt·a_{n+1} — evaluating it at the predictor ṽ
# leaves the acceleration loop of the consistent-traction exchange unresolved
# and is violently unstable. Because the dashpot force is linear in the
# unknown acceleration, the implicit treatment closes on the interface rows:
#   (M + γΔt·Z·W)|_Γ a_{n+1} = M·a_expl|_Γ ,
# where a_expl is the plain explicit update (whose right-hand side already
# contains the dashpot at ṽ and the Robin spring at the known u_{n+1}).
# The interior stays matrix-free explicit; only this small SPD system (one
# per interface, reused for the three components) is solved per evaluation,
# which is Liu et al.'s IMEX-Newmark structure. At the Schwarz fixed point
# both sides then satisfy the transmission condition with synchronized
# t_{n+1} velocities, exactly as in the implicit-implicit case.
# Interface DOFs held by Dirichlet BCs keep their prescribed acceleration
# and enter the free rows' right-hand side as data.
function imex_interface_acceleration!(
    a_new::Vector{Float64},
    integrator::CentralDifference,
    model::SolidMechanics,
    solver::Solver,
)
    γΔt = integrator.γ * integrator.time_step
    free = model.free_dofs
    for bc in model.boundary_conditions
        (bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition && bc.adjoint_pairing) || continue
        size(bc.square_projector, 1) > 0 || continue
        Z = bc.impedance
        W = bc.square_projector
        gmap = bc.global_from_local_map
        for comp in 1:3
            dofs = [3 * (g - 1) + comp for g in gmap]
            m = model.lumped_mass[dofs]
            r = m .* a_new[dofs]
            f_mask = free[dofs]
            if all(f_mask)
                a_new[dofs] = (Diagonal(m) + (γΔt * Z) .* W) \ r
            else
                ff = findall(f_mask)
                fx = findall(!, f_mask)
                rf = r[ff] .- (γΔt * Z) .* (W[ff, fx] * a_new[dofs[fx]])
                a_new[dofs[ff]] = (Diagonal(m[ff]) + (γΔt * Z) .* W[ff, ff]) \ rf
            end
        end
    end
    return nothing
end

function pair_bc(bc::SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition, bc_index::Int64)
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = coupled_subsim_of(bc).model
    coupled_bcs = coupled_model.boundary_conditions
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
        end
    end
    return nothing
end

function compute_impedance_schwarz_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition
)
    if dst_bc.adjoint_pairing
        compute_paired_impedance_schwarz_projectors!(dst_model, dst_bc)
        return nothing
    end
    compute_dirichlet_projector(dst_model, dst_bc)
    compute_neumann_projector(dst_model, dst_bc)
    dst_bc.square_projector = get_square_projection_matrix(dst_model, dst_bc)
    return nothing
end

# Adjoint (variationally paired) construction: derive BOTH sides' transfer
# operators from one shared cross-mass matrix B_mn = ∫_Γ φ¹_m φ²_n dS, so
#   Π₁ = W₁⁻¹ B,   Π₂ = W₂⁻¹ Bᵀ,   N₁ = Π₂ᵀ,   N₂ = Π₁ᵀ,
# which is the discrete adjoint condition W₁Π₁ = (W₂Π₂)ᵀ = B of the
# energy-stable DG/mortar couplings. With the shared pair impedance
# Z = 2 Z₁Z₂/(Z₁+Z₂) (harmonic mean; the Riemann-solver weighting, reducing to
# the one-sided value for identical materials) and a shared Robin α, the
# interface power of the dashpot channel telescopes to
# -Z ∫_Γ [[u̇ʰ]]² dS ≤ 0 and the Robin channel becomes a conservative
# interface spring. B is integrated over the facets of the finer side, where
# the piecewise-smooth integrand is resolved best; the coarse-side partition
# of unity and force-conservation certificates are then quadrature-exact up
# to the nonconformity of the interface itself (logged below).
function compute_paired_impedance_schwarz_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition
)
    # The pair is processed once; the partner's pass finds its operators set.
    if size(dst_bc.dirichlet_projector, 1) > 0
        return nothing
    end
    src_sim = coupled_subsim_of(dst_bc)
    src_model = get_fom_model(src_sim)
    src_bc = src_model.boundary_conditions[dst_bc.coupled_bc_index]
    if !(src_bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition) || !src_bc.adjoint_pairing
        norma_abort(
            "`adjoint pairing: true` must be set on BOTH sides of a Schwarz " *
            "impedance nonoverlap pair (missing on the side coupled to '$(dst_bc.name)').",
        )
    end
    W1 = get_square_projection_matrix(dst_model, dst_bc)
    W2 = get_square_projection_matrix(src_model, src_bc)
    n1 = size(W1, 1)
    n2 = size(W2, 1)
    # Integrate B over the finer side's facets (rows = this side, cols = partner).
    B = if n1 >= n2
        get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    else
        Matrix(transpose(get_rectangular_projection_matrix(src_model, src_bc, dst_model, dst_bc)))
    end
    P1 = W1 \ B
    P2 = W2 \ Matrix(transpose(B))
    dst_bc.square_projector = W1
    dst_bc.dirichlet_projector = P1
    dst_bc.neumann_projector = Matrix(transpose(P2))
    src_bc.square_projector = W2
    src_bc.dirichlet_projector = P2
    src_bc.neumann_projector = Matrix(transpose(P1))
    # Shared pair impedance and Robin parameter. A zero pair impedance is the
    # pure-Robin opt-in (`impedance scale: 0`); guard the harmonic mean.
    Z1 = dst_bc.impedance
    Z2 = src_bc.impedance
    Z_pair = Z1 + Z2 > 0.0 ? 2.0 * Z1 * Z2 / (Z1 + Z2) : 0.0
    dst_bc.impedance = src_bc.impedance = Z_pair
    α1 = dst_bc.robin_parameter
    α2 = src_bc.robin_parameter
    if !isapprox(α1, α2; rtol=1.0e-12, atol=0.0)
        norma_abort(
            "Adjoint pairing requires ONE Robin parameter per interface, but " *
            "the sides '$(dst_bc.name)' and '$(src_bc.name)' specify " *
            "$(α1) and $(α2). Set the same `robin parameter` on both sides " *
            "(the Robin spring is conservative only when the two sides' " *
            "interface forces pair through the same coefficient; the converged " *
            "solution does not depend on its value, only the Schwarz " *
            "convergence rate does). To use per-side Robin parameters with the " *
            "legacy per-side transfer instead, set `adjoint pairing: false` " *
            "on both sides of the interface.",
        )
    end
    # The consistent (D'Alembert) traction exchanged under pairing couples the
    # partner's acceleration into this side's force. Implicit (Newmark)
    # subdomains resolve that algebraic loop within the Schwarz iteration;
    # explicit (central difference) subdomains resolve it with the IMEX
    # interface treatment (imex_interface_acceleration!): the interface rows
    # of the acceleration update are solved implicitly so the dashpot acts on
    # the end-of-step velocity, while the interior stays matrix-free explicit.
    for sim_k in (coupled_subsim_of(dst_bc), self_subsim_of(dst_bc))
        if sim_k.integrator isa CentralDifference
            norma_log(
                0,
                :info,
                "Adjoint pairing with an explicit (central difference) subdomain: " *
                "IMEX interface treatment active (implicit interface rows in the " *
                "acceleration update).",
            )
        end
    end
    # Quadrature-quality certificates (exact on conforming interfaces).
    pu_error = max(
        maximum(abs.(P1 * ones(n2) .- 1.0)),
        maximum(abs.(P2 * ones(n1) .- 1.0)),
    )
    norma_logf(
        0,
        :info,
        "Adjoint-paired impedance interface '%s'/'%s': Z = %.4e, α = %.3e, " *
        "partition-of-unity error %.2e.",
        dst_bc.name,
        src_bc.name,
        Z_pair,
        dst_bc.robin_parameter,
        pu_error,
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Impedance overlap Schwarz: absorbing condition on overlap boundaries
# ---------------------------------------------------------------------------
#
# Instead of overwriting DOFs (DBC-DBC), applies an impedance-matching
# (zeroth-order optimized-Schwarz Robin) force built from the partner's
# traction and kinematics, strongly interpolated at this subdomain's overlap
# boundary:
#   boundary_force += W * (t_partner + Z * u̇_partner + α * u_partner)
# with the matching self terms W * (Z * u̇ + α * u) assembled on the
# internal-force side (build_impedance_schwarz_force), so the converged
# transmission condition is
#   t - t_partner + Z * (u̇ - u̇_partner) + α * (u - u_partner) = 0,
# which the monodomain solution satisfies identically. Omitting t_partner
# turns the condition into a relative dashpot t = Z(u̇_partner - u̇) + ...,
# which the monodomain solution fails by exactly t: the converged coupled
# solution then carries a spurious interface velocity jump ≈ t/Z whose power
# drain ~|t|²/Z acts as permanent interface damping.
# DOFs remain free — waves pass through instead of reflecting.

function _get_impedance_scale(bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition)
    scales = bc.impedance_scale
    if length(scales) == 1
        return scales[1]
    end
    parent_sim = bc.parent
    step = max(1, parent_sim.controller.stop)
    return scales[min(step, length(scales))]
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition)
    scale = _get_impedance_scale(bc)
    Z_p = bc.impedance * scale
    Z_s = bc.impedance_shear * scale
    α = bc.robin_parameter
    W = bc.square_projector

    coupled_model_obj = coupled_subsim_of(bc).model
    coupled_solid =
        coupled_model_obj isa SolidMechanics ? coupled_model_obj : coupled_model_obj.fom_model

    # Partner traction t_partner = σ_partner · n at this subdomain's Schwarz
    # boundary. That boundary is an interior surface of the partner mesh, so
    # the partner's assembled internal force cannot supply the traction there
    # (interior rows carry the inertia residual, not σ·n). On node-aligned
    # interfaces the weak traction W·t_p comes from the consistent-traction patch
    # (exact discrete traction; see build_consistent_traction_patch!). Otherwise
    # interpolate the partner's nodal-recovered stress with the same
    # shape-function data used for u̇ and u below; stress recovery is
    # force-enabled for coupled models when this BC type is present (see
    # SolidMechanics() in model.jl).
    use_consistent_traction = bc.traction_patch !== nothing
    if use_consistent_traction && size(bc.traction_patch.transfer, 1) > 0
        # Stage 2c (non-aligned interface): the whole partner right-hand side
        # comes through the single-operator characteristic transfer; the
        # kinematic projector is bypassed for the RHS.
        traction_term, impedance_term, robin_term =
            characteristic_partner_terms(coupled_solid, bc, Z_p, Z_s, α)
        # Representable dashpot: restrict the impedance channel to the
        # transferable subspace. The terms here are weak (W-weighted)
        # vectors; since Π is W-orthogonal (W Π = Πᵀ W), the weak form of
        # the filtered field is the right-multiplication f Π per component
        # (fᵀ ← Πᵀ fᵀ). The matching self term is filtered identically in
        # build_impedance_schwarz_force.
        if size(bc.representable_projector, 1) > 0
            impedance_term = impedance_term * bc.representable_projector
        end
        partner_rhs = traction_term .+ impedance_term .+ robin_term
        for comp in 1:3
            for (i_local, i_global) in enumerate(bc.global_from_local_map)
                dof_i = 3 * (i_global - 1) + comp
                model.boundary_force[dof_i] += partner_rhs[comp, i_local]
            end
        end
        return nothing
    end
    coupled_nodal_stress = Matrix{Float64}(undef, 0, 0)
    weak_traction = Matrix{Float64}(undef, 0, 0)
    if use_consistent_traction
        weak_traction = consistent_partner_traction(coupled_solid, bc)
    else
        recover_stress!(coupled_solid)
        coupled_nodal_stress = if size(coupled_solid.recovered_stress, 1) > 0
            coupled_solid.recovered_stress
        elseif size(coupled_solid.consistent_recovered_stress, 1) > 0
            coupled_solid.consistent_recovered_stress
        else
            norma_abort(
                "Schwarz impedance overlap requires nodal stress recovery on the coupled " *
                "subdomain to evaluate the partner traction. Enable `nodal recovery:` with " *
                "`stress: true` in the coupled subdomain's `model:` block.",
            )
        end
    end
    # Outward unit normals of this subdomain's Schwarz boundary, indexed like
    # global_from_local_map (= unique(side_set_node_indices), same ordering).
    normals = compute_normal(model.mesh, bc.side_set_id, model)

    global_from_local_map = bc.global_from_local_map

    # Partner traction, displacement, and velocity at each overlap boundary
    # node, by variational projection or pointwise interpolation.
    partner_velo, partner_disp, partner_trac =
        impedance_partner_fields(bc, coupled_solid, coupled_nodal_stress, normals)
    num_dst_nodes = size(partner_velo, 2)
    # Representable dashpot: project the partner velocity onto the
    # transferable subspace before the impedance product, matching the
    # filtered self term in build_impedance_schwarz_force, so the dashpot
    # acts on Π(u̇_p − u̇) — the component of the jump both trace spaces can
    # represent, which vanishes at Schwarz convergence. (Under variational
    # transfer u̇_p already lies in the span and only the self side changes.)
    dashpot_velo =
        size(bc.representable_projector, 1) > 0 ? partner_velo * bc.representable_projector' : partner_velo
    partner_Zvdot = zeros(3, num_dst_nodes)
    for i in 1:num_dst_nodes
        # P/S-split tensor impedance: Z_p on the normal velocity component,
        # Z_s on the tangential components (Lysmer-Kuhlemeyer split).
        n = normals[:, i]
        v = dashpot_velo[:, i]
        vn = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
        partner_Zvdot[:, i] = (Z_p * vn) .* n .+ Z_s .* (v .- vn .* n)
    end

    # RHS (partner contribution only): f = W * (t_partner + Z·u̇_partner + α * u_partner)
    # The self terms (W*(Z·u̇_self + α*u_self)) are handled by
    # apply_impedance_bcs_internal_force! called at each Newton iteration.
    # With the consistent-traction patch, W·t_partner is available directly as
    # the weak traction (no separate W application).
    for comp in 1:3
        alpha_u = α * partner_disp[comp, :]
        rhs = W * (partner_Zvdot[comp, :] + alpha_u)
        rhs += use_consistent_traction ? weak_traction[comp, :] : W * partner_trac[comp, :]
        for (i_local, i_global) in enumerate(global_from_local_map)
            dof_i = 3 * (i_global - 1) + comp
            model.boundary_force[dof_i] += rhs[i_local]
        end
    end
    # DOFs remain free — no model.free_dofs modification
end

function build_impedance_overlap_schwarz_stiffness(model::SolidMechanics, integrator::TimeIntegrator)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_io = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition || continue
        scale = _get_impedance_scale(bc)
        Z_p = bc.impedance * scale
        Z_s = bc.impedance_shear * scale
        α = bc.robin_parameter
        W = bc.square_projector
        c_v = if integrator isa Newmark
            integrator.γ / (integrator.β * integrator.time_step)
        else
            0.0
        end
        global_from_local_map = bc.global_from_local_map
        normals = compute_normal(model.mesh, bc.side_set_id, model)
        # Tangent of the self term Σ_j W_ij (T(n_j) u̇_j + α u_j) with the
        # tensor impedance T(n) = Z_p n⊗n + Z_s (I - n⊗n) applied at the
        # source node j: 3x3 block W_ij (c_v T(n_j) + α I). With the
        # representable dashpot the velocity entering the impedance term is
        # Π u̇, so the dashpot part moves into the Π-weighted block below and
        # only the α-spring remains here.
        representable = size(bc.representable_projector, 1) > 0
        for (j_local, j_global) in enumerate(global_from_local_map)
            n = normals[:, j_local]
            T = zeros(3, 3)
            for a in 1:3, b in 1:3
                dash = representable ? 0.0 : c_v * ((Z_p - Z_s) * n[a] * n[b] + (a == b ? Z_s : 0.0))
                T[a, b] = dash + (a == b ? α : 0.0)
            end
            for (i_local, i_global) in enumerate(global_from_local_map)
                w_ij = W[i_local, j_local]
                for a in 1:3, b in 1:3
                    T[a, b] == 0.0 && continue
                    dof_i = 3 * (i_global - 1) + a
                    dof_j = 3 * (j_global - 1) + b
                    K_io[dof_i, dof_j] += w_ij * T[a, b]
                end
            end
        end
        # Filter-weighted dashpot tangents, c_v Σ_l W_il T(n_l) F_lj blocks:
        # F = content_filter for the content-aware absorption term, and
        # F = representable_projector for the filtered self-impedance term.
        for F in (bc.content_filter, bc.representable_projector)
            (size(F, 1) > 0 && c_v != 0.0) || continue
            num_nodes = length(global_from_local_map)
            for (l_local, _) in enumerate(global_from_local_map)
                n = normals[:, l_local]
                Tl = zeros(3, 3)
                for a in 1:3, b in 1:3
                    Tl[a, b] = c_v * ((Z_p - Z_s) * n[a] * n[b] + (a == b ? Z_s : 0.0))
                end
                for i_local in 1:num_nodes, j_local in 1:num_nodes
                    coeff = W[i_local, l_local] * F[l_local, j_local]
                    coeff == 0.0 && continue
                    i_global = global_from_local_map[i_local]
                    j_global = global_from_local_map[j_local]
                    for a in 1:3, b in 1:3
                        Tl[a, b] == 0.0 && continue
                        dof_i = 3 * (i_global - 1) + a
                        dof_j = 3 * (j_global - 1) + b
                        K_io[dof_i, dof_j] += coeff * Tl[a, b]
                    end
                end
            end
        end
    end
    return K_io
end

function pair_bc(bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition, bc_index::Int64)
    # Overlap BCs don't pair — they use strong interpolation from the partner's interior
    return nothing
end

function compute_impedance_overlap_schwarz_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
)
    dst_bc.square_projector = get_square_projection_matrix(dst_model, dst_bc)
    build_consistent_traction_patch!(dst_model, dst_bc)
    if dst_bc.transfer_mode == "variational" || dst_bc.content_absorption || dst_bc.representable_dashpot
        # Variational projector: P = W \\ L is the L2(Γ)-orthogonal projection
        # of the partner fields onto this side's trace space — non-expansive
        # in L2 for any quadrature, and constants transfer exactly (partner
        # partition of unity gives L·1 = W·1). On a node-aligned interface
        # L = W row-permuted and P reduces to the identity selection. Also
        # built (without being used for the transfer) when the content-aware
        # absorption or the representable dashpot needs its column span.
        coupled_solid = get_fom_model(coupled_subsim_of(dst_bc))
        L = get_overlap_rectangular_projection_matrix(
            dst_model, dst_bc, coupled_solid, dst_bc.coupled_block_name, dst_bc.search_tolerance;
            subdivisions=dst_bc.transfer_subdivisions,
        )
        dst_bc.variational_projector = dst_bc.square_projector \ L
    else
        dst_bc.variational_projector = Matrix{Float64}(undef, 0, 0)
    end
    if dst_bc.content_absorption || dst_bc.representable_dashpot
        # W-orthogonal projection Π onto the transferable space at this
        # boundary — the span of the variational projector's columns (one
        # column per partner node with support on the interface). What Π
        # keeps can be represented by the partner and crosses the interface;
        # what I - Π keeps cannot. P·1 = 1 puts constants in the span, and on
        # node-aligned interfaces the span is everything (Π = I).
        # Two independent consumers:
        #   content_filter = I - Π    (content-aware absorption: LK-dissipate
        #                              the non-transferable self content)
        #   representable_projector = Π  (representable dashpot: restrict the
        #                              impedance term to the transferable jump
        #                              so it vanishes at Schwarz convergence)
        P = dst_bc.variational_projector
        W = dst_bc.square_projector
        cols = [j for j in 1:size(P, 2) if maximum(abs.(P[:, j])) > 1.0e-12]
        B = P[:, cols]
        G = B' * W * B
        Π = B * (pinv(G; rtol=1.0e-10) * (B' * W))
        dst_bc.content_filter =
            dst_bc.content_absorption ? Matrix{Float64}(I, size(Π)...) - Π : Matrix{Float64}(undef, 0, 0)
        dst_bc.representable_projector = dst_bc.representable_dashpot ? Π : Matrix{Float64}(undef, 0, 0)
    else
        dst_bc.content_filter = Matrix{Float64}(undef, 0, 0)
        dst_bc.representable_projector = Matrix{Float64}(undef, 0, 0)
    end
    return nothing
end

# Partner kinematic fields (and, when nodal stress is supplied, the partner
# traction) sampled at this BC's boundary nodes: variational L2 projection when
# the projector is present, pointwise interpolation otherwise. Voigt order
# of the nodal stress: xx, yy, zz, yz, xz, xy. Used by apply_bc_detail and
# by the interface-work instrumentation so that the reported powers reflect
# the transfer the BC actually applies.
function impedance_partner_fields(
    bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition,
    coupled_solid::SolidMechanics,
    coupled_nodal_stress::Matrix{Float64},
    normals::Matrix{Float64},
)
    num_dst_nodes = length(bc.global_from_local_map)
    partner_velo = zeros(3, num_dst_nodes)
    partner_disp = zeros(3, num_dst_nodes)
    partner_trac = zeros(3, num_dst_nodes)
    want_traction = size(coupled_nodal_stress, 1) > 0
    traction_from_stress! = (trac, σ, n, i) -> begin
        trac[1, i] = σ[1] * n[1] + σ[6] * n[2] + σ[5] * n[3]
        trac[2, i] = σ[6] * n[1] + σ[2] * n[2] + σ[4] * n[3]
        trac[3, i] = σ[5] * n[1] + σ[4] * n[2] + σ[3] * n[3]
    end
    if bc.transfer_mode == "variational"
        P = bc.variational_projector
        partner_velo .= coupled_solid.velocity * P'
        partner_disp .= coupled_solid.displacement * P'
        if want_traction
            projected_stress = coupled_nodal_stress * P'
            for i in 1:num_dst_nodes
                traction_from_stress!(partner_trac, projected_stress[:, i], normals[:, i], i)
            end
        end
    else
        unique_node_indices = unique(bc.side_set_node_indices)
        for i in eachindex(unique_node_indices)
            coupled_node_indices = bc.coupled_nodes_indices[i]
            N = bc.interpolation_function_values[i]
            for comp in 1:3
                partner_velo[comp, i] = sum(coupled_solid.velocity[comp, coupled_node_indices] .* N)
                partner_disp[comp, i] = sum(coupled_solid.displacement[comp, coupled_node_indices] .* N)
            end
            if want_traction
                σ = coupled_nodal_stress[:, coupled_node_indices] * N
                traction_from_stress!(partner_trac, σ, normals[:, i], i)
            end
        end
    end
    return partner_velo, partner_disp, partner_trac
end

# Variationally consistent partner traction (i.e. consistent nodal reactions;
# the solid-mechanics counterpart of the consistent-boundary-flux method of
# Wheeler and of Hughes et al.): on a node-aligned interface, the weak partner traction
# ∫_Γ t_p·φ_i is exactly the partial assembly of the partner's discrete
# momentum residual (internal force + inertia; body force is not modeled)
# over the partner elements on the EXTERIOR side of Γ — the side this
# subdomain's outward normal points toward, i.e. the elements this subdomain
# is missing. With that definition the monodomain solution satisfies this
# subdomain's coupled equations identically, so the transmission condition
# transmits ALL discrete content, including the mesh-scale part that nodal
# stress recovery misrepresents (measured as O(h) interface dissipation:
# 3.9% → 1.9% per halving of h on the conforming cantilever benchmark).
function build_consistent_traction_patch!(
    model::SolidMechanics, bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
)
    bc.traction_patch = nothing
    bc.partner_traction_mode == "recovered stress" && return nothing
    coupled_subsim = coupled_subsim_of(bc)
    coupled_solid = get_fom_model(coupled_subsim)
    unique_node_indices = unique(bc.side_set_node_indices)
    # Node-aligned check: every boundary node must coincide with a partner
    # node (its interpolation reduces to a single unit weight).
    dst_of_coupled = Dict{Int64,Int64}()
    conforming = true
    for i in eachindex(unique_node_indices)
        N = bc.interpolation_function_values[i]
        a = argmax(N)
        if N[a] < 1.0 - 1.0e-8
            conforming = false
            break
        end
        dst_of_coupled[bc.coupled_nodes_indices[i][a]] = i
    end
    if conforming == false
        # Non-conforming: the consistent traction is available through the
        # offset partner facet surface, which requires the variational
        # transfer machinery (stage 2a). With pointwise transfer, fall back
        # to recovered stress (legacy behavior) unless explicitly requested.
        if bc.transfer_mode == "variational"
            build_offset_traction_patch!(model, bc)
        elseif bc.partner_traction_mode == "consistent traction"
            norma_abort(
                "`partner traction: consistent traction` on a non-aligned interface " *
                "(side set \"$(bc.name)\") requires `transfer: variational`.",
            )
        end
        return nothing
    end
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    element_nodes = Vector{Vector{Int64}}()
    element_block = Vector{Int64}()
    element_index = Vector{Int64}()
    element_rows = Vector{Vector{Tuple{Int64,Int64}}}()
    blocks = Exodus.read_sets(coupled_solid.mesh, Block)
    block_element_type = Vector{ElementType}(undef, length(blocks))
    for (block_index, block) in enumerate(blocks)
        element_type_string = Exodus.read_block_parameters(coupled_solid.mesh, block.id)[1]
        block_element_type[block_index] = element_type_from_string(element_type_string)
        connectivity = get_block_connectivity(coupled_solid.mesh, block.id)
        num_block_elements, num_element_nodes = size(connectivity)
        for e in 1:num_block_elements
            nodes = [connectivity[(e - 1) * num_element_nodes + n] for n in 1:num_element_nodes]
            rows = Vector{Tuple{Int64,Int64}}()
            for (a, g) in enumerate(nodes)
                i = get(dst_of_coupled, g, 0)
                i > 0 && push!(rows, (a, i))
            end
            isempty(rows) && continue
            centroid = vec(sum(coupled_solid.reference[:, nodes]; dims=2)) ./ num_element_nodes
            num_exterior = 0
            num_interior = 0
            for (_, i) in rows
                offset = centroid - model.reference[:, unique_node_indices[i]]
                if dot(offset, normals[:, i]) > 0.0
                    num_exterior += 1
                else
                    num_interior += 1
                end
            end
            if num_exterior > 0 && num_interior > 0
                if bc.partner_traction_mode == "consistent traction"
                    norma_abort(
                        "Ambiguous interface-side classification of a partner element for " *
                        "the consistent traction of side set \"$(bc.name)\".",
                    )
                end
                bc.traction_patch = nothing
                return nothing  # auto: geometry too irregular, fall back
            end
            num_exterior > 0 || continue
            push!(element_nodes, nodes)
            push!(element_block, block_index)
            push!(element_index, e)
            push!(element_rows, rows)
        end
    end
    lumped_mass = coupled_subsim.integrator isa CentralDifference
    bc.traction_patch = ConsistentTractionPatch(
        element_nodes,
        element_block,
        element_index,
        element_rows,
        block_element_type,
        lumped_mass,
        length(unique_node_indices),
        Matrix{Float64}(undef, 0, 0),
        Matrix{Float64}(undef, 0, 0),
        Int64[],
        Matrix{Float64}(undef, 3, 0),
    )
    norma_log(
        0,
        :schwarz,
        "Impedance overlap side set \"$(bc.name)\": node-aligned interface, using " *
        "consistent partner traction ($(length(element_nodes)) partner elements).",
    )
    return nothing
end

# Exodus HEX8 local face connectivity (node orderings; orientation is
# irrelevant here — faces are matched between elements by sorted node sets
# and used only for surface quadrature).
const _hex8_faces = (
    (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (1, 5, 8, 4), (1, 4, 3, 2), (5, 6, 7, 8)
)

# Closest point on a bilinear quadrilateral (columns of X, 3x4) to point p:
# Gauss-Newton on the tangent-orthogonality conditions, clamped to the
# reference square. Returns the reference coordinates and the distance.
function _closest_point_quad4(X::Matrix{Float64}, p::AbstractVector{Float64})
    ξ = zeros(2)
    for _ in 1:20
        N, dN, _ = interpolate(QUAD4, ξ)
        x = X * Vector(N)
        J = X * Matrix(dN)'
        r = J' * (x - p)
        Δ = (J' * J) \ r
        ξ -= Δ
        norm(Δ) < 1.0e-12 && break
    end
    ξ = clamp.(ξ, -1.0, 1.0)
    N, _, _ = interpolate(QUAD4, ξ)
    return ξ, norm(X * Vector(N) - p)
end

# Stage 2a of the variational-transfer design: consistent partner traction
# across a NON-conforming interface. Γ_k cuts through partner elements, so
# the one-sided partial assembly is performed on the offset partner facet
# surface Γ̃ — the boundary between the partner elements exterior to Γ_k
# (centroid on the side this BC's outward normal points toward) and the
# rest, a staircase within one partner element of Γ_k. The weak traction there
# is the exact discrete traction of the partner scheme; it is converted to
# a traction field with W̃⁻¹ and carried to this side's trace space with the
# surface-to-surface variational operator R, precomposed as transfer = R W̃⁻¹.
# The surface offset introduces an O(h_partner) modeling error that acts on
# the exact discrete traction instead of recovered stress.
function build_offset_traction_patch!(
    model::SolidMechanics, bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
)
    coupled_subsim = coupled_subsim_of(bc)
    coupled_solid = get_fom_model(coupled_subsim)
    unique_node_indices = unique(bc.side_set_node_indices)
    num_dst_nodes = length(unique_node_indices)
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    own_points = model.reference[:, unique_node_indices]

    # Classify every partner element by the side of Γ_k its centroid lies on,
    # using the outward normal at the nearest boundary node of this side.
    struct_record = Vector{Tuple{Int64,Int64,Vector{Int64},Bool}}()  # (block, elem, nodes, exterior)
    blocks = Exodus.read_sets(coupled_solid.mesh, Block)
    block_element_type = Vector{ElementType}(undef, length(blocks))
    for (block_index, block) in enumerate(blocks)
        element_type_string = Exodus.read_block_parameters(coupled_solid.mesh, block.id)[1]
        block_element_type[block_index] = element_type_from_string(element_type_string)
        connectivity = get_block_connectivity(coupled_solid.mesh, block.id)
        num_block_elements, num_element_nodes = size(connectivity)
        if num_element_nodes != 8
            if bc.partner_traction_mode == "consistent traction"
                norma_abort(
                    "Consistent traction across a non-aligned interface currently requires " *
                    "hexahedral partner elements (side set \"$(bc.name)\").",
                )
            end
            return nothing  # auto: fall back to recovered stress
        end
        for e in 1:num_block_elements
            nodes = [connectivity[(e - 1) * num_element_nodes + n] for n in 1:num_element_nodes]
            centroid = vec(sum(coupled_solid.reference[:, nodes]; dims=2)) ./ num_element_nodes
            best_i = 1
            best_d2 = Inf
            for i in 1:num_dst_nodes
                d2 = sum((centroid .- own_points[:, i]) .^ 2)
                if d2 < best_d2
                    best_d2 = d2
                    best_i = i
                end
            end
            exterior = dot(centroid - own_points[:, best_i], normals[:, best_i]) > 0.0
            push!(struct_record, (block_index, e, nodes, exterior))
        end
    end

    # Γ̃ faces: partner element faces shared by an exterior and an interior
    # element, found by matching sorted node quadruples.
    face_owner = Dict{NTuple{4,Int64},Vector{Tuple{Bool,NTuple{4,Int64}}}}()
    for (_, _, nodes, exterior) in struct_record
        for f in _hex8_faces
            fn = (nodes[f[1]], nodes[f[2]], nodes[f[3]], nodes[f[4]])
            key = NTuple{4,Int64}(sort(collect(fn)))
            push!(get!(face_owner, key, Vector{Tuple{Bool,NTuple{4,Int64}}}()), (exterior, fn))
        end
    end
    tilde_faces = Vector{NTuple{4,Int64}}()
    for (_, owners) in face_owner
        if length(owners) == 2 && owners[1][1] != owners[2][1]
            push!(tilde_faces, owners[1][2])
        end
    end
    if isempty(tilde_faces)
        if bc.partner_traction_mode == "consistent traction"
            norma_abort(
                "Could not construct the offset partner facet surface for side set " *
                "\"$(bc.name)\" (no exterior/interior partner element boundary found).",
            )
        end
        return nothing
    end
    tilde_index = Dict{Int64,Int64}()
    for fn in tilde_faces, g in fn
        haskey(tilde_index, g) || (tilde_index[g] = length(tilde_index) + 1)
    end
    num_tilde = length(tilde_index)

    # Surface mass matrix W̃ on Γ̃ (2x2 Gauss per quadrilateral face).
    N_q, dN_q, w_q, _ = isoparametric(QUAD4, 4)
    W̃ = zeros(num_tilde, num_tilde)
    for fn in tilde_faces
        face_nodes = collect(fn)
        Xf = coupled_solid.reference[:, face_nodes]
        rows = [tilde_index[g] for g in face_nodes]
        for point in 1:4
            Np = Vector(N_q[:, point])
            dNp = Matrix(dN_q[:, :, point])
            dXdξ = dNp * Xf'
            j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
            W̃[rows, rows] += Np * Np' * j * w_q[point]
        end
    end

    # Surface-to-surface variational operator R = ∫_Γ N_k Ñᵀ dΓ, assembled
    # with this BC's (optionally subdivided) facet quadrature; at each
    # quadrature point the closest Γ̃ face provides the source values.
    local_from_global_map = bc.local_from_global_map
    R = zeros(num_dst_nodes, num_tilde)
    coords = model.reference
    side_set_node_index = 1
    m = bc.transfer_subdivisions
    g = 1.0 / sqrt(3.0)
    for num_side in bc.num_nodes_sides
        side_nodes = bc.side_set_node_indices[side_set_node_index:(side_set_node_index + num_side - 1)]
        side_set_node_index += num_side
        num_side == 4 || norma_abort(
            "Consistent traction across a non-aligned interface requires quadrilateral " *
            "boundary facets (side set \"$(bc.name)\").",
        )
        Xs = coords[:, side_nodes]
        dst_rows = [local_from_global_map[g_] for g_ in side_nodes]
        for i in 1:m, j_ in 1:m, (η₁, η₂) in ((-g, -g), (g, -g), (g, g), (-g, g))
            ξ_sub = [-1.0 + (2 * i - 1) / m + η₁ / m, -1.0 + (2 * j_ - 1) / m + η₂ / m]
            Np, dNp, _ = interpolate(QUAD4, ξ_sub)
            Np = Vector(Np)
            dNp = Matrix(dNp)
            dXdξ = dNp * Xs'
            jac = norm(cross(dXdξ[1, :], dXdξ[2, :]))
            x_qp = Xs * Np
            best = (Inf, zeros(2), tilde_faces[1])
            for fn in tilde_faces
                Xf = coupled_solid.reference[:, collect(fn)]
                ξ_f, dist = _closest_point_quad4(Xf, x_qp)
                dist < best[1] && (best = (dist, ξ_f, fn))
            end
            Ñ, _, _ = interpolate(QUAD4, best[2])
            src_rows = [tilde_index[g_] for g_ in best[3]]
            R[dst_rows, src_rows] += Np * Vector(Ñ)' * jac / m^2
        end
    end
    transfer = Matrix((W̃ \ R')')

    # One-sided patch: exterior partner elements touching any Γ̃ node.
    element_nodes = Vector{Vector{Int64}}()
    element_block = Vector{Int64}()
    element_index = Vector{Int64}()
    element_rows = Vector{Vector{Tuple{Int64,Int64}}}()
    for (block_index, e, nodes, exterior) in struct_record
        exterior || continue
        rows = Vector{Tuple{Int64,Int64}}()
        for (a, g_) in enumerate(nodes)
            idx = get(tilde_index, g_, 0)
            idx > 0 && push!(rows, (a, idx))
        end
        isempty(rows) && continue
        push!(element_nodes, nodes)
        push!(element_block, block_index)
        push!(element_index, e)
        push!(element_rows, rows)
    end
    lumped_mass = coupled_subsim.integrator isa CentralDifference
    tilde_nodes = Vector{Int64}(undef, num_tilde)
    for (g_, idx) in tilde_index
        tilde_nodes[idx] = g_
    end
    tilde_normals = zeros(3, num_tilde)
    for idx in 1:num_tilde
        p = coupled_solid.reference[:, tilde_nodes[idx]]
        best_i = 1
        best_d2 = Inf
        for i in 1:num_dst_nodes
            d2 = sum((p .- own_points[:, i]) .^ 2)
            if d2 < best_d2
                best_d2 = d2
                best_i = i
            end
        end
        tilde_normals[:, idx] = normals[:, best_i]
    end
    bc.traction_patch = ConsistentTractionPatch(
        element_nodes, element_block, element_index, element_rows, block_element_type,
        lumped_mass, num_tilde, transfer, W̃, tilde_nodes, tilde_normals,
    )
    norma_log(
        0,
        :schwarz,
        "Impedance overlap side set \"$(bc.name)\": non-aligned interface, using " *
        "consistent partner traction via the offset facet surface " *
        "($(length(element_nodes)) patch elements, $(length(tilde_faces)) offset faces).",
    )
    return nothing
end

# Weak partner traction W·t_p (3 × num boundary nodes) by partial assembly
# over the consistent-traction patch: -(internal force + inertia) of the exterior-side
# partner elements, restricted to the interface rows. The element kernels
# mirror evaluate() in model.jl (total Lagrangian, first Piola against
# reference gradients); the inertia uses the partner's mass discretization
# (consistent for implicit, lumped for explicit).
function consistent_partner_traction(
    coupled_solid::SolidMechanics, bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
)
    patch = bc.traction_patch
    weak_traction = zeros(3, patch.num_targets)
    for element in eachindex(patch.element_nodes)
        nodes = patch.element_nodes[element]
        block_index = patch.element_block[element]
        block_element_index = patch.element_index[element]
        material = coupled_solid.materials[block_index]
        density = material.ρ
        element_type = patch.block_element_type[block_index]
        num_points = coupled_solid.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        num_element_nodes = length(nodes)
        X = coupled_solid.reference[:, nodes]
        x = X + coupled_solid.displacement[:, nodes]
        acceleration = coupled_solid.acceleration[:, nodes]
        λ = zeros(3, num_element_nodes)
        reduced_mass = patch.lumped_mass ? nothing : zeros(num_element_nodes, num_element_nodes)
        lumped = patch.lumped_mass ? zeros(num_element_nodes) : nothing
        for point in 1:num_points
            Np = N[:, point]
            dNdξ = dN[:, :, point]
            dXdξ = SMatrix{3,3,Float64,9}(dNdξ * X')
            dNdX = dXdξ \ dNdξ
            F = SMatrix{3,3,Float64,9}(x * dNdX')
            local P
            if material isa Elastic
                _, P, _ = constitutive(material, F; need_tangent=false)
            else
                state = coupled_solid.state_old[block_index][block_element_index][point]
                _, P, _, _ = constitutive(material, F, state; need_tangent=false)
            end
            dvol = det(dXdξ) * ip_weights[point]
            λ .+= (P * dNdX) .* dvol
            if patch.lumped_mass
                lumped .+= (density * dvol * sum(Np)) .* Np
            else
                reduced_mass .+= (density * dvol) .* (Np * Np')
            end
        end
        λ .+= patch.lumped_mass ? acceleration .* lumped' : acceleration * reduced_mass
        for (a, i) in patch.element_rows[element]
            weak_traction[1, i] -= λ[1, a]
            weak_traction[2, i] -= λ[2, a]
            weak_traction[3, i] -= λ[3, a]
        end
    end
    # Non-conforming (offset-surface) patch: carry the weak traction from the
    # offset surface Γ̃ to this side's trace space.
    if size(patch.transfer, 1) > 0
        weak_traction = weak_traction * patch.transfer'
    end
    return weak_traction
end

# Stage 2c: the partner's complete Robin datum, assembled on the offset
# surface Γ̃ from the partner's OWN nodal trace values (no interpolation)
# and carried over by the single operator T = R W̃⁻¹. Returned as the three
# channels (traction, impedance, Robin), each weak on this side's trace
# space, so callers can sum them for the right-hand side and the
# instrumentation can report them separately. Transferring the combined
# datum through one operator preserves the characteristic cancellation
# between traction and velocity at mesh scale, which separate transfers
# through different operators destroy.
function characteristic_partner_terms(
    coupled_solid::SolidMechanics,
    bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition,
    Z_p::Float64,
    Z_s::Float64,
    α::Float64,
)
    patch = bc.traction_patch
    traction_term = consistent_partner_traction(coupled_solid, bc)
    num_tilde = patch.num_targets
    zv = zeros(3, num_tilde)
    au = zeros(3, num_tilde)
    for m in 1:num_tilde
        n = patch.tilde_normals[:, m]
        g_ = patch.tilde_nodes[m]
        v = coupled_solid.velocity[:, g_]
        u = coupled_solid.displacement[:, g_]
        vn = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
        zv[:, m] = (Z_p * vn) .* n .+ Z_s .* (v .- vn .* n)
        au[:, m] = α .* u
    end
    impedance_term = (zv * patch.tilde_mass) * patch.transfer'
    robin_term = (au * patch.tilde_mass) * patch.transfer'
    return traction_term, impedance_term, robin_term
end

# Compute the self-impedance internal force: W*(Z*u̇_self + α*u_self)
# Returns a separate force vector rather than modifying model.internal_force,
# because get_dst_force reads model.internal_force for traction transfer and
# must see only the elastic internal force.
function build_impedance_schwarz_force(model::SolidMechanics)
    num_dofs = 3 * size(model.reference, 2)
    f = zeros(num_dofs)
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
            # P/S-split tensor impedance, matching the partner terms in
            # apply_bc_detail (same scale, same tensor, same W weighting).
            scale = _get_impedance_scale(bc)
            Z_p = bc.impedance * scale
            Z_s = bc.impedance_shear * scale
            α = bc.robin_parameter
            W = bc.square_projector
            global_from_local_map = bc.global_from_local_map
            num_nodes = length(global_from_local_map)
            normals = compute_normal(model.mesh, bc.side_set_id, model)
            velo = model.velocity[:, global_from_local_map]
            # Representable dashpot: the self velocity entering the impedance
            # term is projected onto the transferable subspace, matching the
            # filtered partner term in apply_bc_detail so the dashpot acts on
            # Π(u̇_p − u̇). The α·u spring term stays unfiltered.
            dashpot_velo =
                size(bc.representable_projector, 1) > 0 ? velo * bc.representable_projector' : velo
            # Content-aware absorption: LK-dissipate the non-transferable
            # component (I - Π) u̇ of this side's boundary velocity.
            filtered = size(bc.content_filter, 1) > 0 ? velo * bc.content_filter' : zeros(3, 0)
            Zv = zeros(3, num_nodes)
            for (i_local, i_global) in enumerate(global_from_local_map)
                n = normals[:, i_local]
                v = dashpot_velo[:, i_local]
                u = model.displacement[:, i_global]
                vn = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
                Zv[:, i_local] = (Z_p * vn) .* n .+ Z_s .* (v .- vn .* n) .+ α .* u
                if size(filtered, 2) > 0
                    vf = filtered[:, i_local]
                    vfn = n[1] * vf[1] + n[2] * vf[2] + n[3] * vf[3]
                    Zv[:, i_local] .+= (Z_p * vfn) .* n .+ Z_s .* (vf .- vfn .* n)
                end
            end
            for comp in 1:3
                f_imp = W * Zv[comp, :]
                for (i_local, i_global) in enumerate(global_from_local_map)
                    dof_i = 3 * (i_global - 1) + comp
                    f[dof_i] += f_imp[i_local]
                end
            end
        elseif bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition
            Z = bc.impedance
            α = bc.robin_parameter
            W = bc.square_projector
            global_from_local_map = bc.global_from_local_map
            num_nodes = length(global_from_local_map)
            self_velo = zeros(num_nodes)
            self_disp = zeros(num_nodes)
            for comp in 1:3
                for (i_local, i_global) in enumerate(global_from_local_map)
                    self_velo[i_local] = model.velocity[comp, i_global]
                    self_disp[i_local] = model.displacement[comp, i_global]
                end
                f_imp = W * (Z * self_velo + α * self_disp)
                for (i_local, i_global) in enumerate(global_from_local_map)
                    dof_i = 3 * (i_global - 1) + comp
                    f[dof_i] += f_imp[i_local]
                end
            end
        end
    end
    return f
end

# --- Interface work instrumentation (diagnosis; enabled via env var) --------
#
# When NORMA_IMPEDANCE_WORK_CSV is set, a row is appended per impedance-overlap
# BC after each converged Schwarz stop, decomposing the instantaneous interface
# power into partner-traction, dashpot, and Robin parts:
#
#   P = Σ_i u̇_i · [ W (t_p + T(n)(u̇_p − u̇) + α (u_p − u)) ]_i
#
# The exact (monodomain) solution satisfies u̇_p = u̇ and u_p = u on the
# overlap boundary, so the dashpot and Robin powers vanish identically for it:
# their measured work is purely spurious interface dissipation. The RMS
# velocity jump |u̇_p − u̇| and traction mismatch |t_p − t_own| (both from
# recovered stress) are recorded to attribute the residual: recovery error
# shows up as traction mismatch, Schwarz iteration lag and time-level error as
# velocity jump that shrinks (or does not) with the Schwarz tolerance.

function report_impedance_interface_work(sim)
    csv_path = get(ENV, "NORMA_IMPEDANCE_WORK_CSV", "")
    isempty(csv_path) && return nothing
    if isfile(csv_path) == false
        open(csv_path, "w") do io
            println(
                io,
                "time,subdomain,side_set,P_total,P_traction,P_dashpot,P_robin," *
                "velocity_jump_rms,traction_mismatch_rms,P_absorb",
            )
        end
    end
    t = sim.controller.time
    for (subsim_index, subsim) in enumerate(sim.subsims)
        model = subsim.model
        model isa SolidMechanics || continue
        for bc in model.boundary_conditions
            bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition || continue
            scale = _get_impedance_scale(bc)
            Z_p = bc.impedance * scale
            Z_s = bc.impedance_shear * scale
            α = bc.robin_parameter
            W = bc.square_projector
            coupled_model_obj = coupled_subsim_of(bc).model
            coupled_solid =
                coupled_model_obj isa SolidMechanics ? coupled_model_obj : coupled_model_obj.fom_model
            recover_stress!(coupled_solid)
            recover_stress!(model)
            nodal_stress(m) =
                size(m.recovered_stress, 1) > 0 ? m.recovered_stress : m.consistent_recovered_stress
            coupled_nodal_stress = nodal_stress(coupled_solid)
            own_nodal_stress = nodal_stress(model)
            normals = compute_normal(model.mesh, bc.side_set_id, model)
            unique_node_indices = unique(bc.side_set_node_indices)
            num_nodes = length(unique_node_indices)
            partner_velo, partner_disp, trac_p =
                impedance_partner_fields(bc, coupled_solid, coupled_nodal_stress, normals)
            # Mirror the representable-dashpot filter applied by the BC: the
            # dashpot channel uses the Π-projected velocities. The reported
            # velocity-jump RMS stays RAW so it keeps measuring the full
            # (including unrepresentable) interface jump.
            representable = size(bc.representable_projector, 1) > 0
            self_velo = model.velocity[:, unique_node_indices]
            dash_partner = representable ? partner_velo * bc.representable_projector' : partner_velo
            dash_self = representable ? self_velo * bc.representable_projector' : self_velo
            dash = zeros(3, num_nodes)
            robin = zeros(3, num_nodes)
            velo = zeros(3, num_nodes)
            vjump2 = 0.0
            tmis2 = 0.0
            for (i, node_index) in enumerate(unique_node_indices)
                n = normals[:, i]
                v = model.velocity[:, node_index]
                u = model.displacement[:, node_index]
                σ_o = own_nodal_stress[:, node_index]
                t_o = [
                    σ_o[1] * n[1] + σ_o[6] * n[2] + σ_o[5] * n[3],
                    σ_o[6] * n[1] + σ_o[2] * n[2] + σ_o[4] * n[3],
                    σ_o[5] * n[1] + σ_o[4] * n[2] + σ_o[3] * n[3],
                ]
                dv = dash_partner[:, i] .- dash_self[:, i]
                dvn = n[1] * dv[1] + n[2] * dv[2] + n[3] * dv[3]
                dash[:, i] = (Z_p * dvn) .* n .+ Z_s .* (dv .- dvn .* n)
                robin[:, i] = α .* (partner_disp[:, i] .- u)
                velo[:, i] = v
                dv_raw = partner_velo[:, i] .- v
                vjump2 += dot(dv_raw, dv_raw)
                tmis2 += dot(trac_p[:, i] - t_o, trac_p[:, i] - t_o)
            end
            P_trac = 0.0
            P_dash = 0.0
            P_robin = 0.0
            # Channel powers from the terms the BC actually applies. With the
            # stage-2c characteristic transfer, each channel of the partner
            # datum passes through the same operator; the self parts of the
            # dashpot and Robin channels are subtracted at this side's nodes.
            # The recovered-stress traction mismatch below stays informational
            # in all modes.
            patch = bc.traction_patch
            if patch !== nothing && size(patch.transfer, 1) > 0
                traction_term, impedance_term, robin_term =
                    characteristic_partner_terms(coupled_solid, bc, Z_p, Z_s, α)
                if representable
                    impedance_term = impedance_term * bc.representable_projector
                end
                self_zv = zeros(3, num_nodes)
                self_au = zeros(3, num_nodes)
                for (i, node_index) in enumerate(unique_node_indices)
                    n = normals[:, i]
                    v = dash_self[:, i]
                    vn = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
                    self_zv[:, i] = (Z_p * vn) .* n .+ Z_s .* (v .- vn .* n)
                    self_au[:, i] = α .* model.displacement[:, node_index]
                end
                for comp in 1:3
                    P_trac += dot(velo[comp, :], traction_term[comp, :])
                    P_dash += dot(velo[comp, :], impedance_term[comp, :] - W * self_zv[comp, :])
                    P_robin += dot(velo[comp, :], robin_term[comp, :] - W * self_au[comp, :])
                end
            else
                if patch !== nothing
                    weak_traction = consistent_partner_traction(coupled_solid, bc)
                    for comp in 1:3
                        P_trac += dot(velo[comp, :], weak_traction[comp, :])
                    end
                else
                    for comp in 1:3
                        P_trac += dot(velo[comp, :], W * trac_p[comp, :])
                    end
                end
                for comp in 1:3
                    P_dash += dot(velo[comp, :], W * dash[comp, :])
                    P_robin += dot(velo[comp, :], W * robin[comp, :])
                end
            end
            # Content-aware absorption drain (self term; ≤ 0 by construction).
            P_absorb = 0.0
            if size(bc.content_filter, 1) > 0
                filtered = velo * bc.content_filter'
                zvf = zeros(3, num_nodes)
                for i in 1:num_nodes
                    n = normals[:, i]
                    vf = filtered[:, i]
                    vfn = n[1] * vf[1] + n[2] * vf[2] + n[3] * vf[3]
                    zvf[:, i] = (Z_p * vfn) .* n .+ Z_s .* (vf .- vfn .* n)
                end
                for comp in 1:3
                    P_absorb -= dot(velo[comp, :], W * zvf[comp, :])
                end
            end
            P_total = P_trac + P_dash + P_robin + P_absorb
            vjump_rms = sqrt(vjump2 / num_nodes)
            tmis_rms = sqrt(tmis2 / num_nodes)
            open(csv_path, "a") do io
                println(
                    io,
                    "$t,$subsim_index,$(bc.side_set_id),$P_total,$P_trac,$P_dash,$P_robin," *
                    "$vjump_rms,$tmis_rms,$P_absorb",
                )
            end
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------

function find_point_in_mesh(point::Vector{Float64}, model::SolidMechanics, block_id::Int, tol::Float64)
    mesh = model.mesh
    element_type_string = Exodus.read_block_parameters(mesh, Int32(block_id))[1]
    element_type = element_type_from_string(element_type_string)
    element_block_connectivity = get_block_connectivity(mesh, block_id)
    num_block_elements, num_element_nodes = size(element_block_connectivity)
    node_indices = Vector{Int64}()
    found = false
    ξ = zeros(length(point))
    for block_element_index in 1:num_block_elements
        connectivity_indices =
            ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
        node_indices = element_block_connectivity[connectivity_indices]
        element_ref_pos = model.reference[:, node_indices]
        ξ, found = is_inside(element_type, element_ref_pos, point, tol)
        if found == true
            break
        end
    end
    return node_indices, ξ, found
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

        model.displacement[:, node_index] = (coupled_reference[:, coupled_node_indices] + coupled_displacement[:, coupled_node_indices]) * N - model.reference[:, node_index]
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

function coupling_weak_dbc(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    nodal_curr, _, nodal_velo, nodal_acce = get_dst_curr_disp_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    for (i_local, i_global) in enumerate(global_from_local_map)
        @inbounds model.displacement[:, i_global] = nodal_curr[:, i_local] - model.reference[:, i_global]
        @inbounds model.velocity[:, i_global] = nodal_velo[:, i_local]
        @inbounds model.acceleration[:, i_global] = nodal_acce[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        model.free_dofs[global_range] .= false
    end
end

function coupling_weak_nbc(model::SolidMechanics, bc::SolidMechanicsNonOverlapSchwarzBoundaryCondition)
    nodal_force = get_dst_force(bc)
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

function set_internal_force!(model::SolidMechanics, force)
    return model.internal_force = force
end

# Expand interface-sized field (3 × n_interface_nodes) to full-DOF vector (num_dofs)
function _expand_to_full_dofs(field_iface::Matrix{Float64}, global_from_local_map, num_dofs::Int)
    full = zeros(num_dofs)
    num_nodes = num_dofs ÷ 3
    full_3xN = reshape(full, (3, num_nodes))
    for (i_local, i_global) in enumerate(global_from_local_map)
        @inbounds full_3xN[:, i_global] = field_iface[:, i_local]
    end
    return reshape(full_3xN, num_dofs)
end

function apply_bc(model::Model, bc::SolidMechanicsSchwarzBoundaryCondition)
    parent_sim = bc.parent
    controller = parent_sim.controller

    # Skip application if contact is inactive
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact
        return nothing
    end

    coupled_subsim     = coupled_subsim_of(bc)
    coupled_integrator = coupled_subsim.integrator
    coupled_model      = coupled_subsim.model

    # Save current state (copy data, not reference — aliased integrators share memory with model)
    # For coupled RomModel, these are reduced states
    saved_disp = copy(coupled_integrator.displacement)
    saved_velo = copy(coupled_integrator.velocity)
    saved_acce = copy(coupled_integrator.acceleration)

    # Even for RomModel, this is full-dimensional
    saved_∂Ω_f = get_internal_force(coupled_model)

    # Fetch interpolation inputs
    time = model.time
    coupled_index = bc.coupled_handle.id
    coupled_num_dofs = length(coupled_model.free_dofs)

    time_hist = controller.time_hist[coupled_index]
    disp_hist = controller.disp_hist[coupled_index]
    velo_hist = controller.velo_hist[coupled_index]
    acce_hist = controller.acce_hist[coupled_index]
    ∂Ω_f_hist = controller.∂Ω_f_hist[coupled_index]

    # Interpolate or use fallback
    use_predictor = bc isa SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition &&
    controller.use_interface_predictor &&
    controller.iteration_number == 0 &&
    !isempty(controller.predictor_disp[coupled_index])
    if !isempty(time_hist)
        # Piecewise-linear interpolation of the partner trajectory, in both
        # directions of a subcycled pair (Gravouil-Combescure interpolated
        # continuity). This is the measured stable optimum of the temporal
        # transfer operators tried for the subcycled adjoint pairing on the
        # cantilever benchmark (doc/notes/schwarz-interface-energy):
        # replacements for the coarser side's exchange that are exactly
        # work-conjugate to the finer side's interpolation -- the trapezoidal
        # window average, the de-lagged recursion F = 2 avg - F_prev, and the
        # least-squares endpoint fit (all three Schwarz-compatible forms of
        # the Prakash-Hjelmstad linear-in-time interface traction) -- either
        # add half-window-lag dissipation (average) or diverge (recursion:
        # gain two on the partner-state channel doubles the fixed-point
        # iteration's contraction factor; endpoint fit: zero-lag filtering is
        # extrapolatory and pumps energy over long horizons). Exact
        # interface-work telescoping under subcycling requires solving the
        # interface problem directly, as Prakash-Hjelmstad do, not a
        # partitioned exchange.
        interp_disp = interpolate(time_hist, disp_hist, time)
        interp_velo = interpolate(time_hist, velo_hist, time)
        interp_acce = interpolate(time_hist, acce_hist, time)
        interp_∂Ω_f = interpolate(time_hist, ∂Ω_f_hist, time)
    elseif use_predictor
        interp_disp = controller.predictor_disp[coupled_index]
        interp_velo = controller.predictor_velo[coupled_index]
        interp_acce = controller.predictor_acce[coupled_index]
        interp_∂Ω_f = !isempty(controller.predictor_∂Ω_f[coupled_index]) ?
            controller.predictor_∂Ω_f[coupled_index] : controller.stop_∂Ω_f[coupled_index]
    elseif isempty(time_hist) && !isempty(controller.stop_disp[coupled_index])
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


    # Assign interpolated force
    set_internal_force!(coupled_model, interp_∂Ω_f)

    # Apply relaxed update if needed
    if is_swappable_dn_schwarz(bc)
        iter = controller.iteration_number
        # Per-substep-time relaxation slot (see relaxation_slot!): the previous
        # iterate must come from the previous Schwarz sweep at the SAME time. A
        # fresh slot (every slot on iteration 0) has no previous iterate and
        # falls back to the interpolated partner state, as before.
        slot_k = relaxation_slot!(controller, coupled_index, time)
        ensure_slot!(controller.lambda_disp[coupled_index], slot_k)
        ensure_slot!(controller.lambda_velo[coupled_index], slot_k)
        ensure_slot!(controller.lambda_acce[coupled_index], slot_k)
        λ_u_stored = controller.lambda_disp[coupled_index][slot_k]
        λ_u_prev = isempty(λ_u_stored) ? interp_disp : λ_u_stored
        λ_v_stored = controller.lambda_velo[coupled_index][slot_k]
        λ_v_prev = isempty(λ_v_stored) ? interp_velo : λ_v_stored
        λ_a_stored = controller.lambda_acce[coupled_index][slot_k]
        λ_a_prev = isempty(λ_a_stored) ? interp_acce : λ_a_stored

        if controller.relaxation_method === :aitken_secant
            # Relax the single d-form interface unknown (displacement); recover
            # velocity and acceleration consistently from it (see functions above).
            θ = relaxation_aitken_secant_theta!(controller, coupled_index, slot_k, iter, interp_disp, λ_u_prev)
            controller.lambda_disp[coupled_index][slot_k] = θ * interp_disp + (1 - θ) * λ_u_prev
            recover_interface_kinematics!(controller, coupled_index, slot_k, coupled_integrator, interp_velo, interp_acce)
        else
            θ = relaxation_aitken_recursive_theta!(controller, coupled_index, slot_k, iter, interp_disp, λ_u_prev)

            controller.lambda_disp[coupled_index][slot_k] = θ * interp_disp + (1 - θ) * λ_u_prev
            controller.lambda_velo[coupled_index][slot_k] = θ * interp_velo + (1 - θ) * λ_v_prev
            controller.lambda_acce[coupled_index][slot_k] = θ * interp_acce + (1 - θ) * λ_a_prev
        end

        coupled_integrator.displacement .= controller.lambda_disp[coupled_index][slot_k]
        coupled_integrator.velocity .= controller.lambda_velo[coupled_index][slot_k]
        coupled_integrator.acceleration .= controller.lambda_acce[coupled_index][slot_k]
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

function contact_weak_dbc(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    nodal_curr, _, nodal_velo, nodal_acce = get_dst_curr_disp_velo_acce(bc)
    global_from_local_map = bc.global_from_local_map
    normals = compute_normal(model.mesh, bc.side_set_id, model)
    for (i_local, i_global) in enumerate(global_from_local_map)
        normal = normals[:, i_local]
        global_range = (3 * (i_global - 1) + 1):(3 * i_global)
        if bc.friction_type == 0
            @inbounds model.displacement[:, i_global] = transfer_normal_component(
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
            @inbounds model.acceleration[:, i_global] = nodal_velo[:, i_local]
            model.free_dofs[global_range] .= false
        else
            norma_abort("Unknown or not implemented friction type.")
        end
        # Update the rotation matrix
        axis = SVector{3,Float64}(-normalize(normal))
        bc.rotation_matrix = compute_rotation_matrix(axis)
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

function compute_neumann_projector(dst_model::Model, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
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
    L = get_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
    src_lumped = H * ones(size(H, 2))
    conservation_error = maximum(abs.(vec(sum(L; dims=1)) - src_lumped)) / maximum(abs.(src_lumped))
    if conservation_error > 1.0e-10
        L = get_src_integrated_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
    end
    dst_bc.neumann_projector = L * (H \ I)
    return nothing
end

function compute_dirichlet_projector(dst_model::Model, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
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
    L = get_src_integrated_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
    P = (W \ I) * L
    pu_error = maximum(abs.(P * ones(size(P, 2)) .- 1.0))
    if pu_error > 1.0e-10
        L = get_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
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
    src_global_force = get_internal_force(src_model)
    # Adjoint-paired impedance interfaces exchange the dynamically consistent
    # (D'Alembert) reaction M·a + f_int - f_body instead of the static
    # internal force alone. In dynamics the two sides' static reactions are
    # NOT equal and opposite — they differ by the interface nodes' inertia,
    # each side carrying its own tributary mass — so transferring -f_int
    # makes the coupled fixed point inconsistent with the monodomain
    # discretization by O(m_Γ·a), which the cantilever benchmark measures as
    # a steady interface energy drain (-9.6%/ms on CONFORMING meshes). With
    # the inertia included, the two sides' interface rows sum to the
    # monodomain row exactly, so the monodomain trajectory is the fixed
    # point and the transmission channels vanish at Schwarz convergence.
    # (The overlap variant achieves the same through its consistent-traction
    # element patch; here the source's own assembled operators suffice.)
    if dst_bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition && dst_bc.adjoint_pairing
        src_fom = get_fom_model(src_sim)
        a = vec(src_fom.acceleration)
        inertial_force = if size(src_fom.mass, 1) == length(a)
            src_fom.mass * a
        elseif length(src_fom.lumped_mass) == length(a)
            src_fom.lumped_mass .* a
        else
            # Pre-initialization exchange: no mass assembled yet, and the
            # acceleration is still zero — the correction vanishes anyway.
            zeros(length(a))
        end
        src_global_force = src_global_force + inertial_force
        if length(src_fom.body_force) == length(src_global_force)
            src_global_force = src_global_force - src_fom.body_force
        end
    end
    src_force = -extract_local_vector(src_bc, src_global_force, 3)
    neumann_projector = dst_bc.neumann_projector
    num_dst_nodes = size(neumann_projector, 1)
    dst_force = zeros(3 * num_dst_nodes)
    dst_force[1:3:end] = neumann_projector * src_force[1:3:end]
    dst_force[2:3:end] = neumann_projector * src_force[2:3:end]
    dst_force[3:3:end] = neumann_projector * src_force[3:3:end]
    return dst_force
end

function get_dst_curr_disp_velo_acce(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = coupled_subsim_of(dst_bc)
    src_model = src_sim.model
    src_fom = src_model isa RomModel ? src_model.fom_model : src_model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
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
        type_str = T === NodeSet ? "node set" : T === SideSet ? "side set" : "block"
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

function _create_bcs(subsim::SingleDomainSimulation)
    boundary_conditions = Vector{BoundaryCondition}()
    params = subsim.params
    if haskey(params, "boundary conditions") == false
        return boundary_conditions
    end
    input_mesh = params["input_mesh"]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) in bc_params
        for bc_setting_params in bc_type_params
            if bc_type == "Dirichlet"
                # Same "Dirichlet" syntax for both: a "node set" entry applies on
                # a node set, a "side set" entry applies on the nodes of a side set.
                if haskey(bc_setting_params, "side set")
                    boundary_condition = SolidMechanicsSideSetDirichletBoundaryCondition(input_mesh, bc_setting_params)
                elseif haskey(bc_setting_params, "node set")
                    boundary_condition = SolidMechanicsDirichletBoundaryCondition(input_mesh, bc_setting_params)
                else
                    norma_abort("A Dirichlet boundary condition requires either a \"node set\" or a \"side set\" entry.")
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
            elseif bc_type == "Schwarz contact"
                sim = subsim.parent
                sim.controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SolidMechanicsContactSchwarzBoundaryCondition(
                    subsim, coupled_subsim, input_mesh, bc_setting_params
                )
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz overlap" || bc_type == "Schwarz DN nonoverlap" || bc_type == "Schwarz RR nonoverlap" || bc_type == "Schwarz impedance nonoverlap" || bc_type == "Schwarz impedance overlap"
                sim = subsim.parent
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SMCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "OpInf Schwarz overlap"
                sim = subsim.parent
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subsim = sim.subsims[sim.handle_by_name[coupled_subsim_name].id]
                boundary_condition = SMOpInfCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
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

function pair_bc(bc::SolidMechanicsSchwarzBoundaryCondition, bc_index::Int64)
    if bc isa SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
        return nothing
    end
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = coupled_subsim_of(bc).model
    coupled_bcs = coupled_model.boundary_conditions
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name
            if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition &&
               coupled_bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                if bc.is_dirichlet == coupled_bc.is_dirichlet
                    norma_abort("Nonoverlap Schwarz BCs must specify different default BC types (Dirichlet vs Neumann).")
                end
            end
            coupled_bc.is_dirichlet = !bc.is_dirichlet
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
        end
    end
    return nothing
end

