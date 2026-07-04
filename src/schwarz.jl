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

# Returns the relaxation factor θ applied to interp_disp for this Schwarz
# iterate. Fixed mode returns the user-configured constant; Aitken-recursive
# mode uses Irons–Tuck with the previous residual stored on the controller.
function relaxation_aitken_recursive_theta!(
    controller::MultiDomainTimeController,
    slot::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    if controller.relaxation_method !== :aitken_recursive
        return controller.relaxation_parameter
    end
    aitken_N0 = controller.aitken_N0
    if iter < aitken_N0
        controller.aitken_theta_disp[slot] = controller.relaxation_parameter
        controller.aitken_prev_residual_disp[slot] = Float64[]
        norma_logf(1, :schwarz, "Aitken-recursive θ[slot=%d, iter=%d] = %.4f", slot, iter, 1.0)
        return 1.0
    end
    residual = interp_disp .- lambda_prev
    prev_residual = controller.aitken_prev_residual_disp[slot]
    θ_prev = controller.aitken_theta_disp[slot]
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
    controller.aitken_prev_residual_disp[slot] = residual
    controller.aitken_theta_disp[slot] = θ
    norma_logf(1, :schwarz, "Aitken-recursive θ[slot=%d, iter=%d] = %.4f", slot, iter, θ)
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
    slot::Int,
    iter::Int,
    interp_disp::AbstractVector{Float64},
    lambda_prev::AbstractVector{Float64},
)
    residual = interp_disp .- lambda_prev                 # r^(n) = T(g^(n)) - g^(n)
    prev_residual = controller.aitken_prev_residual_disp[slot]   # r^(n-1)
    prev_lambda = controller.aitken_prev_lambda_disp[slot]       # g^(n-1)
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
    controller.aitken_prev_residual_disp[slot] = residual
    controller.aitken_prev_lambda_disp[slot] = copy(lambda_prev)
    norma_logf(1, :schwarz, "Aitken-secant θ[slot=%d, iter=%d] = %.4f", slot, iter, θ)
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
function recover_interface_kinematics!(controller, slot, integrator, interp_velo, interp_acce)
    if integrator isa Newmark
        Δt = integrator.time_step
        β = integrator.β
        γ = integrator.γ
        u = controller.lambda_disp[slot]
        a = (u .- integrator.disp_pre) ./ (β * Δt * Δt)
        controller.lambda_acce[slot] = a
        controller.lambda_velo[slot] = integrator.velo_pre .+ (γ * Δt) .* a
    else
        controller.lambda_velo[slot] = interp_velo
        controller.lambda_acce[slot] = interp_acce
    end
    return nothing
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsImpedanceSchwarzBoundaryCondition)
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
        if (iter == 0)
            n = length(model.boundary_force)
            g = zeros(n)
        else
            g = controller.lambda_disp[coupled_index]
        end
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
                relaxation_aitken_secant_theta!(controller, coupled_index, iter, candidate, g) :
                relaxation_aitken_recursive_theta!(controller, coupled_index, iter, candidate, g)
        end
        controller.lambda_disp[coupled_index] = copy(model.boundary_force)
        for comp in 1:3
            Z_W_vdot = Z * (W * dst_velo[comp, :])
            α_W_u = α * (W * dst_disp[comp, :])
            for (i_local, i_global) in enumerate(global_from_local_map)
                dof_i = 3 * (i_global - 1) + comp
                rhs_i = neumann_force[3 * (i_local - 1) + comp] + Z_W_vdot[i_local] + α_W_u[i_local]
                controller.lambda_disp[coupled_index][dof_i] += (1 - theta) * g[dof_i] + theta * rhs_i
                model.boundary_force[dof_i] = controller.lambda_disp[coupled_index][dof_i]
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
        bc isa SolidMechanicsImpedanceSchwarzBoundaryCondition || continue
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

function pair_bc(bc::SolidMechanicsImpedanceSchwarzBoundaryCondition, bc_index::Int64)
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
    dst_model::SolidMechanics, dst_bc::SolidMechanicsImpedanceSchwarzBoundaryCondition
)
    compute_dirichlet_projector(dst_model, dst_bc)
    compute_neumann_projector(dst_model, dst_bc)
    dst_bc.square_projector = get_square_projection_matrix(dst_model, dst_bc)
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
    partner_Zvdot = zeros(3, num_dst_nodes)
    for i in 1:num_dst_nodes
        # P/S-split tensor impedance: Z_p on the normal velocity component,
        # Z_s on the tangential components (Lysmer-Kuhlemeyer split).
        n = normals[:, i]
        v = partner_velo[:, i]
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
        # source node j: 3x3 block W_ij (c_v T(n_j) + α I).
        for (j_local, j_global) in enumerate(global_from_local_map)
            n = normals[:, j_local]
            T = zeros(3, 3)
            for a in 1:3, b in 1:3
                T[a, b] = c_v * ((Z_p - Z_s) * n[a] * n[b] + (a == b ? Z_s : 0.0)) + (a == b ? α : 0.0)
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
    if dst_bc.transfer_mode == "variational"
        # Variational transfer: P = W \\ L is the L2(Γ)-orthogonal projection of
        # the partner fields onto this side's trace space — non-expansive in
        # L2 for any quadrature, and constants transfer exactly (partner
        # partition of unity gives L·1 = W·1). On a node-aligned interface
        # L = W row-permuted and P reduces to the identity selection.
        coupled_solid = get_fom_model(coupled_subsim_of(dst_bc))
        L = get_overlap_rectangular_projection_matrix(
            dst_model, dst_bc, coupled_solid, dst_bc.coupled_block_name, dst_bc.search_tolerance;
            subdivisions=dst_bc.transfer_subdivisions,
        )
        dst_bc.variational_projector = dst_bc.square_projector \ L
    else
        dst_bc.variational_projector = Matrix{Float64}(undef, 0, 0)
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
    if size(bc.variational_projector, 1) > 0
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
        if bc.partner_traction_mode == "consistent traction"
            norma_abort(
                "`partner traction: consistent traction` requires a node-aligned (conforming) " *
                "interface, but boundary nodes of side set \"$(bc.name)\" do not coincide " *
                "with partner mesh nodes. Use \"auto\" or \"recovered stress\".",
            )
        end
        return nothing  # auto: fall back to recovered stress
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
        element_nodes, element_block, element_index, element_rows, block_element_type, lumped_mass
    )
    norma_log(
        0,
        :schwarz,
        "Impedance overlap side set \"$(bc.name)\": node-aligned interface, using " *
        "consistent partner traction ($(length(element_nodes)) partner elements).",
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
    num_dst_nodes = length(bc.global_from_local_map)
    weak_traction = zeros(3, num_dst_nodes)
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
    return weak_traction
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
            Zv = zeros(3, num_nodes)
            for (i_local, i_global) in enumerate(global_from_local_map)
                n = normals[:, i_local]
                v = model.velocity[:, i_global]
                u = model.displacement[:, i_global]
                vn = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
                Zv[:, i_local] = (Z_p * vn) .* n .+ Z_s .* (v .- vn .* n) .+ α .* u
            end
            for comp in 1:3
                f_imp = W * Zv[comp, :]
                for (i_local, i_global) in enumerate(global_from_local_map)
                    dof_i = 3 * (i_global - 1) + comp
                    f[dof_i] += f_imp[i_local]
                end
            end
        elseif bc isa SolidMechanicsImpedanceSchwarzBoundaryCondition
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
                "velocity_jump_rms,traction_mismatch_rms",
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
                dv = partner_velo[:, i] .- v
                dvn = n[1] * dv[1] + n[2] * dv[2] + n[3] * dv[3]
                dash[:, i] = (Z_p * dvn) .* n .+ Z_s .* (dv .- dvn .* n)
                robin[:, i] = α .* (partner_disp[:, i] .- u)
                velo[:, i] = v
                vjump2 += dot(dv, dv)
                tmis2 += dot(trac_p[:, i] - t_o, trac_p[:, i] - t_o)
            end
            P_trac = 0.0
            P_dash = 0.0
            P_robin = 0.0
            # Traction power from the traction the BC actually applies: the
            # consistent-traction weak traction when the patch is active, the
            # recovered-stress interpolant otherwise. The recovered-stress
            # traction mismatch below stays informational in both modes.
            if bc.traction_patch !== nothing
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
            P_total = P_trac + P_dash + P_robin
            vjump_rms = sqrt(vjump2 / num_nodes)
            tmis_rms = sqrt(tmis2 / num_nodes)
            open(csv_path, "a") do io
                println(
                    io,
                    "$t,$subsim_index,$(bc.side_set_id),$P_total,$P_trac,$P_dash,$P_robin," *
                    "$vjump_rms,$tmis_rms",
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

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsRobinSchwarzBoundaryCondition)
    α = bc.robin_parameter
    W = bc.square_projector
    parent_sim = bc.parent
    num_domains = parent_sim.num_domains
    controller = parent_sim.controller
    iter = controller.iteration_number
    coupled_index = bc.coupled_handle.id
    this_index = bc.self_handle.id
    # Neumann part of Robin RHS: -t_src projected (= get_dst_force which negates internal_force)
    neumann_force = get_dst_force(bc)
    # Displacement part of Robin RHS: α * W * u_src_projected
    # Compute source displacement (current - reference) and project to destination
    src_sim = coupled_subsim_of(bc)
    src_model = get_fom_model(src_sim)
    src_bc_index = bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
    num_src_nodes = length(src_global_from_local_map)
    src_disp = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_disp[:, i_local] = src_model.displacement[:, i_global]
    end
    dirichlet_projector = bc.dirichlet_projector
    num_dst_nodes = size(dirichlet_projector, 1)
    dst_disp = zeros(3, num_dst_nodes)
    for i in 1:3
        dst_disp[i, :] = dirichlet_projector * src_disp[i, :]
    end
    global_from_local_map = bc.global_from_local_map
    #Get relaxation parameter from input file
    theta = controller.relaxation_parameter
    if (this_index < coupled_index) #right subdomain
      for comp in 1:3
          alpha_W_u = α * (W * dst_disp[comp, :])
          for (i_local, i_global) in enumerate(global_from_local_map)
              dof_i = 3 * (i_global - 1) + comp
              model.boundary_force[dof_i] += neumann_force[3 * (i_local - 1) + comp] + alpha_W_u[i_local]
          end
      end
    else #left subdomain (relaxation applied on this side)
      #Set g and lambda to 0 for iter = 0
      if (iter == 0)
        n = length(model.boundary_force)
        g = zeros(n)
      else
        #this plays role of g_1.  hijack lambda_disp to store it.
        #In particular, we set g to past lambda_disp
        g = controller.lambda_disp[coupled_index]
      end
      # Optional dynamic Aitken relaxation factor. The fixed-point iterate is the
      # Robin RHS stored in lambda_disp; its unrelaxed (theta = 1) candidate is
      # T(g) = boundary_force + coupling term, from which the residual r = T(g) - g
      # is formed. The fixed-theta blend below is otherwise unchanged.
      if controller.relaxation_method === :aitken_recursive || controller.relaxation_method === :aitken_secant
          candidate = copy(model.boundary_force)
          for comp in 1:3
              alpha_W_u = α * (W * dst_disp[comp, :])
              for (i_local, i_global) in enumerate(global_from_local_map)
                  dof_i = 3 * (i_global - 1) + comp
                  candidate[dof_i] += neumann_force[3 * (i_local - 1) + comp] + alpha_W_u[i_local]
              end
          end
          theta = controller.relaxation_method === :aitken_secant ?
              relaxation_aitken_secant_theta!(controller, coupled_index, iter, candidate, g) :
              relaxation_aitken_recursive_theta!(controller, coupled_index, iter, candidate, g)
      end
      #initialize lambda_disp  = model.boundary_force
      controller.lambda_disp[coupled_index] = copy(model.boundary_force)
      for comp in 1:3
          alpha_W_u = α * (W * dst_disp[comp, :])
          for (i_local, i_global) in enumerate(global_from_local_map)
              dof_i = 3 * (i_global - 1) + comp
              controller.lambda_disp[coupled_index][dof_i] += (1 - theta) * g[dof_i] + theta * (neumann_force[3 * (i_local - 1) + comp] + alpha_W_u[i_local])
              model.boundary_force[dof_i] = controller.lambda_disp[coupled_index][dof_i]
          end
      end
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
    use_predictor = (
        bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition ||
        bc isa SolidMechanicsRobinSchwarzBoundaryCondition
    ) &&
    controller.use_interface_predictor &&
    controller.iteration_number == 0 &&
    !isempty(controller.predictor_disp[coupled_index])
    if !isempty(time_hist)
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
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition || bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
        iter = controller.iteration_number
        λ_u_prev = iter < 1 ? interp_disp : controller.lambda_disp[coupled_index]
        λ_v_prev = iter < 1 ? interp_velo : controller.lambda_velo[coupled_index]
        λ_a_prev = iter < 1 ? interp_acce : controller.lambda_acce[coupled_index]

        if controller.relaxation_method === :aitken_secant
            # Relax the single d-form interface unknown (displacement); recover
            # velocity and acceleration consistently from it (see functions above).
            θ = relaxation_aitken_secant_theta!(controller, coupled_index, iter, interp_disp, λ_u_prev)
            controller.lambda_disp[coupled_index] = θ * interp_disp + (1 - θ) * λ_u_prev
            recover_interface_kinematics!(controller, coupled_index, coupled_integrator, interp_velo, interp_acce)
        else
            θ = relaxation_aitken_recursive_theta!(controller, coupled_index, iter, interp_disp, λ_u_prev)

            controller.lambda_disp[coupled_index] = θ * interp_disp + (1 - θ) * λ_u_prev
            controller.lambda_velo[coupled_index] = θ * interp_velo + (1 - θ) * λ_v_prev
            controller.lambda_acce[coupled_index] = θ * interp_acce + (1 - θ) * λ_a_prev
        end

        coupled_integrator.displacement .= controller.lambda_disp[coupled_index]
        coupled_integrator.velocity .= controller.lambda_velo[coupled_index]
        coupled_integrator.acceleration .= controller.lambda_acce[coupled_index]
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
    L = get_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
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
    L = get_rectangular_projection_matrix(dst_fom, dst_bc, src_fom, src_bc)
    dst_bc.dirichlet_projector = (W \ I) * L
    return nothing
end

function get_dst_force(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = coupled_subsim_of(dst_bc)
    src_model = src_sim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_force = get_internal_force(src_model)
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
    if bc isa SolidMechanicsOverlapSchwarzBoundaryCondition || bc isa SolidMechanicsOpInfOverlapSchwarzBoundaryCondition
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

function pair_bc(bc::SolidMechanicsRobinSchwarzBoundaryCondition, bc_index::Int64)
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

function compute_robin_schwarz_projectors!(
    dst_model::SolidMechanics, dst_bc::SolidMechanicsRobinSchwarzBoundaryCondition
)
    compute_dirichlet_projector(dst_model, dst_bc)
    compute_neumann_projector(dst_model, dst_bc)
    dst_bc.square_projector = get_square_projection_matrix(dst_model, dst_bc)
    return nothing
end

function build_robin_schwarz_stiffness(model::SolidMechanics)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_rs = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinSchwarzBoundaryCondition || continue
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
