# Schwarz coupling boundary conditions for multi-domain simulations.
#
# Includes: overlap, non-overlap (Dirichlet-Neumann), Robin-Robin,
# and contact Schwarz BCs. Projectors, time history interpolation,
# and interface force/displacement transfer.

# ---------------------------------------------------------------------------
# Impedance-matching Robin Schwarz: t + Z u̇ + α W u = g
# ---------------------------------------------------------------------------
#
# g = -t_src_projected + Z * u̇_src_projected + α * W * u_src_projected
#
# The impedance term Z * u̇ absorbs outgoing waves at the interface,
# preventing reflections that cause energy growth with mixed integrators.

get_fom_model(sim::Simulation) = sim.model isa RomModel ? sim.model.fom_model : sim.model

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsImpedanceSchwarzBoundaryCondition)
    Z = bc.impedance
    α = bc.robin_parameter
    W = bc.square_projector
    parent_sim = bc.coupled_subsim.params["parent_simulation"]
    controller = parent_sim.controller
    iter = controller.iteration_number
    coupled_subsim = bc.coupled_subsim
    coupled_index = parent_sim.subsim_name_index_map[coupled_subsim.name]
    this_subsim = bc.subsim
    this_index = parent_sim.subsim_name_index_map[this_subsim.name]

    # Neumann part: -t_src projected
    neumann_force = get_dst_force(bc)

    # Source displacement and velocity, projected to destination
    src_sim = bc.coupled_subsim
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
    coupled_model = bc.coupled_subsim.model
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
# Instead of overwriting DOFs (DBC-DBC), applies impedance-matching force:
#   boundary_force += -t_partner + Z * u̇_partner  (strong interpolation)
# DOFs remain free — waves pass through instead of reflecting.

function _get_impedance_scale(bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition)
    scales = bc.impedance_scale
    if length(scales) == 1
        return scales[1]
    end
    parent_sim = bc.subsim.params["parent_simulation"]
    step = max(1, parent_sim.controller.stop)
    return scales[min(step, length(scales))]
end

function apply_bc_detail(model::SolidMechanics, bc::SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition)
    Z = bc.impedance * _get_impedance_scale(bc)
    α = bc.robin_parameter
    W = bc.square_projector

    get_coupled_field = if bc.coupled_subsim.model isa SolidMechanics
        (field -> getfield(bc.coupled_subsim.model, field))
    else
        (field -> getfield(bc.coupled_subsim.model.fom_model, field))
    end

    coupled_displacement = get_coupled_field(:displacement)
    coupled_reference = get_coupled_field(:reference)
    coupled_velocity  = get_coupled_field(:velocity)
    coupled_internal  = get_coupled_field(:internal_force)

    global_from_local_map = bc.global_from_local_map
    unique_node_indices = unique(bc.side_set_node_indices)

    # Interpolate partner displacement and velocity at each overlap boundary node
    num_dst_nodes = length(unique_node_indices)
    partner_velo = zeros(3, num_dst_nodes)
    partner_disp = zeros(3, num_dst_nodes)

    for (i, node_index) in enumerate(unique_node_indices)
        coupled_node_indices = bc.coupled_nodes_indices[i]
        N = bc.interpolation_function_values[i]
        for comp in 1:3
            partner_velo[comp, i] = sum(coupled_velocity[comp, coupled_node_indices] .* N)
            partner_disp[comp, i] = sum(coupled_displacement[comp, coupled_node_indices] .* N)
        end
    end

    # RHS (partner contribution only): f = W * (Z * u̇_partner + α * u_partner)
    # The self terms (Z*W*u̇_self + α*W*u_self) are handled by
    # apply_impedance_bcs_internal_force! called at each Newton iteration.
    for comp in 1:3
        Z_vdot = Z * partner_velo[comp, :]
        alpha_u = α * partner_disp[comp, :]
        rhs = W * (Z_vdot + alpha_u)
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
        Z = bc.impedance * _get_impedance_scale(bc)
        α = bc.robin_parameter
        W = bc.square_projector
        c_v = if integrator isa Newmark
            integrator.γ / (integrator.β * integrator.time_step)
        else
            0.0
        end
        coeff = Z * c_v + α
        global_from_local_map = bc.global_from_local_map
        for (i_local, i_global) in enumerate(global_from_local_map)
            for (j_local, j_global) in enumerate(global_from_local_map)
                w_ij = coeff * W[i_local, j_local]
                for comp in 1:3
                    dof_i = 3 * (i_global - 1) + comp
                    dof_j = 3 * (j_global - 1) + comp
                    K_io[dof_i, dof_j] += w_ij
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
    return nothing
end

# Compute the self-impedance internal force: W*(Z*u̇_self + α*u_self)
# Returns a separate force vector rather than modifying model.internal_force,
# because get_dst_force reads model.internal_force for traction transfer and
# must see only the elastic internal force.
function build_impedance_schwarz_force(model::SolidMechanics)
    num_dofs = 3 * size(model.reference, 2)
    f = zeros(num_dofs)
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsImpedanceSchwarzBoundaryCondition ||
           bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
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
    parent_sim = bc.coupled_subsim.params["parent_simulation"]
    num_domains = parent_sim.num_domains
    controller = parent_sim.controller
    iter = controller.iteration_number
    coupled_subsim = bc.coupled_subsim
    coupled_index = parent_sim.subsim_name_index_map[coupled_subsim.name]
    this_subsim = bc.subsim
    this_index = parent_sim.subsim_name_index_map[this_subsim.name]
    # Neumann part of Robin RHS: -t_src projected (= get_dst_force which negates internal_force)
    neumann_force = get_dst_force(bc)
    # Displacement part of Robin RHS: α * W * u_src_projected
    # Compute source displacement (current - reference) and project to destination
    src_sim = bc.coupled_subsim
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
    else #left subdomain
      #Set g and lambda to 0 for iter = 0
      if (iter == 0)
        n = length(model.boundary_force)
        g = zeros(n)
      else  
        #this plays role of g_1.  hijack lambda_disp to store it.
        #In particular, we set g to past lambda_disp
        g = controller.lambda_disp[coupled_index]  
      end
      #println("IKT g norm = ", norm(g)) 
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
    get_coupled_field = if bc.coupled_subsim.model isa SolidMechanics
        (field -> getfield(bc.coupled_subsim.model, field))
    else
        (field -> getfield(bc.coupled_subsim.model.fom_model, field))
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
    src_model = if bc.coupled_subsim.model isa SolidMechanics
        bc.coupled_subsim.model
    else
        bc.coupled_subsim.model.fom_model
    end

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
    parent_sim = bc.coupled_subsim.params["parent_simulation"]
    controller = parent_sim.controller

    # Skip application if contact is inactive
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition && !controller.active_contact
        return nothing
    end

    coupled_subsim = bc.coupled_subsim
    integrator = coupled_subsim.integrator
    coupled_model = coupled_subsim.model

    # Save current state (copy data, not reference — aliased integrators share memory with model)
    saved_disp = copy(integrator.displacement)
    saved_velo = copy(integrator.velocity)
    saved_acce = copy(integrator.acceleration)
    saved_∂Ω_f = get_internal_force(coupled_model)

    # Fetch interpolation inputs
    time = model.time
    coupled_index = parent_sim.subsim_name_index_map[coupled_subsim.name]
    num_dofs = length(coupled_subsim.model.free_dofs)

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
        interp_disp = zeros(num_dofs)
        interp_velo = zeros(num_dofs)
        interp_acce = zeros(num_dofs)
        interp_∂Ω_f = zeros(num_dofs)
    end


    # Assign interpolated force
    set_internal_force!(coupled_model, interp_∂Ω_f)

    # Apply relaxed update if needed
    if bc isa SolidMechanicsContactSchwarzBoundaryCondition || bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
        θ = controller.relaxation_parameter

        iter = controller.iteration_number
        λ_u_prev = iter < 1 ? interp_disp : controller.lambda_disp[coupled_index]
        λ_v_prev = iter < 1 ? interp_velo : controller.lambda_velo[coupled_index]
        λ_a_prev = iter < 1 ? interp_acce : controller.lambda_acce[coupled_index]

        # Classical relaxation
        controller.lambda_disp[coupled_index] = θ * interp_disp + (1 - θ) * λ_u_prev
        controller.lambda_velo[coupled_index] = θ * interp_velo + (1 - θ) * λ_v_prev
        controller.lambda_acce[coupled_index] = θ * interp_acce + (1 - θ) * λ_a_prev

        integrator.displacement .= controller.lambda_disp[coupled_index]
        integrator.velocity .= controller.lambda_velo[coupled_index]
        integrator.acceleration .= controller.lambda_acce[coupled_index]
    else
        integrator.displacement .= interp_disp
        integrator.velocity .= interp_velo
        integrator.acceleration .= interp_acce
    end

    # For ROM coupled subdomains, reconstruct FOM displacement from reduced state before reading it
    if coupled_subsim.model isa RomModel
        reconstruct_fom_fields!(integrator, coupled_subsim.solver, coupled_subsim.model)
    end
    apply_bc_detail(model, bc)

    # Restore previous state (in-place to keep alias intact)
    integrator.displacement .= saved_disp
    integrator.velocity .= saved_velo
    integrator.acceleration .= saved_acce
    set_internal_force!(coupled_model, saved_∂Ω_f)
    if coupled_subsim.model isa RomModel
        reconstruct_fom_fields!(integrator, coupled_subsim.solver, coupled_subsim.model)
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

function compute_neumann_projector(dst_model::SolidMechanics, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_model = dst_bc.coupled_subsim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    H = get_square_projection_matrix(src_model, src_bc)
    L = get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_bc.neumann_projector = L * (H \ I)
    return nothing
end

function compute_dirichlet_projector(dst_model::SolidMechanics, dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_model = dst_bc.coupled_subsim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    W = get_square_projection_matrix(dst_model, dst_bc)
    L = get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_bc.dirichlet_projector = (W \ I) * L
    return nothing
end

function get_dst_force(dst_bc::SolidMechanicsSchwarzBoundaryCondition)
    src_sim = dst_bc.coupled_subsim
    src_model = src_sim.model
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_force = src_model.internal_force
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
    src_sim = dst_bc.coupled_subsim
    src_model = get_fom_model(src_sim)
    src_bc_index = dst_bc.coupled_bc_index
    src_bc = src_model.boundary_conditions[src_bc_index]
    src_global_from_local_map = src_bc.global_from_local_map
    num_src_nodes = length(src_global_from_local_map)
    src_curr = zeros(3, num_src_nodes)
    src_refe = zeros(3, num_src_nodes)
    src_velo = zeros(3, num_src_nodes)
    src_acce = zeros(3, num_src_nodes)
    for (i_local, i_global) in enumerate(src_global_from_local_map)
        src_curr[:, i_local] = src_model.reference[:, i_global] + src_model.displacement[:, i_global]
        src_refe[:, i_local] = src_model.reference[:, i_global]
        src_velo[:, i_local] = src_model.velocity[:, i_global]
        src_acce[:, i_local] = src_model.acceleration[:, i_global]
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

function extract_value(value::Real)
    return Float64(value)
end

function extract_value(symbol::Num)
    return Float64(symbol.val)
end

function create_bcs(params::Parameters)
    boundary_conditions = Vector{BoundaryCondition}()
    if haskey(params, "boundary conditions") == false
        return boundary_conditions
    end
    input_mesh = params["input_mesh"]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) in bc_params
        for bc_setting_params in bc_type_params
            if bc_type == "Dirichlet"
                boundary_condition = SolidMechanicsDirichletBoundaryCondition(input_mesh, bc_setting_params)
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
                sim = params["parent_simulation"]
                sim.controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SolidMechanicsContactSchwarzBoundaryCondition(
                    coupled_subsim, input_mesh, bc_setting_params
                )
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz overlap" || bc_type == "Schwarz DN nonoverlap" || bc_type == "Schwarz RR nonoverlap" || bc_type == "Schwarz impedance nonoverlap" || bc_type == "Schwarz impedance overlap"
                sim = params["parent_simulation"]
                subsim_name = params["name"]
                subdomain_index = sim.subsim_name_index_map[subsim_name]
                subsim = sim.subsims[subdomain_index]
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMCouplingSchwarzBC(subsim, coupled_subsim, input_mesh, bc_type, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "OpInf Schwarz overlap"
                sim = params["parent_simulation"]
                subsim_name = params["name"]
                subdomain_index = sim.subsim_name_index_map[subsim_name]
                subsim = sim.subsims[subdomain_index]
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
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
            # Parse and expand derivatives only once
            if ic_type == "displacement"
                disp_num = eval(Meta.parse(expression))
                velo_num = expand_derivatives(D(disp_num))
            elseif ic_type == "velocity"
                disp_num = nothing  # Not needed
                velo_num = eval(Meta.parse(expression))
            else
                norma_abort(
                    "Invalid initial condition type: '$ic_type'. Supported types are: displacement or velocity."
                )
            end
            for node_index in node_set_node_indices
                values = Dict(
                    t => model.time,
                    x => model.reference[1, node_index],
                    y => model.reference[2, node_index],
                    z => model.reference[3, node_index],
                )
                disp_val = disp_num === nothing ? 0.0 : extract_value(substitute(disp_num, values))
                velo_val = velo_num === nothing ? 0.0 : extract_value(substitute(velo_num, values))
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
    coupled_model = bc.coupled_subsim.model
    coupled_bcs = coupled_model.boundary_conditions
    for (coupled_bc_index, coupled_bc) in enumerate(coupled_bcs)
        if coupled_bc_name == coupled_bc.name
            coupled_bc.is_dirichlet = !bc.is_dirichlet
            bc.coupled_bc_index = coupled_bc_index
            coupled_bc.coupled_bc_index = bc_index
        end
    end
    return nothing
end

function pair_bc(bc::SolidMechanicsRobinSchwarzBoundaryCondition, bc_index::Int64)
    coupled_bc_name = bc.coupled_bc_name
    coupled_model = bc.coupled_subsim.model
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
