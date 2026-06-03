# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Geometry-only scalar lumped mass for L2 projection of QP fields to nodes:
#   m_i = Σ_e Σ_q N_i(ξ_q) w_q |J(ξ_q)|
# Density-free by design — the same projection works regardless of material
# density per block, which matters for multi-material problems.
function build_recovery_mass_lumped(
    input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::Vector{Float64}
    n_nodes = size(reference, 2)
    m = zeros(Float64, n_nodes)
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = reference[:, node_indices]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                @inbounds for i in 1:num_element_nodes
                    m[node_indices[i]] += Np[i] * dvol
                end
            end
        end
    end
    return m
end

function build_recovery_mass_consistent(
    input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::SparseMatrixCSC{Float64,Int64}
    n_nodes = size(reference, 2)
    blocks = Exodus.read_sets(input_mesh, Block)
    nnz_estimate = 0
    for block in blocks
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        _, num_element_nodes = get_element_dim_nodes(element_type)
        num_block_elements = size(get_block_connectivity(input_mesh, block_id), 1)
        nnz_estimate += num_block_elements * num_element_nodes * num_element_nodes
    end
    rows = Vector{Int64}(undef, nnz_estimate)
    cols = Vector{Int64}(undef, nnz_estimate)
    vals = Vector{Float64}(undef, nnz_estimate)
    idx = 0
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        M_e = Matrix{Float64}(undef, num_element_nodes, num_element_nodes)
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = reference[:, node_indices]
            fill!(M_e, 0.0)
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                @inbounds for i in 1:num_element_nodes
                    NiJxW = Np[i] * dvol
                    for j in 1:num_element_nodes
                        M_e[i, j] += NiJxW * Np[j]
                    end
                end
            end
            @inbounds for i in 1:num_element_nodes
                gi = node_indices[i]
                for j in 1:num_element_nodes
                    idx += 1
                    rows[idx] = gi
                    cols[idx] = node_indices[j]
                    vals[idx] = M_e[i, j]
                end
            end
        end
    end
    return sparse(rows, cols, vals, n_nodes, n_nodes)
end

# Lumped recovery operator.  Match Carina: invert only positive lumped-mass
# entries and leave the rest at zero, so orphan nodes (or sign cancellation in
# higher-order quadrature) yield a zero recovered value there rather than
# aborting the run.
function _build_lumped_recovery(
    input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::LumpedRecovery
    m = build_recovery_mass_lumped(input_mesh, reference, num_int_pts)
    inv_m = zeros(Float64, length(m))
    @inbounds for i in eachindex(m)
        m[i] > 0.0 && (inv_m[i] = 1.0 / m[i])
    end
    return LumpedRecovery(inv_m)
end

# Consistent recovery operator: cache a Cholesky factorization of the SPD
# consistent mass matrix for per-component back-substitution at output time.
function _build_consistent_recovery(
    input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::ConsistentRecovery
    M = build_recovery_mass_consistent(input_mesh, reference, num_int_pts)
    return ConsistentRecovery(M, cholesky(Symmetric(M)))
end

function build_recovery_data(
    kind::Symbol, input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::AbstractRecoveryData
    if kind === :none
        return NoRecovery()
    elseif kind === :lumped
        return _build_lumped_recovery(input_mesh, reference, num_int_pts)
    elseif kind === :consistent
        return _build_consistent_recovery(input_mesh, reference, num_int_pts)
    elseif kind === :both
        return BothRecovery(
            _build_lumped_recovery(input_mesh, reference, num_int_pts),
            _build_consistent_recovery(input_mesh, reference, num_int_pts),
        )
    else
        norma_abort("Unknown nodal recovery method: '$kind' (expected :none, :lumped, :consistent, or :both)")
    end
end

# Project per-QP stress (model.stress[b][e][q], length-6 Voigt: xx, yy, zz, yz, xz, xy)
# onto nodal fields.  For BothRecovery the L2 RHS is assembled once and then
# projected with each mass independently:
#   • model.lumped_recovered_stress    (6 × n_nodes)
#   • model.consistent_recovered_stress (6 × n_nodes)
# For single-type recovery the result is stored in model.recovered_stress.
function recover_stress!(model::SolidMechanics)
    rec = model.recovery_data
    rec isa NoRecovery && return model
    if rec isa BothRecovery
        size(model.lumped_recovered_stress, 1) == 0 && return model
        # Assemble the shared L2 RHS into the lumped buffer first.
        lumped_nodal = model.lumped_recovered_stress
        fill!(lumped_nodal, 0.0)
        _assemble_l2_rhs_stress!(lumped_nodal, model)
        # Copy RHS to the consistent buffer before applying inverse masses.
        consistent_nodal = model.consistent_recovered_stress
        consistent_nodal .= lumped_nodal
        # Apply each mass projection independently.
        _apply_inverse_mass!(lumped_nodal, rec.lumped)
        _apply_inverse_mass!(consistent_nodal, rec.consistent)
    else
        size(model.recovered_stress, 1) == 0 && return model
        nodal = model.recovered_stress
        fill!(nodal, 0.0)
        _assemble_l2_rhs_stress!(nodal, model)
        _apply_inverse_mass!(nodal, rec)
    end
    return model
end

# Project per-QP internal variables (model.state[b][e][q], block-local indexing
# into the material's IV name list) onto nodal fields.  For BothRecovery the
# L2 RHS is assembled once and projected with each mass independently:
#   • model.lumped_recovered_internal_variables    (n_iv × n_nodes)
#   • model.consistent_recovered_internal_variables (n_iv × n_nodes)
# For single-type recovery the result is stored in model.recovered_internal_variables.
function recover_internal_variables!(model::SolidMechanics, all_iv_names::Vector{String})
    rec = model.recovery_data
    rec isa NoRecovery && return model
    isempty(all_iv_names) && return model
    if rec isa BothRecovery
        size(model.lumped_recovered_internal_variables, 1) == 0 && return model
        lumped_nodal = model.lumped_recovered_internal_variables
        fill!(lumped_nodal, 0.0)
        _assemble_l2_rhs_internal_variables!(lumped_nodal, model, all_iv_names)
        consistent_nodal = model.consistent_recovered_internal_variables
        consistent_nodal .= lumped_nodal
        _apply_inverse_mass!(lumped_nodal, rec.lumped)
        _apply_inverse_mass!(consistent_nodal, rec.consistent)
    else
        size(model.recovered_internal_variables, 1) == 0 && return model
        nodal = model.recovered_internal_variables
        fill!(nodal, 0.0)
        _assemble_l2_rhs_internal_variables!(nodal, model, all_iv_names)
        _apply_inverse_mass!(nodal, rec)
    end
    return model
end

# Project per-QP von Mises stress onto a nodal field stored in
# model.recovered_von_mises (1 × n_nodes).  σ_vm is computed at each QP from
# the deviatoric Cauchy stress and projected directly — equivalent to running
# σ_vm = sqrt(3/2 sᵢⱼsᵢⱼ) at the QP level.  Note this is *not* equal to the
# σ_vm derived from the nodal-projected stress tensor (vM is nonlinear in σ).
function recover_von_mises_stress!(model::SolidMechanics)
    rec = model.recovery_data
    rec isa NoRecovery && return model
    if rec isa BothRecovery
        size(model.lumped_recovered_von_mises, 1) == 0 && return model
        lumped_nodal = model.lumped_recovered_von_mises
        fill!(lumped_nodal, 0.0)
        _assemble_l2_rhs_von_mises!(lumped_nodal, model)
        consistent_nodal = model.consistent_recovered_von_mises
        consistent_nodal .= lumped_nodal
        _apply_inverse_mass!(lumped_nodal, rec.lumped)
        _apply_inverse_mass!(consistent_nodal, rec.consistent)
    else
        size(model.recovered_von_mises, 1) == 0 && return model
        nodal = model.recovered_von_mises
        fill!(nodal, 0.0)
        _assemble_l2_rhs_von_mises!(nodal, model)
        _apply_inverse_mass!(nodal, rec)
    end
    return model
end

# Project per-QP deformation gradient F = I + ∇u onto a nodal field stored in
# model.recovered_F (9 × n_nodes), column-major: F_11, F_21, F_31, F_12, ...
# Computed on the fly from displacement (no per-QP storage of F).
function recover_deformation_gradient!(model::SolidMechanics)
    rec = model.recovery_data
    rec isa NoRecovery && return model
    if rec isa BothRecovery
        size(model.lumped_recovered_F, 1) == 0 && return model
        lumped_nodal = model.lumped_recovered_F
        fill!(lumped_nodal, 0.0)
        _assemble_l2_rhs_deformation_gradient!(lumped_nodal, model)
        consistent_nodal = model.consistent_recovered_F
        consistent_nodal .= lumped_nodal
        _apply_inverse_mass!(lumped_nodal, rec.lumped)
        _apply_inverse_mass!(consistent_nodal, rec.consistent)
    else
        size(model.recovered_F, 1) == 0 && return model
        nodal = model.recovered_F
        fill!(nodal, 0.0)
        _assemble_l2_rhs_deformation_gradient!(nodal, model)
        _apply_inverse_mass!(nodal, rec)
    end
    return model
end

function _assemble_l2_rhs_stress!(nodal::Matrix{Float64}, model::SolidMechanics)
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        block_stress = model.stress[block_index]
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = model.reference[:, node_indices]
            element_stress = block_stress[block_element_index]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                qp_stress = element_stress[point]
                @inbounds for i in 1:num_element_nodes
                    n = node_indices[i]
                    NiJxW = Np[i] * dvol
                    for c in 1:6
                        nodal[c, n] += NiJxW * qp_stress[c]
                    end
                end
            end
        end
    end
    return nothing
end

function _assemble_l2_rhs_internal_variables!(
    nodal::Matrix{Float64}, model::SolidMechanics, all_iv_names::Vector{String}
)
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    n_iv = length(all_iv_names)
    for (block_index, block) in enumerate(blocks)
        mat_iv_names = internal_variable_names(model.materials[block_index])
        # Map each global IV index to its local position in this block's
        # material (0 if absent).  Computed once per block.
        block_iv_local = Vector{Int}(undef, n_iv)
        @inbounds for k in 1:n_iv
            idx = findfirst(==(all_iv_names[k]), mat_iv_names)
            block_iv_local[k] = idx === nothing ? 0 : idx
        end
        any(>(0), block_iv_local) || continue
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        block_state = model.state[block_index]
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = model.reference[:, node_indices]
            element_state = block_state[block_element_index]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                qp_state = element_state[point]
                @inbounds for i in 1:num_element_nodes
                    n = node_indices[i]
                    NiJxW = Np[i] * dvol
                    for k in 1:n_iv
                        local_idx = block_iv_local[k]
                        local_idx == 0 && continue
                        nodal[k, n] += NiJxW * qp_state[local_idx]
                    end
                end
            end
        end
    end
    return nothing
end

@inline function von_mises_from_voigt(σ)
    return sqrt(
        0.5 * ((σ[1] - σ[2])^2 + (σ[2] - σ[3])^2 + (σ[3] - σ[1])^2) +
        3.0 * (σ[4]^2 + σ[5]^2 + σ[6]^2),
    )
end

function _assemble_l2_rhs_von_mises!(nodal::Matrix{Float64}, model::SolidMechanics)
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        block_stress = model.stress[block_index]
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = model.reference[:, node_indices]
            element_stress = block_stress[block_element_index]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                σ = element_stress[point]   # Voigt: xx, yy, zz, yz, xz, xy
                σ_vm = von_mises_from_voigt(σ)
                @inbounds for i in 1:num_element_nodes
                    nodal[1, node_indices[i]] += Np[i] * dvol * σ_vm
                end
            end
        end
    end
    return nothing
end

function _assemble_l2_rhs_deformation_gradient!(nodal::Matrix{Float64}, model::SolidMechanics)
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_reference_position = model.reference[:, node_indices]
            element_current_position = element_reference_position + model.displacement[:, node_indices]
            for point in 1:num_points
                Np = N[:, point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference_position')
                dvol = det(dXdξ) * ip_weights[point]
                # F[i,k] = ∂y_i/∂X_k = sum_a y_i(a) ∂N_a/∂X_k; matches src/model.jl
                dNdX = dXdξ \ dNdξ
                F = SMatrix{3,3,Float64,9}(element_current_position * dNdX')
                @inbounds for i in 1:num_element_nodes
                    n = node_indices[i]
                    NiJxW = Np[i] * dvol
                    # Column-major flattening: F[1,1], F[2,1], F[3,1], F[1,2], ...
                    for col in 1:3
                        base = 3 * (col - 1)
                        for row in 1:3
                            nodal[base + row, n] += NiJxW * F[row, col]
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function _apply_inverse_mass!(nodal::Matrix{Float64}, rec::LumpedRecovery)
    inv_m = rec.inv_m
    n_nodes = size(nodal, 2)
    @inbounds for n in 1:n_nodes
        w = inv_m[n]
        for c in 1:size(nodal, 1)
            nodal[c, n] *= w
        end
    end
    return nothing
end

function _apply_inverse_mass!(nodal::Matrix{Float64}, rec::ConsistentRecovery)
    @views for c in 1:size(nodal, 1)
        nodal[c, :] = rec.factor \ nodal[c, :]
    end
    return nothing
end
