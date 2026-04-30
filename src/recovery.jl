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

build_recovery_mass_lumped(model::SolidMechanics) =
    build_recovery_mass_lumped(model.mesh, model.reference, model.num_int_pts)

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

build_recovery_mass_consistent(model::SolidMechanics) =
    build_recovery_mass_consistent(model.mesh, model.reference, model.num_int_pts)

function build_recovery_data(
    kind::Symbol, input_mesh::ExodusDatabase, reference::Matrix{Float64}, num_int_pts::Vector{Int}
)::AbstractRecoveryData
    if kind === :none
        return NoRecovery()
    elseif kind === :lumped
        m = build_recovery_mass_lumped(input_mesh, reference, num_int_pts)
        # Match Carina: invert only positive entries; leave the rest at zero so
        # orphan nodes (or sign cancellation in higher-order quadrature) yield a
        # zero recovered value at those nodes rather than aborting the run.
        n = length(m)
        inv_m = zeros(Float64, n)
        @inbounds for i in 1:n
            if m[i] > 0.0
                inv_m[i] = 1.0 / m[i]
            end
        end
        return LumpedRecovery(inv_m)
    elseif kind === :consistent
        M = build_recovery_mass_consistent(input_mesh, reference, num_int_pts)
        factor = cholesky(Symmetric(M))
        return ConsistentRecovery(M, factor)
    else
        norma_abort("Unknown stress recovery kind: '$kind' (expected :none, :lumped, or :consistent)")
    end
end

# Project per-QP stress (model.stress[b][e][q], length-6 Voigt: xx, yy, zz, yz, xz, xy)
# onto a nodal field stored in model.recovered_stress (6 × n_nodes).
function recover_stress!(model::SolidMechanics)
    rec = model.recovery_data
    rec isa NoRecovery && return model
    nodal = model.recovered_stress
    fill!(nodal, 0.0)
    _assemble_l2_rhs_stress!(nodal, model)
    _apply_inverse_mass!(nodal, rec)
    return model
end

# Project per-QP internal variables (model.state[b][e][q], block-local indexing
# into the material's IV name list) onto a nodal field stored in
# model.recovered_internal_variables (n_iv × n_nodes), where n_iv is the union
# of IV names across all materials.  Blocks whose material lacks a given IV
# contribute zero to that component's RHS but still contribute their geometric
# tributary to the (shared) projection mass.
function recover_internal_variables!(model::SolidMechanics, all_iv_names::Vector{String})
    rec = model.recovery_data
    rec isa NoRecovery && return model
    isempty(all_iv_names) && return model
    nodal = model.recovered_internal_variables
    fill!(nodal, 0.0)
    _assemble_l2_rhs_internal_variables!(nodal, model, all_iv_names)
    _apply_inverse_mass!(nodal, rec)
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
