# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# L2-project a nodal field from src_model onto dst_model under the assumption
# that src and dst cover the same physical geometry (element decomposition
# and order may differ).  At each destination QP, the host source element is
# located via find_point_in_mesh (looping src blocks), the source field is
# interpolated at the recovered parametric coordinates, and the result is
# L2-projected onto the destination basis using dst_model.recovery_data.
function transfer_field(
    src_model::SolidMechanics,
    src_nodal_field::AbstractMatrix{Float64},
    dst_model::SolidMechanics;
    src_blocks::Union{Nothing,Vector{Int}}=nothing,
    tol::Float64=1.0e-6,
)::Matrix{Float64}
    if dst_model.recovery_data isa NoRecovery
        norma_abort(
            "transfer_field: destination model must have 'stress recovery' set to 'lumped' or 'consistent'",
        )
    end

    C = size(src_nodal_field, 1)
    n_dst_nodes = size(dst_model.reference, 2)
    out = zeros(C, n_dst_nodes)

    src_block_ids = src_blocks === nothing ?
        [Int(b.id) for b in Exodus.read_sets(src_model.mesh, Block)] : src_blocks

    src_element_types = Dict{Int,ElementType}()
    for bid in src_block_ids
        et_str = Exodus.read_block_parameters(src_model.mesh, Int32(bid))[1]
        src_element_types[bid] = element_type_from_string(et_str)
    end

    dst_mesh = dst_model.mesh
    dst_blocks = Exodus.read_sets(dst_mesh, Block)
    point_buf = Vector{Float64}(undef, 3)
    for (dst_block_index, dst_block) in enumerate(dst_blocks)
        dst_block_id = dst_block.id
        dst_et_str = Exodus.read_block_parameters(dst_mesh, dst_block_id)[1]
        dst_et = element_type_from_string(dst_et_str)
        dst_npts = dst_model.num_int_pts[dst_block_index]
        N_dst, dN_dst, dst_w = isoparametric(dst_et, dst_npts)
        dst_conn = get_block_connectivity(dst_mesh, dst_block_id)
        n_dst_elem, n_dst_nodes_per_elem = size(dst_conn)
        for e in 1:n_dst_elem
            conn_idx = ((e - 1) * n_dst_nodes_per_elem + 1):(e * n_dst_nodes_per_elem)
            node_indices = dst_conn[conn_idx]
            elem_ref_pos = dst_model.reference[:, node_indices]
            for q in 1:dst_npts
                Np = N_dst[:, q]
                dNdξ = dN_dst[:, :, q]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * elem_ref_pos')
                dvol = det(dXdξ) * dst_w[q]
                point_buf[1] = elem_ref_pos[1, :] ⋅ Np
                point_buf[2] = elem_ref_pos[2, :] ⋅ Np
                point_buf[3] = elem_ref_pos[3, :] ⋅ Np
                src_node_indices = Int64[]
                src_ξ = zeros(3)
                found = false
                src_bid_found = 0
                for bid in src_block_ids
                    src_node_indices, src_ξ, found = find_point_in_mesh(point_buf, src_model, bid, tol)
                    if found
                        src_bid_found = bid
                        break
                    end
                end
                if !found
                    norma_abortf(
                        "transfer_field: destination QP at (%.4e, %.4e, %.4e) not found in any source block",
                        point_buf[1], point_buf[2], point_buf[3],
                    )
                end
                src_et = src_element_types[src_bid_found]
                N_src, _, _ = interpolate(src_et, src_ξ)
                v = @view(src_nodal_field[:, src_node_indices]) * N_src
                @inbounds for i in 1:n_dst_nodes_per_elem
                    n = node_indices[i]
                    NiJxW = Np[i] * dvol
                    for c in 1:C
                        out[c, n] += NiJxW * v[c]
                    end
                end
            end
        end
    end

    _apply_inverse_mass!(out, dst_model.recovery_data)
    return out
end

# Transfer per-QP internal variables between meshes covering the same
# physical geometry.  IVs common to both src and dst materials are L2-
# projected from src QPs to src nodes (via src.recovery_data), L2-
# transferred to dst nodes (via dst.recovery_data), and then sampled at
# each dst QP and written into dst.state[b][e][q][k].  IVs in dst not in
# src remain at their initial values (set in the model constructor); IVs
# in src not in dst are silently dropped.
function _transfer_qp_internal_variables!(dst::SolidMechanics, src::SolidMechanics)
    src_names = collect_internal_variable_names(src.materials)
    dst_names = collect_internal_variable_names(dst.materials)
    common_set = intersect(Set(src_names), Set(dst_names))
    isempty(common_set) && return nothing
    common_list = [n for n in src_names if n in common_set]   # preserve src order

    src_iv_nodal = zeros(length(common_list), size(src.reference, 2))
    _assemble_l2_rhs_internal_variables!(src_iv_nodal, src, common_list)
    _apply_inverse_mass!(src_iv_nodal, src.recovery_data)

    dst_iv_nodal = transfer_field(src, src_iv_nodal, dst)

    _sample_nodal_at_dst_qps!(dst, dst_iv_nodal, common_list)

    # Mirror state into state_old so the next evaluate's read of state_old
    # reflects the transferred values.
    dst.state_old = deepcopy(dst.state)
    return nothing
end

function _sample_nodal_at_dst_qps!(
    dst::SolidMechanics,
    nodal_iv::AbstractMatrix{Float64},
    common_list::Vector{String},
)
    input_mesh = dst.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = dst.num_int_pts[block_index]
        N, _, _ = isoparametric(element_type, num_points)

        mat_iv_names = internal_variable_names(dst.materials[block_index])
        # block_iv_local[k] is the local IV index in this block's material
        # for the k-th name in common_list (0 if absent).
        block_iv_local = [something(findfirst(==(n), mat_iv_names), 0) for n in common_list]
        any(>(0), block_iv_local) || continue

        element_block_connectivity = get_block_connectivity(input_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        block_state = dst.state[block_index]
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = element_block_connectivity[connectivity_indices]
            element_state = block_state[block_element_index]
            for point in 1:num_points
                Np = N[:, point]
                qp_state = element_state[point]
                @inbounds for k in eachindex(common_list)
                    local_idx = block_iv_local[k]
                    local_idx == 0 && continue
                    val = 0.0
                    for a in 1:num_element_nodes
                        val += Np[a] * nodal_iv[k, node_indices[a]]
                    end
                    qp_state[local_idx] = val
                end
            end
        end
    end
    return nothing
end
