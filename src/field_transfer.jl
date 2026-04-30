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
    for dst_block in dst_blocks
        dst_block_id = dst_block.id
        dst_et_str = Exodus.read_block_parameters(dst_mesh, dst_block_id)[1]
        dst_et = element_type_from_string(dst_et_str)
        dst_npts = default_num_int_pts(dst_et)
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
