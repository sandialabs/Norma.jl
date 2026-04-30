# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Geometry-only scalar lumped mass for L2 projection of QP fields to nodes:
#   m_i = Σ_e Σ_q N_i(ξ_q) w_q |J(ξ_q)|
# Density-free by design — the same projection works regardless of material
# density per block, which matters for multi-material problems.
function build_recovery_mass_lumped(input_mesh::ExodusDatabase, reference::Matrix{Float64})::Vector{Float64}
    n_nodes = size(reference, 2)
    m = zeros(Float64, n_nodes)
    blocks = Exodus.read_sets(input_mesh, Block)
    for block in blocks
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
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

build_recovery_mass_lumped(model::SolidMechanics) = build_recovery_mass_lumped(model.mesh, model.reference)

function build_recovery_data(kind::Symbol, input_mesh::ExodusDatabase, reference::Matrix{Float64})::AbstractRecoveryData
    if kind === :none
        return NoRecovery()
    elseif kind === :lumped
        m = build_recovery_mass_lumped(input_mesh, reference)
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
    else
        norma_abort("Unknown stress recovery kind: '$kind' (expected :none or :lumped)")
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

function _assemble_l2_rhs_stress!(nodal::Matrix{Float64}, model::SolidMechanics)
    input_mesh = model.mesh
    blocks = Exodus.read_sets(input_mesh, Block)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(input_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
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
