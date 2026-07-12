# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# ---------------------------------------------------------------------------
# Arlequin (distance-ratio) partition-of-unity energy for overlapping Schwarz.
#
# A multidomain decomposition with overlapping subdomains meshes the physical
# overlap region twice, so naively summing per-subdomain energies double-counts
# it.  This module builds a smooth partition of unity across the overlap and
# reports blended kinetic / stored / total energy that is free of the double
# count and reduces to the monodomain energy at Schwarz convergence.
#
# For a point x in the overlap of subdomain k with partner j, let d_own(x) be
# the distance to subdomain k's own Schwarz boundary Gamma_k and d_partner(x)
# the distance to the partner's Schwarz boundary Gamma_j.  The weight
#
#     w_k(x) = d_own(x) / (d_own(x) + d_partner(x))
#
# is 0 at Gamma_k (this subdomain is driven there, so hand the energy to the
# partner), 1 at Gamma_j (the partner is driven there), and w_k + w_j = 1
# everywhere.  Outside the overlap w_k = 1 (this subdomain's exclusive region).
# Distances are the clamped closest-point distance to a side set, measured on
# the reference configuration, so the weights are a fixed function of the
# undeformed geometry and are computed once and cached.
#
# A vanishing d_own + d_partner means the two Schwarz boundaries touch (the
# overlap pinches to zero width), which indicates a bad decomposition; such a
# point is rejected outright rather than regularized.
# ---------------------------------------------------------------------------

# Reject the weight when the overlap band pinches below this fraction of the
# local element size: the two Schwarz boundaries effectively coincide there.
const OVERLAP_DEGENERACY_TOL = 1.0e-3

# Per-quadrature-point blending weights for one subdomain, laid out like
# SolidMechanics.stored_energy but resolved per integration point:
# weights[block_index][element_index][point].
struct ArlequinWeights
    weights::Vector{Vector{Vector{Float64}}}
end

# Reference-configuration weights never change during a run, so memoize them by
# model identity.  Prototype-scoped: this avoids threading a new field through
# SolidMechanics and its constructors.
const ARLEQUIN_WEIGHT_CACHE = IdDict{SolidMechanics,ArlequinWeights}()

# One overlapping-Schwarz coupling of a subdomain: this subdomain's own Schwarz
# side set, the partner model, and the partner's Schwarz side set.
struct OverlapPartner
    own_side_set_id::Int64
    partner_model::SolidMechanics
    partner_side_set_id::Int64
    tol::Float64
end

# Clamped closest-point distance from `point` to side set `side_set_id` of
# `model`, on the reference (undeformed) configuration.  Mirrors the facet scan
# of project_point_to_containing_facet but on model.reference and returning only
# the distance; falls back to the nearest-node distance if every facet's Newton
# projection diverges, so the result is always finite.
function distance_to_side_set_reference(point::Vector{Float64}, model::SolidMechanics, side_set_id::Integer)
    mesh = model.mesh
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    coords = model.reference
    best_facet = Inf
    best_node = Inf
    ss_node_index = 1
    for num_nodes_side in num_nodes_sides
        face_node_indices = side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)]
        ss_node_index += num_nodes_side
        face_nodes = coords[:, face_node_indices]
        best_node = min(best_node, get_minimum_distance_to_nodes(face_nodes, point))
        # strict=false reports a diverging projection as an infinite distance
        # instead of aborting, so a badly shaped facet is simply skipped.
        _, _, distance, _ = closest_point_projection(face_nodes, point; strict=false)
        isfinite(distance) && (best_facet = min(best_facet, abs(distance)))
    end
    return isfinite(best_facet) ? best_facet : best_node
end

# True when `point` lies inside any element of `model`, evaluated on the
# reference configuration (find_point_in_mesh tests model.reference).  This is
# the physical overlap membership test and is independent of mesh conformity.
function point_in_model_reference(point::Vector{Float64}, model::SolidMechanics, tol::Float64)
    for block in Exodus.read_sets(model.mesh, Block)
        _, _, found = find_point_in_mesh(point, model, Int64(block.id), tol)
        found && return true
    end
    return false
end

# The overlapping-Schwarz couplings of a subdomain.  Empty when the subdomain
# has no overlap coupling (single domain, or non-overlap / contact Schwarz, all
# of which partition the domain with a measure-zero interface and need no
# correction), in which case every weight is 1.
function overlap_partners(subsim::SingleDomainSimulation)
    model = subsim.model
    model isa SolidMechanics ||
        norma_abort("Blended overlap energy is only implemented for SolidMechanics subdomains.")
    partners = OverlapPartner[]
    for bc in model.boundary_conditions
        # Every overlapping coupling (DBC or impedance) subtypes this abstract, so
        # new overlap variants are picked up automatically; non-overlap and
        # contact Schwarz partition across a measure-zero interface and keep unit
        # weights.
        bc isa SolidMechanicsOverlapCouplingSchwarzBoundaryCondition || continue
        partner_subsim = coupled_subsim_of(bc)
        partner_model = partner_subsim.model
        partner_model isa SolidMechanics || norma_abort(
            "Blended overlap energy requires a full-order (SolidMechanics) partner; subdomain " *
            "\"$(partner_subsim.name)\" is not supported.",
        )
        partner_side_set_id = partner_schwarz_side_set(bc, partner_model)
        push!(partners, OverlapPartner(bc.side_set_id, partner_model, partner_side_set_id, bc.search_tolerance))
    end
    return partners
end

# The partner's Schwarz boundary side set is carried by the partner's overlap BC
# that couples back to this subdomain (its coupled_handle points at us).
function partner_schwarz_side_set(
    bc::SolidMechanicsOverlapCouplingSchwarzBoundaryCondition, partner_model::SolidMechanics
)
    for partner_bc in partner_model.boundary_conditions
        if partner_bc isa SolidMechanicsOverlapCouplingSchwarzBoundaryCondition &&
            partner_bc.coupled_handle.id == bc.self_handle.id
            return partner_bc.side_set_id
        end
    end
    return norma_abort("Could not locate the partner Schwarz boundary for overlap BC \"$(bc.name)\".")
end

# Arlequin partition-of-unity weight at reference point `x` for the subdomain
# whose partners are `partners`.  `h` is the local element size, used only to
# scale the degeneracy tolerance; `label` names the subdomain for diagnostics.
function arlequin_weight(
    x::Vector{Float64}, model::SolidMechanics, partners::Vector{OverlapPartner}, h::Float64, label::String
)
    isempty(partners) && return 1.0
    # Distance to this subdomain's own Schwarz boundary (nearest of its overlap
    # side sets).
    d_own = minimum(distance_to_side_set_reference(x, model, p.own_side_set_id) for p in partners)
    # Sum distances to the Schwarz boundaries of the partners whose interiors
    # contain x (the physical overlap this point sits in).
    d_partner = 0.0
    in_overlap = false
    for p in partners
        if point_in_model_reference(x, p.partner_model, p.tol)
            in_overlap = true
            d_partner += distance_to_side_set_reference(x, p.partner_model, p.partner_side_set_id)
        end
    end
    in_overlap || return 1.0  # exclusive region of this subdomain
    denom = d_own + d_partner
    if denom < OVERLAP_DEGENERACY_TOL * h
        norma_abort(
            "Degenerate overlap in subdomain \"$label\" at reference point " *
            "$(round.(x; digits=6)): the two Schwarz boundaries are within $denom (< " *
            "$(OVERLAP_DEGENERACY_TOL * h)) of each other, so the overlap has near-zero width " *
            "here. This indicates a bad domain decomposition — widen the overlap or correct " *
            "the Schwarz side sets.",
        )
    end
    return d_own / denom
end

# Bounding-box diagonal of an element's reference nodes, a cheap characteristic
# length for the degeneracy tolerance.
function element_characteristic_length(element_reference::AbstractMatrix{Float64})
    return norm([maximum(element_reference[i, :]) - minimum(element_reference[i, :]) for i in 1:3])
end

# Compute (and do not cache) the per-quadrature-point Arlequin weights for a
# subdomain.  Aborts on a degenerate overlap.
function compute_arlequin_weights(subsim::SingleDomainSimulation)
    model = subsim.model
    model isa SolidMechanics ||
        norma_abort("Blended overlap energy is only implemented for SolidMechanics subdomains.")
    model.mesh_smoothing &&
        norma_abort("Blended overlap energy does not support mesh-smoothing (EMS) subdomains.")
    partners = overlap_partners(subsim)
    mesh = model.mesh
    blocks = Exodus.read_sets(mesh, Block)
    weights = Vector{Vector{Vector{Float64}}}()
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, _, _ = isoparametric(element_type, num_points)
        conn = get_block_connectivity(mesh, block_id)
        num_block_elements, num_element_nodes = size(conn)
        block_weights = Vector{Vector{Float64}}()
        for element_index in 1:num_block_elements
            # evaluate() recovers element node lists by LINEAR indexing into the
            # reshaped connectivity (element-major), not by row slicing.
            connectivity_indices =
                ((element_index - 1) * num_element_nodes + 1):(element_index * num_element_nodes)
            node_indices = conn[connectivity_indices]
            element_reference = model.reference[:, node_indices]
            h = element_characteristic_length(element_reference)
            element_weights = Vector{Float64}(undef, num_points)
            for point in 1:num_points
                x = Vector{Float64}(element_reference * N[:, point])
                element_weights[point] = arlequin_weight(x, model, partners, h, subsim.name)
            end
            push!(block_weights, element_weights)
        end
        push!(weights, block_weights)
    end
    return ArlequinWeights(weights)
end

# Cached per-quadrature-point Arlequin weights for a subdomain.
function get_arlequin_weights(subsim::SingleDomainSimulation)
    return get!(() -> compute_arlequin_weights(subsim), ARLEQUIN_WEIGHT_CACHE, subsim.model)
end

# Strain-energy density W(F) at a single point, recomputed exactly as evaluate()
# does: the same constitutive call, branching on the material family and passing
# the old internal state for history-dependent (inelastic) materials.  Only W is
# needed here, so the tangent is skipped and the returned stress/new state are
# discarded (this is a read-only re-integration, it must not advance state).
function strain_energy_density(material::Material, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64})
    if material isa Elastic
        W, _, _ = constitutive(material, F; need_tangent=false)
    else
        W, _, _, _ = constitutive(material, F, state_old; need_tangent=false)
    end
    return W
end

# Blended stored (strain) and kinetic energy of one subdomain, using the cached
# Arlequin weights.  Both energies are re-integrated per quadrature point with
# the weight applied inside the integrand — exact regardless of mesh conformity:
#
#     stored = Σ_qp w(x_qp) W(F_qp) dV_qp,   kinetic = Σ_qp w(x_qp) ½ρ|v_qp|² dV_qp.
#
# The strain-energy density W is recomputed from the deformation gradient with
# the same kinematics and constitutive model as evaluate(), so with unit weights
# stored reduces to model.strain_energy and kinetic to 0.5 vᵀ M v (consistent
# mass) to machine precision.
function blended_subdomain_energy(subsim::SingleDomainSimulation)
    model = subsim.model
    aw = get_arlequin_weights(subsim)
    mesh = model.mesh
    blocks = Exodus.read_sets(mesh, Block)
    dynamic = is_dynamic(subsim.integrator)
    # Central difference advances the LUMPED-mass system with leapfrog
    # (staggered-velocity) structure, so the discrete energy it conserves
    # (linear problems) is the staggered product with the lumped mass,
    #   ½ v_{n-½}ᵀ M_L v_{n+½} + SE(u_n)
    #     = ½ v_nᵀ M_L v_n − (Δt²/8) a_nᵀ M_L a_n + SE(u_n).
    # Evaluating the consistent-mass integer-step kinetic energy ½ρ|v_n|²
    # instead misprices the short-wavelength content of sharp wavefronts twice
    # over (consistent vs lumped Rayleigh quotients differ at high wavenumber,
    # and the integer-step velocity aliases the staggered content) and reports
    # an O(10%) spurious dip/oscillation on the cantilever release benchmark.
    # For central-difference subdomains the kinetic term therefore accumulates
    # the row-sum-lumped nodal masses of the SAME blended quadrature (exactly
    # add_lumped_mass!'s lumping, weight applied inside) against the nodal
    # staggered product |v_i|² − (Δt/2)²|a_i|². Implicit (consistent-mass)
    # integrators keep the consistent-mass integer-step evaluation.
    integrator = subsim.integrator
    half_dt = 0.0
    if integrator isa CentralDifference
        nominal_dt = integrator.minimum_time_step == integrator.maximum_time_step ?
            integrator.maximum_time_step : integrator.time_step
        half_dt = 0.5 * nominal_dt
    end
    stored = 0.0
    kinetic = 0.0
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        material = model.materials[block_index]
        density = material.ρ
        element_type_string = Exodus.read_block_parameters(mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        N, dN, ip_weights = isoparametric(element_type, num_points)
        conn = get_block_connectivity(mesh, block_id)
        num_block_elements, num_element_nodes = size(conn)
        for element_index in 1:num_block_elements
            connectivity_indices =
                ((element_index - 1) * num_element_nodes + 1):(element_index * num_element_nodes)
            node_indices = conn[connectivity_indices]
            element_reference = model.reference[:, node_indices]
            element_current = element_reference + model.displacement[:, node_indices]
            element_velocity = model.velocity[:, node_indices]
            element_acceleration = half_dt > 0.0 ? model.acceleration[:, node_indices] : nothing
            element_weights = aw.weights[block_index][element_index]
            for point in 1:num_points
                w = element_weights[point]
                dNdξ = dN[:, :, point]
                dXdξ = SMatrix{3,3,Float64,9}(dNdξ * element_reference')
                dvol = det(dXdξ) * ip_weights[point]
                # Deformation gradient and strain-energy density, exactly as
                # evaluate() computes them.
                dNdX = dXdξ \ dNdξ
                F = SMatrix{3,3,Float64,9}(element_current * dNdX')
                state_old = model.state_old[block_index][element_index][point]
                W = strain_energy_density(material, F, state_old)
                stored += w * W * dvol
                if dynamic
                    if half_dt > 0.0
                        s = sum(N[:, point])
                        for i in 1:num_element_nodes
                            m_i = w * density * dvol * N[i, point] * s
                            v_i = @views element_velocity[:, i]
                            a_i = @views element_acceleration[:, i]
                            kinetic +=
                                0.5 * m_i * (dot(v_i, v_i) - half_dt * half_dt * dot(a_i, a_i))
                        end
                    else
                        v = element_velocity * N[:, point]
                        kinetic += w * 0.5 * density * dot(v, v) * dvol
                    end
                end
            end
        end
    end
    return stored, kinetic
end

# Blended (stored, kinetic, total) energy of a multidomain simulation, summed
# over subdomains with the overlap double count removed by the Arlequin weights.
function total_blended_energy(sim::MultiDomainSimulation)
    stored = 0.0
    kinetic = 0.0
    for subsim in sim.subsims
        subsim_stored, subsim_kinetic = blended_subdomain_energy(subsim)
        stored += subsim_stored
        kinetic += subsim_kinetic
    end
    return stored, kinetic, stored + kinetic
end

# Append the current blended energy to <sim.name>-energy.csv (header written on
# the first stop).
function write_blended_energy_csv(sim::MultiDomainSimulation)
    stored, kinetic, total = total_blended_energy(sim)
    filename = sim.name * "-energy.csv"
    stop = sim.controller.stop
    time = sim.controller.time
    open(filename, stop == 0 ? "w" : "a") do io
        stop == 0 && println(io, "time,stored_energy,kinetic_energy,total_energy")
        println(io, join((time, stored, kinetic, total), ","))
    end
    return nothing
end
