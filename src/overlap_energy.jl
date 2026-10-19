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

# Distance from `point` to side set `side_set_id` of `model` on the reference
# (undeformed) configuration.
function distance_to_side_set_reference(point::Vector{Float64}, model::SolidMechanics, side_set_id::Integer)
    facets = get_facet_grid(model, side_set_id)
    p = SVector{3,Float64}(point)
    return box_grid_nearest(
        facets.grid, p, facets.box_min, facets.box_max, facet -> clamped_facet_distance(facets.facet_coords[facet], p)
    )
end

# True when `point` lies inside any element of `model`, evaluated on the
# reference configuration (is_inside tests the reference nodes).  This is the
# physical overlap membership test and is independent of mesh conformity.
function point_in_model_reference(point::Vector{Float64}, model::SolidMechanics, tol::Float64)
    elements = get_element_grid(model)
    for item in box_grid_items(elements.grid, SVector{3,Float64}(point))
        block = model.blocks[elements.block_index[item]]
        nodes = model.reference[:, view(block.connectivity, :, elements.element_index[item])]
        _, found = is_inside(block.element_type, nodes, point, tol)
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
    model isa SolidMechanics || norma_abort("Blended overlap energy is only implemented for SolidMechanics subdomains.")
    partners = OverlapPartner[]
    for bc in model.boundary_conditions
        # Every overlapping coupling subtypes this abstract, so
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
    w = arlequin_weight_or_nan(x, model, partners, h)
    isnan(w) && abort_degenerate_overlap(x, model, partners, h, label)
    return w
end

# Constant weights of the overlap matching coupling: a point inside n partners
# has the weight 1/(1 + n), so that the weights of the subdomains sharing a
# point sum to one; outside every partner the weight is one. These are the
# weights of the assembly of such a subdomain (ASSEMBLY_WEIGHT_CACHE), and
# the blended energy and the rate invariant use the same.
function constant_overlap_weight(x::Vector{Float64}, partners::Vector{OverlapPartner})
    inside = count(p -> point_in_model_reference(x, p.partner_model, p.tol), partners)
    return 1.0 / (1.0 + inside)
end

# True when the subdomain couples by the overlap matching coupling, whose
# weights are the constant ones.
function uses_constant_overlap_weights(subsim::SingleDomainSimulation)
    return any(bc isa SolidMechanicsOverlapMatchingSchwarzBoundaryCondition for bc in subsim.model.boundary_conditions)
end

# Store the per-quadrature-point weights of a subdomain's assembly: the
# constant overlap weights of the overlap matching coupling. Every assembled
# quantity is then weighted, so the matrices are rebuilt.
function set_assembly_weights!(subsim::SingleDomainSimulation)
    model = subsim.model
    ASSEMBLY_WEIGHT_CACHE[model] = get_arlequin_weights(subsim).weights
    model.compute_mass = model.compute_lumped_mass = model.compute_stiffness = true
    return nothing
end

# The weight, or NaN where the overlap is degenerate; the threaded setup uses
# this form and reports the degeneracy serially afterwards.
function arlequin_weight_or_nan(x::Vector{Float64}, model::SolidMechanics, partners::Vector{OverlapPartner}, h::Float64)
    isempty(partners) && return 1.0
    # Sum distances to the Schwarz boundaries of the partners whose interiors
    # contain x (the physical overlap this point sits in). Membership is
    # tested first, so points of the exclusive region need no distances.
    d_partner = 0.0
    in_overlap = false
    for p in partners
        if point_in_model_reference(x, p.partner_model, p.tol)
            in_overlap = true
            d_partner += distance_to_side_set_reference(x, p.partner_model, p.partner_side_set_id)
        end
    end
    in_overlap || return 1.0  # exclusive region of this subdomain
    # Distance to this subdomain's own Schwarz boundary (nearest of its overlap
    # side sets).
    d_own = minimum(distance_to_side_set_reference(x, model, p.own_side_set_id) for p in partners)
    denom = d_own + d_partner
    denom < OVERLAP_DEGENERACY_TOL * h && return NaN
    return d_own / denom
end

function abort_degenerate_overlap(
    x::Vector{Float64}, model::SolidMechanics, partners::Vector{OverlapPartner}, h::Float64, label::String
)
    d_own = minimum(distance_to_side_set_reference(x, model, p.own_side_set_id) for p in partners)
    d_partner = sum(if point_in_model_reference(x, p.partner_model, p.tol)
        distance_to_side_set_reference(x, p.partner_model, p.partner_side_set_id)
    else
        0.0
    end for p in partners)
    denom = d_own + d_partner
    return norma_abort(
        "Degenerate overlap in subdomain \"$label\" at reference point " *
        "$(round.(x; digits=6)): the two Schwarz boundaries are within $denom (< " *
        "$(OVERLAP_DEGENERACY_TOL * h)) of each other, so the overlap has near-zero width " *
        "here. This indicates a bad domain decomposition: widen the overlap or correct " *
        "the Schwarz side sets.",
    )
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
    model isa SolidMechanics || norma_abort("Blended overlap energy is only implemented for SolidMechanics subdomains.")
    model.mesh_smoothing && norma_abort("Blended overlap energy does not support mesh-smoothing (EMS) subdomains.")
    partners = overlap_partners(subsim)
    # Build the search grids serially before the threaded classification: the
    # builders read the Exodus files and fill the caches, neither of which is
    # safe from several threads at once.
    for p in partners
        get_element_grid(p.partner_model)
        get_facet_grid(p.partner_model, p.partner_side_set_id)
        get_facet_grid(model, p.own_side_set_id)
    end
    constant = uses_constant_overlap_weights(subsim)
    weights = Vector{Vector{Vector{Float64}}}(undef, length(model.blocks))
    for (block_index, block) in enumerate(model.blocks)
        weights[block_index] = compute_block_arlequin_weights(model, block, partners, constant)
    end
    # Report a degenerate overlap serially, with the full diagnostic.
    for (block_index, block) in enumerate(model.blocks), e in 1:block.num_elements
        element_weights = weights[block_index][e]
        any(isnan, element_weights) || continue
        point = findfirst(isnan, element_weights)
        element_reference = gather_nodal(model.reference, view(block.connectivity, :, e), block.N)
        x = Vector{Float64}(element_reference * block.N[:, point])
        abort_degenerate_overlap(x, model, partners, element_characteristic_length(element_reference), subsim.name)
    end
    return ArlequinWeights(weights)
end

# Function barrier on the block's static shape functions. The grids are read
# only, so the elements are classified in parallel.
function compute_block_arlequin_weights(
    model::SolidMechanics, block::ElementBlockData, partners::Vector{OverlapPartner}, constant::Bool=false
)
    N = block.N
    conn = block.connectivity
    num_points = block.num_points
    block_weights = Vector{Vector{Float64}}(undef, block.num_elements)
    @threads for e in 1:block.num_elements
        element_reference = gather_nodal(model.reference, view(conn, :, e), N)
        h = element_characteristic_length(element_reference)
        element_weights = Vector{Float64}(undef, num_points)
        for point in 1:num_points
            x = Vector{Float64}(element_reference * N[:, point])
            element_weights[point] =
                constant ? constant_overlap_weight(x, partners) : arlequin_weight_or_nan(x, model, partners, h)
        end
        block_weights[e] = element_weights
    end
    return block_weights
end

# Cached per-quadrature-point Arlequin weights for a subdomain.
function get_arlequin_weights(subsim::SingleDomainSimulation)
    return get!(() -> compute_arlequin_weights(subsim), ARLEQUIN_WEIGHT_CACHE, subsim.model)
end

# Strain-energy density W(F) at a single point, recomputed exactly as evaluate()
# does: the same constitutive update, with the lattice orientation of the
# element and the old internal state for history-dependent (inelastic)
# materials.  Only W is needed here, so the tangent is skipped and the returned
# stress/new state are discarded (this is a read-only re-integration, it must
# not advance state).  The energy is written at stops, after the steps have
# been committed, so the time increment is zero: a rate-dependent material
# then returns the energy of its committed state.
function strain_energy_density(
    material::Material, F::SMatrix{3,3,Float64,9}, Q::SMatrix{3,3,Float64,9}, state_old::Vector{Float64}
)
    W, _, _, _ = constitutive_update(material, F, Q, state_old, 0.0; need_tangent=false)
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
# mass) to machine precision. For central difference the kinetic energy is the
# product ½ v_{n-½}ᵀ M_L v_{n+½} of the half-step velocities, ½ vᵀ M_L v -
# (Δt²/8) aᵀ M_L a, which the scheme conserves exactly on a linear problem. It
# can be negative: for a body released at rest it is -(Δt²/8) aᵀ M_L a at t = 0.
function blended_subdomain_energy(subsim::SingleDomainSimulation)
    model = subsim.model
    aw = get_arlequin_weights(subsim)
    dynamic = is_dynamic(subsim.integrator)
    integrator = subsim.integrator
    half_dt = 0.0
    if integrator isa CentralDifference
        nominal_dt = if integrator.minimum_time_step == integrator.maximum_time_step
            integrator.maximum_time_step
        else
            integrator.time_step
        end
        half_dt = 0.5 * nominal_dt
    end
    stored_tl = zeros(maxthreadid())
    kinetic_tl = zeros(maxthreadid())
    for (block_index, block) in enumerate(model.blocks)
        material = model.materials[block_index]
        blended_block_energy!(
            stored_tl,
            kinetic_tl,
            model,
            block,
            block_index,
            material,
            material.ρ,
            aw.weights[block_index],
            dynamic,
            half_dt,
        )
    end
    return sum(stored_tl), sum(kinetic_tl)
end

# Function barrier on the block's static shape functions; the elements are
# integrated in parallel into per-thread accumulators without allocation.
function blended_block_energy!(
    stored_tl::Vector{Float64},
    kinetic_tl::Vector{Float64},
    model::SolidMechanics,
    block::ElementBlockData,
    block_index::Int64,
    material::Material,
    density::Float64,
    block_weights::Vector{Vector{Float64}},
    dynamic::Bool,
    half_dt::Float64,
)
    N = block.N
    dN = block.dN
    ip_weights = block.weights
    conn = block.connectivity
    num_points = block.num_points
    num_element_nodes = block.num_nodes_per_element
    @threads for e in 1:block.num_elements
        node_indices = view(conn, :, e)
        element_reference = gather_nodal(model.reference, node_indices, N)
        element_current = element_reference + gather_nodal(model.displacement, node_indices, N)
        element_velocity = gather_nodal(model.velocity, node_indices, N)
        element_acceleration = gather_nodal(model.acceleration, node_indices, N)
        element_weights = block_weights[e]
        stored = 0.0
        kinetic = 0.0
        for point in 1:num_points
            w = element_weights[point]
            dNdξ = dN[:, :, point]
            dXdξ = dNdξ * element_reference'
            dvol = det(dXdξ) * ip_weights[point]
            # Deformation gradient and strain-energy density, exactly as
            # evaluate() computes them.
            dNdX = dXdξ \ dNdξ
            F = element_current * dNdX'
            state_old = model.state_old[block_index][e][point]
            Q = element_lattice_orientation(model, material, block_index, e)
            W = strain_energy_density(material, F, Q, state_old)
            stored += w * W * dvol
            if dynamic
                Np = N[:, point]
                if half_dt > 0.0
                    s = sum(Np)
                    for i in 1:num_element_nodes
                        m_i = w * density * dvol * Np[i] * s
                        v2 = element_velocity[1, i]^2 + element_velocity[2, i]^2 + element_velocity[3, i]^2
                        a2 = element_acceleration[1, i]^2 + element_acceleration[2, i]^2 + element_acceleration[3, i]^2
                        kinetic += 0.5 * m_i * (v2 - half_dt * half_dt * a2)
                    end
                else
                    v = element_velocity * Np
                    kinetic += w * 0.5 * density * dot(v, v) * dvol
                end
            end
        end
        t = threadid()
        stored_tl[t] += stored
        kinetic_tl[t] += kinetic
    end
    return nothing
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

# ---------------------------------------------------------------------------
# Rate invariant: the invariant of the time-differentiated system, the energy
# functional evaluated on the velocity and the acceleration (energy per time squared).
#
# Differentiating the linear equation of motion M ü + K u = f in time gives
# M a' + K v = f', whose Newmark discrete energy (Prakash and Hjelmstad 2004,
# Eqs. (56)-(58), (71)) is
#
#     Ẽ = ½ aᵀ A a + ½ vᵀ K v,   A = M + (Δt²/2)(2β - γ) K.
#
# Continuity of the interface velocity (Gravouil and Combescure 2001) conserves
# the sum of Ẽ over the subdomains at equal time steps, while continuity of
# the displacement with trapezoidal Newmark (β = 1/4, γ = 1/2) conserves the
# physical energy ½ vᵀ M v + SE(u). With SE(x) the strain energy evaluated with
# the nodal field x in place of the displacement, ½ xᵀ K x = SE(x) for a linear
# elastic material, whose strain energy is the quadratic form of the small
# strain of the displacement gradient (strain_energy for Linear_Elastic in
# constitutive.jl), so
#
#     Ẽ = ½ aᵀ M a + (Δt²/2)(2β - γ) SE(a) + SE(v),
#
# which is ½ aᵀ M a + SE(v) for trapezoidal Newmark with the consistent mass and
# ½ aᵀ M_L a - (Δt²/4) SE(a) + SE(v) for central difference (β = 0, γ = 1/2)
# with the lumped mass M_L. Ẽ is evaluated only for subdomains whose materials
# are all linear elastic and whose integrator is Newmark without HHT-α or
# central difference; it is NaN otherwise.
# ---------------------------------------------------------------------------

# True when Ẽ is defined for the subdomain.
function rate_invariant_applies(subsim::SingleDomainSimulation)
    model = subsim.model
    model isa SolidMechanics || return false
    model.mesh_smoothing && return false
    all(m -> m isa Linear_Elastic, model.materials) || return false
    integrator = subsim.integrator
    integrator isa CentralDifference && return true
    return integrator isa Newmark && integrator.hht_alpha == 0.0
end

# SE(v), SE(a), and the consistent-mass term ½ aᵀ M a of a linear elastic
# subdomain, integrated with the element quadrature of evaluate().
function quadratic_energies(model::SolidMechanics, weights::Union{Nothing,Vector{Vector{Vector{Float64}}}}=nothing)
    se_v_tl = zeros(maxthreadid())
    se_a_tl = zeros(maxthreadid())
    kin_a_tl = zeros(maxthreadid())
    for (block_index, block) in enumerate(model.blocks)
        material = model.materials[block_index]
        block_weights = weights === nothing ? nothing : weights[block_index]
        quadratic_block_energies!(se_v_tl, se_a_tl, kin_a_tl, model, block, material, block_weights)
    end
    return sum(se_v_tl), sum(se_a_tl), sum(kin_a_tl)
end

# Function barrier on the block's static shape functions.
function quadratic_block_energies!(
    se_v_tl::Vector{Float64},
    se_a_tl::Vector{Float64},
    kin_a_tl::Vector{Float64},
    model::SolidMechanics,
    block::ElementBlockData,
    material::Linear_Elastic,
    block_weights::Union{Nothing,Vector{Vector{Float64}}}=nothing,
)
    N = block.N
    dN = block.dN
    ip_weights = block.weights
    conn = block.connectivity
    num_points = block.num_points
    density = material.ρ
    @threads for e in 1:block.num_elements
        node_indices = view(conn, :, e)
        element_reference = gather_nodal(model.reference, node_indices, N)
        element_velocity = gather_nodal(model.velocity, node_indices, N)
        element_acceleration = gather_nodal(model.acceleration, node_indices, N)
        se_v = 0.0
        se_a = 0.0
        kin_a = 0.0
        for point in 1:num_points
            dNdξ = dN[:, :, point]
            dXdξ = dNdξ * element_reference'
            dvol = det(dXdξ) * ip_weights[point]
            block_weights === nothing || (dvol *= block_weights[e][point])
            dNdX = dXdξ \ dNdξ
            F_v = SMatrix{3,3,Float64,9}(I) + element_velocity * dNdX'
            F_a = SMatrix{3,3,Float64,9}(I) + element_acceleration * dNdX'
            se_v += strain_energy(material, F_v) * dvol
            se_a += strain_energy(material, F_a) * dvol
            a = element_acceleration * N[:, point]
            kin_a += 0.5 * density * dot(a, a) * dvol
        end
        t = threadid()
        se_v_tl[t] += se_v
        se_a_tl[t] += se_a
        kin_a_tl[t] += kin_a
    end
    return nothing
end

# Mask of the degrees of freedom of the nodes in any Schwarz side set.
function schwarz_interface_dofs(model::SolidMechanics)
    mask = falses(length(model.free_dofs))
    for bc in model.boundary_conditions
        bc isa SolidMechanicsSchwarzBoundaryCondition || continue
        hasproperty(bc, :side_set_node_indices) || continue
        for node in bc.side_set_node_indices
            mask[(3 * node - 2):(3 * node)] .= true
        end
    end
    return mask
end

# Ẽ of one subdomain (NaN where it does not apply) and, for central
# difference, the staggered kinetic energy ½ vᵀ M_L v - (Δt²/8) aᵀ M_L a split
# into the rows of the Schwarz interface nodes and the remaining rows (NaN for
# other integrators). The lumped mass is diagonal, so the split is exact.
function subdomain_rate_invariant(subsim::SingleDomainSimulation)
    rate_invariant_applies(subsim) || return NaN, NaN, NaN
    model = subsim.model
    integrator = subsim.integrator
    Δt = integrator.time_step
    β, γ = integrator isa CentralDifference ? (0.0, integrator.γ) : (integrator.β, integrator.γ)
    # The rate invariant of an overlap matching subdomain is blended with the
    # constant weights of its assembly, so that the overlap is counted once as
    # in the coupled problem; its lumped mass carries the same weights.
    weights = uses_constant_overlap_weights(subsim) ? get_arlequin_weights(subsim).weights : nothing
    se_v, se_a, kin_a_consistent = quadratic_energies(model, weights)
    staggered_interface = NaN
    staggered_interior = NaN
    if integrator isa CentralDifference
        m = model.lumped_mass
        a = vec(model.acceleration)
        v = vec(model.velocity)
        length(m) == length(a) || return NaN, NaN, NaN
        kin_a = 0.5 * dot(m, a .* a)
        rows = 0.5 .* m .* (v .* v .- 0.25 * Δt * Δt .* a .* a)
        interface = schwarz_interface_dofs(model)
        staggered_interface = sum(rows[interface])
        staggered_interior = sum(rows[.!interface])
    else
        kin_a = kin_a_consistent
    end
    invariant = kin_a + 0.5 * Δt * Δt * (2.0 * β - γ) * se_a + se_v
    return invariant, staggered_interface, staggered_interior
end

# Interface displacement and velocity jumps of every DN pair at the end of the
# stop (see dn_interface_residuals in schwarz.jl), relative and as root mean
# square values over the interface, with the column names.
function dn_jump_columns(sim::MultiDomainSimulation)
    names = String[]
    values = Float64[]
    for (bc, r) in dn_interface_residuals(sim)
        pair = r.dirichlet_name * "_" * r.neumann_name
        push!(names, "displacement_jump_" * pair, "velocity_jump_" * pair)
        push!(values, r.displacement_jump, r.velocity_jump)
        push!(names, "displacement_jump_rms_" * pair, "velocity_jump_rms_" * pair)
        push!(values, r.displacement_jump_rms, r.velocity_jump_rms)
        # Interface impulse residual of a constrained nonoverlap pair over the stop.
        if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition && is_constrained_dn(bc)
            net, relative = is_direct_dn(bc) ? (fill(NaN, 3), NaN) : dn_impulse_residual(sim, bc)
            for (comp, label) in enumerate(("x", "y", "z"))
                push!(names, "impulse_residual_" * label * "_" * pair)
                push!(values, net[comp])
            end
            push!(names, "impulse_residual_relative_" * pair)
            push!(values, relative)
        end
    end
    return names, values
end

# Append the current blended energy to <sim.name>-energy.csv (header written on
# the first stop). The first four columns are the blended physical energy; the
# columns after them are the rate invariant Ẽ (total and per
# subdomain), the staggered kinetic energy of the central difference subdomains
# on the Schwarz interface rows and on the remaining rows, and the one-sided
# displacement and velocity jumps of each Dirichlet-Neumann pair.
function write_blended_energy_csv(sim::MultiDomainSimulation)
    stored, kinetic, total = total_blended_energy(sim)
    invariant = Float64[]
    staggered_interface = 0.0
    staggered_interior = 0.0
    any_explicit = false
    for subsim in sim.subsims
        invariant_k, interface_k, interior_k = subdomain_rate_invariant(subsim)
        push!(invariant, invariant_k)
        if subsim.integrator isa CentralDifference
            any_explicit = true
            staggered_interface += interface_k
            staggered_interior += interior_k
        end
    end
    if !any_explicit
        staggered_interface = staggered_interior = NaN
    end
    jump_names, jump_values = dn_jump_columns(sim)
    filename = sim.name * "-energy.csv"
    stop = sim.controller.stop
    time = sim.controller.time
    open(filename, stop == 0 ? "w" : "a") do io
        if stop == 0
            header = ["time", "stored_energy", "kinetic_energy", "total_energy", "rate_invariant_total"]
            append!(header, ["rate_invariant_" * subsim.name for subsim in sim.subsims])
            append!(header, ["staggered_kinetic_interface", "staggered_kinetic_interior"])
            append!(header, jump_names)
            println(io, join(header, ","))
        end
        row = Any[time, stored, kinetic, total, sum(invariant)]
        append!(row, invariant)
        append!(row, [staggered_interface, staggered_interior])
        append!(row, jump_values)
        return println(io, join(row, ","))
    end
    return nothing
end
