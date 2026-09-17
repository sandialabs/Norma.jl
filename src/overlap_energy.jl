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

# ---------------------------------------------------------------------------
# Spatial search on the reference configuration.
#
# The weights need, for every quadrature point of a subdomain, whether the
# point lies inside the partner mesh and its distances to the two Schwarz
# side sets. Scanning every partner element and every facet per point made
# the setup quadratic (hours on meshes of a few tens of thousands of
# elements), so elements and facets are binned once into uniform grids of
# cells by their reference bounding boxes: a point-in-mesh query tests only
# the elements binned in the point's cell, and a distance query visits the
# cells in rings around the point until no unvisited cell can hold a closer
# facet.
# ---------------------------------------------------------------------------

struct BoxGrid
    lo::SVector{3,Float64}
    cell::Float64
    dims::SVector{3,Int64}
    items::Vector{Vector{Int64}}  # item indices binned per cell, column-major over dims
end

const BOX_GRID_MAX_CELLS_PER_AXIS = 256
const EMPTY_ITEMS = Int64[]

@inline function box_grid_cell(grid::BoxGrid, p::SVector{3,Float64})
    return SVector{3,Int64}(floor(Int64, (p[i] - grid.lo[i]) / grid.cell) + 1 for i in 1:3)
end

@inline function box_grid_linear(grid::BoxGrid, c::SVector{3,Int64})
    return c[1] + grid.dims[1] * ((c[2] - 1) + grid.dims[2] * (c[3] - 1))
end

@inline function box_grid_contains(grid::BoxGrid, c::SVector{3,Int64})
    return 1 <= c[1] <= grid.dims[1] && 1 <= c[2] <= grid.dims[2] && 1 <= c[3] <= grid.dims[3]
end

function BoxGrid(box_min::Vector{SVector{3,Float64}}, box_max::Vector{SVector{3,Float64}})
    num_items = length(box_min)
    num_items == 0 && return BoxGrid(zero(SVector{3,Float64}), 1.0, SVector{3,Int64}(1, 1, 1), [Int64[]])
    lo = reduce((a, b) -> min.(a, b), box_min)
    hi = reduce((a, b) -> max.(a, b), box_max)
    extent = hi - lo
    mean_diagonal = sum(norm(box_max[i] - box_min[i]) for i in 1:num_items) / num_items
    # A cell about twice the mean item diagonal keeps the items per cell small
    # while an item overlaps only a few cells; the cap bounds the memory.
    cell = max(2.0 * mean_diagonal, maximum(extent) / BOX_GRID_MAX_CELLS_PER_AXIS, 1.0e-12)
    dims = SVector{3,Int64}(max(1, ceil(Int64, extent[i] / cell)) for i in 1:3)
    # Shift the origin so that boxes on the lower faces fall inside cell 1.
    lo = lo - SVector{3,Float64}(1.0e-6 * cell, 1.0e-6 * cell, 1.0e-6 * cell)
    grid = BoxGrid(lo, cell, dims, [Int64[] for _ in 1:prod(dims)])
    for item in 1:num_items
        c_lo = max.(box_grid_cell(grid, box_min[item]), 1)
        c_hi = min.(box_grid_cell(grid, box_max[item]), dims)
        for k in c_lo[3]:c_hi[3], j in c_lo[2]:c_hi[2], i in c_lo[1]:c_hi[1]
            push!(grid.items[box_grid_linear(grid, SVector{3,Int64}(i, j, k))], item)
        end
    end
    return grid
end

# Items whose boxes overlap the cell containing `p` (empty outside the grid).
function box_grid_items(grid::BoxGrid, p::SVector{3,Float64})
    c = box_grid_cell(grid, p)
    return box_grid_contains(grid, c) ? grid.items[box_grid_linear(grid, c)] : EMPTY_ITEMS
end

# Distance from `p` to the box [lo, hi], zero inside: a lower bound on the
# distance to anything the box contains.
@inline function box_distance(p::SVector{3,Float64}, lo::SVector{3,Float64}, hi::SVector{3,Float64})
    d = max.(lo - p, p - hi, 0.0)
    return norm(d)
end

# Minimum over the items of `distance(item)`, visiting the cells in rings of
# growing Chebyshev index distance around the cell of `p` and stopping once
# every unvisited cell is farther from `p` than the best distance found. An
# item is evaluated only when its box is closer than the best distance so
# far, and an item binned in several cells is simply examined more than once.
function box_grid_nearest(
    grid::BoxGrid, p::SVector{3,Float64}, box_min::Vector{SVector{3,Float64}}, box_max::Vector{SVector{3,Float64}}, distance::F
) where {F}
    c = box_grid_cell(grid, p)
    best = Inf
    max_ring = 0
    for i in 1:3
        max_ring = max(max_ring, c[i] - 1, grid.dims[i] - c[i])
    end
    for ring in 0:max_ring
        ring > 0 && best <= (ring - 1) * grid.cell && break
        for k in (c[3] - ring):(c[3] + ring), j in (c[2] - ring):(c[2] + ring), i in (c[1] - ring):(c[1] + ring)
            max(abs(i - c[1]), abs(j - c[2]), abs(k - c[3])) == ring || continue
            cc = SVector{3,Int64}(i, j, k)
            box_grid_contains(grid, cc) || continue
            for item in grid.items[box_grid_linear(grid, cc)]
                box_distance(p, box_min[item], box_max[item]) < best || continue
                best = min(best, distance(item))
            end
        end
    end
    return best
end

# Elements of a model binned by their reference bounding boxes.
struct ElementGrid
    grid::BoxGrid
    block_index::Vector{Int64}
    element_index::Vector{Int64}
end

function build_element_grid(model::SolidMechanics)
    box_min = SVector{3,Float64}[]
    box_max = SVector{3,Float64}[]
    block_index = Int64[]
    element_index = Int64[]
    for (b, block) in enumerate(model.blocks)
        conn = block.connectivity
        for e in 1:block.num_elements
            nodes = view(model.reference, :, view(conn, :, e))
            push!(box_min, SVector{3,Float64}(minimum(view(nodes, i, :)) for i in 1:3))
            push!(box_max, SVector{3,Float64}(maximum(view(nodes, i, :)) for i in 1:3))
            push!(block_index, b)
            push!(element_index, e)
        end
    end
    return ElementGrid(BoxGrid(box_min, box_max), block_index, element_index)
end

# Facets of a side set with their reference coordinates, binned by bounding box.
struct FacetGrid
    grid::BoxGrid
    facet_coords::Vector{Matrix{Float64}}
    box_min::Vector{SVector{3,Float64}}
    box_max::Vector{SVector{3,Float64}}
end

function build_facet_grid(model::SolidMechanics, side_set_id::Integer)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(model.mesh, side_set_id)
    coords = model.reference
    num_facets = length(num_nodes_sides)
    facet_coords = Vector{Matrix{Float64}}(undef, num_facets)
    box_min = Vector{SVector{3,Float64}}(undef, num_facets)
    box_max = Vector{SVector{3,Float64}}(undef, num_facets)
    ss_node_index = 1
    for (facet, num_nodes_side) in enumerate(num_nodes_sides)
        indices = Int64.(side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)])
        ss_node_index += num_nodes_side
        nodes = coords[:, indices]
        facet_coords[facet] = nodes
        box_min[facet] = SVector{3,Float64}(minimum(nodes[i, :]) for i in 1:3)
        box_max[facet] = SVector{3,Float64}(maximum(nodes[i, :]) for i in 1:3)
    end
    return FacetGrid(BoxGrid(box_min, box_max), facet_coords, box_min, box_max)
end

# The grids depend only on the reference configuration, so they are built once
# per model (and per side set) and memoized by model identity, like the weights.
const ELEMENT_GRID_CACHE = IdDict{SolidMechanics,ElementGrid}()
const FACET_GRID_CACHE = IdDict{SolidMechanics,Dict{Int64,FacetGrid}}()

function get_element_grid(model::SolidMechanics)
    return get!(() -> build_element_grid(model), ELEMENT_GRID_CACHE, model)
end

function get_facet_grid(model::SolidMechanics, side_set_id::Integer)
    grids = get!(() -> Dict{Int64,FacetGrid}(), FACET_GRID_CACHE, model)
    return get!(() -> build_facet_grid(model, side_set_id), grids, Int64(side_set_id))
end

# Distance from `p` to a facet, clamped to the facet: the Newton projection of
# closest_point_projection extends the facet surface beyond its edges, and on a
# side set with corners (the caps and the lateral surface of a cylinder, say)
# the extension of a far facet can pass arbitrarily close to a point. A
# projection that lands outside the parametric domain is clamped to it, and
# the result is bounded by the distances to the facet's edges, which are also
# what a diverging projection falls back to.
function clamped_facet_distance(nodes::Matrix{Float64}, p::SVector{3,Float64})
    num_nodes = size(nodes, 2)
    _, ξ, distance, _ = project_onto_facet(nodes, p)
    best = Inf
    if isfinite(distance)
        tol = 1.0e-8
        if num_nodes == 4
            if abs(ξ[1]) <= 1.0 + tol && abs(ξ[2]) <= 1.0 + tol
                return abs(distance)
            end
            ξc = clamp.(ξ, -1.0, 1.0)
            Nc, _, _ = interpolate(QUAD4, ξc)
            best = norm(SMatrix{3,4,Float64,12}(nodes) * Nc - p)
        elseif ξ[1] >= -tol && ξ[2] >= -tol && ξ[1] + ξ[2] <= 1.0 + tol
            return abs(distance)
        end
    end
    for a in 1:num_nodes
        b = a == num_nodes ? 1 : a + 1
        xa = SVector{3,Float64}(nodes[1, a], nodes[2, a], nodes[3, a])
        xb = SVector{3,Float64}(nodes[1, b], nodes[2, b], nodes[3, b])
        ab = xb - xa
        t = clamp(dot(p - xa, ab) / max(dot(ab, ab), 1.0e-300), 0.0, 1.0)
        best = min(best, norm(xa + t * ab - p))
    end
    return best
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
    w = arlequin_weight_or_nan(x, model, partners, h)
    isnan(w) && abort_degenerate_overlap(x, model, partners, h, label)
    return w
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
    d_partner = sum(
        point_in_model_reference(x, p.partner_model, p.tol) ?
        distance_to_side_set_reference(x, p.partner_model, p.partner_side_set_id) : 0.0 for p in partners
    )
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
    model isa SolidMechanics ||
        norma_abort("Blended overlap energy is only implemented for SolidMechanics subdomains.")
    model.mesh_smoothing &&
        norma_abort("Blended overlap energy does not support mesh-smoothing (EMS) subdomains.")
    partners = overlap_partners(subsim)
    # Build the search grids serially before the threaded classification: the
    # builders read the Exodus files and fill the caches, neither of which is
    # safe from several threads at once.
    for p in partners
        get_element_grid(p.partner_model)
        get_facet_grid(p.partner_model, p.partner_side_set_id)
        get_facet_grid(model, p.own_side_set_id)
    end
    weights = Vector{Vector{Vector{Float64}}}(undef, length(model.blocks))
    for (block_index, block) in enumerate(model.blocks)
        weights[block_index] = compute_block_arlequin_weights(model, block, partners)
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
function compute_block_arlequin_weights(model::SolidMechanics, block::ElementBlockData, partners::Vector{OverlapPartner})
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
            element_weights[point] = arlequin_weight_or_nan(x, model, partners, h)
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
    dynamic = is_dynamic(subsim.integrator)
    integrator = subsim.integrator
    half_dt = 0.0
    if integrator isa CentralDifference
        nominal_dt = integrator.minimum_time_step == integrator.maximum_time_step ?
            integrator.maximum_time_step : integrator.time_step
        half_dt = 0.5 * nominal_dt
    end
    stored_tl = zeros(maxthreadid())
    kinetic_tl = zeros(maxthreadid())
    for (block_index, block) in enumerate(model.blocks)
        material = model.materials[block_index]
        blended_block_energy!(
            stored_tl, kinetic_tl, model, block, block_index, material, material.ρ, aw.weights[block_index], dynamic, half_dt
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
            W = strain_energy_density(material, F, state_old)
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
