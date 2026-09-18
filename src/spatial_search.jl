# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# ---------------------------------------------------------------------------
# Spatial search on the reference configuration.
#
# Locating points in a partner mesh (every overlap coupling, the blended
# energy weights) and measuring distances to a side set (the blended energy
# weights) by scanning every element or facet per point made those setups
# quadratic: hours on meshes of a few tens of thousands of elements. Elements
# and facets are binned once into uniform grids of cells by their reference
# bounding boxes. A point-in-mesh query tests only the elements binned in the
# point's cell, and a distance query visits the cells in rings around the
# point until no unvisited cell can hold a closer facet. The grids depend
# only on the reference configuration, which does not change during a run,
# and are memoized by model identity.
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

function BoxGrid(box_min::Vector{SVector{3,Float64}}, box_max::Vector{SVector{3,Float64}}; cell::Float64=0.0)
    num_items = length(box_min)
    num_items == 0 && return BoxGrid(zero(SVector{3,Float64}), 1.0, SVector{3,Int64}(1, 1, 1), [Int64[]])
    lo = reduce((a, b) -> min.(a, b), box_min)
    hi = reduce((a, b) -> max.(a, b), box_max)
    extent = hi - lo
    if cell <= 0.0
        mean_diagonal = sum(norm(box_max[i] - box_min[i]) for i in 1:num_items) / num_items
        # A cell about twice the mean item diagonal keeps the items per cell
        # small while an item overlaps only a few cells. Point sets have no
        # diagonal, so they get the mean point spacing instead.
        if mean_diagonal > 0.0
            cell = 2.0 * mean_diagonal
        else
            # Mean spacing over the directions the set actually spans: a flat
            # set (an interface) must not get a cell of its vanishing thickness,
            # or a far query would walk thousands of empty rings.
            spanned = [extent[i] for i in 1:3 if extent[i] > 1.0e-6 * maximum(extent)]
            cell = 2.0 * prod(spanned)^(1.0 / length(spanned)) / num_items^(1.0 / length(spanned))
        end
    end
    # The cap bounds the memory.
    cell = max(cell, maximum(extent) / BOX_GRID_MAX_CELLS_PER_AXIS, 1.0e-12)
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
    return box_grid_nearest_item(grid, p, box_min, box_max, distance)[1]
end

# As box_grid_nearest, also returning the index of the nearest item (0 if the
# grid is empty).
function box_grid_nearest_item(
    grid::BoxGrid, p::SVector{3,Float64}, box_min::Vector{SVector{3,Float64}}, box_max::Vector{SVector{3,Float64}}, distance::F
) where {F}
    c = box_grid_cell(grid, p)
    best = Inf
    best_item = 0
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
                d = distance(item)
                if d < best
                    best = d
                    best_item = item
                end
            end
        end
    end
    return best, best_item
end

# A grid over a set of points, for nearest-point queries.
struct PointGrid
    grid::BoxGrid
    points::Vector{SVector{3,Float64}}
end

function PointGrid(points::Vector{SVector{3,Float64}})
    return PointGrid(BoxGrid(points, points), points)
end

function nearest_point_index(pg::PointGrid, p::SVector{3,Float64})
    return box_grid_nearest_item(pg.grid, p, pg.points, pg.points, i -> norm(pg.points[i] - p))[2]
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
            lo = SVector{3,Float64}(minimum(view(nodes, i, :)) for i in 1:3)
            hi = SVector{3,Float64}(maximum(view(nodes, i, :)) for i in 1:3)
            # is_inside accepts points up to a tenth of the element's extent
            # outside its bounding box, so the bins must cover that margin.
            pad = 0.1 * (hi - lo)
            push!(box_min, lo - pad)
            push!(box_max, hi + pad)
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

