# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Converts a MicroStructPy mesh (trimesh.txt or rastermesh.txt) to an Exodus
# mesh for Norma, with a lattice orientation drawn uniformly over SO(3) for
# every grain:
#
#   julia --project=<Norma.jl> microstructpy_to_exodus.jl <mesh.txt> <mesh.g> [--tet10] [--blocks] [--seed N]
#
# The Exodus mesh has
#   - one block named polycrystal, or with --blocks one block per grain,
#     named grain_<g>;
#   - TETRA4 elements, TETRA10 with --tet10 (straight edges), or HEX8 for a
#     voxel mesh;
#   - node sets nsx-, nsx+, nsy-, nsy+, nsz-, nsz+ and side sets ssx-, ...,
#     ssz+ on the faces of the bounding box, and a node set of one node at
#     each corner of the box, named after its faces, nsx-y-z- to nsx+y+z+;
#   - element variables at time 0: grain (from 1), and phi1, Phi and phi2,
#     the Bunge Euler angles in degrees of the lattice orientation of the
#     element's grain, which a material reads with
#     `orientation: {mesh variables: [phi1, Phi, phi2]}`.
#
# A periodic MicroStructPy mesh lists the pairs of nodes on the opposite faces
# of every periodic axis. The converter checks that the mesh conforms across
# those faces: each pair is a translation by the length of the box, the pairs
# cover both faces one to one, and the sides of the elements on one face map
# onto those on the other. For a periodic axis, the node set of the upper face
# lists the images of the nodes of the lower face in order: the i-th node of
# nsx+ is the image of the i-th node of nsx-. Norma's `Periodic` boundary
# condition pairs the nodes of the two sets.

using Exodus
using LinearAlgebra
using Norma
using Printf
using Random

# Sections of a MicroStructPy text mesh: name => lines.
function read_microstructpy(path::String)
    sections = Dict{String,Vector{String}}()
    current = ""
    for line in eachline(path)
        if !isempty(line) && !isspace(line[1])
            current = strip(split(line, ':')[1])
            sections[current] = String[]
        elseif !isempty(strip(line))
            push!(sections[current], strip(line))
        end
    end
    points = reduce(hcat, [parse.(Float64, split(l, ',')) for l in sections["Mesh Points"]])
    elements = reduce(hcat, [parse.(Int, split(l, ',')) .+ 1 for l in sections["Mesh Elements"]])
    grains = [parse(Int, l) + 1 for l in sections["Element Attributes"]]
    # The node pairs of a periodic mesh, as lines "axis, lower, upper"
    # numbered from 0, gathered by axis into lower => upper maps.
    periodic = Dict{Int,Dict{Int,Int}}()
    flags = haskey(sections, "Periodic Axes") ? parse.(Int, split(only(sections["Periodic Axes"]), ',')) : Int[]
    for line in get(sections, "Periodic Nodes", String[])
        axis, lower, upper = parse.(Int, split(line, ','))
        get!(periodic, axis + 1, Dict{Int,Int}())[lower + 1] = upper + 1
    end
    for (axis, flag) in enumerate(flags)
        flag == 1 && !haskey(periodic, axis) && error("Periodic axis $axis has no node pairs.")
    end
    return points, elements, grains, periodic
end

signed_volume(X, e) = det(hcat(X[:, e[2]] - X[:, e[1]], X[:, e[3]] - X[:, e[1]], X[:, e[4]] - X[:, e[1]]))

# Mid-edge nodes of TETRA10 in Exodus order: (1,2), (2,3), (3,1), (1,4), (2,4), (3,4).
const TETRA10_EDGES = ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))

# The mid-edge nodes are added at the midpoints of the edges; also returns
# the mid-edge node of each edge, by its corner nodes in increasing order.
function quadratic_tetrahedra(points::Matrix{Float64}, elements::Matrix{Int})
    positions = [points[:, i] for i in axes(points, 2)]
    midnode = Dict{Tuple{Int,Int},Int}()
    connectivity = zeros(Int, 10, size(elements, 2))
    for e in axes(elements, 2)
        corners = elements[:, e]
        connectivity[1:4, e] = corners
        for (k, (a, b)) in enumerate(TETRA10_EDGES)
            edge = minmax(corners[a], corners[b])
            connectivity[4 + k, e] = get!(midnode, edge) do
                push!(positions, (positions[edge[1]] + positions[edge[2]]) / 2)
                return length(positions)
            end
        end
    end
    return reduce(hcat, positions), connectivity, midnode
end

# The mid-edge nodes of an edge on a lower face pair with those of its image.
function pair_midnodes!(images::Dict{Int,Int}, midnode::Dict{Tuple{Int,Int},Int})
    for ((a, b), node) in midnode
        haskey(images, a) && haskey(images, b) || continue
        image = get(midnode, minmax(images[a], images[b]), 0)
        image == 0 && error("The image of edge ($a, $b) on a periodic face is not an edge of the mesh.")
        images[node] = image
    end
    return images
end

# Checks that a periodic mesh conforms across the faces of `axis`: every node
# of the lower face pairs with a node of the upper face at its translate by
# the length of the box, one to one, and the sides of the elements on the
# lower face map onto those on the upper face.
function check_periodic_faces(
    points,
    connectivity,
    element_type,
    axis,
    images::Dict{Int,Int},
    lower::Vector{Int},
    upper::Vector{Int},
    period,
    tolerance,
)
    sort(collect(keys(images))) == lower ||
        error("The node pairs of axis $axis do not cover the nodes of its lower face.")
    sort(collect(values(images))) == upper ||
        error("The node pairs of axis $axis do not cover the nodes of its upper face one to one.")
    shift = [j == axis ? period : 0.0 for j in 1:3]
    for (a, b) in images
        norm(points[:, b] - points[:, a] - shift) <= tolerance ||
            error("Nodes $a and $b of axis $axis are not periodic images.")
    end
    lower_set, upper_set = Set(lower), Set(upper)
    lower_sides, upper_sides = Set{Set{Int}}(), Set{Set{Int}}()
    for e in axes(connectivity, 2), side_nodes in SIDES[element_type]
        corners = connectivity[collect(side_nodes), e]
        all(in(lower_set), corners) && push!(lower_sides, Set(corners))
        all(in(upper_set), corners) && push!(upper_sides, Set(corners))
    end
    Set(Set(images[n] for n in side) for side in lower_sides) == upper_sides ||
        error("The element sides on the faces of axis $axis do not match: the mesh does not conform.")
    return length(lower_sides)
end

# Corner nodes of the sides of each element type, in Exodus numbering.
const SIDES = Dict(
    "TETRA4" => ((1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2)),
    "TETRA10" => ((1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2)),
    "HEX8" => ((1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (1, 5, 8, 4), (1, 4, 3, 2), (5, 6, 7, 8)),
)

const FACES = (("x-", 1, :min), ("x+", 1, :max), ("y-", 2, :min), ("y+", 2, :max), ("z-", 3, :min), ("z+", 3, :max))

function main(args)
    length(args) >= 2 || error("usage: microstructpy_to_exodus.jl <mesh.txt> <mesh.g> [--tet10] [--blocks] [--seed N]")
    input, output = args[1], args[2]
    tet10 = "--tet10" in args
    per_grain_blocks = "--blocks" in args
    seed = (i = findfirst(==("--seed"), args)) === nothing ? 1 : parse(Int, args[i + 1])
    points, elements, grains, periodic = read_microstructpy(input)
    num_corners = size(elements, 1)
    element_type = if num_corners == 8
        "HEX8"
    elseif num_corners == 4
        "TETRA4"
    else
        error("Unknown element with $num_corners nodes.")
    end
    if element_type == "TETRA4"
        for e in axes(elements, 2)
            if signed_volume(points, elements[:, e]) < 0
                elements[2, e], elements[3, e] = elements[3, e], elements[2, e]
            end
        end
    else
        for e in axes(elements, 2)
            hex = elements[:, e]
            signed_volume(points, (hex[1], hex[2], hex[4], hex[5])) > 0 ||
                error("Voxel $e is not in Exodus HEX8 order.")
        end
    end
    connectivity = elements
    if tet10
        element_type == "TETRA4" || error("--tet10 applies to tetrahedral meshes.")
        points, connectivity, midnode = quadratic_tetrahedra(points, elements)
        element_type = "TETRA10"
        for images in values(periodic)
            pair_midnodes!(images, midnode)
        end
    end
    num_nodes = size(points, 2)
    num_grains = maximum(grains)
    rng = Xoshiro(seed)
    orientations = [Norma.uniform_random_rotation(rng) for _ in 1:num_grains]

    # Blocks: all elements, or the elements of each grain.
    blocks = per_grain_blocks ? [findall(==(g), grains) for g in 1:num_grains] : [collect(eachindex(grains))]
    filter!(!isempty, blocks)
    block_names = per_grain_blocks ? ["grain_$(grains[b[1]])" for b in blocks] : ["polycrystal"]
    global_index = zeros(Int, length(grains))
    next = 0
    for block in blocks, e in block
        next += 1
        global_index[e] = next
    end

    # Faces of the bounding box.
    lo = minimum(points; dims=2)
    hi = maximum(points; dims=2)
    tolerance = 1.0e-9 * maximum(hi - lo)
    on_face(node, axis, side) = abs(points[axis, node] - (side == :min ? lo[axis] : hi[axis])) <= tolerance
    node_sets = [findall(n -> on_face(n, axis, side), 1:num_nodes) for (_, axis, side) in FACES]
    # The upper face of a periodic axis lists the images of the lower face.
    for axis in sort(collect(keys(periodic)))
        images = periodic[axis]
        lower, upper = node_sets[2 * axis - 1], node_sets[2 * axis]
        num_sides = check_periodic_faces(
            points, connectivity, element_type, axis, images, lower, upper, hi[axis] - lo[axis], tolerance
        )
        node_sets[2 * axis] = [images[n] for n in lower]
        @printf("periodic along %s: %d node pairs, %d element sides per face\n", "xyz"[axis], length(lower), num_sides)
    end
    # The node at each corner of the box, if there is one.
    corner_sets = Tuple{String,Vector{Int}}[]
    face_sign(side) = side == :min ? "-" : "+"
    for sz in (:min, :max), sy in (:min, :max), sx in (:min, :max)
        nodes = findall(n -> on_face(n, 1, sx) && on_face(n, 2, sy) && on_face(n, 3, sz), 1:num_nodes)
        isempty(nodes) && continue
        push!(corner_sets, ("nsx" * face_sign(sx) * "y" * face_sign(sy) * "z" * face_sign(sz), nodes))
    end
    side_sets = [(Int[], Int[]) for _ in FACES]
    for e in axes(connectivity, 2), (k, side_nodes) in enumerate(SIDES[element_type])
        corners = connectivity[collect(side_nodes), e]
        for (f, (_, axis, side)) in enumerate(FACES)
            if all(on_face(n, axis, side) for n in corners)
                push!(side_sets[f][1], global_index[e])
                push!(side_sets[f][2], k)
            end
        end
    end

    init = Exodus.Initialization{Int32}(
        Int32(3),
        Int32(num_nodes),
        Int32(length(grains)),
        Int32(length(blocks)),
        Int32(length(FACES) + length(corner_sets)),
        Int32(length(FACES)),
    )
    isfile(output) && rm(output; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(output, "w", init)
    Exodus.write_coordinates(exo, points)
    for (b, block) in enumerate(blocks)
        Exodus.write_block(exo, b, element_type, Matrix{Int32}(connectivity[:, block]))
        Exodus.write_name(exo, Block, b, block_names[b])
    end
    for (f, (name, _, _)) in enumerate(FACES)
        node_set = Exodus.NodeSet(Int32(f), Vector{Int32}(node_sets[f]))
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, "ns" * name)
        side_set = Exodus.SideSet{Int32,Vector{Int32}}(
            Int32(f), Vector{Int32}(side_sets[f][1]), Vector{Int32}(side_sets[f][2]), Int32[], Int32[]
        )
        Exodus.write_set(exo, side_set)
        Exodus.write_name(exo, side_set, "ss" * name)
    end
    for (c, (name, nodes)) in enumerate(corner_sets)
        node_set = Exodus.NodeSet(Int32(length(FACES) + c), Vector{Int32}(nodes))
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, name)
    end
    angles = [rad2deg.(collect(Norma.bunge_angles(Q))) for Q in orientations]
    names = ["grain", "phi1", "Phi", "phi2"]
    Exodus.write_number_of_variables(exo, ElementVariable, length(names))
    Exodus.write_names(exo, ElementVariable, names)
    Exodus.write_time(exo, 1, 0.0)
    for (b, block) in enumerate(blocks)
        Exodus.write_values(exo, ElementVariable, 1, b, "grain", Float64.(grains[block]))
        for (k, name) in enumerate(names[2:end])
            Exodus.write_values(exo, ElementVariable, 1, b, name, [angles[grains[e]][k] for e in block])
        end
    end
    Exodus.close(exo)
    @printf(
        "%s: %d %s elements, %d nodes, %d grains, %d block(s)\n",
        output,
        length(grains),
        element_type,
        num_nodes,
        num_grains,
        length(blocks)
    )
    return nothing
end

# Run as a script; test/polycrystal.jl and test/periodic.jl include the file and
# call main.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
