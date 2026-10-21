# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Splits every HEX8 of LCM's RubiksCube.g into 6 TETRA10 elements, for the
# composite tet comparison with LCM. It is not part of the test suite:
#
#   julia --project=. test/lcm-reference/make_rubiks_cube_tet10.jl <RubiksCube.g> <RubiksCubeTet10.g>
#
# Each cube is cut along its diagonal from its lowest to its highest corner
# into the 6 Kuhn tetrahedra, the same way in every cube, so the faces match
# between neighbors. The midside nodes are at the middle of the straight
# edges and are shared. The blocks keep their ids and get one 6-tet set
# each; the node sets keep their ids and gain the midside nodes on their
# faces. Names are LCM's defaults, block_<id> and nodelist_<id>.

using Exodus
using LinearAlgebra

function main(input, output)
    exo = ExodusDatabase(input, "r")
    X = read_coordinates(exo)
    block_ids = read_ids(exo, Block)
    node_set_ids = read_ids(exo, NodeSet)
    node_sets = Dict(id => Int.(Exodus.read_node_set_nodes(exo, id)) for id in node_set_ids)
    hexes = Dict{Int,Matrix{Int}}()
    for id in block_ids
        _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(exo, id)
        hexes[id] = reshape(Int.(Exodus.read_block_connectivity(exo, id, num_elements * num_nodes)), num_nodes, :)
    end
    close(exo)
    points = [X[:, n] for n in axes(X, 2)]
    midside = Dict{Tuple{Int,Int},Int}()
    function midside_node(a, b)
        key = minmax(a, b)
        return get!(midside, key) do
            push!(points, (points[a] + points[b]) / 2)
            return length(points)
        end
    end
    # Kuhn tetrahedra of the cube with corners c[(i, j, k)], i, j, k ∈ {0, 1}:
    # the paths from (0,0,0) to (1,1,1) along the edges, one per order of axes.
    orders = [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]
    tets = Dict{Int,Vector{Vector{Int}}}()
    for id in block_ids
        tets[id] = Vector{Int}[]
        for e in axes(hexes[id], 2)
            nodes = hexes[id][:, e]
            center = sum(X[:, nodes]; dims=2) / 8
            corner = Dict(Tuple(Int.(X[:, n] .> center[:])) => n for n in nodes)
            for order in orders
                path = [(0, 0, 0)]
                for axis in order
                    next = collect(path[end])
                    next[axis] = 1
                    push!(path, Tuple(next))
                end
                v = [corner[p] for p in path]
                volume = det(hcat((X[:, v[k]] - X[:, v[1]] for k in 2:4)...))
                volume < 0 && ((v[3], v[4]) = (v[4], v[3]))
                push!(
                    tets[id],
                    vcat(v, [midside_node(v[a], v[b]) for (a, b) in ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))]),
                )
            end
        end
    end
    P = hcat(points...)
    # A node set is a face of the box: it gains every midside node whose two
    # ends are in it.
    for (id, nodes) in node_sets
        members = Set(nodes)
        for ((a, b), m) in midside
            a in members && b in members && push!(nodes, m)
        end
        sort!(nodes)
    end
    init = Exodus.Initialization{Int32}(
        Int32(3),
        Int32(size(P, 2)),
        Int32(sum(length, values(tets))),
        Int32(length(block_ids)),
        Int32(length(node_set_ids)),
        Int32(0),
    )
    isfile(output) && rm(output; force=true)
    out = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(output, "w", init)
    Exodus.write_coordinates(out, P)
    for id in block_ids
        Exodus.write_block(out, Int(id), "TETRA10", Matrix{Int32}(hcat(tets[id]...)))
        Exodus.write_name(out, Block, Int(id), "block_$id")
    end
    for id in node_set_ids
        node_set = Exodus.NodeSet(Int32(id), Vector{Int32}(node_sets[id]))
        Exodus.write_set(out, node_set)
        Exodus.write_name(out, node_set, "nodelist_$id")
    end
    Exodus.close(out)
    println(
        output, ": ", size(P, 2), " nodes, ", sum(length, values(tets)), " TETRA10 in ", length(block_ids), " blocks"
    )
    return nothing
end

length(ARGS) == 2 || error("usage: make_rubiks_cube_tet10.jl <RubiksCube.g> <RubiksCubeTet10.g>")
main(ARGS[1], ARGS[2])
