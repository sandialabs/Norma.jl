# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Writes the mesh of make_cylinder.py as Exodus: one block (block_1), node
# sets nodelist_1 (z = 0), nodelist_2 (z = 2) and nodelist_3 (every node),
# the names LCM gives them by default. Elements are oriented positively and
# their edge nodes checked against the Exodus TETRA10 order.
using Exodus, DelimitedFiles, LinearAlgebra
x = Matrix(readdlm("coords.txt")');
conn = Matrix(Int32.(readdlm("conn.txt", Int)'))
# orient and check edge nodes
for e in 1:size(conn, 2)
    X = x[:, conn[:, e]]
    if det(hcat(X[:, 2]-X[:, 1], X[:, 3]-X[:, 1], X[:, 4]-X[:, 1])) < 0
        conn[:, e] = conn[[1, 3, 2, 4, 7, 6, 5, 8, 10, 9], e]
    end
    X = x[:, conn[:, e]]
    for (m, (a, b)) in enumerate(((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)))
        norm(X[:, 4 + m] - (X[:, a]+X[:, b])/2) < 0.2 * norm(X[:, a]-X[:, b]) || error("edge order")
    end
end
H = maximum(x[3, :])
bottom = Int32.(findall(abs.(x[3, :]) .< 1e-9));
top = Int32.(findall(abs.(x[3, :] .- H) .< 1e-9))
init = Exodus.Initialization{Int32}(Int32(3), Int32(size(x, 2)), Int32(size(conn, 2)), Int32(1), Int32(3), Int32(0))
exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(ARGS[1], "w", init)
Exodus.write_coordinates(exo, x)
Exodus.write_block(exo, 1, "TETRA10", conn);
Exodus.write_name(exo, Block, 1, "block_1")
for (i, (nodes, name)) in enumerate(((bottom, "nodelist_1"), (top, "nodelist_2"), (Int32.(1:size(x, 2)), "nodelist_3")))
    s = Exodus.NodeSet(Int32(i), nodes)
    Exodus.write_set(exo, s)
    Exodus.write_name(exo, s, name)
end
close(exo)
# curvature check: max distance of edge nodes from chord midpoints
dev = maximum(
    norm(x[:, conn[4 + m, e]] - (x[:, conn[a, e]] + x[:, conn[b, e]]) / 2) for e in 1:size(conn, 2) for
    (m, (a, b)) in enumerate(((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)))
)
println("bottom ", length(bottom), " top ", length(top), " max midside offset ", dev)
