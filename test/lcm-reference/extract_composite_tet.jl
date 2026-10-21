# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Reference data of LCM's composite tetrahedron (`Use Composite Tet 10: true`)
# for test/composite-tet10-lcm.jl, written to lcm-reference/composite-tet:
#
#   julia --project=. test/lcm-reference/extract_composite_tet.jl <LCM checkout> <runs>
#
# <runs> holds LCM's outputs of the inputs in composite-tet/inputs:
# clamped/clamped_ct.e and clamped/clamped_ct_expl.e (the ClampedSDBC1DCt
# test, written at every step) and cylinder/cylinder.e. ProjectTet10 comes
# from its gold file in the LCM checkout.
#
# Nodal files have one row per (step, node), with the reference coordinates
# and, for ProjectTet10, LCM's nodal projection of the linear field of its
# test; element files one row per (step, element, point), with the Cauchy
# stress and, when LCM wrote it, the deformation gradient, both row-major. The
# first line of each names its source: the gold file and the last LCM commit
# that changed it, or the LCM output and the commit of the LCM checkout, which
# the runs are of.

using Exodus
using Printf

const OUTPUT = joinpath(@__DIR__, "composite-tet")
const MESHES = joinpath(@__DIR__, "meshes")

format_value(v::Real) = string(Float64(v))

function write_nodal(path, exo, steps, nodes, prefix, source; extra=String[])
    times = Exodus.read_times(exo)
    X = read_coordinates(exo)
    open(path, "w") do io
        println(io, "# ", source)
        println(io, join(["step,time,node,X,Y,Z,ux,uy,uz"; extra], ','))
        for step in steps
            u = [read_values(exo, NodalVariable, step, 1, prefix * c) for c in ("x", "y", "z")]
            v = [read_values(exo, NodalVariable, step, 1, name) for name in extra]
            for n in nodes
                println(
                    io,
                    join(
                        [
                            string(step),
                            format_value(times[step]),
                            string(n),
                            format_value.(X[:, n])...,
                            (format_value(u[c][n]) for c in 1:3)...,
                            (format_value(w[n]) for w in v)...,
                        ],
                        ',',
                    ),
                )
            end
        end
    end
    return nothing
end

# LCM names the component (i, j) of a tensor at point q of a tetrahedron
# <field>_<j>_<k>, with k = 3 (q - 1) + i written with two digits: the
# components of a point are in column-major order.
function write_element(path, exo, steps, num_points, fields, source)
    names = Exodus.read_names(exo, ElementVariable)
    times = Exodus.read_times(exo)
    block_id = Exodus.read_sets(exo, Block)[1].id
    header = ["step", "time", "element", "point"]
    for field in fields, k in 1:9
        push!(header, "$(field)_$k")
    end
    open(path, "w") do io
        println(io, "# ", source)
        println(io, join(header, ','))
        for step in steps
            values = Dict(name => read_values(exo, ElementVariable, step, block_id, name) for name in names)
            num_elements = length(first(values)[2])
            for e in 1:num_elements, q in 1:num_points
                row = [string(step), format_value(times[step]), string(e), string(q)]
                for field in fields, i in 1:3, j in 1:3
                    push!(row, format_value(values[@sprintf("%s_%d_%02d", field, j, 3 * (q - 1) + i)][e]))
                end
                println(io, join(row, ','))
            end
        end
    end
    return nothing
end

# LCM's ClampedSDBC1DCt mesh with the node set every node is in, which the
# Gaussian initial displacement of the case is applied to.
function write_clamped_mesh(source, target)
    src = ExodusDatabase(source, "r")
    coords = read_coordinates(src)
    block = Exodus.read_sets(src, Block)[1]
    _, num_elements, num_element_nodes = Exodus.read_block_parameters(src, block.id)
    connectivity = reshape(
        Exodus.read_block_connectivity(src, block.id, num_elements * num_element_nodes), :, num_elements
    )
    node_sets = Exodus.read_sets(src, NodeSet)
    node_set_names = Exodus.read_names(src, NodeSet)
    close(src)
    num_nodes = size(coords, 2)
    init = Exodus.Initialization{Int32}(
        Int32(3), Int32(num_nodes), Int32(num_elements), Int32(1), Int32(length(node_sets) + 1), Int32(0)
    )
    rm(target; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(target, "w", init)
    Exodus.write_coordinates(exo, coords)
    Exodus.write_block(exo, 1, "TETRA10", Matrix{Int32}(connectivity))
    Exodus.write_name(exo, Block, 1, "Block0")
    for (k, (set, name)) in enumerate(zip(node_sets, node_set_names))
        node_set = Exodus.NodeSet(Int32(k), Int32.(set.nodes))
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, name)
    end
    node_set = Exodus.NodeSet(Int32(length(node_sets) + 1), Int32.(1:num_nodes))
    Exodus.write_set(exo, node_set)
    Exodus.write_name(exo, node_set, "NodeSetAll")
    close(exo)
    return nothing
end

function main(args)
    length(args) == 2 || error("usage: extract_composite_tet.jl <LCM checkout> <runs>")
    lcm_root, runs = args
    mkpath(OUTPUT)
    # The last commit that touched a gold file, and the checked-out commit,
    # which the runs are of.
    gold_source(relative) = "LCM gold file: $relative (LCM $(readchomp(`git -C $lcm_root log -1 --abbrev=7 --format=%h -- $relative`)))"
    run_source(output) = "output $output of an LCM run (LCM $(readchomp(`git -C $lcm_root rev-parse --short=7 HEAD`)))"

    project = joinpath(lcm_root, "tests", "LCM", "ProjectTet10")
    # Its node sets have no names; LCM calls them nodelist_<id>.
    mesh = joinpath(MESHES, "cube-single-tet10.g")
    cp(joinpath(project, "cube-single-tet10.g"), mesh; force=true)
    exo = ExodusDatabase(mesh, "rw")
    write_names(exo, NodeSet, ["nodelist_$id" for id in read_ids(exo, NodeSet)])
    close(exo)
    exo = ExodusDatabase(joinpath(project, "project_tet10_gold.e"), "r")
    source = gold_source("tests/LCM/ProjectTet10/project_tet10_gold.e")
    write_nodal(
        joinpath(OUTPUT, "ProjectTet10-nodal.csv"),
        exo,
        [1],
        1:size(read_coordinates(exo), 2),
        "solution_",
        source;
        extra=["proj_nodal_linear"],
    )
    write_element(joinpath(OUTPUT, "ProjectTet10.csv"), exo, [1], 4, ["Cauchy_Stress"], source)
    close(exo)

    clamped = joinpath(lcm_root, "tests", "LCM", "DynamicsTempus", "ClampedSDBC1DCt_AnalyticMass")
    write_clamped_mesh(joinpath(clamped, "clamped_ct.g"), joinpath(MESHES, "clamped-ct.g"))
    # The nodes on the edge x = y = 0 of the bar, 2001 of its 18009.
    for (file, case, steps) in (
        ("clamped_ct.e", "ClampedCt_Implicit", [1, 2, 4, 6, 11, 16]),
        ("clamped_ct_expl.e", "ClampedCt_Explicit", [1, 51, 101, 151]),
    )
        exo = ExodusDatabase(joinpath(runs, "clamped", file), "r")
        X = read_coordinates(exo)
        edge = [n for n in 1:size(X, 2) if abs(X[1, n]) < 1e-12 && abs(X[2, n]) < 1e-12]
        write_nodal(joinpath(OUTPUT, "$case-nodal.csv"), exo, steps, edge, "disp_", run_source("clamped/$file"))
        close(exo)
    end

    exo = ExodusDatabase(joinpath(runs, "cylinder", "cylinder.e"), "r")
    num_steps = length(Exodus.read_times(exo))
    steps = [1, (num_steps + 1) ÷ 2, num_steps]
    source = run_source("cylinder/cylinder.e")
    write_nodal(
        joinpath(OUTPUT, "Cylinder-nodal.csv"), exo, steps, 1:size(read_coordinates(exo), 2), "displacement_", source
    )
    write_element(joinpath(OUTPUT, "Cylinder.csv"), exo, [num_steps], 4, ["Cauchy_Stress", "F"], source)
    close(exo)
    return nothing
end

main(ARGS)
