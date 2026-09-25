# Meshes of the dynamic stability study.
#
# The beam is a brick, meshed here directly as a structured grid of HEX8
# elements, the same grid that Cubit produces for `create brick` with a
# uniform size: x along the length from `x0`, y across the height centered
# on zero, z across the width from zero.  Each mesh carries the node sets
# nsx- (the face x = x0) and nsall, and the side sets ssx- and ssx+ (the
# faces x = x0 and x = x0 + length).
#
# The cylinders come from the Cubit journal of
# examples/jmp/concentric-cylinders.  Level 1 copies the meshes committed
# there; a higher level plays the journal back with the variable
# `refinement` set to the level, which needs Cubit (the executable `cubit`
# on PATH, or its path in the environment variable CUBIT).

using Exodus
using Norma

const CYLINDER_EXAMPLE = normpath(joinpath(@__DIR__, "..", "..", "examples", "jmp", "concentric-cylinders"))
const CYLINDER_MESHES = [
    "monolithic/cylinders.g", "nonoverlap/inner.g", "nonoverlap/outer.g", "overlap/inner.g", "overlap/outer.g"
]

# Number of elements along a length for a target size, as Cubit rounds it.
intervals(length, h) = max(1, round(Int, length / h))

function write_brick(file::AbstractString, x0::Real, len::Real, h::Real, block_name::AbstractString)
    nx = intervals(len, h)
    ny = intervals(BEAM.height, h)
    nz = intervals(BEAM.width, h)
    node(i, j, k) = 1 + i + (nx + 1) * (j + (ny + 1) * k)
    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    coordinates = zeros(3, num_nodes)
    for k in 0:nz, j in 0:ny, i in 0:nx
        coordinates[:, node(i, j, k)] = [x0 + len * i / nx, BEAM.height * (j / ny - 0.5), BEAM.width * k / nz]
    end
    num_elements = nx * ny * nz
    connectivity = zeros(Int32, 8, num_elements)
    element(i, j, k) = 1 + i + nx * (j + ny * k)
    for k in 0:(nz - 1), j in 0:(ny - 1), i in 0:(nx - 1)
        connectivity[:, element(i, j, k)] = [
            node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
            node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1), node(i, j + 1, k + 1),
        ]
    end
    init = Exodus.Initialization{Int32}(Int32(3), Int32(num_nodes), Int32(num_elements), Int32(1), Int32(2), Int32(2))
    rm(file; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(file, "w", init)
    Exodus.write_coordinates(exo, coordinates)
    Exodus.write_block(exo, 1, "HEX8", connectivity)
    Exodus.write_name(exo, Block, 1, block_name)
    x_minus = Int32[node(0, j, k) for k in 0:nz for j in 0:ny]
    for (id, name, nodes) in ((1, "nsx-", x_minus), (2, "nsall", Int32.(1:num_nodes)))
        node_set = Exodus.NodeSet(Int32(id), nodes)
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, name)
    end
    # Exodus numbers the faces of a HEX8 with side 4 at ξ = -1 (x-) and side 2
    # at ξ = +1 (x+).
    first_elements = Int32[element(0, j, k) for k in 0:(nz - 1) for j in 0:(ny - 1)]
    last_elements = Int32[element(nx - 1, j, k) for k in 0:(nz - 1) for j in 0:(ny - 1)]
    for (id, name, elements, side) in ((1, "ssx-", first_elements, 4), (2, "ssx+", last_elements, 2))
        sides = fill(Int32(side), length(elements))
        side_set = Exodus.SideSet{Int32,Vector{Int32}}(Int32(id), elements, sides, Int32[], Int32[])
        Exodus.write_set(exo, side_set)
        Exodus.write_name(exo, side_set, name)
    end
    Exodus.close(exo)
    return (nodes=num_nodes, elements=num_elements)
end

# The beam meshes of a case, written into `dir` under the names of its
# subdomains: beam.g for the reference; clamped.g and free.g for a coupled
# case.  The free part keeps the size h / level and the clamped part is
# coarsened by the mesh ratio.
function write_beam_meshes(c::Case, dir::AbstractString)
    h = BEAM.h / c.level
    if c.coupling == "mono"
        write_brick(joinpath(dir, "beam.g"), 0.0, BEAM.length, h, "beam")
        return nothing
    end
    if startswith(c.coupling, "ov")
        clamped_end = BEAM.split + BEAM.overlap / 2
        free_start = BEAM.split - BEAM.overlap / 2
    else
        clamped_end = free_start = BEAM.split
    end
    write_brick(joinpath(dir, "clamped.g"), 0.0, clamped_end, h / c.ratio, "clamped")
    write_brick(joinpath(dir, "free.g"), free_start, BEAM.length - free_start, h, "free")
    return nothing
end

# The five cylinder meshes at a refinement level, in `root`/cyl-L<level>,
# made once and reused by every case of that level.
function cylinder_mesh_dir(root::AbstractString, level::Int)
    dir = joinpath(root, "cyl-L$level")
    all(isfile(joinpath(dir, m)) for m in CYLINDER_MESHES) && return dir
    for sub in ("monolithic", "nonoverlap", "overlap")
        mkpath(joinpath(dir, sub))
    end
    if level == 1
        for m in CYLINDER_MESHES
            cp(joinpath(CYLINDER_EXAMPLE, m), joinpath(dir, m); force=true)
        end
        return dir
    end
    cubit = get(ENV, "CUBIT", "cubit")
    Sys.which(cubit) === nothing && !isfile(cubit) &&
        error("Cubit is needed for the cylinder meshes at level $level: put cubit on PATH or set CUBIT")
    wrapper = joinpath(dir, "refinement.jou")
    journal = joinpath(CYLINDER_EXAMPLE, "concentric-cylinders.jou")
    write(wrapper, "\${refinement = $level}\nplayback \"$journal\"\n")
    println("Meshing the cylinders at level $level with Cubit; this takes a while")
    cd(dir) do
        command = `$cubit -batch -nographics -nojournal -noecho refinement.jou`
        run(pipeline(command; stdout="cubit.log", stderr="cubit.log"))
    end
    for m in CYLINDER_MESHES
        isfile(joinpath(dir, m)) || error("Cubit did not write $m; see $(joinpath(dir, "cubit.log"))")
    end
    return dir
end

# The cylinder meshes of a case, linked into its directory under the names
# of its subdomains.
function link_cylinder_meshes(c::Case, dir::AbstractString, mesh_root::AbstractString)
    source = cylinder_mesh_dir(mesh_root, c.level)
    files = if c.coupling == "mono"
        ["monolithic/cylinders.g" => "cylinders.g"]
    else
        sub = startswith(c.coupling, "ov") ? "overlap" : "nonoverlap"
        ["$sub/inner.g" => "inner.g", "$sub/outer.g" => "outer.g"]
    end
    for (from, to) in files
        target = joinpath(dir, to)
        (islink(target) || isfile(target)) && rm(target)
        symlink(relpath(joinpath(source, from), dir), target)
    end
    return nothing
end
