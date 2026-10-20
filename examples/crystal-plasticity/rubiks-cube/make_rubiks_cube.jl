# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Builds the meshes of LCM's RubiksCube polycrystal (tests/LCM/
# CrystalPlasticity/RubiksCube): the cube [-0.5, 0.5]³ split into 3 × 3 × 3
# grains with LCM's 27 random lattice orientations, each grain meshed with
# n × n × n HEX8 elements:
#
#   julia --project=<Norma.jl> make_rubiks_cube.jl <LCM checkout> [n ...]
#
# writes rubiks-cube-<n>.g for each n (default 1 and 2). Every mesh has one
# block per grain, named block_<k> as in LCM, node sets nsx-, ..., nsz+ and
# side sets ssx-, ..., ssz+ on the faces, and element variables phi1, Phi
# and phi2: the Bunge Euler angles in degrees of the element's grain.

using Exodus
using LinearAlgebra
using Norma
using Printf
using StaticArrays
using YAML

const LCM_RUBIKS_CUBE = "tests/LCM/CrystalPlasticity/RubiksCube"

# Lattice orientation of each LCM block, and the 3 × 3 × 3 cell it occupies.
function lcm_grains(lcm_root::String)
    mesh = ExodusDatabase(joinpath(lcm_root, LCM_RUBIKS_CUBE, "RubiksCube.g"), "r")
    coords = read_coordinates(mesh)
    cells = Dict{Int,NTuple{3,Int}}()
    for block_id in read_ids(mesh, Block)
        _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(mesh, block_id)
        nodes = Exodus.read_block_connectivity(mesh, block_id, num_elements * num_nodes)
        center = vec(sum(coords[:, nodes]; dims=2) / length(nodes))
        cells[Int(block_id)] = Tuple(clamp.(floor.(Int, (center .+ 0.5) .* 3) .+ 1, 1, 3))
    end
    close(mesh)
    materials = YAML.load_file(joinpath(lcm_root, LCM_RUBIKS_CUBE, "RubiksCube_Implicit_Materials.yaml"))["LCM"]
    orientations = Dict{Int,SMatrix{3,3,Float64,9}}()
    for (name, block) in materials["ElementBlocks"]
        elasticity = materials["Materials"][block["material"]]["Crystal Elasticity"]
        basis = [Float64.(elasticity["Basis Vector $k"]) for k in 1:3]
        orientations[parse(Int, split(name, '_')[end])] = Norma.orientation_from_basis_vectors(basis)
    end
    return cells, orientations
end

const SIDES = ((1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (1, 5, 8, 4), (1, 4, 3, 2), (5, 6, 7, 8))
const FACES = (("x-", 1, -0.5), ("x+", 1, 0.5), ("y-", 2, -0.5), ("y+", 2, 0.5), ("z-", 3, -0.5), ("z+", 3, 0.5))

function write_rubiks_cube(path::String, n::Int, cells, orientations)
    m = 3n                          # elements per side
    node(i, j, k) = 1 + i + (m + 1) * (j + (m + 1) * k)
    points = [(-0.5 + i / m, -0.5 + j / m, -0.5 + k / m) for k in 0:m for j in 0:m for i in 0:m]
    coordinates = reduce(hcat, [collect(p) for p in points])
    block_ids = sort(collect(keys(cells)))
    blocks = Dict(b => Vector{Int}[] for b in block_ids)
    for (b, (ci, cj, ck)) in cells
        for k in ((ck - 1) * n):(ck * n - 1), j in ((cj - 1) * n):(cj * n - 1), i in ((ci - 1) * n):(ci * n - 1)
            push!(
                blocks[b],
                [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i + 1, j + 1, k + 1),
                    node(i, j + 1, k + 1),
                ],
            )
        end
    end
    num_elements = sum(length(blocks[b]) for b in block_ids)
    init = Exodus.Initialization{Int32}(
        Int32(3), Int32(size(coordinates, 2)), Int32(num_elements), Int32(length(block_ids)), Int32(6), Int32(6)
    )
    isfile(path) && rm(path; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(path, "w", init)
    Exodus.write_coordinates(exo, coordinates)
    side_sets = [(Int32[], Int32[]) for _ in FACES]
    element = 0
    for b in block_ids
        Exodus.write_block(exo, b, "HEX8", Matrix{Int32}(reduce(hcat, blocks[b])))
        Exodus.write_name(exo, Block, b, "block_$b")
        for hex in blocks[b]
            element += 1
            for (s, side) in enumerate(SIDES), (f, (_, axis, value)) in enumerate(FACES)
                if all(abs(coordinates[axis, hex[c]] - value) < 1.0e-12 for c in side)
                    push!(side_sets[f][1], element)
                    push!(side_sets[f][2], s)
                end
            end
        end
    end
    for (f, (name, axis, value)) in enumerate(FACES)
        nodes = Int32.(findall(p -> abs(coordinates[axis, p] - value) < 1.0e-12, axes(coordinates, 2)))
        node_set = Exodus.NodeSet(Int32(f), nodes)
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, "ns" * name)
        side_set = Exodus.SideSet{Int32,Vector{Int32}}(Int32(f), side_sets[f][1], side_sets[f][2], Int32[], Int32[])
        Exodus.write_set(exo, side_set)
        Exodus.write_name(exo, side_set, "ss" * name)
    end
    names = ["phi1", "Phi", "phi2"]
    Exodus.write_number_of_variables(exo, ElementVariable, 3)
    Exodus.write_names(exo, ElementVariable, names)
    Exodus.write_time(exo, 1, 0.0)
    for b in block_ids
        angles = rad2deg.(collect(Norma.bunge_angles(orientations[b])))
        for (k, name) in enumerate(names)
            Exodus.write_values(exo, ElementVariable, 1, b, name, fill(angles[k], length(blocks[b])))
        end
    end
    Exodus.close(exo)
    @printf("%s: %d grains, %d HEX8 elements\n", path, length(block_ids), num_elements)
    return nothing
end

function main(args)
    isempty(args) && error("usage: make_rubiks_cube.jl <LCM checkout> [n ...]")
    cells, orientations = lcm_grains(args[1])
    length(unique(values(cells))) == 27 || error("The LCM blocks do not fill the 3 × 3 × 3 cells.")
    levels = length(args) > 1 ? parse.(Int, args[2:end]) : [1, 2]
    for n in levels
        write_rubiks_cube(joinpath(@__DIR__, "rubiks-cube-$n.g"), n, cells, orientations)
    end
    return nothing
end

main(ARGS)
