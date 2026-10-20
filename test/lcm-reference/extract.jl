# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Extracts the reference ("gold") outputs of the elastic tests of LCM's
# crystal plasticity into CSV files that the Norma tests read. It is not part
# of the test suite; run it only to regenerate the fixtures:
#
#   julia --project=. test/lcm-reference/extract.jl <LCM checkout> [output dir]
#
# Each case gives <case>.csv, the element variables, one row per (step,
# block, element, point); when every integration point of an element carries
# the same values, the points are collapsed into one row with point = 0. The
# first line of the file names its source: the gold file and the last LCM
# commit that changed it.
#
# Tensor components are numbered row-major for a 3 × 3 tensor, so that
# Cauchy_Stress_k is σ[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1] in every file. The 2D
# cases only carry the in-plane components 1, 2, 4 and 5.
#
# LCM names element variables differently for each layout:
#   :hex8   <field>_<k>_<q>  tensor component k at integration point q, and
#           <field>_<q>      scalar fields such as gamma_<s> or eqps
#   :quad4  a 2 × 2 tensor at 4 points is flattened into <field>_<a>_<b> with
#           flat index f = 2 (b - 1) + (a - 1), point f ÷ 4 + 1 and row-major
#           2 × 2 component f % 4 + 1; scalar fields are <field>_<q>

using Exodus
using LinearAlgebra
using Printf
using Statistics

const CRYSTAL_PLASTICITY = "tests/LCM/CrystalPlasticity"

const CASES = [
    ("ElasticShear2d", :quad4, "$CRYSTAL_PLASTICITY/SingleElement2d/ElasticShear2d.gold.e"),
    ("ElasticShear3d", :hex8, "$CRYSTAL_PLASTICITY/SingleElement3d/ElasticShear3d.gold.e"),
    ("ElasticShearRotated3d", :hex8, "$CRYSTAL_PLASTICITY/SingleElement3d/ElasticShearRotated3d.gold.e"),
]

# The LCM material file of each case, copied to materials/<case>.yaml.
const MATERIAL_FILES = Dict(
    "ElasticShear2d" => "$CRYSTAL_PLASTICITY/SingleElement2d/CrystalElasticity2d_Material.yaml",
    "ElasticShear3d" => "$CRYSTAL_PLASTICITY/SingleElement3d/CrystalElasticity3d_Material.yaml",
    "ElasticShearRotated3d" => "$CRYSTAL_PLASTICITY/SingleElement3d/CrystalElasticityRotated3d_Material.yaml",
)

# Meshes of the cases that the FEM tests run, copied to meshes/ with the
# block and node set names LCM gives them by default (block_<id>,
# nodelist_<id>), which the Norma inputs refer to.
const MESHES = Dict("RubiksCube.g" => "$CRYSTAL_PLASTICITY/RubiksCube/RubiksCube.g")

function copy_named_mesh(source::String, target::String)
    input = ExodusDatabase(source, "r")
    Exodus.copy(input, target)
    close(input)
    mesh = ExodusDatabase(target, "rw")
    write_names(mesh, Block, ["block_$id" for id in read_ids(mesh, Block)])
    write_names(mesh, NodeSet, ["nodelist_$id" for id in read_ids(mesh, NodeSet)])
    close(mesh)
    return nothing
end

# Deformation gradient of every element at every step, the affine fit of the
# nodal displacements over the element's nodes, which is exact for the
# homogeneous deformations of these cases. Keyed by block id; each entry is
# num_elements × 9 × num_steps, row-major. Nothing for meshes without
# displacements.
function nodal_deformation_gradients(exo, block_ids, num_steps)
    nodal_names = read_names(exo, NodalVariable)
    coords = read_coordinates(exo)
    dim = size(coords, 1)
    # Older gold files number the displacement components instead of naming them.
    components = if "displacement_x" in nodal_names
        ["displacement_x", "displacement_y", "displacement_z"][1:dim]
    elseif "displacement_1" in nodal_names
        ["displacement_1", "displacement_2", "displacement_3"][1:dim]
    else
        return nothing
    end
    displacements = [
        reduce(vcat, [permutedims(Float64.(read_values(exo, NodalVariable, step, c))) for c in components]) for
        step in 1:num_steps
    ]
    gradients = Dict{Int,Array{Float64,3}}()
    for block_id in block_ids
        _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(exo, block_id)
        connectivity = reshape(
            Exodus.read_block_connectivity(exo, block_id, num_elements * num_nodes), num_nodes, num_elements
        )
        F = zeros(num_elements, 9, num_steps)
        for e in 1:num_elements, step in 1:num_steps
            nodes = connectivity[:, e]
            X = coords[:, nodes]
            x = X + displacements[step][:, nodes]
            Xc = X .- mean(X; dims=2)
            xc = x .- mean(x; dims=2)
            F3 = Matrix{Float64}(I, 3, 3)
            F3[1:dim, 1:dim] = (xc * Xc') / (Xc * Xc')
            F[e, :, step] = vec(permutedims(F3))
        end
        gradients[Int(block_id)] = F
    end
    return gradients
end

# Variables that describe the discretization rather than the material response.
const SKIPPED_VARIABLES = ["Weights"]

# Trailing integers of a name: "Cauchy_Stress_3_8" -> ("Cauchy_Stress", [3, 8]).
function trailing_integers(name::AbstractString)
    parts = split(name, '_')
    integers = Int[]
    while length(parts) > 1
        value = tryparse(Int, parts[end])
        value === nothing && break
        pushfirst!(integers, value)
        pop!(parts)
    end
    return join(parts, '_'), integers
end

# Maps an LCM element-variable name to (column name, integration point).
function column_and_point(name::AbstractString, layout::Symbol, tensor_fields::Set{String})
    stem, integers = trailing_integers(name)
    if stem in tensor_fields && length(integers) >= 2
        a, b = integers[end - 1], integers[end]
        prefix = length(integers) > 2 ? stem * "_" * join(integers[1:(end - 2)], '_') : stem
        if layout == :quad4
            flat = 2 * (b - 1) + (a - 1)
            component_2d = flat % 4 + 1
            row, col = (component_2d - 1) ÷ 2 + 1, (component_2d - 1) % 2 + 1
            return @sprintf("%s_%d", prefix, 3 * (row - 1) + col), flat ÷ 4 + 1
        end
        return @sprintf("%s_%d", prefix, a), b
    end
    isempty(integers) && return String(name), 0
    head, _ = rsplit(name, '_'; limit=2)
    return String(head), integers[end]
end

# Fields whose names carry a tensor component before the point index, i.e.
# <field>_<a>_<b> with a = 1..9 and b = 1..8 (:hex8) or a = 1..2 and b = 1..8
# (:quad4). Per-slip-system scalars <field>_<s>_<q> also have two trailing
# integers, but in 2D their point index only runs to 4.
function tensor_fields_of(names::Vector{String}, layout::Symbol)
    first_indices = Dict{String,Set{Int}}()
    second_indices = Dict{String,Set{Int}}()
    for name in names
        stem, integers = trailing_integers(name)
        length(integers) == 2 || continue
        push!(get!(first_indices, stem, Set{Int}()), integers[1])
        push!(get!(second_indices, stem, Set{Int}()), integers[2])
    end
    expected_first = layout == :quad4 ? Set(1:2) : Set(1:9)
    expected_second = Set(1:8)
    return Set(
        stem for
        stem in keys(first_indices) if first_indices[stem] == expected_first && second_indices[stem] == expected_second
    )
end

format_value(v::Real) = isfinite(v) ? string(Float64(v)) : "NaN"

function extract_case(name, layout, gold, output_dir, origin, lcm_commit)
    exo = ExodusDatabase(gold, "r")
    times = read_times(exo)
    num_steps = length(times)
    block_ids = read_ids(exo, Block)
    element_names = read_names(exo, ElementVariable)
    tensor_fields = tensor_fields_of(element_names, layout)
    columns = String[]
    points = Int[]
    source = Dict{Tuple{String,Int},String}()
    for var_name in element_names
        column, point = column_and_point(var_name, layout, tensor_fields)
        column in SKIPPED_VARIABLES && continue
        column in columns || push!(columns, column)
        point in points || push!(points, point)
        source[(column, point)] = var_name
    end
    sort!(points)

    # values[block][(column, point)] is a num_elements × num_steps matrix
    values = Dict{Int,Dict{Tuple{String,Int},Matrix{Float64}}}()
    for block_id in block_ids
        values[block_id] = Dict(
            key => reduce(
                hcat, [Float64.(read_values(exo, ElementVariable, step, block_id, var)) for step in 1:num_steps]
            ) for (key, var) in source
        )
    end

    gradients = nodal_deformation_gradients(exo, block_ids, num_steps)
    if gradients !== nothing
        for k in 1:9
            push!(columns, "F_nodal_$k")
            for block_id in block_ids, point in points
                values[block_id][("F_nodal_$k", point)] = gradients[Int(block_id)][:, k, :]
            end
        end
    end

    # Collapse the integration points of an element when they all agree.
    max_spread = 0.0
    for block_id in block_ids, column in columns
        reference = get(values[block_id], (column, points[1]), nothing)
        reference === nothing && continue
        scale = max(1.0, maximum(abs, reference))
        for point in points[2:end]
            other = get(values[block_id], (column, point), nothing)
            other === nothing && continue
            max_spread = max(max_spread, maximum(abs, other .- reference) / scale)
        end
    end
    collapse = length(points) > 1 && max_spread <= 1.0e-9
    row_points = collapse ? [0] : points

    open(joinpath(output_dir, name * ".csv"), "w") do io
        println(io, "# ", origin, " (LCM ", lcm_commit, ")")
        @printf(
            io,
            "# %d integration point(s) per element; max relative spread between points %.1e",
            length(points),
            max_spread
        )
        println(io, collapse ? " (collapsed into point 0)" : "")
        block_names = read_names(exo, Block)
        println(io, "# block names: ", join(("$id => $name" for (id, name) in zip(block_ids, block_names)), ", "))
        gradients === nothing || println(io, "# F_nodal_k: affine fit of the nodal displacements over each element")
        println(io, join(vcat(["step", "time", "block", "element", "point"], columns), ","))
        for step in 1:num_steps, block_id in block_ids
            num_elements = size(first(values[block_id])[2], 1)
            for element in 1:num_elements, point in row_points
                source_point = point == 0 ? points[1] : point
                row = [string(step), format_value(times[step]), string(block_id), string(element), string(point)]
                for column in columns
                    matrix = get(values[block_id], (column, source_point), nothing)
                    push!(row, matrix === nothing ? "NaN" : format_value(matrix[element, step]))
                end
                println(io, join(row, ","))
            end
        end
    end

    close(exo)
    @printf(
        "%-40s %-6s %3d steps %d block(s) %d point(s)%s %3d columns, spread %.1e\n",
        name,
        layout,
        num_steps,
        length(block_ids),
        length(points),
        collapse ? " collapsed" : "          ",
        length(columns),
        max_spread
    )
    return nothing
end

function main(args)
    length(args) in (1, 2) || error("usage: extract.jl <LCM checkout> [output dir]")
    lcm_root = abspath(args[1])
    output_dir = length(args) == 2 ? abspath(args[2]) : @__DIR__
    mkpath(joinpath(output_dir, "materials"))
    mkpath(joinpath(output_dir, "meshes"))
    for (name, relative) in MESHES
        copy_named_mesh(joinpath(lcm_root, relative), joinpath(output_dir, "meshes", name))
    end
    for (name, layout, relative) in CASES
        cp(joinpath(lcm_root, MATERIAL_FILES[name]), joinpath(output_dir, "materials", name * ".yaml"); force=true)
        # The last commit that touched the gold file, not HEAD: a local branch of
        # the LCM checkout may carry unrelated changes.
        lcm_commit = readchomp(`git -C $lcm_root log -1 --format=%h -- $relative`)
        origin = "LCM gold file: " * relative
        extract_case(name, layout, joinpath(lcm_root, relative), output_dir, origin, lcm_commit)
    end
    return nothing
end

main(ARGS)
