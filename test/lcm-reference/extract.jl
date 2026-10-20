# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Extracts the reference ("gold") outputs of LCM's crystal plasticity and
# anisotropic hyperelasticity tests into CSV files that the Norma tests read.
# It is not part of the test suite; run it only to regenerate the fixtures:
#
#   julia --project=. test/lcm-reference/extract.jl <LCM checkout> [output dir] [--runs <dir>] [--every-step]
#
# LCM's polycrystal, RubiksCube, has no gold file. Its cases (RUN_CASES) are
# read from the outputs of LCM runs in the folder given with --runs, together
# with the material files those runs used; without --runs they are skipped.
# They keep only the steps and elements their tests read, unless
# --every-step is given.
#
# Each case gives two files:
#   <case>.csv        element variables, one row per (step, block, element,
#                     point); when every integration point of an element
#                     carries the same values, the points are collapsed into
#                     one row with point = 0
#   <case>-nodal.csv  nodal coordinates and variables, one row per
#                     (step, node); written only for the run cases whose
#                     test reads it
#
# The first line of each file names its source: the gold file and the last
# LCM commit that changed it, or the output of an LCM run and the LCM commit
# of the run.
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
#   :tet10  a 3 × 3 tensor at 4 points (LCM's composite tet at cubature
#           degree 2) is <field>_<a>_<b> with b = 01..12: the component
#           (row (b - 1) % 3 + 1, column a) at point (b - 1) ÷ 3 + 1; scalar
#           fields are <field>_<q>
#   :point  one material point (Material Point Simulator), no point suffix

using Exodus
using LinearAlgebra
using Printf
using Statistics

const CRYSTAL_PLASTICITY = "tests/LCM/CrystalPlasticity"
const MATERIAL_POINT = "tests/LCM/MaterialPointSimulator"

const CASES = [
    ("ElasticShear2d", :quad4, "$CRYSTAL_PLASTICITY/SingleElement2d/ElasticShear2d.gold.e"),
    ("PlasticShear2d", :quad4, "$CRYSTAL_PLASTICITY/SingleElement2d/PlasticShear2d.gold.e"),
    ("ElasticShear3d", :hex8, "$CRYSTAL_PLASTICITY/SingleElement3d/ElasticShear3d.gold.e"),
    ("ElasticShearRotated3d", :hex8, "$CRYSTAL_PLASTICITY/SingleElement3d/ElasticShearRotated3d.gold.e"),
    ("SingleSlip_Explicit", :hex8, "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlip/SingleSlip_Explicit.gold.e"),
    ("SingleSlip_Implicit", :hex8, "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlip/SingleSlip_Implicit.gold.e"),
    ("SingleSlipHard_Explicit", :hex8, "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipHard/SingleSlipHard_Explicit.gold.e"),
    ("SingleSlipHard_Implicit", :hex8, "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipHard/SingleSlipHard_Implicit.gold.e"),
    ("SingleSlipSaturation", :hex8, "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipSaturation/SingleSlipSaturation.gold.e"),
    ("MultiSlipPlaneHard_Implicit", :hex8, "$CRYSTAL_PLASTICITY/MultiSlipHard/MultiSlipPlaneHard_Implicit.gold.e"),
    (
        "MultiSlipPlaneHard_Implicit_Active_Sets",
        :hex8,
        "$CRYSTAL_PLASTICITY/MultiSlipHard/MultiSlipPlaneHard_Implicit_Active_Sets.gold.e",
    ),
    ("MinisolverStep_Newton", :hex8, "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_Newton.gold.e"),
    (
        "MinisolverStep_NewtonLineSearch",
        :hex8,
        "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_NewtonLineSearch.gold.e",
    ),
    ("MinisolverStep_TrustRegion", :hex8, "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_TrustRegion.gold.e"),
    (
        "MinisolverStep_ConjugateGradient",
        :hex8,
        "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_ConjugateGradient.gold.e",
    ),
    (
        "MinisolverStep_LineSearchRegularized",
        :hex8,
        "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_LineSearchRegularized.gold.e",
    ),
    ("ThermallyActivatedSlip", :hex8, "$CRYSTAL_PLASTICITY/ThermallyActivatedSlip/ThermallyActivatedSlip.gold.e"),
    ("MultiFamily", :hex8, "$CRYSTAL_PLASTICITY/MultiFamily/MultiFamily.gold.e"),
    (
        "DislocationDensityHardening",
        :hex8,
        "$CRYSTAL_PLASTICITY/DislocationDensityHardening/DislocationDensityHardening.gold.e",
    ),
    (
        "DislocationDensityTwoFamilies",
        :hex8,
        "$CRYSTAL_PLASTICITY/DislocationDensityTwoFamilies/DislocationDensityTwoFamilies.gold.e",
    ),
    (
        "QuadSlipDislocationDensityTraction",
        :hex8,
        "$CRYSTAL_PLASTICITY/QuadSlipDislocationDensityTraction/QuadSlipDislocationDensityTraction.gold.e",
    ),
    ("CP-uniaxial", :point, "$MATERIAL_POINT/CrystalPlasticity_MPS/CP-uniaxial.gold.exo"),
    ("AHD-uniaxial", :point, "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-uniaxial.gold.exo"),
    ("AHD-shear", :point, "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-shear.gold.exo"),
    ("AHD-hydrostatic", :point, "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-hydrostatic.gold.exo"),
]

# The LCM material file of each case, copied to materials/<case>.yaml. The
# Material Point Simulator cases keep their material in the input file.
const MATERIAL_FILES = Dict(
    "ElasticShear2d" => "$CRYSTAL_PLASTICITY/SingleElement2d/CrystalElasticity2d_Material.yaml",
    "PlasticShear2d" => "$CRYSTAL_PLASTICITY/SingleElement2d/CrystalPlasticity2d_Material.yaml",
    "ElasticShear3d" => "$CRYSTAL_PLASTICITY/SingleElement3d/CrystalElasticity3d_Material.yaml",
    "ElasticShearRotated3d" => "$CRYSTAL_PLASTICITY/SingleElement3d/CrystalElasticityRotated3d_Material.yaml",
    "SingleSlip_Explicit" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlip/SingleSlip_Explicit_Material.yaml",
    "SingleSlip_Implicit" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlip/SingleSlip_Implicit_Material.yaml",
    "SingleSlipHard_Explicit" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipHard/SingleSlipHard_Material_Explicit.yaml",
    "SingleSlipHard_Implicit" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipHard/SingleSlipHard_Material_Implicit.yaml",
    "SingleSlipSaturation" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlipSaturation/SingleSlipSaturation_Material.yaml",
    "MultiSlipPlaneHard_Implicit" => "$CRYSTAL_PLASTICITY/MultiSlipHard/MultiSlipPlaneHard_Material_Implicit.yaml",
    "MultiSlipPlaneHard_Implicit_Active_Sets" => "$CRYSTAL_PLASTICITY/MultiSlipHard/MultiSlipPlaneHard_Material_Implicit_Active_Sets.yaml",
    "MinisolverStep_Newton" => "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_Newton_Material.yaml",
    "MinisolverStep_NewtonLineSearch" => "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_NewtonLineSearch_Material.yaml",
    "MinisolverStep_TrustRegion" => "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_TrustRegion_Material.yaml",
    "MinisolverStep_ConjugateGradient" => "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_ConjugateGradient_Material.yaml",
    "MinisolverStep_LineSearchRegularized" => "$CRYSTAL_PLASTICITY/MinisolverStep/MinisolverStep_LineSearchRegularized_Material.yaml",
    "ThermallyActivatedSlip" => "$CRYSTAL_PLASTICITY/ThermallyActivatedSlip/ThermallyActivatedSlip_Material.yaml",
    "MultiFamily" => "$CRYSTAL_PLASTICITY/MultiFamily/MultiFamily_Material.yaml",
    "DislocationDensityHardening" => "$CRYSTAL_PLASTICITY/DislocationDensityHardening/DislocationDensityHardening_Material.yaml",
    "DislocationDensityTwoFamilies" => "$CRYSTAL_PLASTICITY/DislocationDensityTwoFamilies/DislocationDensityTwoFamilies_Material.yaml",
    "QuadSlipDislocationDensityTraction" => "$CRYSTAL_PLASTICITY/QuadSlipDislocationDensityTraction/QuadSlipDislocationDensity_Material.yaml",
    "CP-uniaxial" => "$MATERIAL_POINT/CrystalPlasticity_MPS/CP-uniaxial.yaml",
    "AHD-uniaxial" => "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-uniaxial.yaml",
    "AHD-shear" => "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-shear.yaml",
    "AHD-hydrostatic" => "$MATERIAL_POINT/AnisotropicHyperelasticDamage/AHD-hydrostatic.yaml",
)

# Meshes of the cases that the FEM tests run, copied to meshes/ with the
# block and node set names LCM gives them by default (block_<id>,
# nodelist_<id>), which the Norma inputs refer to.
const MESHES = Dict(
    "SingleSlip.g" => "$CRYSTAL_PLASTICITY/SingleSlip/SingleSlip.g",
    "MultiSlipPlaneHard.g" => "$CRYSTAL_PLASTICITY/MultiSlipHard/MultiSlipPlaneHard.g",
    "RubiksCube.g" => "$CRYSTAL_PLASTICITY/RubiksCube/RubiksCube.g",
)

# Every step, and every stride-th step with the last.
every_step(num_steps) = collect(1:num_steps)
stride(k) = num_steps -> sort!(unique(vcat(collect(1:k:num_steps), [num_steps])))

# Cases from LCM runs rather than gold files, both files in the --runs
# folder, the inputs of the runs in inputs/. Each keeps only what
# test/lcm-rubiks-cube.jl reads:
#   steps  the steps written, a function of the number of steps; a case that
#          keeps only some also gets <case>-times.csv, the time of every LCM
#          step, which a run must take
#   face   if a number, only the elements with a node on the face x = face
#   nodal  whether <case>-nodal.csv is written
# The cases are RubiksCube_Implicit.yaml as LCM ships it, of which the test
# reads the axial reaction at the last step; the same with its volume
# averaging switched off (Weighted Volume Average J and Volume Average
# Pressure false in every block), whose every fourth step includes those of
# the largest deviations; and the latter on RubiksCubeTet10.g
# (make_rubiks_cube_tet10.jl) with LCM's composite tet (Use Composite Tet 10
# in every block, cubature degree 2). Their deformation is not homogeneous,
# so no F_nodal is written.
const RUN_CASES = [
    (
        name="RubiksCube",
        layout=:hex8,
        output="RubiksCube_Implicit.e",
        material="RubiksCube_Implicit_Materials.yaml",
        steps=num_steps -> [num_steps],
        face=0.5,
        nodal=false,
    ),
    (
        name="RubiksCube_NoAverage",
        layout=:hex8,
        output="RubiksCube_Implicit_NoAvg.e",
        material="RubiksCube_Implicit_NoAvg_Materials.yaml",
        steps=stride(4),
        face=nothing,
        nodal=true,
    ),
    (
        name="RubiksCubeTet10",
        layout=:tet10,
        output="RubiksCubeTet10.e",
        material="RubiksCubeTet10_Materials.yaml",
        steps=stride(15),
        face=nothing,
        nodal=true,
    ),
]

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
# displacements (the Material Point Simulator cases, which write F).
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
    layout == :point && return String(name), 1
    stem, integers = trailing_integers(name)
    if stem in tensor_fields && length(integers) >= 2
        a, b = integers[end - 1], integers[end]
        prefix = length(integers) > 2 ? stem * "_" * join(integers[1:(end - 2)], '_') : stem
        if layout == :quad4
            flat = 2 * (b - 1) + (a - 1)
            component_2d = flat % 4 + 1
            row, col = (component_2d - 1) ÷ 2 + 1, (component_2d - 1) % 2 + 1
            return @sprintf("%s_%d", prefix, 3 * (row - 1) + col), flat ÷ 4 + 1
        elseif layout == :tet10
            row, col = (b - 1) % 3 + 1, a
            return @sprintf("%s_%d", prefix, 3 * (row - 1) + col), (b - 1) ÷ 3 + 1
        end
        return @sprintf("%s_%d", prefix, a), b
    end
    isempty(integers) && return String(name), 0
    head, _ = rsplit(name, '_'; limit=2)
    return String(head), integers[end]
end

# Fields whose names carry a tensor component before the point index, i.e.
# <field>_<a>_<b> with a = 1..9 and b = 1..8 (:hex8), a = 1..3 and b = 1..12
# (:tet10) or a = 1..2 and b = 1..8
# (:quad4). Per-slip-system scalars <field>_<s>_<q> also have two trailing
# integers, but in 2D their point index only runs to 4.
function tensor_fields_of(names::Vector{String}, layout::Symbol)
    layout == :point && return Set{String}()
    first_indices = Dict{String,Set{Int}}()
    second_indices = Dict{String,Set{Int}}()
    for name in names
        stem, integers = trailing_integers(name)
        length(integers) == 2 || continue
        push!(get!(first_indices, stem, Set{Int}()), integers[1])
        push!(get!(second_indices, stem, Set{Int}()), integers[2])
    end
    expected_first = if layout == :quad4
        Set(1:2)
    elseif layout == :tet10
        Set(1:3)
    else
        Set(1:9)
    end
    expected_second = layout == :tet10 ? Set(1:12) : Set(1:8)
    return Set(
        stem for
        stem in keys(first_indices) if first_indices[stem] == expected_first && second_indices[stem] == expected_second
    )
end

format_value(v::Real) = isfinite(v) ? string(Float64(v)) : "NaN"

# The elements of each block with a node on the face x = `face`.
function elements_on_face(exo, block_ids, face)
    coords = read_coordinates(exo)
    elements = Dict{Int,Vector{Int}}()
    for block_id in block_ids
        _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(exo, block_id)
        connectivity = reshape(
            Exodus.read_block_connectivity(exo, block_id, num_elements * num_nodes), num_nodes, num_elements
        )
        elements[Int(block_id)] = [
            e for e in 1:num_elements if any(abs.(coords[1, connectivity[:, e]] .- face) .< 1.0e-12)
        ]
    end
    return elements
end

function extract_case(
    name, layout, gold, output_dir, origin, lcm_commit; affine_fit=true, steps=every_step, face=nothing, nodal=true
)
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

    gradients = affine_fit ? nodal_deformation_gradients(exo, block_ids, num_steps) : nothing
    written_steps = steps(num_steps)
    kept_elements = face === nothing ? nothing : elements_on_face(exo, block_ids, face)
    if written_steps != 1:num_steps
        open(joinpath(output_dir, name * "-times.csv"), "w") do io
            println(io, "# ", origin, " (LCM ", lcm_commit, "); the time of every step")
            println(io, "step,time")
            return foreach(step -> println(io, step, ",", format_value(times[step])), 1:num_steps)
        end
    end
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
        if written_steps != 1:num_steps || kept_elements !== nothing
            print(
                io,
                "# kept: step",
                length(written_steps) == 1 ? " " : "s ",
                join(written_steps, ", "),
                " of ",
                num_steps,
            )
            println(io, kept_elements === nothing ? "" : "; the elements with a node on x = $face")
        end
        println(io, join(vcat(["step", "time", "block", "element", "point"], columns), ","))
        for step in written_steps, block_id in block_ids
            num_elements = size(first(values[block_id])[2], 1)
            elements = kept_elements === nothing ? (1:num_elements) : kept_elements[Int(block_id)]
            for element in elements, point in row_points
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

    nodal && write_nodal_file(exo, joinpath(output_dir, name * "-nodal.csv"), origin, lcm_commit, written_steps)
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

function write_nodal_file(exo, path, origin, lcm_commit, written_steps)
    times = read_times(exo)
    nodal_names = read_names(exo, NodalVariable)
    coords = read_coordinates(exo)
    num_dims, num_nodes = size(coords)
    nodal = Dict(
        n => Dict(s => Float64.(read_values(exo, NodalVariable, s, n)) for s in written_steps) for n in nodal_names
    )
    open(path, "w") do io
        println(io, "# ", origin, " (LCM ", lcm_commit, ")")
        println(io, join(vcat(["step", "time", "node"], ["x", "y", "z"][1:num_dims], nodal_names), ","))
        for step in written_steps, node in 1:num_nodes
            row = [string(step), format_value(times[step]), string(node)]
            append!(row, [format_value(coords[d, node]) for d in 1:num_dims])
            append!(row, [format_value(nodal[n][step][node]) for n in nodal_names])
            println(io, join(row, ","))
        end
    end
    return nothing
end

function main(args)
    args = copy(args)
    every = "--every-step" in args
    filter!(!=("--every-step"), args)
    runs_dir = nothing
    k = findfirst(==("--runs"), args)
    if k !== nothing
        k < length(args) || error("--runs needs a folder")
        runs_dir = abspath(args[k + 1])
        deleteat!(args, k:(k + 1))
    end
    length(args) in (1, 2) || error("usage: extract.jl <LCM checkout> [output dir] [--runs <dir>] [--every-step]")
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
        extract_case(name, layout, joinpath(lcm_root, relative), output_dir, origin, lcm_commit; nodal=false)
    end
    runs_dir === nothing && return nothing
    # The runs are of the checked-out LCM commit.
    lcm_commit = readchomp(`git -C $lcm_root rev-parse --short=7 HEAD`)
    for case in RUN_CASES
        cp(joinpath(runs_dir, case.material), joinpath(output_dir, "materials", case.name * ".yaml"); force=true)
        origin = "output " * case.output * " of an LCM run"
        gold = joinpath(runs_dir, case.output)
        if every
            extract_case(case.name, case.layout, gold, output_dir, origin, lcm_commit; affine_fit=false)
        else
            (; steps, face, nodal) = case
            extract_case(
                case.name, case.layout, gold, output_dir, origin, lcm_commit; affine_fit=false, steps, face, nodal
            )
        end
    end
    return nothing
end

main(ARGS)
