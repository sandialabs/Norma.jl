# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
function average_components(V::Vector{Float64})
    x = mean(V[1:3:end])
    y = mean(V[2:3:end])
    z = mean(V[3:3:end])
    return [x y z]
end

function maximum_components(V::Vector{Float64})
    x = maximum(V[1:3:end])
    y = maximum(V[2:3:end])
    z = maximum(V[3:3:end])
    return [x y z]
end

function minimum_components(V::Vector{Float64})
    x = minimum(V[1:3:end])
    y = minimum(V[2:3:end])
    z = minimum(V[3:3:end])
    return [x y z]
end

function average_components(stress::Vector)
    xx = 0.0
    yy = 0.0
    zz = 0.0
    yz = 0.0
    xz = 0.0
    xy = 0.0
    num_stress = 0
    for block_index in 1:length(stress)
        for block_element_index in 1:length(stress[block_index])
            for point in 1:length(stress[block_index][block_element_index])
                xx += stress[block_index][block_element_index][point][1]
                yy += stress[block_index][block_element_index][point][2]
                zz += stress[block_index][block_element_index][point][3]
                yz += stress[block_index][block_element_index][point][4]
                xz += stress[block_index][block_element_index][point][5]
                xy += stress[block_index][block_element_index][point][6]
                num_stress += 1
            end
        end
    end
    return [xx yy zz yz xz xy] ./ num_stress
end

using DelimitedFiles
using StaticArrays

# Rows of an LCM reference case (test/lcm-reference/README.md), each a Dict
# from column name to value, in file order.
function lcm_reference(case::AbstractString)
    path = joinpath(@__DIR__, "lcm-reference", case * ".csv")
    data, header = readdlm(path, ',', Float64, '\n'; header=true, comments=true)
    names = String.(vec(header))
    return [Dict(names[j] => data[r, j] for j in eachindex(names)) for r in axes(data, 1)]
end

# The rows of an LCM reference case at its last step.
function lcm_final_rows(case::AbstractString)
    rows = lcm_reference(case)
    last_step = maximum(row["step"] for row in rows)
    return filter(row -> row["step"] == last_step, rows)
end

# Row-major 3 × 3 tensor `field` of an LCM reference row; components missing
# from the file (the out-of-plane ones of the 2D cases) are NaN.
function lcm_tensor(row::AbstractDict, field::AbstractString)
    return SMatrix{3,3,Float64,9}(get(row, "$(field)_$(3 * (i - 1) + j)", NaN) for i in 1:3, j in 1:3)
end

# Cauchy stress of a cubic St. Venant-Kirchhoff crystal with orientation Q,
# written out index by index, independently of anisotropy.jl.
function cubic_svk_cauchy(C11, C12, C44, Q, F)
    δ(a, b) = a == b ? 1.0 : 0.0
    CL = [
        C12 * δ(i, j) * δ(k, l) +
        C44 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) +
        (C11 - C12 - 2C44) * δ(i, j) * δ(j, k) * δ(k, l) for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    ]
    C = [
        sum(Q[i, p] * Q[j, q] * Q[k, r] * Q[l, s] * CL[p, q, r, s] for p in 1:3, q in 1:3, r in 1:3, s in 1:3) for
        i in 1:3, j in 1:3, k in 1:3, l in 1:3
    ]
    E = 0.5 * (F' * F - I)
    S = [sum(C[i, j, k, l] * E[k, l] for k in 1:3, l in 1:3) for i in 1:3, j in 1:3]
    return F * S * F' / det(F)
end

using Exodus
using Symbolics
@variables t, x, y, z

function create_force(expression::String, mesh::ExodusDatabase, side_set_id::Int64, time::Float64)
    force_num = eval(Meta.parse(expression))
    force_fun = eval(build_function(force_num, [t, x, y, z]; expression=Val(false)))
    coords = read_coordinates(mesh)
    num_nodes = size(coords)[2]
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    local_from_global_map = Norma.get_side_set_local_from_global_map(mesh, side_set_id)
    num_nodes = length(local_from_global_map)
    force = zeros(num_nodes)
    ss_node_index = 1
    for side in num_nodes_sides
        side_nodes = side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = coords[:, side_nodes]
        nodal_force_component = Norma.get_side_set_nodal_forces(side_coordinates, force_fun, time)
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        force[local_indices] += nodal_force_component
        ss_node_index += side
    end
    return force
end

function create_displacement(expression::String, mesh::ExodusDatabase, side_set_id::Int64, time::Float64)
    disp_num = eval(Meta.parse(expression))
    disp_fun = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
    coords = read_coordinates(mesh)
    num_nodes = size(coords)[2]
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    local_from_global_map = Norma.get_side_set_local_from_global_map(mesh, side_set_id)
    num_nodes = length(local_from_global_map)
    displacements = zeros(num_nodes)
    ss_node_index = 1
    for side in num_nodes_sides
        side_nodes = side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        for side_node in side_nodes
            nc = coords[:, side_node]
            disp_val = Float64(disp_fun((time, nc[1], nc[2], nc[3])))
            local_index = get.(Ref(local_from_global_map), side_node, 0)
            displacements[local_index] = disp_val
        end
        ss_node_index += side
    end
    return displacements
end

# Polycrystals (polycrystal.jl and periodic.jl).
using Random
using Statistics

const POLYCRYSTAL_EXAMPLES = "../examples/crystal-plasticity"

# Copper constants of LCM's SingleElement tests, and the FCC crystal of
# LCM's RubiksCube.
function copper_elastic()
    return Dict{String,Any}(
        "model" => "anisotropic linear elastic",
        "symmetry" => "cubic",
        "C11" => 168.4,
        "C12" => 121.4,
        "C44" => 75.4,
        "density" => 1.0,
    )
end
function rubiks_copper()
    return Dict{String,Any}(
        "model" => "crystal plasticity",
        "symmetry" => "cubic",
        "C11" => 204600.0,
        "C12" => 137700.0,
        "C44" => 126200.0,
        "density" => 1.0,
        "integration scheme" => "implicit",
        "slip families" => [
            Dict(
                "slip systems" => "fcc",
                "flow rule" => Dict("type" => "power law", "reference slip rate" => 1.0, "rate exponent" => 20.0),
                "hardening law" => Dict(
                    "type" => "linear minus recovery",
                    "hardening modulus" => 355.0,
                    "recovery modulus" => 2.9,
                    "initial hardening state" => 122.0,
                ),
            ),
        ],
    )
end
mesh_orientation() = Dict("mesh variables" => ["phi1", "Phi", "phi2"])

# The input of a quasi-static Norma run of `mesh` with `materials` (block name
# => material parameters) and `boundary_conditions`, the dictionary of the
# input file or a list of Dirichlet conditions (node set, component,
# function). With `adaptive`, a step whose local or global solve fails is
# retried with half the time step, down to a fiftieth, as LCM's RubiksCube
# does. The Exodus output is written every `exodus_interval` if it is
# positive. `predictor` is that of the quasi-static time integrator.
function polycrystal_parameters(
    mesh,
    materials,
    boundary_conditions;
    final_time=1.0,
    time_step=1.0,
    linear_solver="cg",
    adaptive=false,
    exodus_interval=0.0,
    predictor="constant",
)
    if boundary_conditions isa AbstractVector
        boundary_conditions = Dict(
            "Dirichlet" => [
                Dict{String,Any}("node set" => n, "component" => c, "function" => f) for
                (n, c, f) in boundary_conditions
            ],
        )
    end
    blocks = Dict{String,Any}()
    material = Dict{String,Any}("blocks" => blocks)
    for (block, params) in materials
        label = "material_" * block
        blocks[block] = label
        material[label] = params
    end
    name = splitext(mesh)[1]
    params = Norma.Parameters(
        "name" => name * ".yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => name * ".e",
        "Exodus output interval" => exodus_interval,
        "CSV output interval" => 0,
        "model" => Dict("type" => "solid mechanics", "material" => material),
        "time integrator" => merge(
            Dict{String,Any}(
                "type" => "quasi static",
                "initial time" => 0.0,
                "final time" => final_time,
                "time step" => time_step,
                "predictor" => predictor,
            ),
            if adaptive
                Dict{String,Any}(
                    "minimum time step" => time_step / 50,
                    "maximum time step" => time_step,
                    "decrease factor" => 0.5,
                    "increase factor" => 2.0,
                )
            else
                Dict{String,Any}()
            end,
        ),
        "boundary conditions" => boundary_conditions,
        "solver" => Dict(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => linear_solver,
            "minimum iterations" => 1,
            "maximum iterations" => 20,
            "absolute tolerance" => 1.0e-10,
            "relative tolerance" => 1.0e-12,
            "linear solver absolute tolerance" => 1.0e-14,
            "linear solver relative tolerance" => 1.0e-12,
        ),
    )
    return params
end

# The run of polycrystal_parameters, whose Exodus output is kept only if it
# was written.
function run_polycrystal(mesh, materials, boundary_conditions; exodus_interval=0.0, options...)
    sim = Norma.run(
        polycrystal_parameters(mesh, materials, boundary_conditions; exodus_interval=exodus_interval, options...)
    )
    exodus_interval > 0.0 || rm(splitext(mesh)[1] * ".e"; force=true)
    return sim
end

# Volume of every element, block by block.
function element_volumes(model)
    volumes = Vector{Float64}[]
    for block in model.blocks
        block_volumes = zeros(block.num_elements)
        for e in 1:block.num_elements
            X = Norma.gather_nodal(model.reference, view(block.connectivity, :, e), block.N)
            for p in 1:block.num_points
                block_volumes[e] += det(block.dN[:, :, p] * X') * block.weights[p]
            end
        end
        push!(volumes, block_volumes)
    end
    return volumes
end

# Fourth-order compliance of a Voigt stiffness with engineering shears.
function compliance_tensor(C6::SMatrix{6,6,Float64,36})
    S6 = inv(C6)
    factor(a) = a <= 3 ? 1.0 : 0.5
    S = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        a, b = Norma.voigt_index(i, j), Norma.voigt_index(k, l)
        S[i, j, k, l] = S6[a, b] * factor(a) * factor(b)
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(S)
end

function compliance_to_voigt(S::SArray{Tuple{3,3,3,3},Float64,4,81})
    factor(a) = a <= 3 ? 1.0 : 2.0
    return SMatrix{6,6,Float64,36}(
        S[Norma.VOIGT_PAIRS[a]..., Norma.VOIGT_PAIRS[b]...] * factor(a) * factor(b) for a in 1:6, b in 1:6
    )
end

# Voigt and Reuss averages, as Voigt matrices, of the lattice elasticity CC
# rotated into each element, weighted by volume.
function voigt_reuss(model, CC)
    volumes = element_volumes(model)
    total = sum(sum, volumes)
    S = compliance_tensor(Norma.tensor_to_voigt(CC))
    C_voigt = zero(SArray{Tuple{3,3,3,3},Float64,4,81})
    S_reuss = zero(SArray{Tuple{3,3,3,3},Float64,4,81})
    for (b, block_volumes) in enumerate(volumes), (e, v) in enumerate(block_volumes)
        Q = model.lattice_orientation[b][e]
        C_voigt += (v / total) * Norma.rotate_tensor4(Q, CC)
        S_reuss += (v / total) * Norma.rotate_tensor4(Q, S)
    end
    return Norma.tensor_to_voigt(C_voigt), compliance_to_voigt(S_reuss), total
end

# Axial reaction on the nodes of a node set with a prescribed x displacement.
function axial_reaction(model, node_set)
    bc = only(
        bc for bc in model.boundary_conditions if
        bc isa Norma.SolidMechanicsDirichletBoundaryCondition && bc.name == node_set && bc.offset == 1
    )
    return sum(model.internal_force[3 * (n - 1) + 1] for n in bc.node_set_node_indices)
end

function mesh_block_names(mesh)
    database = Exodus.ExodusDatabase(mesh, "r")
    names = Exodus.read_names(database, Block)
    close(database)
    return names
end

const FACE_SETS = ("nsx-", "nsx+", "nsy-", "nsy+", "nsz-", "nsz+")
# Uniaxial stress along x on rollers: x-, y- and z- faces fixed normally.
function uniaxial_rollers(displacement)
    return [("nsx-", "x", "0.0"), ("nsy-", "y", "0.0"), ("nsz-", "z", "0.0"), ("nsx+", "x", displacement)]
end
# Uniaxial strain ε along x imposed on the whole boundary: u = (ε x, 0, 0).
function affine_uniaxial_strain(ε)
    return vcat([(s, "x", "$ε * x") for s in FACE_SETS], [(s, c, "0.0") for s in FACE_SETS for c in ("y", "z")])
end

using YAML

# The abort message of a restart of the cuboid restart example
# (examples/ahead/single/cuboid/dynamic-restart, index 2) with `material` on
# its one block, or nothing when the restart is accepted. With `schwarz`, the
# restart check of two subdomains that use the material, which runs before
# any model is built.
function restart_abort_message(material::Norma.Parameters; schwarz::Bool=false)
    folder = "../examples/ahead/single/cuboid/dynamic-restart"
    cp(joinpath(folder, "cuboid-restart-in.e"), "material-restart-in.e"; force=true)
    blocks = Norma.Parameters("blocks" => Norma.Parameters("cuboid" => "solid"), "solid" => material)
    Norma.NORMA_TEST_MODE[] = true
    try
        if schwarz
            YAML.write_file(
                "material-restart-domain.yaml",
                Dict(
                    "input mesh file" => "material-restart-in.e",
                    "model" => Dict("type" => "solid mechanics", "material" => blocks),
                ),
            )
            params = Norma.Parameters(
                "domains" => ["material-restart-domain.yaml", "material-restart-domain.yaml"],
                "restart" => Norma.Parameters("index" => 2),
                "final time" => 1.0,
            )
            Norma.process_multidomain_restart!(params)
        else
            params = YAML.load_file(joinpath(folder, "cuboid.yaml"); dicttype=Norma.Parameters)
            params["name"] = "material-restart"
            params["input mesh file"] = "material-restart-in.e"
            params["output mesh file"] = "material-restart-out.e"
            params["model"]["material"] = blocks
            sim = Norma.create_simulation(params)
            Exodus.close(sim.params["input_mesh"])
            Exodus.close(sim.params["output_mesh"])
        end
        return nothing
    catch error
        error isa Norma.NormaAbortException || rethrow()
        return error.msg
    finally
        Norma.NORMA_TEST_MODE[] = false
        for file in ("material-restart-in.e", "material-restart-out.e", "material-restart-domain.yaml")
            rm(file; force=true)
        end
    end
end
