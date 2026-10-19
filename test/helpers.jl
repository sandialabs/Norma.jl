# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
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

function average_components(stress::Vector{Vector{Vector{Vector{Float64}}}})
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
