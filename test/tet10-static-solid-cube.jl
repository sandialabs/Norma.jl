# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using YAML

@testset "Tet10 Static Solid Cube" begin
    cp("../examples/element-types/tet10/cube/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/element-types/tet10/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    min_disp_x = minimum(integrator.displacement[1:3:end])
    min_disp_y = minimum(integrator.displacement[2:3:end])
    max_disp_z = maximum(integrator.displacement[3:3:end])
    avg_stress = average_components(model.stress)
    @test min_disp_x ≈ -0.25 rtol = 1.0e-06
    @test min_disp_y ≈ -0.25 rtol = 1.0e-06
    @test max_disp_z ≈ 1.00 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

# The cube cut into two blocks at z = 0, with the 4-point rule on the lower
# block and the default 5-point rule on the upper one. The Exodus output
# declares the element variables of all blocks at once, and defines on each
# block those of its own points only: on the lower block, stress_zz_5 and
# von_mises_stress_5 do not exist, where with every variable defined on every
# block they read as zeros.
@testset "Tet10 Static Solid Cube Output Of Blocks With Different Rules" begin
    source = Exodus.ExodusDatabase("../examples/element-types/tet10/cube/cube.g", "r")
    X = Exodus.read_coordinates(source)
    block = only(Exodus.read_ids(source, Block))
    _, num_elements, num_element_nodes = Exodus.read_block_parameters(source, block)
    connectivity = reshape(
        Exodus.read_block_connectivity(source, block, num_elements * num_element_nodes), num_element_nodes, num_elements
    )
    node_sets = [
        (name, Exodus.read_node_set_nodes(source, id)) for
        (name, id) in zip(Exodus.read_names(source, NodeSet), Exodus.read_ids(source, NodeSet))
    ]
    close(source)
    lower = [sum(X[3, connectivity[1:4, e]]) < 0.0 for e in 1:num_elements]
    init = Exodus.Initialization{Int32}(
        Int32(3), Int32(size(X, 2)), Int32(num_elements), Int32(2), Int32(length(node_sets)), Int32(0)
    )
    rm("cube-blocks.g"; force=true)
    mesh = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}("cube-blocks.g", "w", init)
    Exodus.write_coordinates(mesh, X)
    for (id, (name, elements)) in enumerate((("lower", lower), ("upper", .!lower)))
        Exodus.write_block(mesh, id, "TETRA10", Matrix{Int32}(connectivity[:, elements]))
        Exodus.write_name(mesh, Block, id, name)
    end
    for (id, (name, nodes)) in enumerate(node_sets)
        node_set = Exodus.NodeSet(Int32(id), Int32.(nodes))
        Exodus.write_set(mesh, node_set)
        Exodus.write_name(mesh, node_set, name)
    end
    close(mesh)
    params = YAML.load_file("../examples/element-types/tet10/cube/cube.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cube-blocks.yaml"
    params["input mesh file"] = "cube-blocks.g"
    params["output mesh file"] = "cube-blocks.e"
    params["Exodus output interval"] = 1.0
    params["model"]["material"]["blocks"] = Dict{String,Any}("lower" => "elastic", "upper" => "elastic")
    params["model"]["num integration points"] = Dict{String,Any}("lower" => 4)
    params["time integrator"]["time step"] = 1.0
    Norma.run(params)
    output = Exodus.ExodusDatabase("cube-blocks.e", "r")
    try
        ids = Dict(zip(Exodus.read_names(output, Block), Exodus.read_ids(output, Block)))
        # Uniaxial stress of 1e9 along z, exact with both rules.
        stress(block, variable) = Exodus.read_values(output, ElementVariable, 2, ids[block], variable)
        @test maximum(abs.(stress("lower", "stress_zz_4") .- 1.0e9)) < 1.0e-6 * 1.0e9
        @test maximum(abs.(stress("upper", "stress_zz_5") .- 1.0e9)) < 1.0e-6 * 1.0e9
        for variable in ("stress_zz_5", "von_mises_stress_5")
            @test_throws Exodus.ExodusError stress("lower", variable)
        end
        @test length(stress("lower", "stored_energy")) == count(lower)
    finally
        close(output)
        rm("cube-blocks.g"; force=true)
        rm("cube-blocks.e"; force=true)
    end
end
