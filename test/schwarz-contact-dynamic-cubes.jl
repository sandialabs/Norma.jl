# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz Contact Implicit Cubes Tied" begin
    cp("../examples/contact/implicit-dynamic/friction-cubes/cubes.yaml", "cubes.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/friction-cubes/cube-1.yaml", "cube-1.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/friction-cubes/cube-2.yaml", "cube-2.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/friction-cubes/cube-1.g", "cube-1.g"; force=true)
    cp("../examples/contact/implicit-dynamic/friction-cubes/cube-2.g", "cube-2.g"; force=true)
    input_file = "cubes.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = -1.0e-06
    params["time step"] = 1e-6
    params["final time"] = 4e-6
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_1 = subsims[1].model
    model_2 = subsims[2].model

    rm("cubes.yaml"; force=true)
    rm("cube-1.yaml"; force=true)
    rm("cube-2.yaml"; force=true)
    rm("cube-1.g"; force=true)
    rm("cube-2.g"; force=true)
    rm("cube-1.e"; force=true)
    rm("cube-2.e"; force=true)

    # Enforce that the side set displacements are the same between the two blocks
    x1 = model_1.reference[1, :] + model_1.displacement[1, :]
    y1 = model_1.reference[2, :] + model_1.displacement[2, :]
    z1 = model_1.reference[3, :] + model_1.displacement[3, :]
    x2 = model_2.reference[1, :] + model_2.displacement[1, :]
    y2 = model_2.reference[2, :] + model_2.displacement[2, :]
    z2 = model_2.reference[3, :] + model_2.displacement[3, :]

    face_pairs = Dict(19 => 1, 20 => 2, 21 => 3, 22 => 4, 23 => 9, 24 => 10, 25 => 13, 26 => 14, 27 => 17)

    for (idx1, idx2) in face_pairs
        coordinate_1 = [x1[idx1], y1[idx1], z1[idx1]]
        coordinate_2 = [x2[idx2], y2[idx2], z2[idx2]]
        @test coordinate_1[1] ≈ coordinate_2[1] atol = 1e-6
        @test coordinate_1[2] ≈ coordinate_2[2] atol = 1e-6
        @test coordinate_1[3] ≈ coordinate_2[3] atol = 1e-6
    end
end

@testset "Schwarz Contact Explicit Cubes Tied" begin
    cp("../examples/contact/explicit-dynamic/friction-cubes/cubes.yaml", "cubes.yaml"; force=true)
    cp("../examples/contact/explicit-dynamic/friction-cubes/cube-1.yaml", "cube-1.yaml"; force=true)
    cp("../examples/contact/explicit-dynamic/friction-cubes/cube-2.yaml", "cube-2.yaml"; force=true)
    cp("../examples/contact/explicit-dynamic/friction-cubes/cube-1.g", "cube-1.g"; force=true)
    cp("../examples/contact/explicit-dynamic/friction-cubes/cube-2.g", "cube-2.g"; force=true)
    input_file = "cubes.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = -1.0e-06
    params["time step"] = 1e-6
    params["final time"] = 4e-6
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_1 = subsims[1].model
    model_2 = subsims[2].model

    rm("cubes.yaml"; force=true)
    rm("cube-1.yaml"; force=true)
    rm("cube-2.yaml"; force=true)
    rm("cube-1.g"; force=true)
    rm("cube-2.g"; force=true)
    rm("cube-1.e"; force=true)
    rm("cube-2.e"; force=true)

    # Enforce that the side set displacements are the same between the two blocks
    x1 = model_1.reference[1, :] + model_1.displacement[1, :]
    y1 = model_1.reference[2, :] + model_1.displacement[2, :]
    z1 = model_1.reference[3, :] + model_1.displacement[3, :]
    x2 = model_2.reference[1, :] + model_2.displacement[1, :]
    y2 = model_2.reference[2, :] + model_2.displacement[2, :]
    z2 = model_2.reference[3, :] + model_2.displacement[3, :]

    face_pairs = Dict(19 => 1, 20 => 2, 21 => 3, 22 => 4, 23 => 9, 24 => 10, 25 => 13, 26 => 14, 27 => 17)

    for (idx1, idx2) in face_pairs
        coordinate_1 = [x1[idx1], y1[idx1], z1[idx1]]
        coordinate_2 = [x2[idx2], y2[idx2], z2[idx2]]
        @test coordinate_1 ≈ coordinate_2 atol = 5e-3
    end
end

