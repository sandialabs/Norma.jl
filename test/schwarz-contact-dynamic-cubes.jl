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

    rm("cubes.yaml")
    rm("cube-1.yaml")
    rm("cube-2.yaml")
    rm("cube-1.g")
    rm("cube-2.g")
    rm("cube-1.e")
    rm("cube-2.e")

    # Enforce that the side set displacements are the same between the two blocks
    x1 = model_1.current[1, :]
    y1 = model_1.current[2, :]
    z1 = model_1.current[3, :]
    x2 = model_2.current[1, :]
    y2 = model_2.current[2, :]
    z2 = model_2.current[3, :]

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

    rm("cubes.yaml")
    rm("cube-1.yaml")
    rm("cube-2.yaml")
    rm("cube-1.g")
    rm("cube-2.g")
    rm("cube-1.e")
    rm("cube-2.e")

    # Enforce that the side set displacements are the same between the two blocks
    x1 = model_1.current[1, :]
    y1 = model_1.current[2, :]
    z1 = model_1.current[3, :]
    x2 = model_2.current[1, :]
    y2 = model_2.current[2, :]
    z2 = model_2.current[3, :]

    face_pairs = Dict(19 => 1, 20 => 2, 21 => 3, 22 => 4, 23 => 9, 24 => 10, 25 => 13, 26 => 14, 27 => 17)

    for (idx1, idx2) in face_pairs
        coordinate_1 = [x1[idx1], y1[idx1], z1[idx1]]
        coordinate_2 = [x2[idx2], y2[idx2], z2[idx2]]
        @test coordinate_1 ≈ coordinate_2 atol = 5e-3
    end
end

@testset "Schwarz Contact Inclined Explicit Cubes" begin
    model_fine = nothing
    model_coarse = nothing

    angles = [0.0, 22.5, 45, 90]
    for (i, angle_deg) in enumerate(angles)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cubes-test$i.yaml", "cubes-test$i.yaml"; force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-1.yaml", "cube-test$i-1.yaml"; force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-2.yaml", "cube-test$i-2.yaml"; force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-1.g", "cube-test$i-1.g"; force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-2.g", "cube-test$i-2.g"; force=true)
        input_file = "cubes-test$i.yaml"
        params = YAML.load_file(input_file; dicttype=Norma.Parameters)
        params["initial time"] = -1.0e-06
        params["final time"] = 1.0e-5
        params["name"] = input_file
        sim = Norma.run(params)
        subsim_temp = sim.subsims
        model_fine_temp = subsim_temp[1].model.current
        model_coarse_temp = subsim_temp[2].model.current

        rm("cubes-test$i.yaml")
        rm("cube-test$i-1.yaml")
        rm("cube-test$i-2.yaml")
        rm("cube-test$i-1.g")
        rm("cube-test$i-2.g")
        rm("cube-test$i-1.e")
        rm("cube-test$i-2.e")

        if (i == 1)
            model_fine = model_fine_temp
            model_coarse = model_coarse_temp
        else
            # Rotate these displacements about z
            angle = angle_deg * π / 180
            c = cos(angle)
            s = sin(angle)
            local_rotation_matrix = [c s 0; -s c 0; 0 0 1]
            model_fine_rotated = zeros((3, 68))
            model_coarse_rotated = zeros((3, 27))
            for i in range(1, 68)
                base = (i - 1) * (3) + 1
                model_fine_rotated[:, i] = local_rotation_matrix * model_fine_temp[:, i]
            end
            for i in range(1, 27)
                base = (i - 1) * (3) + 1
                model_coarse_rotated[:, i] = local_rotation_matrix * model_coarse_temp[:, i]
            end

            @test model_fine_rotated ≈ model_fine rtol = 1e-5
            @test model_coarse_rotated ≈ model_coarse rtol = 1e-5
        end
    end
end

@testset "Schwarz Contact Inclined Implicit Cubes" begin
    model_fine = nothing
    model_coarse = nothing

    angles = [0.0, 22.5, 45, 90]
    for (i, angle_deg) in enumerate(angles)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cubes-test$i.yaml", "cubes-test$i.yaml"; force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-1.yaml", "cube-test$i-1.yaml"; force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-2.yaml", "cube-test$i-2.yaml"; force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-1.g", "cube-test$i-1.g"; force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-2.g", "cube-test$i-2.g"; force=true)
        input_file = "cubes-test$i.yaml"
        params = YAML.load_file(input_file; dicttype=Norma.Parameters)
        params["initial time"] = -1.0e-06
        params["time step"] = 1e-6
        params["final time"] = 4e-6
        params["name"] = input_file
        sim = Norma.run(params)
        subsim_temp = sim.subsims
        model_fine_temp = subsim_temp[1].model.current
        model_coarse_temp = subsim_temp[2].model.current

        rm("cubes-test$i.yaml")
        rm("cube-test$i-1.yaml")
        rm("cube-test$i-2.yaml")
        rm("cube-test$i-1.g")
        rm("cube-test$i-2.g")
        rm("cube-test$i-1.e")
        rm("cube-test$i-2.e")

        if (i == 1)
            model_fine = model_fine_temp
            model_coarse = model_coarse_temp
        else
            # Rotate these displacements about z
            angle = angle_deg * π / 180
            c = cos(angle)
            s = sin(angle)
            local_rotation_matrix = [c s 0; -s c 0; 0 0 1]
            model_fine_rotated = zeros((3, 68))
            model_coarse_rotated = zeros((3, 27))
            for i in range(1, 68)
                base = (i - 1) * (3) + 1
                model_fine_rotated[:, i] = local_rotation_matrix * model_fine_temp[:, i]
            end
            for i in range(1, 27)
                base = (i - 1) * (3) + 1
                model_coarse_rotated[:, i] = local_rotation_matrix * model_coarse_temp[:, i]
            end

            @test model_fine_rotated ≈ model_fine rtol = 1e-5
            @test model_coarse_rotated ≈ model_coarse rtol = 1e-5
        end
    end
end
