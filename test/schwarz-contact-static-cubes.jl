# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz Contact Static Cubes Hex8" begin
    cp("../examples/contact/static/cubes/cubes.yaml", "cubes.yaml"; force=true)
    cp("../examples/contact/static/cubes/cube-1.yaml", "cube-1.yaml"; force=true)
    cp("../examples/contact/static/cubes/cube-2.yaml", "cube-2.yaml"; force=true)
    cp("../examples/contact/static/cubes/cube-1.g", "cube-1.g"; force=true)
    cp("../examples/contact/static/cubes/cube-2.g", "cube-2.g"; force=true)
    input_file = "cubes.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = -2.0
    params["final time"] = 0.0
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cubes.yaml")
    rm("cube-1.yaml")
    rm("cube-2.yaml")
    rm("cube-1.g")
    rm("cube-2.g")
    rm("cube-1.e")
    rm("cube-2.e")
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    max_disp_y_fine = maximum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    max_disp_x_coarse = maximum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    max_disp_y_coarse = maximum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    max_disp_z_coarse = maximum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ 0.0 atol = 2.0e-05
    @test max_disp_y_fine ≈ 0.025 rtol = 2.0e-04
    @test max_disp_z_fine ≈ 0.025 rtol = 2.0e-04
    @test max_disp_x_coarse ≈ 0.0 atol = 2.0e-05
    @test max_disp_y_coarse ≈ 0.025 rtol = 2.0e-04
    @test max_disp_z_coarse ≈ 0.025 rtol = 2.0e-04
    @test avg_stress_fine[1] ≈ -1.0e08 rtol = 2.0e-04
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_fine[3] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_coarse[1] ≈ -1.0e08 rtol = 2.0e-04
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_coarse[3] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-04
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-04
end

@testset "Schwarz Contact Static Cubes Tied" begin
    cp("../examples/contact/static/friction-cubes/cubes.yaml", "cubes.yaml"; force=true)
    cp("../examples/contact/static/friction-cubes/cube-1.yaml", "cube-1.yaml"; force=true)
    cp("../examples/contact/static/friction-cubes/cube-2.yaml", "cube-2.yaml"; force=true)
    cp("../examples/contact/static/friction-cubes/cube-1.g", "cube-1.g"; force=true)
    cp("../examples/contact/static/friction-cubes/cube-2.g", "cube-2.g"; force=true)
    input_file = "cubes.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = -2.0
    params["final time"] = 0.0
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

@testset "Schwarz Contact Inclined Static Cubes" begin
    model_fine = nothing
    model_coarse = nothing

    angles = [0.0, 22.5, 45, 90]
    for (i, angle_deg) in enumerate(angles)
        cp("../examples/contact/static/inclined-cubes/cubes-test$i.yaml", "cubes-test$i.yaml"; force=true)
        cp("../examples/contact/static/inclined-cubes/cube-test$i-1.yaml", "cube-test$i-1.yaml"; force=true)
        cp("../examples/contact/static/inclined-cubes/cube-test$i-2.yaml", "cube-test$i-2.yaml"; force=true)
        cp("../examples/contact/static/inclined-cubes/cube-test$i-1.g", "cube-test$i-1.g"; force=true)
        cp("../examples/contact/static/inclined-cubes/cube-test$i-2.g", "cube-test$i-2.g"; force=true)
        input_file = "cubes-test$i.yaml"
        params = YAML.load_file(input_file; dicttype=Norma.Parameters)
        params["initial time"] = -2.0
        params["final time"] = 0.0
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
            model_fine_rotated = zeros((3, 64))
            model_coarse_rotated = zeros((3, 27))
            for i in range(1, 64)
                base = (i - 1) * (3) + 1
                model_fine_rotated[:, i] = local_rotation_matrix * model_fine_temp[:, i]
            end
            for i in range(1, 27)
                base = (i - 1) * (3) + 1
                model_coarse_rotated[:, i] = local_rotation_matrix * model_coarse_temp[:, i]
            end

            @test model_fine_rotated ≈ model_fine rtol = 2e-4
            @test model_coarse_rotated ≈ model_coarse rtol = 2e-4
        end
    end
end
