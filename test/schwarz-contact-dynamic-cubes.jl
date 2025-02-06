# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "schwarz-contact-inclined-explicit-cubes" begin
    
    model_fine = nothing
    model_coarse = nothing

    angles = [ 0.0, 22.5, 45, 90 ]
    for (i, angle_deg) in enumerate(angles)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cubes-test$i.yaml", "cubes-test$i.yaml", force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-1.yaml", "cube-test$i-1.yaml", force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-2.yaml", "cube-test$i-2.yaml", force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-1.g", "cube-test$i-1.g", force=true)
        cp("../examples/contact/explicit-dynamic/inclined-cubes/cube-test$i-2.g", "cube-test$i-2.g", force=true)
        input_file = "cubes-test$i.yaml"
        params = YAML.load_file(input_file; dicttype=Dict{String,Any})
        params["initial time"] = -1.0e-06
        params["final time"] = 1.0e-5
        sim = Norma.run(params, input_file)
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

        if ( i == 1 )
            model_fine = model_fine_temp
            model_coarse = model_coarse_temp
        else
            # Rotate these displacements about z
            angle = angle_deg * π / 180
            c = cos(angle)
            s = sin(angle)
            local_rotation_matrix = [ c s 0; -s c 0 ; 0 0 1]
            model_fine_rotated = zeros((3,68))
            model_coarse_rotated = zeros((3,27))
            for i in range(1,68)
                base = (i-1)*(3) + 1
                model_fine_rotated[:, i] = local_rotation_matrix * model_fine_temp[:, i] 
            end
            for i in range(1,27)
                base = (i-1)*(3) + 1
                model_coarse_rotated[:, i] = local_rotation_matrix * model_coarse_temp[:, i] 
            end

            @test model_fine_rotated ≈ model_fine rtol=1e-5
            @test model_coarse_rotated ≈ model_coarse rtol=1e-5
        end
    end

end

@testset "schwarz-contact-inclined-implicit-cubes" begin
    
    model_fine = nothing
    model_coarse = nothing

    angles = [ 0.0, 22.5, 45, 90 ]
    for (i, angle_deg) in enumerate(angles)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cubes-test$i.yaml", "cubes-test$i.yaml", force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-1.yaml", "cube-test$i-1.yaml", force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-2.yaml", "cube-test$i-2.yaml", force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-1.g", "cube-test$i-1.g", force=true)
        cp("../examples/contact/implicit-dynamic/inclined-cubes/cube-test$i-2.g", "cube-test$i-2.g", force=true)
        input_file = "cubes-test$i.yaml"
        params = YAML.load_file(input_file; dicttype=Dict{String,Any})
        params["initial time"] = -1.0e-06
        params["final time"] = 1.0e-3
        sim = Norma.run(params, input_file)
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

        if ( i == 1 )
            model_fine = model_fine_temp
            model_coarse = model_coarse_temp
        else
            # Rotate these displacements about z
            angle = angle_deg * π / 180
            c = cos(angle)
            s = sin(angle)
            local_rotation_matrix = [ c s 0; -s c 0 ; 0 0 1]
            model_fine_rotated = zeros((3,68))
            model_coarse_rotated = zeros((3,27))
            for i in range(1,68)
                base = (i-1)*(3) + 1
                model_fine_rotated[:, i] = local_rotation_matrix * model_fine_temp[:, i] 
            end
            for i in range(1,27)
                base = (i-1)*(3) + 1
                model_coarse_rotated[:, i] = local_rotation_matrix * model_coarse_temp[:, i] 
            end

            @test model_fine_rotated ≈ model_fine rtol=2e-4
            @test model_coarse_rotated ≈ model_coarse rtol=2e-4
        end
    end

end