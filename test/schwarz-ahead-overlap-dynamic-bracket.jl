# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Bracket TET4-TET4 Dissipative Newmark" begin
    cp("../examples/ahead/overlap/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-1.yaml", "bracket-1.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-2.yaml", "bracket-2.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-1.g", "../bracket-1.g"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-2.g", "../bracket-2.g"; force=true)
    input_file = "bracket.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-6
    params["final time"] = 5.0e-5
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_bracket1 = subsims[1].model
    model_bracket2 = subsims[2].model

    rm("bracket.yaml")
    rm("bracket-1.yaml")
    rm("bracket-2.yaml")
    rm("../bracket-1.g")
    rm("../bracket-2.g")
    rm("bracket-1.e")
    rm("bracket-2.e")

    min_disp_x_bracket1 = minimum(model_bracket1.current[1, :] - model_bracket1.reference[1, :])
    min_disp_y_bracket1 = minimum(model_bracket1.current[2, :] - model_bracket1.reference[2, :])
    max_disp_z_bracket1 = maximum(model_bracket1.current[3, :] - model_bracket1.reference[3, :])
    min_disp_x_bracket2 = minimum(model_bracket2.current[1, :] - model_bracket2.reference[1, :])
    min_disp_y_bracket2 = minimum(model_bracket2.current[2, :] - model_bracket2.reference[2, :])
    max_disp_z_bracket2 = maximum(model_bracket2.current[3, :] - model_bracket2.reference[3, :])
    avg_stress_bracket1 = average_components(model_bracket1.stress)
    avg_stress_bracket2 = average_components(model_bracket2.stress)

    println("min_disp_x_bracket1 = ", min_disp_x_bracket1, "\n")
    println("min_disp_y_bracket1 = ", min_disp_y_bracket1, "\n")
    println("max_disp_z_bracket1 = ", max_disp_z_bracket1, "\n")
    println("min_disp_x_bracket2 = ", min_disp_x_bracket2, "\n")
    println("min_disp_y_bracket2 = ", min_disp_y_bracket2, "\n")
    println("max_disp_z_bracket2 = ", max_disp_z_bracket2, "\n")
    println("avg_stress_bracket1 = ", avg_stress_bracket1, "\n")
    println("avg_stress_bracket2 = ", avg_stress_bracket2, "\n")

    @test min_disp_x_bracket1 ≈ -3.239534728147805e-5 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -5.390040635191995e-5 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.0001136313123976402 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.00011223576082500242 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -5.617541938109094e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.0008779212667680046 atol = 1e-12
    @test avg_stress_bracket1 ≈ [-72712.82232539689 -90192.50915687135 4480.231591848708 27863.45785162852 -1.0987341235848952e6 -35974.56740305576] atol = 1e-12
    @test avg_stress_bracket2 ≈ [1.4716110552352623e6 -3498.5059655966143 -380718.1494447633 -166878.95669528542 -4.988667087022611e6 -145718.2020229697] atol = 1e-12
    @test sim.controller.schwarz_iters ≈ [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3] atol = 0

end

@testset "Schwarz AHeaD Overlap Dynamic Bracket HEX8-TET4 Dissipative Newmark" begin
    cp("../examples/ahead/overlap/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-1.yaml", "bracket-1.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-2.yaml", "bracket-2.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-1-hex.g", "../bracket-1.g"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-2-tet.g", "../bracket-2.g"; force=true)
    input_file = "bracket.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-6
    params["final time"] = 5.0e-5
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_bracket1 = subsims[1].model
    model_bracket2 = subsims[2].model

    rm("bracket.yaml")
    rm("bracket-1.yaml")
    rm("bracket-2.yaml")
    rm("../bracket-1.g")
    rm("../bracket-2.g")
    rm("bracket-1.e")
    rm("bracket-2.e")

    min_disp_x_bracket1 = minimum(model_bracket1.current[1, :] - model_bracket1.reference[1, :])
    min_disp_y_bracket1 = minimum(model_bracket1.current[2, :] - model_bracket1.reference[2, :])
    max_disp_z_bracket1 = maximum(model_bracket1.current[3, :] - model_bracket1.reference[3, :])
    min_disp_x_bracket2 = minimum(model_bracket2.current[1, :] - model_bracket2.reference[1, :])
    min_disp_y_bracket2 = minimum(model_bracket2.current[2, :] - model_bracket2.reference[2, :])
    max_disp_z_bracket2 = maximum(model_bracket2.current[3, :] - model_bracket2.reference[3, :])
    avg_stress_bracket1 = average_components(model_bracket1.stress)
    avg_stress_bracket2 = average_components(model_bracket2.stress)

    println("min_disp_x_bracket1 = ", min_disp_x_bracket1, "\n")
    println("min_disp_y_bracket1 = ", min_disp_y_bracket1, "\n")
    println("max_disp_z_bracket1 = ", max_disp_z_bracket1, "\n")
    println("min_disp_x_bracket2 = ", min_disp_x_bracket2, "\n")
    println("min_disp_y_bracket2 = ", min_disp_y_bracket2, "\n")
    println("max_disp_z_bracket2 = ", max_disp_z_bracket2, "\n")
    println("avg_stress_bracket1 = ", avg_stress_bracket1, "\n")
    println("avg_stress_bracket2 = ", avg_stress_bracket2, "\n")

    @test min_disp_x_bracket1 ≈ -3.255687122828699e-5 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -1.7910917248203928e-5 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.00012620821252524398 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.00017930987076132665 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -4.052985881622395e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.0016696791571882948 atol = 1e-12
    @test avg_stress_bracket1 ≈ [5261.147871837818 1.8229072090992345e6 165837.7136189107 256680.3604909731 245289.15366988638 -147308.3450816948] atol = 1e-12
    @test avg_stress_bracket2 ≈ [4.982364543425214e6 2.8123926560396124e6 -405418.24914701295 1.3230848216495877e6 607114.2790286107 2.233382099904088e6] atol = 1e-12
    @test sim.controller.schwarz_iters ≈ [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6] atol = 0

end

