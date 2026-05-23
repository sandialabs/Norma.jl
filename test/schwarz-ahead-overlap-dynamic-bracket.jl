# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Bracket HEX8-TET4 Dissipative Newmark" begin
    cp("../examples/ahead/overlap/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-1.yaml", "bracket-1.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-2.yaml", "bracket-2.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-1.g", "../bracket-1.g"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-2.g", "../bracket-2.g"; force=true)
    input_file = "bracket.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-5
    params["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_bracket1 = subsims[1].model
    model_bracket2 = subsims[2].model

    rm("bracket.yaml"; force=true)
    rm("bracket-1.yaml"; force=true)
    rm("bracket-2.yaml"; force=true)
    rm("../bracket-1.g"; force=true)
    rm("../bracket-2.g"; force=true)
    rm("bracket-1.e"; force=true)
    rm("bracket-2.e"; force=true)

    min_disp_x_bracket1 = minimum(model_bracket1.displacement[1, :])
    min_disp_y_bracket1 = minimum(model_bracket1.displacement[2, :])
    max_disp_z_bracket1 = maximum(model_bracket1.displacement[3, :])
    min_disp_x_bracket2 = minimum(model_bracket2.displacement[1, :])
    min_disp_y_bracket2 = minimum(model_bracket2.displacement[2, :])
    max_disp_z_bracket2 = maximum(model_bracket2.displacement[3, :])
    avg_stress_bracket1 = average_components(model_bracket1.stress)
    avg_stress_bracket2 = average_components(model_bracket2.stress)

    @test min_disp_x_bracket1 ≈ -5.230704712576306e-5 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -5.933993412379768e-6 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.00014337338409010306 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.0005079636241489688 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -5.188214899938537e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.007592937380626189 atol = 1e-12
    @test avg_stress_bracket1 ≈
        [-510160.2985190721 13618.32596222584 97405.25762366697 44143.889064639196 1.0434950811102245e7 1.516827315829358e6] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-5.102052720843673e6 -63490.06767812974 3.721911911481932e6 -747550.2645947285 1.7087825185346685e7 434117.5061109137] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        2,
        5,
        7,
        9,
        10,
        10,
        9,
        11,
        11,
        13,
        14,
        14,
        10,
        13,
        14,
        15,
        14,
        13,
        12,
        9,
        9,
        10,
        10,
        9,
        8,
        7,
        5,
        7,
        8,
        9,
        9,
        10,
        10,
        9,
        9,
        8,
        6,
        8,
        10,
        10,
        11,
        11,
        10,
        8,
        8,
        10,
        11,
        11,
        11,
        10,
    ] atol = 0
end

@testset "Schwarz AHeaD Overlap Dynamic Bracket TET4-TET4 Dissipative Newmark" begin
    cp("../examples/ahead/overlap/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-1.yaml", "bracket-1.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/dynamic/bracket-2.yaml", "bracket-2.yaml"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-1-tet.g", "../bracket-1.g"; force=true)
    cp("../examples/ahead/overlap/bracket/bracket-2-tet.g", "../bracket-2.g"; force=true)
    input_file = "bracket.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-5
    params["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_bracket1 = subsims[1].model
    model_bracket2 = subsims[2].model

    rm("bracket.yaml"; force=true)
    rm("bracket-1.yaml"; force=true)
    rm("bracket-2.yaml"; force=true)
    rm("../bracket-1.g"; force=true)
    rm("../bracket-2.g"; force=true)
    rm("bracket-1.e"; force=true)
    rm("bracket-2.e"; force=true)

    min_disp_x_bracket1 = minimum(model_bracket1.displacement[1, :])
    min_disp_y_bracket1 = minimum(model_bracket1.displacement[2, :])
    max_disp_z_bracket1 = maximum(model_bracket1.displacement[3, :])
    min_disp_x_bracket2 = minimum(model_bracket2.displacement[1, :])
    min_disp_y_bracket2 = minimum(model_bracket2.displacement[2, :])
    max_disp_z_bracket2 = maximum(model_bracket2.displacement[3, :])
    avg_stress_bracket1 = average_components(model_bracket1.stress)
    avg_stress_bracket2 = average_components(model_bracket2.stress)

    @test min_disp_x_bracket1 ≈ -0.00013489122345448565 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -6.88670267136593e-5 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.0019608960014893886 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.00017674720316077086 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -7.039442537440377e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.0037945898323606952 atol = 1e-12
    @test avg_stress_bracket1 ≈
        [891484.3827023633 475950.11641397886 -40682.88983479411 38071.02301581592 2.8594511190481563e6 -85643.29050974561] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-587643.9418079571 315211.6345669816 496153.6551374258 9562.532715049952 1.6393813991316264e6 -228194.52047169508] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        2,
        2,
        3,
        3,
        3,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
    ] atol = 0
end
