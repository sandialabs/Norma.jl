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

    @test min_disp_x_bracket1 ≈ -5.230704712576306e-5 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -5.933993412379768e-6 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.00014337338409010306 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.0005079636241489688 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -5.188214899938537e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.007592937380626189 atol = 1e-12
    @test avg_stress_bracket1 ≈
        [-443291.51669701113 19019.88910771844 25134.91272134775 39068.427259961936 1.0434914230313826e7 1.5174866682280425e6] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-4.213214851786042e6 89029.68204894006 2.680554292983087e6 -774307.7048034652 1.735630636595213e7 416093.38218943775] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        3,
        6,
        8,
        10,
        11,
        11,
        10,
        12,
        12,
        14,
        15,
        15,
        11,
        14,
        15,
        16,
        15,
        14,
        13,
        10,
        10,
        11,
        11,
        10,
        9,
        8,
        6,
        8,
        9,
        10,
        10,
        11,
        11,
        10,
        10,
        9,
        7,
        9,
        11,
        11,
        12,
        12,
        11,
        9,
        9,
        11,
        12,
        12,
        12,
        11,
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

    @test min_disp_x_bracket1 ≈ -0.00013489122345448565 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -6.88670267136593e-5 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 0.0019608960014893886 atol = 1e-12
    @test min_disp_x_bracket2 ≈ -0.00017674720316077086 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -7.039442537440377e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.0037945898323606952 atol = 1e-12
    @test avg_stress_bracket1 ≈
        [928393.0804969968 484726.67620943306 -86368.14791471219 39681.35369128554 2.8561263044046783e6 -84999.59564231828] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-498781.1969274577 380458.68521612795 342043.85822372086 18008.589711749457 1.654179334345221e6 -231776.00506447483] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        3,
        3,
        4,
        4,
        4,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        4,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        3,
        3,
    ] atol = 0
end
