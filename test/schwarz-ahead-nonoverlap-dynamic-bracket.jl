# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Bracket TET4-TET4 Dissipative Newmark" begin
    cp("../examples/ahead/nonoverlap/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/bracket/dynamic/bracket-1.yaml", "bracket-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/bracket/dynamic/bracket-2.yaml", "bracket-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/bracket/bracket-1.g", "../bracket-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/bracket/bracket-2.g", "../bracket-2.g"; force=true)
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

    @test min_disp_x_bracket1 ≈ -2.3109772377374713e-5 atol = 1e-8
    @test min_disp_y_bracket1 ≈ -2.623329287039511e-5 atol = 1e-8
    @test max_disp_z_bracket1 ≈ 8.31361801725248e-5 atol = 1e-8
    @test min_disp_x_bracket2 ≈ -0.00010880437118793661 atol = 1e-8
    @test min_disp_y_bracket2 ≈ -4.797628710582998e-5 atol = 1e-8
    @test max_disp_z_bracket2 ≈ 0.0008416981567619049 atol = 1e-8
    @test avg_stress_bracket1 ≈
        [1.3835933504327782e6 140004.03823244444 108241.53414225587 2528.519705772454 -4.432997150138795e6 -81274.16079474054] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-962012.7875065166 256050.96159196054 -32566.10066348576 -381336.3590060364 1.0873933562525716e6 -1.2683959067522327e6] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        1,
        2,
        2,
        1,
        2,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        5,
        5,
        5,
        5,
        5,
        5,
        5,
        5,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
    ] atol = 0
end
