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

    @test min_disp_x_bracket1 ≈ -2.292363909493874e-5 atol = 1e-12
    @test min_disp_y_bracket1 ≈ -2.6772610574016253e-5 atol = 1e-12
    @test max_disp_z_bracket1 ≈ 8.406701507824004e-5 atol = 1e-8
    @test min_disp_x_bracket2 ≈ -0.00010879133317491518 atol = 1e-12
    @test min_disp_y_bracket2 ≈ -4.821150003511687e-5 atol = 1e-12
    @test max_disp_z_bracket2 ≈ 0.0008416546150997427 atol = 1e-8
    @test avg_stress_bracket1 ≈
        [1.4344587897170822e6 147404.28440152534 80028.09147772216 1331.5200561233653 -4.470431261192811e6 -102990.28645337945] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [-887460.5393152214 320509.9388991009 -380020.14019136556 -350857.5055218009 1.0378814655089243e6 -1.285348432558279e6] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        4,
        4,
        4,
        5,
        7,
        8,
        9,
        9,
        9,
        9,
        8,
        8,
        8,
        9,
        10,
        10,
        11,
        11,
        11,
        11,
        11,
        11,
        11,
        11,
        10,
        10,
        10,
        10,
        9,
        9,
        9,
    ] atol = 0
end
