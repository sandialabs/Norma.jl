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

    # Baseline of the relative Schwarz criterion normalized by the solution (not
    # by the positions X + u, under which every early stop accepted its first
    # iterate). With the example's absolute tolerance of 1e-6 m the largest
    # displacement differs by 2.6% from the answer converged to 1e-14 (4.0%
    # under the former criterion).
    @test min_disp_x_bracket1 ≈ -2.4039728277800348e-5 atol = 1e-8
    @test min_disp_y_bracket1 ≈ -2.998071050631785e-5 atol = 1e-8
    @test max_disp_z_bracket1 ≈ 8.24591776516464e-5 atol = 1e-8
    @test min_disp_x_bracket2 ≈ -0.00010887463311789825 atol = 1e-8
    @test min_disp_y_bracket2 ≈ -5.13008163845687e-5 atol = 1e-8
    @test max_disp_z_bracket2 ≈ 0.0008415476201673325 atol = 1e-8
    @test avg_stress_bracket1 ≈
        [652910.6580490598 -12123.223816426462 50154.572220476926 5493.09267826974 -4.39286450144604e6 -146484.54939631541] atol =
        1.0e1
    @test avg_stress_bracket2 ≈
        [759835.7486153797 405964.9277498441 -7906.673323030443 -399843.83671783155 887611.4735772951 235560.99596637356] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 6, 7, 8, 8, 8, 8, 7, 7, 7, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 8, 8, 8] atol = 0
end
