# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Laser Weld HEX8-HEX8 Symmetry BCs" begin
    cp("../examples/ahead/nonoverlap/laser-weld/dynamic/symmetry/laser-weld.yaml", "laser-weld.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/dynamic/symmetry/holder-0.yaml", "holder-0.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/dynamic/symmetry/holder-1.yaml", "holder-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/dynamic/symmetry/gauge.yaml", "gauge.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/holder-0.g", "../../holder-0.g"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/holder-1.g", "../../holder-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/laser-weld/gauge.g", "../../gauge.g"; force=true)
    input_file = "laser-weld.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 0.05
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_holder0 = subsims[1].model
    model_gauge = subsims[2].model
    model_holder1 = subsims[3].model

    rm("laser-weld.yaml"; force=true)
    rm("holder-0.yaml"; force=true)
    rm("holder-1.yaml"; force=true)
    rm("gauge.yaml"; force=true)
    rm("../../holder-0.g"; force=true)
    rm("../../holder-1.g"; force=true)
    rm("../../gauge.g"; force=true)
    rm("holder-0.e"; force=true)
    rm("holder-1.e"; force=true)
    rm("gauge.e"; force=true)

    min_disp_x_holder0 = minimum(model_holder0.displacement[1, :])
    min_disp_y_holder0 = minimum(model_holder0.displacement[2, :])
    max_disp_z_holder0 = maximum(model_holder0.displacement[3, :])
    min_disp_x_gauge = minimum(model_gauge.displacement[1, :])
    min_disp_y_gauge = minimum(model_gauge.displacement[2, :])
    max_disp_z_gauge = maximum(model_gauge.displacement[3, :])
    min_disp_x_holder1 = minimum(model_holder1.displacement[1, :])
    min_disp_y_holder1 = minimum(model_holder1.displacement[2, :])
    max_disp_z_holder1 = maximum(model_holder1.displacement[3, :])
    avg_stress_holder0 = average_components(model_holder0.stress)
    avg_stress_gauge = average_components(model_gauge.stress)
    avg_stress_holder1 = average_components(model_holder1.stress)

    @test min_disp_x_holder0 ≈ -4.851628886395898e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.0 atol = 1e-12
    @test min_disp_x_gauge ≈ -7.731358906070013e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.0002794662877980967 atol = 1e-12
    @test max_disp_z_gauge ≈ 4.334015744528366e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -4.8101884538691166e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.00023579973361378126 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0 atol = 1e-12
    @test avg_stress_holder0 ≈
        [-651.1464725557873 1.7780739370373667e6 -51812.70666114553 -46494.7446878697 1154.5056911923787 1353.1974148960885] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [16539.536519038895 1.841840302350875e6 116828.34012429455 76.90140926572899 29330.37279517582 196.6113248495793] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [-390.26264976813616 1.7780748572807203e6 -51773.82697053141 46424.020268542656 1145.3620331332904 -1462.8716682712588] atol =
        1.0e1
    # Iteration counts and converged values reflect the slot-keyed relaxation
    # state (commit 0deb579): with the controller step (0.01) larger than the
    # adaptively-substepped subdomains, these windowed stops now blend the
    # fixed-θ relaxation per substep-time slot instead of across time.
    @test sim.controller.schwarz_iters ≈ [33, 40, 44, 47, 49] atol = 0
end
