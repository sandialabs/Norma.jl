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

    @test min_disp_x_holder0 ≈ -4.8503353649823056e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.0 atol = 1e-12
    @test min_disp_x_gauge ≈ -7.732050750191575e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.0002794601407257627 atol = 1e-12
    @test max_disp_z_gauge ≈ 4.334628364924738e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -4.810143238079646e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.00023577100314417488 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0 atol = 1e-12
    @test avg_stress_holder0 ≈
        [-639.6261094759388 1.7782199539431625e6 -51807.83120981262 -46504.79388129707 1156.5355292243473 1346.1170342658352] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [16522.192587960613 1.8419942446213954e6 116817.7849691327 88.09734397642102 29332.667642010096 205.35636784946962] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [-384.0829568066664 1.7782249323431996e6 -51768.35215766566 46423.40270612918 1145.662052893047 -1463.950585851161] atol =
        1.0e1
    # Iteration counts and converged values reflect the slot-keyed relaxation
    # state (commit 0deb579): with the controller step (0.01) larger than the
    # adaptively-substepped subdomains, these windowed stops now blend the
    # fixed-θ relaxation per substep-time slot instead of across time. They also
    # reflect relaxing only the Dirichlet side of an interface (commit 38bd028):
    # the two holders couple to the same gauge, so before that gate the second
    # of them re-blended the shared iterate the first had just written.
    @test sim.controller.schwarz_iters ≈ [34, 42, 47, 50, 52] atol = 0
end
