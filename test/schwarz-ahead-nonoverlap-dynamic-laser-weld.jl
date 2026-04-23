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

    @test min_disp_x_holder0 ≈ -4.8291007049908186e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.0 atol = 1e-12
    @test min_disp_x_gauge ≈ -7.768875531261044e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.00027639374248203397 atol = 1e-12
    @test max_disp_z_gauge ≈ 4.43695678256778e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -4.89520944503552e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.0002347720124472319 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0 atol = 1e-12
    @test avg_stress_holder0 ≈
        [-346.98214729874337 1.779760561862967e6 -50233.06035476495 -45951.88333708876 1255.3331344265928 1554.3140211168936] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [16328.8628137128 1.84664104865588e6 121467.20818681682 4584.026400893279 29936.877739450956 69.36019489750036] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [-1054.751995236128 1.7715978058931958e6 -54259.9126842316 50591.113276641314 1620.4225341530955 -1379.8420874727585] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [25, 37, 43, 47, 50] atol = 0
end
