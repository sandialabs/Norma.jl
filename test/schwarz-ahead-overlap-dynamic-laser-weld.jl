# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Laser Weld HEX8-HEX8 Clamped BCs" begin
    cp("../examples/ahead/overlap/laser-weld/dynamic/clamped/laser-weld.yaml", "laser-weld.yaml"; force=true)
    cp("../examples/ahead/overlap/laser-weld/dynamic/clamped/holder-0.yaml", "holder-0.yaml"; force=true)
    cp("../examples/ahead/overlap/laser-weld/dynamic/clamped/holder-1.yaml", "holder-1.yaml"; force=true)
    cp("../examples/ahead/overlap/laser-weld/dynamic/clamped/gauge.yaml", "gauge.yaml"; force=true)
    cp("../examples/ahead/overlap/laser-weld/holder-0.g", "../../holder-0.g"; force=true)
    cp("../examples/ahead/overlap/laser-weld/holder-1.g", "../../holder-1.g"; force=true)
    cp("../examples/ahead/overlap/laser-weld/gauge.g", "../../gauge.g"; force=true)
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

    rm("laser-weld.yaml")
    rm("holder-0.yaml")
    rm("holder-1.yaml")
    rm("gauge.yaml")
    rm("../../holder-0.g")
    rm("../../holder-1.g")
    rm("../../gauge.g")
    rm("holder-0.e")
    rm("holder-1.e")
    rm("gauge.e")

    min_disp_x_holder0 = minimum(model_holder0.current[1, :] - model_holder0.reference[1, :])
    min_disp_y_holder0 = minimum(model_holder0.current[2, :] - model_holder0.reference[2, :])
    max_disp_z_holder0 = maximum(model_holder0.current[3, :] - model_holder0.reference[3, :])
    min_disp_x_gauge = minimum(model_gauge.current[1, :] - model_gauge.reference[1, :])
    min_disp_y_gauge = minimum(model_gauge.current[2, :] - model_gauge.reference[2, :])
    min_disp_z_gauge = minimum(model_gauge.current[3, :] - model_gauge.reference[3, :])
    min_disp_x_holder1 = minimum(model_holder1.current[1, :] - model_holder1.reference[1, :])
    min_disp_y_holder1 = minimum(model_holder1.current[2, :] - model_holder1.reference[2, :])
    max_disp_z_holder1 = maximum(model_holder1.current[3, :] - model_holder1.reference[3, :])
    avg_stress_holder0 = average_components(model_holder0.stress)
    avg_stress_gauge = average_components(model_gauge.stress)
    avg_stress_holder1 = average_components(model_holder1.stress)

    @test min_disp_x_holder0 ≈ -3.0733236697749744e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.00020594161227460822 atol = 1e-12
    @test min_disp_x_gauge ≈ -6.11625728770826e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.0005587147296370887 atol = 1e-12
    @test min_disp_z_gauge ≈ 2.3880604363396563e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -3.0734132231508005e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.00025514717233354745 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0002059887389275507 atol = 1e-12
    @test avg_stress_holder0 ≈
        [56175.3792091251 1.6118582225707504e6 56432.31379334061 -2152.420453309265 3.3326750781270676e-9 1.2065091337850238e-8] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [12641.620130337215 1.6499174233606872e6 28966.287163254172 -179.1767306388687 -1.3967327005127079e-9 -1.079985106175825e-9] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [56175.5776544634 1.611857519176203e6 56430.374623511234 2138.110887978947 5.693916591553716e-9 -8.169333644521733e-8] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [23, 27, 29, 30, 31] atol = 0
end
