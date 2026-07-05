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

    @test min_disp_x_holder0 ≈ -4.851961547122641e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.0 atol = 1e-12
    @test min_disp_x_gauge ≈ -7.731601561910473e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.00027946470502994306 atol = 1e-12
    @test max_disp_z_gauge ≈ 4.334356104108237e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -4.809985907257732e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.00023578919305653514 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0 atol = 1e-12
    @test avg_stress_holder0 ≈
        [-652.6493536915126 1.7781311578399208e6 -51816.284666868196 -46496.526237025195 1154.4198413653457 1353.642505548652] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [16536.400139704365 1.8419002146714842e6 116827.12132240279 80.43536267459744 29330.79334629724 197.50639392421243] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [-385.3424601165005 1.7781315932129258e6 -51769.00760511988 46422.108997527466 1145.1823849946404 -1464.4671474624108] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [26, 37, 43, 46, 49] atol = 0
end
