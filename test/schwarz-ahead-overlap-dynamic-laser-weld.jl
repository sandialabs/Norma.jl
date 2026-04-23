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
    min_disp_z_gauge = minimum(model_gauge.displacement[3, :])
    min_disp_x_holder1 = minimum(model_holder1.displacement[1, :])
    min_disp_y_holder1 = minimum(model_holder1.displacement[2, :])
    max_disp_z_holder1 = maximum(model_holder1.displacement[3, :])
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
        [56184.76411164522 1.6118497760906e6 56431.37537113758 -194.60281369086042 -8.16978323806931e-9 -4.3784983739432925e-8] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [12652.963424860558 1.6507177256895711e6 28154.641539885597 -178.457550278936 -2.1974017066487e-9 -1.2105597220298754e-8] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [56184.9630748282 1.6118490402532257e6 56429.46812580038 179.82423795232404 3.209146749820017e-9 -9.850264177657664e-8] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [22, 26, 28, 29, 30] atol = 0
end
