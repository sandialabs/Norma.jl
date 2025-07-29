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
    max_disp_z_gauge = maximum(model_gauge.current[3, :] - model_gauge.reference[3, :])
    min_disp_x_holder1 = minimum(model_holder1.current[1, :] - model_holder1.reference[1, :])
    min_disp_y_holder1 = minimum(model_holder1.current[2, :] - model_holder1.reference[2, :])
    max_disp_z_holder1 = maximum(model_holder1.current[3, :] - model_holder1.reference[3, :])
    avg_stress_holder0 = average_components(model_holder0.stress)
    avg_stress_gauge = average_components(model_gauge.stress)
    avg_stress_holder1 = average_components(model_holder1.stress)

    @test min_disp_x_holder0 ≈ -4.8524923688300636e-5 atol = 1e-12
    @test min_disp_y_holder0 ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z_holder0 ≈ 0.0 atol = 1e-12
    @test min_disp_x_gauge ≈ -7.768426106180559e-5 atol = 1e-12
    @test min_disp_y_gauge ≈ -0.00027640414578330996 atol = 1e-12
    @test max_disp_z_gauge ≈ 4.4372026740832626e-5 atol = 1e-12
    @test min_disp_x_holder1 ≈ -4.8953646255081584e-5 atol = 1e-12
    @test min_disp_y_holder1 ≈ 0.00023478463555726137 atol = 1e-12
    @test max_disp_z_holder1 ≈ 0.0 atol = 1e-12
    @test avg_stress_holder0 ≈
        [-243.0717661305563 1.7798335903690055e6 -49963.13502708086 -46273.081927730775 1209.6636290360905 1525.1001490585427] atol =
        1.0e1
    @test avg_stress_gauge ≈
        [16219.565177010196 1.845519057684895e6 122146.1419563506 4904.312109770015 29899.525543416865 94.33582939764187] atol =
        1.0e1
    @test avg_stress_holder1 ≈
        [-1056.0441890865773 1.7715466435733493e6 -54239.702779201296 50590.033301702606 1620.6847557192993 -1375.575102437882] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [26, 38, 44, 48, 51] atol = 0
end
