# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Plate HEX8-HEX8 Dissipative Newmark" begin
    cp("../examples/ahead/overlap/plate/dynamic/plate.yaml", "plate.yaml"; force=true)
    cp("../examples/ahead/overlap/plate/dynamic/plate-1.yaml", "plate-1.yaml"; force=true)
    cp("../examples/ahead/overlap/plate/dynamic/plate-2.yaml", "plate-2.yaml"; force=true)
    cp("../examples/ahead/overlap/plate/plate-1.g", "../plate-1.g"; force=true)
    cp("../examples/ahead/overlap/plate/plate-2.g", "../plate-2.g"; force=true)
    input_file = "plate.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-5
    params["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_plate1 = subsims[1].model
    model_plate2 = subsims[2].model

    rm("plate.yaml"; force=true)
    rm("plate-1.yaml"; force=true)
    rm("plate-2.yaml"; force=true)
    rm("../plate-1.g"; force=true)
    rm("../plate-2.g"; force=true)
    rm("plate-1.e"; force=true)
    rm("plate-2.e"; force=true)

    min_disp_x_plate1 = minimum(model_plate1.displacement[1, :])
    min_disp_y_plate1 = minimum(model_plate1.displacement[2, :])
    max_disp_z_plate1 = maximum(model_plate1.displacement[3, :])
    min_disp_x_plate2 = minimum(model_plate2.displacement[1, :])
    min_disp_y_plate2 = minimum(model_plate2.displacement[2, :])
    max_disp_z_plate2 = maximum(model_plate2.displacement[3, :])
    avg_stress_plate1 = average_components(model_plate1.stress)
    avg_stress_plate2 = average_components(model_plate2.stress)

    @test min_disp_x_plate1 ≈ -9.921851879713967e-5 atol = 1e-12
    @test min_disp_y_plate1 ≈ -1.2983244377276493e-5 atol = 1e-12
    @test max_disp_z_plate1 ≈ 0.0006671119244265724 atol = 1e-12
    @test min_disp_x_plate2 ≈ -9.931932093781848e-5 atol = 1e-12
    @test min_disp_y_plate2 ≈ -1.6908937976650718e-5 atol = 1e-12
    @test max_disp_z_plate2 ≈ 0.0011022795117623524 atol = 1e-12
    @test avg_stress_plate1 ≈
        [287501.9620223973 3399.8276901241115 125544.50901440646 2.3343535035681456e-7 2.1548960065833375e7 5.8586988299607465e-6] atol =
        1.0e1
    @test avg_stress_plate2 ≈
        [212642.4726788842 -26447.452635092955 41385.32446019253 1.264027626045087e-6 3.0982580107926014e6 1.074629208233091e-7] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
    ] atol = 0
end
