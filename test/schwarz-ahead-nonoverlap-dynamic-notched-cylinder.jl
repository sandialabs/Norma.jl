# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Notched Cylinder HEX8-HEX8" begin
    cp(
        "../examples/ahead/nonoverlap/notched-cylinder/dynamic/notched-cylinder.yaml",
        "notched-cylinder.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/nonoverlap/notched-cylinder/dynamic/notched-cylinder-1.yaml",
        "notched-cylinder-1.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/nonoverlap/notched-cylinder/dynamic/notched-cylinder-2.yaml",
        "notched-cylinder-2.yaml";
        force=true,
    )
    cp("../examples/ahead/nonoverlap/notched-cylinder/notched-cylinder-1.g", "../notched-cylinder-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/notched-cylinder/notched-cylinder-2.g", "../notched-cylinder-2.g"; force=true)
    input_file = "notched-cylinder.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 0.05
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_1 = subsims[1].model
    model_2 = subsims[2].model

    rm("notched-cylinder.yaml"; force=true)
    rm("notched-cylinder-1.yaml"; force=true)
    rm("notched-cylinder-2.yaml"; force=true)
    rm("../notched-cylinder-1.g"; force=true)
    rm("../notched-cylinder-2.g"; force=true)
    rm("notched-cylinder-1.e"; force=true)
    rm("notched-cylinder-2.e"; force=true)

    min_disp_x_1 = minimum(model_1.displacement[1, :])
    min_disp_y_1 = minimum(model_1.displacement[2, :])
    max_disp_z_1 = maximum(model_1.displacement[3, :])
    min_disp_x_2 = minimum(model_2.displacement[1, :])
    min_disp_y_2 = minimum(model_2.displacement[2, :])
    min_disp_z_2 = minimum(model_2.displacement[3, :])
    avg_stress_1 = average_components(model_1.stress)
    avg_stress_2 = average_components(model_2.stress)

    @test min_disp_x_1 ≈ -4.037470161949048e-6 atol = 1e-12
    @test min_disp_y_1 ≈ -4.0374652743727105e-6 atol = 1e-12
    @test max_disp_z_1 ≈ 2.8125101839471256e-5 atol = 1e-12
    @test min_disp_x_2 ≈ -3.9155928465825706e-6 atol = 1e-12
    @test min_disp_y_2 ≈ -3.915589896008824e-6 atol = 1e-12
    @test min_disp_z_2 ≈ 2.5053561481486768e-5 atol = 1e-12
    @test avg_stress_1 ≈
        [32439.05762287872 32465.61490606161 540278.4092679358 72018.91185980223 72309.93048014899 14488.221920501783] atol =
        1.0e1
    @test avg_stress_2 ≈
        [-20448.838840349512 -20470.797065627095 377569.57314743317 14471.770489157607 14562.262796808678 5106.555522648798] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [7, 9, 10, 10, 11] atol = 0
end
