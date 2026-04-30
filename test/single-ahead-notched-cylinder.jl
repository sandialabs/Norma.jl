# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Notched Cylinder HEX8" begin
    cp("../examples/ahead/single/notched-cylinder/dynamic/notched-cylinder.yaml", "notched-cylinder.yaml"; force=true)
    cp("../examples/ahead/single/notched-cylinder/notched-cylinder.g", "../notched-cylinder.g"; force=true)
    input_file = "notched-cylinder.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 0.01
    params["time integrator"]["final time"] = 0.1
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("notched-cylinder.yaml"; force=true)
    rm("../notched-cylinder.g"; force=true)
    rm("notched-cylinder.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -1.9074596682447376e-5 atol = 1e-12
    @test min_disp_y ≈ -1.907460302130043e-5 atol = 1e-12
    @test max_disp_z ≈ 0.0001566191478555093 atol = 1e-12
    @test avg_stress ≈
        [49709.55780582162 49709.58813454913 2.378064481860813e6 240111.2224446391 239864.70358086814 55329.43423102153] atol =
        1e-6
end

@testset "AHeaD Single Quasistatic Notched Cylinder HEX8" begin
    cp(
        "../examples/ahead/single/notched-cylinder/quasistatic/notched-cylinder.yaml",
        "notched-cylinder.yaml";
        force=true,
    )
    cp("../examples/ahead/single/notched-cylinder/notched-cylinder.g", "../notched-cylinder.g"; force=true)
    input_file = "notched-cylinder.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 0.1
    params["time integrator"]["final time"] = 0.5
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("notched-cylinder.yaml"; force=true)
    rm("../notched-cylinder.g"; force=true)
    rm("notched-cylinder.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.00036607877763784186 atol = 1e-12
    @test min_disp_y ≈ -0.00036607891968151035 atol = 1e-12
    @test max_disp_z ≈ 0.0032000000000000015 atol = 1e-12
    @test avg_stress ≈
        [724971.4756383626 724979.2799407415 4.69334188431597e7 4.575342747126081e6 4.570406479227432e6 1.019033684998046e6] atol =
        1.0e1
end
