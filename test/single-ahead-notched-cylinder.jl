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

    rm("notched-cylinder.yaml")
    rm("../notched-cylinder.g")
    rm("notched-cylinder.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -1.9074596682447376e-5 atol = 1e-12
    @test min_disp_y ≈ -1.907460302130043e-5 atol = 1e-12
    @test max_disp_z ≈ 0.0001566191478555093 atol = 1e-12
    @test avg_stress ≈
        [50251.61286085953 50251.86478138371 2.3769801501592183e6 241114.77959964555 240866.46234339915 55675.444425864385] atol =
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

    rm("notched-cylinder.yaml")
    rm("../notched-cylinder.g")
    rm("notched-cylinder.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.00036607877763784186 atol = 1e-12
    @test min_disp_y ≈ -0.00036607891968151035 atol = 1e-12
    @test max_disp_z ≈ 0.0032000000000000015 atol = 1e-12
    @test avg_stress ≈
        [917864.3278006782 917964.6350355675 4.654754063590306e7 4.923095033584741e6 4.9175174443381205e6 1.1421705828295958e6] atol =
        1.0e1
end
