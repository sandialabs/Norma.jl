# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Notched Cylinder HEX8-HEX8" begin
    cp("../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder.yaml", "notched-cylinder.yaml"; force=true)
    cp(
        "../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder-1.yaml",
        "notched-cylinder-1.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder-2.yaml",
        "notched-cylinder-2.yaml";
        force=true,
    )
    cp("../examples/ahead/overlap/notched-cylinder/notched-cylinder-1.g", "../notched-cylinder-1.g"; force=true)
    cp("../examples/ahead/overlap/notched-cylinder/notched-cylinder-2.g", "../notched-cylinder-2.g"; force=true)
    input_file = "notched-cylinder.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 0.1
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_1 = subsims[1].model
    model_2 = subsims[2].model

    rm("notched-cylinder.yaml")
    rm("notched-cylinder-1.yaml")
    rm("notched-cylinder-2.yaml")
    rm("../notched-cylinder-1.g")
    rm("../notched-cylinder-2.g")
    rm("notched-cylinder-1.e")
    rm("notched-cylinder-2.e")

    min_disp_x_1 = minimum(model_1.current[1, :] - model_1.reference[1, :])
    min_disp_y_1 = minimum(model_1.current[2, :] - model_1.reference[2, :])
    max_disp_z_1 = maximum(model_1.current[3, :] - model_1.reference[3, :])
    min_disp_x_2 = minimum(model_2.current[1, :] - model_2.reference[1, :])
    min_disp_y_2 = minimum(model_2.current[2, :] - model_2.reference[2, :])
    min_disp_z_2 = minimum(model_2.current[3, :] - model_2.reference[3, :])
    avg_stress_1 = average_components(model_1.stress)
    avg_stress_2 = average_components(model_2.stress)

    @test min_disp_x_1 ≈ -1.6007099461888552e-5 atol = 1e-6
    @test min_disp_y_1 ≈ -1.6007073405502337e-5 atol = 1e-6
    @test max_disp_z_1 ≈ 0.00011007952600178977 atol = 1e-6
    @test min_disp_x_2 ≈ -1.6186344441136702e-5 atol = 1e-6
    @test min_disp_y_2 ≈ -1.618628933883204e-5 atol = 1e-6
    @test min_disp_z_2 ≈ 5.1726659805836905e-5 atol = 1e-6
    @test avg_stress_1 ≈
        [119658.24151708037 119759.96550264646 2.105505961105632e6 280369.0355588143 281504.50820556487 57168.0765508701] atol =
        1.0e1
    @test avg_stress_2 ≈
        [-101340.89672529901 -101446.28712096746 1.5977255795638144e6 121800.33493609141 122483.30918545111 37064.478126464295] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [1, 1, 1, 1, 4, 4, 4, 4, 5, 5] atol = 0
end

@testset "Schwarz AHeaD Overlap Dynamic Notched Cylinder TET10-HEX8" begin
    cp("../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder.yaml", "notched-cylinder.yaml"; force=true)
    cp(
        "../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder-1.yaml",
        "notched-cylinder-1.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/notched-cylinder/dynamic/notched-cylinder-2.yaml",
        "notched-cylinder-2.yaml";
        force=true,
    )
    cp("../examples/ahead/overlap/notched-cylinder/notched-cylinder-1-tet10.g", "../notched-cylinder-1.g"; force=true)
    cp("../examples/ahead/overlap/notched-cylinder/notched-cylinder-2.g", "../notched-cylinder-2.g"; force=true)
    input_file = "notched-cylinder.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 0.1
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_1 = subsims[1].model
    model_2 = subsims[2].model

    rm("notched-cylinder.yaml")
    rm("notched-cylinder-1.yaml")
    rm("notched-cylinder-2.yaml")
    rm("../notched-cylinder-1.g")
    rm("../notched-cylinder-2.g")
    rm("notched-cylinder-1.e")
    rm("notched-cylinder-2.e")

    min_disp_x_1 = minimum(model_1.current[1, :] - model_1.reference[1, :])
    min_disp_y_1 = minimum(model_1.current[2, :] - model_1.reference[2, :])
    max_disp_z_1 = maximum(model_1.current[3, :] - model_1.reference[3, :])
    min_disp_x_2 = minimum(model_2.current[1, :] - model_2.reference[1, :])
    min_disp_y_2 = minimum(model_2.current[2, :] - model_2.reference[2, :])
    min_disp_z_2 = minimum(model_2.current[3, :] - model_2.reference[3, :])
    avg_stress_1 = average_components(model_1.stress)
    avg_stress_2 = average_components(model_2.stress)

    @test min_disp_x_1 ≈ -1.6651285885280198e-5 atol = 1e-12
    @test min_disp_y_1 ≈ -1.662806113268689e-5 atol = 1e-12
    @test max_disp_z_1 ≈ 0.00010485387156757009 atol = 1e-12
    @test min_disp_x_2 ≈ -1.6461199941086857e-5 atol = 1e-12
    @test min_disp_y_2 ≈ -1.643006006920439e-5 atol = 1e-12
    @test min_disp_z_2 ≈ 4.895828234158239e-5 atol = 1e-12
    @test avg_stress_1 ≈
        [74920.96247178296 72298.92321097263 1.678913240195309e6 259770.47737702358 271427.7166969666 76184.88555466678] atol =
        1.0e1
    @test avg_stress_2 ≈
        [-89472.31510791164 -89526.66568169543 1.6857645398230806e6 114694.99228486177 114706.49638144106 34816.156918097106] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [1, 1, 1, 1, 1, 3, 3, 3, 3, 3] atol = 0
end
