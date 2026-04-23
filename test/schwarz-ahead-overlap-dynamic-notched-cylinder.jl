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

    @test min_disp_x_1 ≈ -1.6007099461888552e-5 atol = 1e-6
    @test min_disp_y_1 ≈ -1.6007073405502337e-5 atol = 1e-6
    @test max_disp_z_1 ≈ 0.00011007952600178977 atol = 1e-6
    @test min_disp_x_2 ≈ -1.6186344441136702e-5 atol = 1e-6
    @test min_disp_y_2 ≈ -1.618628933883204e-5 atol = 1e-6
    @test min_disp_z_2 ≈ 5.1726659805836905e-5 atol = 1e-6
    @test avg_stress_1 ≈
        [119098.66905646767 119201.64540623137 2.1066238536625747e6 279333.00110998185 280460.3455847958 56808.524330534376] atol =
        1.0e1
    @test avg_stress_2 ≈
        [-101423.59028091333 -101528.52409858494 1.5978905100970943e6 121611.29785378547 122292.56850566593 37011.15803770428] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [0, 0, 0, 0, 3, 3, 3, 3, 4, 4] atol = 0
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

    @test min_disp_x_1 ≈ -1.6651285885280198e-5 atol = 1e-12
    @test min_disp_y_1 ≈ -1.662806113268689e-5 atol = 1e-12
    @test max_disp_z_1 ≈ 0.00010485387156757009 atol = 1e-12
    @test min_disp_x_2 ≈ -1.6461199941086857e-5 atol = 1e-12
    @test min_disp_y_2 ≈ -1.643006006920439e-5 atol = 1e-12
    @test min_disp_z_2 ≈ 4.895828234158239e-5 atol = 1e-12
    @test avg_stress_1 ≈
        [74255.45497241359 71739.39938200016 1.680138271523533e6 258948.1851516863 270573.5517476705 75797.27494781958] atol =
        1.0e1
    @test avg_stress_2 ≈
        [-89548.343641492 -89603.24718234365 1.6859171498567525e6 114497.79715981241 114510.3753457256 34767.594667242905] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [0, 0, 0, 0, 0, 2, 2, 2, 2, 2] atol = 0
end
