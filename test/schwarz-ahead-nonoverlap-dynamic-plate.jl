# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Plate HEX8-HEX8 Dissipative Newmark" begin
    cp("../examples/ahead/nonoverlap/plate/dynamic/plate.yaml", "plate.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/plate/dynamic/plate-1.yaml", "plate-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/plate/dynamic/plate-2.yaml", "plate-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/plate/plate-1.g", "../plate-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/plate/plate-2.g", "../plate-2.g"; force=true)
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

    @test min_disp_x_plate1 ≈ -5.118470446491376e-5 atol = 1e-12
    @test min_disp_y_plate1 ≈ -2.7317646505833904e-5 atol = 1e-12
    @test max_disp_z_plate1 ≈ 0.00019738422607168404 atol = 1e-12
    @test min_disp_x_plate2 ≈ -8.911519006287206e-5 atol = 1e-12
    @test min_disp_y_plate2 ≈ -4.8536073557900016e-5 atol = 1e-12
    @test max_disp_z_plate2 ≈ 0.0006959607671721849 atol = 1e-12
    @test avg_stress_plate1 ≈
        [870477.1328662889 59551.25146841846 28072.00550264897 -1.9753060769289732e-5 3.8743046672905e6 -1.1942523997277021e-5] atol =
        1.0e1
    @test avg_stress_plate2 ≈
        [489422.0518422433 -126742.85120534342 78364.30882320109 7.712837250437588e-6 6.510011260192466e6 -9.027880635888627e-7] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [53, 64, 65, 72, 76, 74, 70, 71, 74, 75, 74, 72, 69, 67, 64, 62, 61, 61, 63, 65, 67, 68, 67, 65, 63, 62, 62, 63, 63, 64, 64, 64, 64, 64, 64, 64, 64, 65, 66, 66, 67, 67, 67, 66, 65, 64, 63, 61, 59, 57] atol = 0
end
