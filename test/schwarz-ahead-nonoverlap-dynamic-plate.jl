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

    rm("plate.yaml")
    rm("plate-1.yaml")
    rm("plate-2.yaml")
    rm("../plate-1.g")
    rm("../plate-2.g")
    rm("plate-1.e")
    rm("plate-2.e")

    min_disp_x_plate1 = minimum(model_plate1.current[1, :] - model_plate1.reference[1, :])
    min_disp_y_plate1 = minimum(model_plate1.current[2, :] - model_plate1.reference[2, :])
    max_disp_z_plate1 = maximum(model_plate1.current[3, :] - model_plate1.reference[3, :])
    min_disp_x_plate2 = minimum(model_plate2.current[1, :] - model_plate2.reference[1, :])
    min_disp_y_plate2 = minimum(model_plate2.current[2, :] - model_plate2.reference[2, :])
    max_disp_z_plate2 = maximum(model_plate2.current[3, :] - model_plate2.reference[3, :])
    avg_stress_plate1 = average_components(model_plate1.stress)
    avg_stress_plate2 = average_components(model_plate2.stress)

    @test min_disp_x_plate1 ≈ -5.2251935992955734e-5 atol = 1e-12
    @test min_disp_y_plate1 ≈ -2.6616926060499257e-5 atol = 1e-12
    @test max_disp_z_plate1 ≈ 0.0001981244917465038 atol = 1e-12
    @test min_disp_x_plate2 ≈ -9.305900913551823e-5 atol = 1e-12
    @test min_disp_y_plate2 ≈ -4.800795014556214e-5 atol = 1e-12
    @test max_disp_z_plate2 ≈ 0.0007216778269555964 atol = 1e-12
    @test avg_stress_plate1 ≈
        [977595.9630637506 85772.43956109945 141163.5317646388 -49663.11299128937 3.137536668989802e6 3493.7573966998607] atol =
        1.0e1
    @test avg_stress_plate2 ≈
        [540844.912808216 -117477.17517570367 -138702.5759558706 30518.477639077442 7.223736578711799e6 12844.149940346517] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        54,
        65,
        66,
        73,
        77,
        74,
        70,
        72,
        75,
        76,
        75,
        73,
        71,
        68,
        65,
        63,
        62,
        62,
        63,
        66,
        68,
        69,
        68,
        66,
        64,
        63,
        63,
        63,
        64,
        64,
        64,
        65,
        65,
        65,
        64,
        64,
        64,
        65,
        66,
        67,
        67,
        67,
        67,
        67,
        66,
        66,
        64,
        63,
        61,
        59,
    ] atol = 0
end
