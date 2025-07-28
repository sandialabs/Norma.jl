# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Torsion HEX8-HEX8" begin
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion.yaml", "torsion.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion-1.yaml", "torsion-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion-2.yaml", "torsion-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/torsion-1.g", "../torsion-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/torsion-2.g", "../torsion-2.g"; force=true)
    input_file = "torsion.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-6
    params["final time"] = 5.0e-5
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_torsion1 = subsims[1].model
    model_torsion2 = subsims[2].model

    rm("torsion.yaml")
    rm("torsion-1.yaml")
    rm("torsion-2.yaml")
    rm("../torsion-1.g")
    rm("../torsion-2.g")
    rm("torsion-1.e")
    rm("torsion-2.e")

    min_disp_x_torsion1 = minimum(model_torsion1.current[1, :] - model_torsion1.reference[1, :])
    min_disp_y_torsion1 = minimum(model_torsion1.current[2, :] - model_torsion1.reference[2, :])
    max_disp_z_torsion1 = maximum(model_torsion1.current[3, :] - model_torsion1.reference[3, :])
    min_disp_x_torsion2 = minimum(model_torsion2.current[1, :] - model_torsion2.reference[1, :])
    min_disp_y_torsion2 = minimum(model_torsion2.current[2, :] - model_torsion2.reference[2, :])
    max_disp_z_torsion2 = maximum(model_torsion2.current[3, :] - model_torsion2.reference[3, :])
    avg_stress_torsion1 = average_components(model_torsion1.stress)
    avg_stress_torsion2 = average_components(model_torsion2.stress)

    @test min_disp_x_torsion1 ≈ -0.005190915825901907 atol = 1e-12
    @test min_disp_y_torsion1 ≈ -0.0051909158259019 atol = 1e-12
    @test max_disp_z_torsion1 ≈ 0.00011253980579573053 atol = 1e-12
    @test min_disp_x_torsion2 ≈ -0.005191500055775754 atol = 1e-12
    @test min_disp_y_torsion2 ≈ -0.005191500055775934 atol = 1e-12
    @test max_disp_z_torsion2 ≈ 2.234108562854109e-5 atol = 1e-12
    @test avg_stress_torsion1 ≈
        [1.7705537665662847e6 1.7706066742399093e6 736180.8008143709 -238.39535563270692 365.22587280095036 -2.557130848843019] atol =
        1.0e1
    @test avg_stress_torsion2 ≈
        [1.4257262240109073e6 1.4257762858988189e6 549384.8122747836 -69.31989121219449 81.6909752368978 -20.893920474565377] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
    ] atol = 0
end
