# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 ROM-FOM with Interface Predictor" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-rom-fom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-rom-fom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-rom-fom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-rom-fom/opinf-operator-1.npz", "opinf-operator-1.npz"; force=true)
    input_file = "cuboids.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 1.0
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_cuboid1 = subsims[1].model
    model_cuboid2 = subsims[2].model

    rm("cuboid.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("opinf-operator-1.npz"; force=true)

    min_disp_x_cuboid1 = minimum(model_cuboid1.fom_model.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.fom_model.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.fom_model.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.fom_model.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    @test min_disp_x_cuboid1 ≈ -0.16666795382949284 atol = 1e-4
    @test min_disp_y_cuboid1 ≈ -0.16666812371062245 atol = 1e-4
    @test max_disp_z_cuboid1 ≈ 0.500000342149204 atol = 1e-4
    @test min_disp_x_cuboid2 ≈ -0.16666666925334647 atol = 1e-4
    @test min_disp_y_cuboid2 ≈ -0.16666669868308254 atol = 1e-4
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-4

    @test avg_stress_cuboid1 ≈
        [-36730.699690253394 -46288.90805986865 6.666400082335727e8 -8434.913477209397 -7178.033853696234 -10275.466479944356] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [-6178.126053395448 -8111.4466235595755 6.666716435235487e8 181.54281905189106 -147.2799569560275 -902.7403261531296] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16] atol = 0
end
