# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 FOM-FOM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-fom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-fom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-fom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
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

    min_disp_x_cuboid1 = minimum(model_cuboid1.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    @test min_disp_x_cuboid1 ≈ -0.16666795382949284 atol = 1e-8
    @test min_disp_y_cuboid1 ≈ -0.16666812371062245 atol = 1e-8
    @test max_disp_z_cuboid1 ≈ 0.5000003140329964 atol = 1e-8
    @test min_disp_x_cuboid2 ≈ -0.16666667992128556 atol = 1e-8
    @test min_disp_y_cuboid2 ≈ -0.16666670983277013 atol = 1e-8
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-8
    @test avg_stress_cuboid1 ≈
        [-502.4644176257619 -591.6834172675153 6.666661928177706e8 -363.78013190731104 -303.79600091161717 -131.7465006700666] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [256.6644325077068 358.6386882842053 6.666658236719099e8 -115.20327389463883 -91.56168693370653 190.9867681440033] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16] atol = 0
end
