# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 FOM-ROM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-linear-elastic-fom-rom/opinf-operator-2.npz", "opinf-operator-2.npz"; force=true)
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
    rm("opinf-operator-2.npz"; force=true)

    min_disp_x_cuboid1 = minimum(model_cuboid1.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.fom_model.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.fom_model.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.fom_model.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.fom_model.stress)

    @test min_disp_x_cuboid1 ≈ -0.16666795382949284 atol = 1e-4
    @test min_disp_y_cuboid1 ≈ -0.16666812371062245 atol = 1e-4
    @test max_disp_z_cuboid1 ≈ 0.500000342149204 atol = 1e-3
    @test min_disp_x_cuboid2 ≈ -0.16666666925334647 atol = 1e-4
    @test min_disp_y_cuboid2 ≈ -0.16666669868308254 atol = 1e-4
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-8

    @test avg_stress_cuboid1 ≈
        [111458.64198295772 134655.65888941273 6.667178354513713e8 60875.15367809931 50309.05183966297 35027.64966104812] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [106278.62096209801 116525.73576283199 6.66711638084274e8 -10569.012522677203 -7721.812259974251 21044.28734074556] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [9, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16] atol = 0
end
