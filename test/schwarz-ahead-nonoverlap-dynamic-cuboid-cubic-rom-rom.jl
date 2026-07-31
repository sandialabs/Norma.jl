# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 Cubic ROM-ROM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-rom/copinf-operator-1.npz", "copinf-operator-1.npz"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-rom/copinf-operator-2.npz", "copinf-operator-2.npz"; force=true)
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
    rm("copinf-operator-1.npz"; force=true)
    rm("copinf-operator-2.npz"; force=true)

    min_disp_x_cuboid1 = minimum(model_cuboid1.fom_model.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.fom_model.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.fom_model.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.fom_model.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.fom_model.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.fom_model.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.fom_model.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.fom_model.stress)

    @test min_disp_x_cuboid1 ≈  -0.11889390491310992 atol = 1e-8
    @test min_disp_y_cuboid1 ≈ -0.11889390491310986 atol = 1e-8
    @test max_disp_z_cuboid1 ≈ 0.4999964535116122 atol = 1e-8
    @test min_disp_x_cuboid2 ≈ -0.11889390491310993 atol = 1e-8
    @test min_disp_y_cuboid2 ≈ -0.11889390491310993 atol = 1e-8
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-8

    @test avg_stress_cuboid1 ≈
        [32732.535671921814 32732.53567193764 5.2106047463055396e8 -5007.607013025704 -5007.6070129977115 -7.008613719871314] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [2452.229964393693 2452.229964406546 5.21067714704261e8 148.41234819314997 148.41234817819614 -22.64427872912446] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ 
       [4, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] atol = 0
end
