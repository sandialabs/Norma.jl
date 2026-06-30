# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 Cubic ROM-FOM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-fom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-fom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-fom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-cubic-rom-fom/copinf-operator-1.npz", "copinf-operator-1.npz"; force=true)
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

    min_disp_x_cuboid1 = minimum(model_cuboid1.fom_model.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.fom_model.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.fom_model.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.fom_model.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    @test min_disp_x_cuboid1 ≈ -0.11864677526086059 atol = 1e-8
    @test min_disp_y_cuboid1 ≈ -0.11949932080200247 atol = 1e-8
    @test max_disp_z_cuboid1 ≈ 0.5003173879483614 atol = 1e-8
    @test min_disp_x_cuboid2 ≈ -0.11882175148082282 atol = 1e-8
    @test min_disp_y_cuboid2 ≈ -0.11913818555963603 atol = 1e-8
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-8
 
    @test avg_stress_cuboid1 ≈
        [155228.88994111048 -280652.1931872213 5.212788326869257e8 90135.00937846146 -42100.498370989386 -481.72595738959785] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [45457.614317274296 -73020.18351931273 5.208529395720143e8 55831.218017986845 -42042.73640307174 -1705.2349833544874] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16] atol = 0
end
