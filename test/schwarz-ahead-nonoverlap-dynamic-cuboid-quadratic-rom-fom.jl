# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 Quadratic ROM-FOM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-fom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-fom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-fom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-fom/qopinf-operator-1.npz", "qopinf-operator-1.npz"; force=true)
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
    rm("qopinf-operator-1.npz"; force=true)

    min_disp_x_cuboid1 = minimum(model_cuboid1.fom_model.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.fom_model.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.fom_model.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.displacement[3, :])
    avg_stress_cuboid1 = average_components(model_cuboid1.fom_model.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    @test min_disp_x_cuboid1 ≈ -0.11889492959831663  atol = 1e-8
    @test min_disp_y_cuboid1 ≈ -0.11889477668721234 atol = 1e-8
    @test max_disp_z_cuboid1 ≈ 0.500000342149204 atol = 1e-4
    @test min_disp_x_cuboid2 ≈ -0.11889410034198968  atol = 1e-8
    @test min_disp_y_cuboid2 ≈ -0.11889392491436034 atol = 1e-8
    @test max_disp_z_cuboid2 ≈ 1.0  atol = 1e-4
 
    @test avg_stress_cuboid1 ≈
        [-1154.1249830647962 -1071.0915433145747 5.210651730230153e8 290.8307281075439 274.5475705453412 108.54281591111055] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [173.63015319499377 212.86507083956633 5.210638738923003e8 79.80109299459701 106.78857759614577 114.1044592159878] atol =
        1.0e1
    @test sim.controller.schwarz_iters ≈ [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 14, 14, 14, 13, 12] atol = 0
end
