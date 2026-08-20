# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Schwarz Nonoverlap Static Cuboid Hex8 Aitken Defaults" begin
    src = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/"
    cp(src * "cuboids-aitken.yaml", "cuboids-aitken.yaml"; force=true)
    cp(src * "cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp(src * "cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp(src * "cuboid-1.g", "cuboid-1.g"; force=true)
    cp(src * "cuboid-2.g", "cuboid-2.g"; force=true)

    sim = Norma.run("cuboids-aitken.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model

    rm("cuboids-aitken.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)

    @test sim.controller.relaxation_method == :aitken_recursive
    @test sim.failed == false

    # Aitken must reach the same physical solution as classical relaxation.
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.5 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.5 rtol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
end

@testset "Schwarz Nonoverlap Static Cuboid Hex8 Aitken Specified Theta0 and N0" begin
    src = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/"
    cp(src * "cuboids-aitken.yaml", "cuboids-aitken.yaml"; force=true)
    cp(src * "cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp(src * "cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp(src * "cuboid-1.g", "cuboid-1.g"; force=true)
    cp(src * "cuboid-2.g", "cuboid-2.g"; force=true)

    write("cuboids-aitken.yaml", replace(read("cuboids-aitken.yaml", String),
    "relaxation: aitken recursive" =>
    "relaxation: aitken recursive\nrelaxation parameter: 0.1\naitken N0 parameter: 2"))
    
    sim = Norma.run("cuboids-aitken.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model

    rm("cuboids-aitken.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)

    @test sim.controller.relaxation_method == :aitken_recursive
    @test sim.failed == false

    # Below N0 the relaxation factor is the input theta, not 1: this run sets
    # `relaxation parameter: 0.1` and `aitken N0 parameter: 2`, so the first two
    # sweeps of every step must relax by 0.1 (issue #218). A key with no
    # relaxation history reaches the same branch those sweeps take.
    fresh_key = (99, 98, 1)
    zero_iterate = zeros(3)
    for iter in 0:1
        @test Norma.relaxation_aitken_recursive_theta!(
            sim.controller, fresh_key, 1, iter, zero_iterate, zero_iterate) == 0.1
    end

    # Aitken must reach the same physical solution as classical relaxation.
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.5 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.5 rtol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
end
