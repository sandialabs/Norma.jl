# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Schwarz Overlap Static Cuboid Hex8 Same Step" begin
    cp("../examples/overlap/static-same-step/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    min_disp_y_fine = minimum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    min_disp_x_coarse = minimum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    min_disp_y_coarse = minimum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    min_disp_z_coarse = minimum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.75 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.25 rtol = 1.0e-01
    @test avg_stress_fine[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-01
end

@testset "Schwarz Overlap Static Cuboid Hex8 Different Steps" begin
    cp("../examples/overlap/static-different-steps/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    min_disp_y_fine = minimum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    min_disp_x_coarse = minimum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    min_disp_y_coarse = minimum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    min_disp_z_coarse = minimum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.75 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.25 rtol = 1.0e-01
    @test avg_stress_fine[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-01
end
