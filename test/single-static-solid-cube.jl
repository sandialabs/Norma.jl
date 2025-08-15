# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Static Solid Cube" begin
    cp("../examples/single/static-solid/cube/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

@testset "Single Static Solid Cube Line Search" begin
    cp("../examples/single/static-solid/cube/cube-line-search.yaml", "cube-line-search.yaml"; force=true)
    cp("../examples/single/static-solid/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-line-search.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-line-search.yaml")
    rm("cube.g")
    rm("cube.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

@testset "Single Static Solid Cube Steepest Descent" begin
    cp("../examples/single/static-solid/cube/cube-steepest-descent.yaml", "cube-steepest-descent.yaml"; force=true)
    cp("../examples/single/static-solid/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-steepest-descent.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-steepest-descent.yaml")
    rm("cube.g")
    rm("cube.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 3.0
    @test avg_stress[2] ≈ 0.0 atol = 3.0
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 3.0
    @test avg_stress[5] ≈ 0.0 atol = 3.0
    @test avg_stress[6] ≈ 0.0 atol = 3.0
end
