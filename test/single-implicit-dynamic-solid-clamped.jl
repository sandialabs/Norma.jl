# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Implicit Dynamic Solid Clamped" begin
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    simulation = Norma.run("clamped.yaml")
    integrator = simulation.integrator
    rm("clamped.yaml")
    rm("clamped.g")
    rm("clamped.e")
    max_disp = maximum_components(integrator.displacement)
    max_velo = maximum_components(integrator.velocity)
    max_acce = maximum_components(integrator.acceleration)
    min_disp = minimum_components(integrator.displacement)
    min_velo = minimum_components(integrator.velocity)
    min_acce = minimum_components(integrator.acceleration)
    @test max_disp[1] ≈ 0.0 atol = 1.0e-06
    @test max_disp[2] ≈ 0.0 atol = 1.0e-06
    @test max_disp[3] ≈ 0.00882497 rtol = 8.0e-05
    @test max_velo[1] ≈ 0.0 atol = 1.0e-06
    @test max_velo[2] ≈ 0.0 atol = 1.0e-06
    @test max_velo[3] ≈ 98.7781 rtol = 5.0e-05
    @test max_acce[1] ≈ 0.0 atol = 1.0e-06
    @test max_acce[2] ≈ 0.0 atol = 1.0e-06
    @test max_acce[3] ≈ 7.95606e6 rtol = 5.0e-04
    @test min_disp[1] ≈ 0.0 atol = 1.0e-06
    @test min_disp[2] ≈ 0.0 atol = 1.0e-06
    @test min_disp[3] ≈ 0.0 atol = 1.0e-06
    @test min_velo[1] ≈ 0.0 atol = 1.0e-06
    @test min_velo[2] ≈ 0.0 atol = 1.0e-06
    @test min_velo[3] ≈ -220.624 rtol = 2.0e-04
    @test min_acce[1] ≈ 0.0 atol = 1.0e-06
    @test min_acce[2] ≈ 0.0 atol = 1.0e-06
    @test min_acce[3] ≈ -1.65468e7 rtol = 9.0e-04
end
