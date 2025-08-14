# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Static Solid Cube with Spatially-Distributed DBC" begin
    cp("../examples/single/static-solid/cube/cube-sd-dbc.yaml", "cube-sd-dbc.yaml"; force=true)
    cp("../examples/single/static-solid/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-sd-dbc.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-sd-dbc.yaml")
    rm("cube.g")
    rm("cube-sd-dbc.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ 0.0 atol = 1.0e-06
    @test avg_disp[2] ≈ -0.15807560137457044 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.3810708848289451 rtol = 1.0e-06
    @test avg_stress[1] ≈ 2.6666666666666663e8 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0666666666666669e9 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end
