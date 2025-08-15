# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Implicit Dynamic Solid Cube Pressure BC Stretch" begin
    cp("../examples/single/implicit-dynamic-solid/pressure-bc-stretch/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/pressure-bc-stretch/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.1250304399154171 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.12503043991541707 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.5000809368243189 rtol = 1.0e-06
    @test avg_stress[1] ≈ -32774.51499193255 rtol = 1.0e-06
    @test avg_stress[2] ≈ -32774.514991994016 rtol = 1.0e-06
    @test avg_stress[3] ≈ 1.00013184461106e9 rtol = 1.0e-03
    @test avg_stress[4] ≈ 1433.1287287250655 rtol = 1.0e-06
    @test avg_stress[5] ≈ 1433.1287287160662 rtol = 1.0e-06
    @test avg_stress[6] ≈ -724.3872738532795 rtol = 1.0e-06
end
