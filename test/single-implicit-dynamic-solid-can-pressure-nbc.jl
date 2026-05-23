# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Implicit Dynamic Solid Can Pressure BC" begin
    cp("../examples/single/implicit-dynamic-solid/pressure-bc-can/can.yaml", "can.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/pressure-bc-can/can.g", "can.g"; force=true)
    simulation = Norma.run("can.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("can.yaml"; force=true)
    rm("can.g"; force=true)
    rm("can.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ 9.042916992588757e-6 rtol = 1.0e-06
    @test avg_disp[2] ≈ 1.1433688342508286e-5 rtol = 1.0e-06
    @test avg_disp[3] ≈ -7.136350840569688e-7 rtol = 1.0e-06
    @test avg_stress[1] ≈ 3.630589420887601e8 rtol = 1.0e-06
    @test avg_stress[2] ≈ 3.747781411025189e8 rtol = 1.0e-06
    @test avg_stress[3] ≈ 1.4693877877755922e8 rtol = 1.0e-06
    @test avg_stress[4] ≈ -146272.58880190144 rtol = 1.0e-06
    @test avg_stress[5] ≈ 366103.93697419297 rtol = 1.0e-06
    @test avg_stress[6] ≈ 1.0667487852609965e6 rtol = 1.0e-06
end
