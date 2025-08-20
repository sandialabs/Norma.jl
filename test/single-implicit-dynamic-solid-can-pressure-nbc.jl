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
    rm("can.yaml")
    rm("can.g")
    rm("can.e")
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ 9.042916992588757e-6 rtol = 1.0e-06
    @test avg_disp[2] ≈ 1.1433688342508286e-5 rtol = 1.0e-06
    @test avg_disp[3] ≈ -7.136350840569688e-7 rtol = 1.0e-06
    @test avg_stress[1] ≈ 3.6281556858938783e8 rtol = 1.0e-06
    @test avg_stress[2] ≈ 3.745392952723529e8 rtol = 1.0e-06
    @test avg_stress[3] ≈ 1.474209981070922e8 rtol = 1.0e-06
    @test avg_stress[4] ≈ -134646.4674700197 rtol = 1.0e-06
    @test avg_stress[5] ≈ 356235.66342965094 rtol = 1.0e-06
    @test avg_stress[6] ≈ 1.060198839489482e6 rtol = 1.0e-06
end
