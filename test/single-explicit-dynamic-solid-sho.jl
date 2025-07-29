# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Explicit Dynamic Solid Sho" begin
    cp("../examples/single/explicit-dynamic-solid/sho/sho.yaml", "sho.yaml"; force=true)
    cp("../examples/single/explicit-dynamic-solid/sho/sho.g", "sho.g"; force=true)
    simulation = Norma.run("sho.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("sho.yaml")
    rm("sho.g")
    rm("sho.e")
    displacement = integrator.displacement[end]
    velocity = integrator.velocity[end]
    acceleration = integrator.acceleration[end]
    @test displacement ≈ 0.0 atol = 2.0e-03
    @test velocity ≈ -1.0 rtol = 2.0e-03
    @test acceleration ≈ 0.0 atol = 2.0e-03
end
