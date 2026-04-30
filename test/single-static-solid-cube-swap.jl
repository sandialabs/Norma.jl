# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Static Solid Cube Mid-Run Swap" begin
    cp("../examples/single/static-solid/cube-swap/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cube-phase2.yaml", "cube-phase2.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cube.g", "cube.g"; force=true)
    sim = Norma.run("cube.yaml")
    rm("cube.yaml"; force=true)
    rm("cube-phase2.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    rm("cube-phase2.e"; force=true)

    # The replacement was substituted in place: name reflects the phase-2 YAML
    # and the plan list (now mirroring the replacement's, which has none) is
    # empty.
    @test sim.name == "cube-phase2"
    @test isempty(sim.swaps)
    # Simulation reached final time without failure.
    @test sim.failed == false
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-09

    # Phase-2 matches phase-1 materials, so the post-swap solution should
    # match the no-swap baseline behavior of the standard cube test.
    avg_disp = average_components(sim.integrator.displacement)
    avg_stress = average_components(sim.model.stress)
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
