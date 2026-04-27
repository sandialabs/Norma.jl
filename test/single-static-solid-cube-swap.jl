# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Static Solid Cube Mid-Run Swap" begin
    cp("../examples/single/static-solid/cube-swap/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cuboid-1-phase2.yaml", "cuboid-1-phase2.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cube.g", "cube.g"; force=true)
    sim = Norma.run("cube.yaml")
    rm("cube.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-1-phase2.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    rm("cube-phase2.e"; force=true)

    # Swap occurred: slot 1 now holds the phase-2 replacement.
    @test sim.subsims[1].name == "cuboid-1-phase2"
    @test sim.swaps[1].applied == true
    # Handle identity preserved through the swap.
    @test sim.subsims[1].handle.id == 1
    @test sim.handle_by_name["cuboid-1"].id == 1
    # Simulation reached final time without failure.
    @test sim.failed == false
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-09

    # Phase-2 matches phase-1 materials, so the post-swap solution should
    # match the no-swap baseline stress of 1e9.
    avg_stress_fine = average_components(sim.subsims[1].model.stress)
    @test avg_stress_fine[3] ≈ 1.0e9 rtol = 1.0e-06
    
    # Check regular test criteria like in non-swap cube test 
    avg_disp = average_components(sim.subsims[1].integrator.displacement)
    avg_stress = average_components(sim.subsims[1].model.stress)
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

