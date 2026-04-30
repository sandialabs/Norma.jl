# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

@testset "Schwarz Overlap Static Cuboid Hex8 Mid-Run Swap" begin
    src = "../examples/overlap/static-same-step/cuboids-swap/"
    cp(src * "cuboids.yaml", "cuboids.yaml"; force=true)
    cp(src * "cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp(src * "cuboid-1-phase2.yaml", "cuboid-1-phase2.yaml"; force=true)
    cp(src * "cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp(src * "cuboid-1.g", "cuboid-1.g"; force=true)
    cp(src * "cuboid-2.g", "cuboid-2.g"; force=true)

    sim = Norma.run("cuboids.yaml")

    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-1-phase2.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-1-phase2.e"; force=true)
    rm("cuboid-2.e"; force=true)

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
    # match the no-swap baseline stress of 5e8.
    avg_stress_fine = average_components(sim.subsims[1].model.stress)
    @test avg_stress_fine[3] ≈ 5.0e8 rtol = 1.0e-06
end
