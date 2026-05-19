# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests for StressRecoverySwapCriterion with direction: coarsen — the swap
# fires when the relative Frobenius-norm difference between the lumped and
# consistent L2-projected stress fields falls BELOW the tolerance.
#
# For a uniform hex8 cube under uniaxial loading the stress field is
# elementally constant; both projection methods recover it exactly, so their
# relative difference is at the level of floating-point rounding (~1e-14),
# far below the 1 % tolerance.  The coarsen criterion therefore fires on the
# first evaluation, and the simulation runs to completion as the phase-2
# (replacement) model.

@testset "Single Static Solid Cube Mid-Run Stress Recovery-based Swap" begin
    cp("../examples/single/static-solid/cube-swap/standard/cube.g", "cube.g"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cube-sr-swap.yaml", "cube-sr-swap.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/standard/cube-sr-swap-phase2.yaml", "cube-sr-swap-phase2.yaml"; force=true)

    sim = Norma.run("cube-sr-swap.yaml")

    rm("cube-sr-swap.yaml";        force=true)
    rm("cube-sr-swap-phase2.yaml"; force=true)
    rm("cube.g";                   force=true)
    rm("cube-sr-swap.e";           force=true)
    rm("cube-sr-swap-phase2.e";    force=true)

    # The swap must have fired: name reflects the phase-2 YAML and the
    # plan list (mirroring the replacement, which has none) is empty.
    @test sim.name == "cube-sr-swap-phase2"
    @test isempty(sim.swaps)

    # Simulation ran to completion without failure.
    @test sim.failed == false
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-09

    # Solution quality: same uniaxial stress state as the vanilla cube test.
    avg_disp   = average_components(sim.integrator.displacement)
    avg_stress = average_components(sim.model.stress)
    @test avg_disp[1]   ≈  -0.125   rtol = 1.0e-06
    @test avg_disp[2]   ≈  -0.125   rtol = 1.0e-06
    @test avg_disp[3]   ≈   0.500   rtol = 1.0e-06
    @test avg_stress[1] ≈   0.0     atol = 1.0e-06
    @test avg_stress[2] ≈   0.0     atol = 1.0e-06
    @test avg_stress[3] ≈   1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈   0.0     atol = 1.0e-06
    @test avg_stress[5] ≈   0.0     atol = 1.0e-06
    @test avg_stress[6] ≈   0.0     atol = 1.0e-06
end

