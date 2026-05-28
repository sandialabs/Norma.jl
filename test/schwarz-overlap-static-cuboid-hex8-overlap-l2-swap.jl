# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests for SchwarzOverlapL2SwapCriterion in both directions.
#
# The criterion monitors the L2 norm of the displacement mismatch at the
# overlap interface between two subdomains.  `compute overlap L2 relative error: disp`
# must be set on the relevant Schwarz overlap BC; the criterion then fires
# when the *maximum* overlap L2 error over all such BCs satisfies the
# direction condition:
#
#   direction: coarsen — fire when the error is BELOW the tolerance
#   direction: refine  — fire when the error EXCEEDS the tolerance
#
# Both tests use the same two-cuboid geometry and materials as the base
# overlap swap example (cuboids-swap/).

# ---------------------------------------------------------------------------
# Coarsen direction
# ---------------------------------------------------------------------------
# At the first time step the displacement field is zero, so the overlap L2
# error is exactly 0.0 — well below the 1.0 tolerance.  The coarsen criterion
# fires immediately before any Schwarz iteration and replaces cuboid-1 with
# cuboid-1-phase2.  Because the replacement carries the same material and
# mesh, the final stress matches the single-step uniaxial reference.

@testset "Schwarz Overlap Static Cuboid Hex8 Overlap-L2-Error Swap (coarsen)" begin
    src = "../examples/overlap/static-same-step/cuboids-swap-overlap-l2/"
    cp(src * "cuboids-coarsen.yaml",   "cuboids-coarsen.yaml";   force=true)
    cp(src * "cuboid-1.yaml",          "cuboid-1.yaml";           force=true)
    cp(src * "cuboid-1-phase2.yaml",   "cuboid-1-phase2.yaml";    force=true)
    cp(src * "cuboid-2.yaml",          "cuboid-2.yaml";           force=true)
    cp(src * "cuboid-1.g",             "cuboid-1.g";              force=true)
    cp(src * "cuboid-2.g",             "cuboid-2.g";              force=true)

    sim = Norma.run("cuboids-coarsen.yaml")

    rm("cuboids-coarsen.yaml";   force=true)
    rm("cuboid-1.yaml";          force=true)
    rm("cuboid-1-phase2.yaml";   force=true)
    rm("cuboid-2.yaml";          force=true)
    rm("cuboid-1.g";             force=true)
    rm("cuboid-2.g";             force=true)
    rm("cuboid-1.e";             force=true)
    rm("cuboid-1-phase2.e";      force=true)
    rm("cuboid-2.e";             force=true)

    # Swap must have fired: slot 1 now holds the phase-2 replacement.
    @test sim.subsims[1].name == "cuboid-1-phase2"
    @test sim.swaps[1].applied == true
    # Handle identity is preserved through the swap.
    @test sim.subsims[1].handle.id == 1
    @test sim.handle_by_name["cuboid-1"].id == 1
    # Simulation reached final time without failure.
    @test sim.failed == false
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-09
    # Same material and mesh as the original, so the uniaxial stress state is
    # unaffected: σ_zz ≈ 5e8 Pa (E = 1e9, applied strain ≈ 0.5).
    avg_stress = average_components(sim.subsims[1].model.stress)
    @test avg_stress[3] ≈ 5.0e8 rtol = 1.0e-06
end

# ---------------------------------------------------------------------------
# Refine direction
# ---------------------------------------------------------------------------
# At the first time step the displacement is zero, so the overlap L2 error is
# 0.0 and the refine criterion does NOT fire (0 is not > 1e-15).  After the
# first Schwarz solve, the nodal displacements are non-zero and the coupling
# mismatch — evaluated at the start of the second time step from the converged
# step-1 state — is a small positive number (on the order of the Schwarz
# residual, well above 1e-15).  The refine criterion fires on step 2 and the
# swap is applied before that step's Schwarz iteration.

@testset "Schwarz Overlap Static Cuboid Hex8 Overlap-L2-Error Swap (refine)" begin
    src = "../examples/overlap/static-same-step/cuboids-swap-overlap-l2/"
    cp(src * "cuboids-refine.yaml",    "cuboids-refine.yaml";     force=true)
    cp(src * "cuboid-1.yaml",          "cuboid-1.yaml";           force=true)
    cp(src * "cuboid-1-phase2.yaml",   "cuboid-1-phase2.yaml";    force=true)
    cp(src * "cuboid-2.yaml",          "cuboid-2.yaml";           force=true)
    cp(src * "cuboid-1.g",             "cuboid-1.g";              force=true)
    cp(src * "cuboid-2.g",             "cuboid-2.g";              force=true)

    sim = Norma.run("cuboids-refine.yaml")

    rm("cuboids-refine.yaml";    force=true)
    rm("cuboid-1.yaml";          force=true)
    rm("cuboid-1-phase2.yaml";   force=true)
    rm("cuboid-2.yaml";          force=true)
    rm("cuboid-1.g";             force=true)
    rm("cuboid-2.g";             force=true)
    rm("cuboid-1.e";             force=true)
    rm("cuboid-1-phase2.e";      force=true)
    rm("cuboid-2.e";             force=true)

    # Swap must have fired on step 2.
    @test sim.subsims[1].name == "cuboid-1-phase2"
    @test sim.swaps[1].applied == true
    # Handle identity is preserved through the swap.
    @test sim.subsims[1].handle.id == 1
    @test sim.handle_by_name["cuboid-1"].id == 1
    # Simulation reached final time without failure.
    @test sim.failed == false
    @test sim.controller.time ≈ 2.0 rtol = 1.0e-09
end
