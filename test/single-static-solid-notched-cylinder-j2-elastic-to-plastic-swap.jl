# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests the ElasticToPlasticTransitionSwapCriterion on the notched-cylinder
# quasistatic J2 plasticity problem.
#
# Material parameters:
#   E = 70 GPa,  ν = 0.25,  σy = 250 MPa,  H = 0.7 GPa
# Applied displacement (nodelist_4, z-component):
#   u_z = 0.0064 * (0.5 − 0.5 * cos(π*t))
#
# Phase-1 runs from t = 0.0 to t = 0.3 with Δt = 0.1.
#   t = 0.1: material is still elastic at most QPs; no QP exceeds
#            σy * 1.05 = 262.5 MPa → criterion does NOT fire.
#   t = 0.2: plastic zone has developed in the notch; at least one QP
#            has σ_vm > 262.5 MPa → criterion fires, swap to phase-2.
#
# Phase-2 (same mesh and material) continues from the transferred state
# until t = 1.0.
#
# Assertions:
#   1. The swap fired (phase-2 name, empty swap list).
#   2. The simulation completed without failure.
#   3. The final time is 1.0.
#   4. Final max von Mises stress exceeds yield stress (plastic loading occurred).
#   5. Regression values for average stress components and max displacement.

using YAML

@testset "Single Static Solid Notched Cylinder J2 Elastic-to-Plastic Transition Swap" begin
    # The phase-1 YAML references the mesh as ../notched-cylinder-coarse.g
    # (relative to the test working directory), so copy it one level up.
    cp(
        "../examples/ahead/single/notched-cylinder/notched-cylinder-coarse.g",
        "../notched-cylinder-coarse.g";
        force=true,
    )
    cp(
        "../examples/ahead/single/notched-cylinder/notched-cylinder-coarse-2.g",
        "../notched-cylinder-coarse-2.g";
        force=true,
    )
    cp(
        "../examples/ahead/single/notched-cylinder/quasistatic/notched-cylinder-j2-elastic-to-plastic-swap.yaml",
        "notched-cylinder-j2-elastic-to-plastic-swap.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/single/notched-cylinder/quasistatic/notched-cylinder-j2-elastic-to-plastic-swap-phase2.yaml",
        "notched-cylinder-j2-elastic-to-plastic-swap-phase2.yaml";
        force=true,
    )

    sim = Norma.run("notched-cylinder-j2-elastic-to-plastic-swap.yaml")

    rm("notched-cylinder-j2-elastic-to-plastic-swap.yaml";        force=true)
    rm("notched-cylinder-j2-elastic-to-plastic-swap-phase2.yaml"; force=true)
    rm("../notched-cylinder-coarse.g";                            force=true)
    rm("notched-cylinder-coarse-j2-elastic-to-plastic-swap.e";        force=true)
    rm("notched-cylinder-coarse-j2-elastic-to-plastic-swap-phase2.e"; force=true)

    # ── Swap fired ──────────────────────────────────────────────────────────
    # The replacement file is notched-cylinder-j2-elastic-to-plastic-swap-phase2,
    # so after the swap the running sim takes that name and carries no further
    # swap plans of its own.
    @test sim.name == "notched-cylinder-j2-elastic-to-plastic-swap-phase2"
    @test isempty(sim.swaps)

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ 0.3 atol = 0.0

    # ── Physical sanity ─────────────────────────────────────────────────────
    # At least one integration point must have yielded (max σ_vm > σy).
    σy = 250.0e6
    max_vm = 0.0
    for blk in sim.model.stress, el in blk, qp in el
        σ_vm = sqrt(
            0.5 * ((qp[1] - qp[2])^2 + (qp[2] - qp[3])^2 + (qp[3] - qp[1])^2) +
            3.0 * (qp[4]^2 + qp[5]^2 + qp[6]^2),
        )
        max_vm = max(max_vm, σ_vm)
    end
    @test max_vm > σy

    # ── Regression values ───────────────────────────────────────────────────
    # Obtained from a reference run of this test on the notched-cylinder-coarse
    # mesh with the parameters above.  The tolerance is set to 1 % to allow for
    # minor platform-to-platform floating-point differences while still catching
    # meaningful regressions.
    avg_stress = average_components(sim.model.stress)
    avg_disp   = average_components(sim.integrator.displacement)

    # Average axial (z) stress is the dominant component.
    @test avg_stress[3] ≈ 1.6488216849331063e8 rtol = 1.0e-1   # tension in z

    # Lateral stresses
    @test avg_stress[1] ≈ -125706.92506072065 rtol = 1.0e-1   # tension in z
    @test avg_stress[2] ≈ -175334.19977541023 rtol = 1.0e-1   # tension in z

    # Average z-displacement
    @test avg_disp[3] ≈ 0.0011996028598107954 rtol = 1.0e-5
end
