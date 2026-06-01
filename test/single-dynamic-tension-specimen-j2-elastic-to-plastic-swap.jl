# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests the ElasticToPlasticTransitionSwapCriterion on the tension-specimen
# dynamic (Newmark) J2 plasticity problem.
#
# Material parameters:
#   E = 70 GPa,  ν = 0.25,  ρ = 1000 kg/m³,  σy = 250 MPa,  H = 0.7 GPa
# Time integrator: Newmark (β = 0.25, γ = 0.5 — trapezoidal rule)
# Applied displacement (top_y node set, y-component):
#   u_y = 0.005 * 0.5 * (1 − cos(π*t))
#
# Phase-1 runs from t = 0 to t = 0.3 with Δt = 0.1.
# The criterion fires once any integration point exceeds σy * 1.05 = 262.5 MPa,
# after which the simulation swaps to phase-2 (same mesh and material).
# Phase-2 continues from the transferred state until t = 1.0.

@testset "Single Dynamic Tension Specimen J2 Elastic-to-Plastic Transition Swap" begin
    # The phase-1 YAML references the mesh as ../tension-specimen-coarse.g
    # (relative to the test working directory), so copy it one level up.
    cp(
        "../examples/ahead/single/tension-specimen/tension-specimen-coarse.g",
        "../tension-specimen-coarse.g";
        force=true,
    )
    cp(
        "../examples/ahead/single/tension-specimen/dynamic/tension-specimen-j2-elastic-to-plastic-swap.yaml",
        "tension-specimen-j2-elastic-to-plastic-swap.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/single/tension-specimen/dynamic/tension-specimen-j2-elastic-to-plastic-swap-phase2.yaml",
        "tension-specimen-j2-elastic-to-plastic-swap-phase2.yaml";
        force=true,
    )

    sim = Norma.run("tension-specimen-j2-elastic-to-plastic-swap.yaml")

    rm("tension-specimen-j2-elastic-to-plastic-swap.yaml";        force=true)
    rm("tension-specimen-j2-elastic-to-plastic-swap-phase2.yaml"; force=true)
    rm("../tension-specimen-coarse.g";                            force=true)
    rm("tension-specimen-j2-elastic-to-plastic-swap.e";           force=true)
    rm("tension-specimen-j2-elastic-to-plastic-swap-phase2.e";    force=true)

    # ── Swap fired ──────────────────────────────────────────────────────────
    @test sim.name == "tension-specimen-j2-elastic-to-plastic-swap-phase2"
    @test isempty(sim.swaps)

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ 0.3 atol = 0

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
    @test max_vm ≈ 3.650555617152344e8 rtol = 1.0e-3

    # Average y-displacement
    avg_disp = average_components(sim.integrator.displacement)
    @test avg_disp[2] ≈ 0.0008299533558821846 rtol = 1.0e-8
end
