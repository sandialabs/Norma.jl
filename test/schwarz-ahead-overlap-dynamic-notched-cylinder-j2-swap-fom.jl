# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests the ElasticToPlasticTransitionSwapCriterion on the AHeaD Schwarz
# overlap two-domain notched-cylinder FOM-FOM problem (dynamic J2 plasticity).
#
# Configuration (examples/ahead/overlap/notched-cylinder/dynamic-j2-swap-fom):
#   Two overlapping domains (notched-cylinder-1 / notched-cylinder-2) both
#   using J2 plasticity with:
#       E  = 70 GPa,  ν = 0.25,  ρ = 1000 kg/m³
#       σy = 250 MPa, H = 0.7 GPa
#   Time integrator: Newmark (β = 0.25, γ = 0.5 — trapezoidal rule)
#   Applied displacement (+Z_top node set, z-component):
#       u_z = 0.0064 * (0.5 − 0.5 * cos(π*t))
#   Phase-1 runs from t = 0.0 to t = 0.3 with Δt = 0.1 (3 time steps).
#
# Swap criterion (tolerance = 0.05):
#   Fires when at least one integration point of notched-cylinder-1 exceeds
#   σy * (1 + 0.05) = 262.5 MPa.  The criterion is expected to trigger
#   during the first few time steps once the notch region yields, replacing
#   the running notched-cylinder-1 subsim with notched-cylinder-1-phase2.
#
# Assertions:
#   1. The swap fired: subsim[1] carries the phase-2 name and the swap plan
#      is marked applied.
#   2. The simulation completed to final time without failure.
#   3. Physical sanity: the plastic zone has developed — at least one
#      integration point in each domain exceeds σy.
#   4. Regression values for displacement and stress at the end of the run.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Notched Cylinder J2 Elastic-to-Plastic Swap FOM" begin
    example_dir = "../examples/ahead/overlap/notched-cylinder/dynamic-j2-swap-fom"
    mesh_dir    = "../examples/ahead/overlap/notched-cylinder"

    # ── Stage input files ──────────────────────────────────────────────────
    # YAML files are copied to the test working directory; mesh files are
    # expected one level up (as specified by "input mesh file: ../…" in the
    # sub-domain YAMLs).
    cp("$example_dir/notched-cylinder.yaml",          "notched-cylinder.yaml";          force=true)
    cp("$example_dir/notched-cylinder-1.yaml",        "notched-cylinder-1.yaml";        force=true)
    cp("$example_dir/notched-cylinder-2.yaml",        "notched-cylinder-2.yaml";        force=true)
    cp("$example_dir/notched-cylinder-1-phase2.yaml", "notched-cylinder-1-phase2.yaml"; force=true)
    cp("$mesh_dir/notched-cylinder-1-coarse.g",       "../notched-cylinder-1-coarse.g"; force=true)
    cp("$mesh_dir/notched-cylinder-2-coarse.g",       "../notched-cylinder-2-coarse.g"; force=true)

    # ── Run ────────────────────────────────────────────────────────────────
    sim = Norma.run("notched-cylinder.yaml")

    # ── Cleanup ────────────────────────────────────────────────────────────
    rm("notched-cylinder.yaml";          force=true)
    rm("notched-cylinder-1.yaml";        force=true)
    rm("notched-cylinder-2.yaml";        force=true)
    rm("notched-cylinder-1-phase2.yaml"; force=true)
    rm("../notched-cylinder-1-coarse.g"; force=true)
    rm("../notched-cylinder-2-coarse.g"; force=true)
    rm("notched-cylinder-1.e";           force=true)
    rm("notched-cylinder-1-phase2.e";    force=true)
    rm("notched-cylinder-2.e";           force=true)

    subsim_1 = sim.subsims[1]
    subsim_2 = sim.subsims[2]
    model_1  = subsim_1.model
    model_2  = subsim_2.model

    # ── Swap fired ─────────────────────────────────────────────────────────
    # After the criterion fires, slot 1 holds the phase-2 replacement.
    @test subsim_1.name == "notched-cylinder-1-phase2"
    @test sim.swaps[1].applied == true
    # The domain handle for the original subsim name must still resolve to
    # slot 1 (the handle is preserved through the swap).
    @test sim.subsims[1].handle.id == 1
    @test sim.handle_by_name["notched-cylinder-1"].id == 1

    # ── Completion ─────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ 0.3 atol = 0.0

    # ── Physical sanity ────────────────────────────────────────────────────
    # At least one integration point in each domain must have yielded.
    σy = 250.0e6
    function max_von_mises(model)
        mv = 0.0
        for blk in model.stress, el in blk, qp in el
            σ_vm = sqrt(
                0.5 * ((qp[1] - qp[2])^2 + (qp[2] - qp[3])^2 + (qp[3] - qp[1])^2) +
                3.0  * (qp[4]^2 + qp[5]^2 + qp[6]^2),
            )
            mv = max(mv, σ_vm)
        end
        return mv
    end
    @test max_von_mises(model_1) > σy
    @test max_von_mises(model_2) > σy

    # ── Regression values ───────────────────────────────────────────────────
    # Reference values obtained from a verified run on the notched-cylinder
    # coarse mesh with the parameters above.  Regenerated 2026-08-28 after the
    # J2 state-commit fix: state and state_old used to alias the same arrays,
    # so plastic state was committed on every residual assembly (every Newton
    # iteration and every Schwarz iteration) instead of on step acceptance;
    # this Schwarz-coupled plastic problem integrated its internal variables
    # from each iteration's trial history and the old values reflect that.  This stiff J2 + Schwarz problem
    # is highly sensitive to floating-point rounding: the converged stresses
    # vary by up to ~3 % across platforms (BLAS/CPU differences amplified
    # through the nonlinear and Schwarz iterations), with the Schwarz-coupled
    # domain-2 components most affected.  Stresses therefore use rtol=5e-2 and
    # displacements rtol=5e-3 — loose enough to be machine-independent yet
    # still tight enough to catch a real regression in the swap physics.
    avg_stress_1 = average_components(model_1.stress)
    avg_stress_2 = average_components(model_2.stress)

    # Domain 1 (notched region)
    @test avg_stress_1[1] ≈ 3.551515280362134e7 rtol = 5.0e-2   # σ_xx
    @test avg_stress_1[2] ≈ 3.577406231880717e7 rtol = 5.0e-2   # σ_yy
    @test avg_stress_1[3] ≈ 2.3412178041665545e8 rtol = 5.0e-2   # σ_zz (axial, dominant)
    @test avg_stress_1[4] ≈ 4.49667130971914e7 rtol = 5.0e-2   # σ_yz
    @test avg_stress_1[5] ≈ 4.4320888713590436e7 rtol = 5.0e-2   # σ_xz
    @test avg_stress_1[6] ≈ 2.2459202815243542e7 rtol = 5.0e-2   # σ_xy

    # Domain 2 (away from notch)
    @test avg_stress_2[1] ≈ -2.069915563053306e7 rtol = 5.0e-2   # σ_xx
    @test avg_stress_2[2] ≈ -2.050951210213198e7 rtol = 5.0e-2   # σ_yy
    @test avg_stress_2[3] ≈ 1.4799049836122042e8 rtol = 5.0e-2   # σ_zz (axial)
    @test avg_stress_2[4] ≈ 1.3364302745664774e7 rtol = 5.0e-2   # σ_yz
    @test avg_stress_2[5] ≈ 1.3137457764035271e7 rtol = 5.0e-2   # σ_xz
    @test avg_stress_2[6] ≈ 8.245908260137818e6 rtol = 5.0e-2   # σ_xy

    # Displacements
    min_disp_x_1 = minimum(model_1.displacement[1, :])
    min_disp_y_1 = minimum(model_1.displacement[2, :])
    max_disp_z_1 = maximum(model_1.displacement[3, :])
    min_disp_x_2 = minimum(model_2.displacement[1, :])
    min_disp_y_2 = minimum(model_2.displacement[2, :])
    min_disp_z_2 = minimum(model_2.displacement[3, :])

    @test min_disp_x_1 ≈ -0.0008581867167312236 rtol = 5.0e-3
    @test min_disp_y_1 ≈ -0.0008545716050116027 rtol = 5.0e-3
    @test max_disp_z_1 ≈ 0.0012634193505906832 rtol = 5.0e-3
    @test min_disp_x_2 ≈ -0.00016698327118298134 rtol = 5.0e-3
    @test min_disp_y_2 ≈ -0.00016694097282405826 rtol = 5.0e-3
    @test min_disp_z_2 ≈ 0.0009764701314749769 rtol = 5.0e-3

    # Schwarz iteration counts should be non-trivial (coupling active)
    @test all(sim.controller.schwarz_iters .>= 0)
    @test any(sim.controller.schwarz_iters .> 0)
end
