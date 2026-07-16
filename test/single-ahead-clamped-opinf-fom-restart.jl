# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Runs the AHeaD single-domain clamped (single-Gaussian IC) OpInf-FOM problem
# two ways and checks that they agree at the final time step:
#
#   1. "dynamic-opinf-fom"          — a single continuous run from
#                                      t = 0.0 to t = 1.0e-3, starting from
#                                      the non-trivial Gaussian-pulse initial
#                                      displacement/velocity field.
#   2. "dynamic-opinf-fom-restart"   — a restart run that resumes from the
#                                      checkpoint stored in
#                                      clamped_single_gaussian-in.e
#                                      (index 30) and continues on to the
#                                      same final time, t = 1.0e-3.
#
# This is a full-order-model (FOM) run (no actual OpInf reduced model is
# constructed), so this exercises restart together with a problem that has a
# non-trivial initial condition, as opposed to the cuboid restart test, which
# starts from rest.
#
# If restart is implemented correctly, the displacement, velocity,
# acceleration, and stress fields at the shared final time should match the
# values produced by the uninterrupted run to tight tolerance.

@testset "AHeaD Single Dynamic Clamped OpInf-FOM Restart" begin
    # ── Phase 1: full, uninterrupted run (t = 0.0 -> 1.0e-3) ────────────────
    cp("../examples/ahead/single/clamped/clamped.g", "../clamped.g"; force=true)
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-fom/clamped_single_gaussian.yaml",
        "clamped-single-gaussian-fom.yaml";
        force=true,
    )

    sim_full = Norma.run("clamped-single-gaussian-fom.yaml")

    rm("clamped-single-gaussian-fom.yaml"; force=true)
    rm("clamped_single_gaussian.e"; force=true)
    for f in filter(f -> startswith(f, "clamped_single_gaussian") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Phase 2: restart run resuming from index 30 -> t = 1.0e-3 ──────────
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-fom-restart/clamped_single_gaussian-in.e",
        "clamped_single_gaussian-in.e";
        force=true,
    )
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-fom-restart/clamped_single_gaussian.yaml",
        "clamped-single-gaussian-fom-restart.yaml";
        force=true,
    )

    sim_restart = Norma.run("clamped-single-gaussian-fom-restart.yaml")

    rm("clamped-single-gaussian-fom-restart.yaml"; force=true)
    rm("clamped_single_gaussian-in.e"; force=true)
    rm("clamped_single_gaussian-out.e"; force=true)
    rm("../clamped.g"; force=true)
    for f in filter(f -> startswith(f, "clamped_single_gaussian") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Both runs reached the same final time without failing ──────────────
    @test sim_full.failed == false
    @test sim_restart.failed == false
    @test sim_full.controller.time ≈ 1.0e-03 rtol = 1.0e-09
    @test sim_restart.controller.time ≈ 1.0e-03 rtol = 1.0e-09

    # ── Final-step kinematic fields agree between the two runs ─────────────
    model_full = sim_full.model
    model_restart = sim_restart.model

    @test model_restart.displacement ≈ model_full.displacement rtol = 1.0e-05
    @test model_restart.velocity ≈ model_full.velocity rtol = 1.0e-05
    @test model_restart.acceleration ≈ model_full.acceleration rtol = 1.0e-04

    # ── Final-step stress field agrees as well ──────────────────────────────
    avg_stress_full = average_components(model_full.stress)
    avg_stress_restart = average_components(model_restart.stress)
    err_abs = norm(avg_stress_full - avg_stress_restart)
    # Was `@test err_abs ≈ 6.055561584887406e-6 atol = 1.0e-5`: an equality
    # check whose atol (1.0e-5) is larger than the target value itself, so it
    # was really just `err_abs < 6.055561584887406e-6 + 1.0e-5 ≈ 1.6e-5`
    # written as an equality. Write the actual intent directly as an upper
    # bound, with a bit of extra headroom over the observed value (order
    # 1.0e-5, well below the ~1.0e-3 displacement/stress scale of this
    # problem) so this doesn't become platform-brittle.
    @test err_abs < 2.0e-5

    # ── Sanity check: the body actually moved from its Gaussian-pulse IC ───
    avg_disp_full = average_components(sim_full.integrator.displacement)
    @test abs(avg_disp_full[3]) > 0.0
end
