# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Runs the AHeaD single-domain clamped (single-Gaussian IC) OpInf-ROM problem
# two ways and checks that they agree at the final time step:
#
#   1. "dynamic-opinf-rom"          — a single continuous run from
#                                      t = 0.0 to t = 1.0e-3, starting from
#                                      the non-trivial Gaussian-pulse initial
#                                      displacement/velocity field, evolved
#                                      with a linear OpInf reduced-order
#                                      model.
#   2. "dynamic-opinf-rom-restart"  — a restart run that resumes from the
#                                      checkpoint stored in
#                                      clamped_single_gaussian-in.e
#                                      (index 30) and continues on to the
#                                      same final time, t = 1.0e-3, using the
#                                      same ROM.
#
# Restart for ROM models is layered on top of FOM restart: the model's
# internal FOM displacement/velocity are restored from the snapshot exactly
# like a standalone FOM restart (SolidMechanics()), and that restored state
# is then projected onto the reduced basis (see
# apply_ics(::Parameters, ::RomModel, ...) in opinf_ics_bcs.jl). If this is
# implemented correctly, the reduced state/velocity — and the full-order
# fields reconstructed from them — should match the values produced by the
# uninterrupted run to tight tolerance.
#
# CSV output is disabled for both runs (params["CSV output interval"] = 0)
# simply to keep the test directory clean; it has no bearing on correctness,
# and Exodus output (which is what actually triggers the ROM's
# reduced-state -> full-order-field reconstruction at the final time step)
# is left untouched.

using YAML

@testset "AHeaD Single Dynamic Clamped OpInf-ROM Restart" begin
    # ── Phase 1: full, uninterrupted ROM run (t = 0.0 -> 1.0e-3) ────────────
    cp("../examples/ahead/single/clamped/clamped.g", "../clamped.g"; force=true)
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-rom/clamped_single_gaussian.yaml",
        "clamped-single-gaussian-rom.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-rom/linear-opinf-operator.npz",
        "linear-opinf-operator.npz";
        force=true,
    )

    # Load into a params dict (instead of running the file directly) so CSV
    # output can be turned off without touching the original example file.
    params_full = YAML.load_file("clamped-single-gaussian-rom.yaml"; dicttype=Norma.Parameters)
    params_full["CSV output interval"] = 0.0
    params_full["name"] = "clamped-single-gaussian-rom"

    sim_full = Norma.run(params_full)

    rm("clamped-single-gaussian-rom.yaml"; force=true)
    rm("clamped_single_gaussian.e"; force=true)
    rm("linear-opinf-operator.npz"; force=true)
    for f in filter(f -> startswith(f, "clamped_single_gaussian") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Phase 2: restart run resuming from index 30 -> t = 1.0e-3 ──────────
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-rom-restart/clamped_single_gaussian-in.e",
        "clamped_single_gaussian-in.e";
        force=true,
    )
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-rom-restart/clamped_single_gaussian.yaml",
        "clamped-single-gaussian-rom-restart.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/single/clamped/dynamic-opinf-rom-restart/linear-opinf-operator.npz",
        "linear-opinf-operator.npz";
        force=true,
    )

    params_restart = YAML.load_file("clamped-single-gaussian-rom-restart.yaml"; dicttype=Norma.Parameters)
    params_restart["CSV output interval"] = 0.0
    params_restart["name"] = "clamped-single-gaussian-rom-restart"

    sim_restart = Norma.run(params_restart)

    rm("clamped-single-gaussian-rom-restart.yaml"; force=true)
    rm("clamped_single_gaussian-in.e"; force=true)
    rm("clamped_single_gaussian-out.e"; force=true)
    rm("linear-opinf-operator.npz"; force=true)
    rm("../clamped.g"; force=true)
    for f in filter(f -> startswith(f, "clamped_single_gaussian") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Both runs reached the same final time without failing ──────────────
    @test sim_full.failed == false
    @test sim_restart.failed == false
    @test sim_full.controller.time ≈ 1.0e-03 rtol = 1.0e-09
    @test sim_restart.controller.time ≈ 1.0e-03 rtol = 1.0e-09

    model_full = sim_full.model
    model_restart = sim_restart.model

    # ── Final-step reconstructed full-order kinematic fields also agree ────
    # (model.fom_model.displacement/velocity/acceleration are refreshed from
    # reduced_state/reduced_velocity at every Exodus output step; both runs
    # share the same output interval and both land exactly on the final
    # time, so these are current as of t = 1.0e-3.)
    @test model_restart.fom_model.displacement ≈ model_full.fom_model.displacement rtol = 1.0e-05
    @test model_restart.fom_model.velocity ≈ model_full.fom_model.velocity rtol = 1.0e-05
    @test model_restart.fom_model.acceleration ≈ model_full.fom_model.acceleration rtol = 1.0e-04

    # ── Final-step stress field agrees as well ──────────────────────────────
    avg_stress_full = average_components(model_full.fom_model.stress)
    avg_stress_restart = average_components(model_restart.fom_model.stress)
    err_abs = norm(avg_stress_full - avg_stress_restart)
    @test err_abs ≈ 6.682167652892721e-7 atol = 1.0e-6

    # ── Sanity check: the body actually moved from its Gaussian-pulse IC ───
    avg_disp_full = average_components(vec(model_full.fom_model.displacement))
    @test abs(avg_disp_full[3]) > 0.0
end
