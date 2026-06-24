# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests a 3-subdomain overlapping-Schwarz dynamic simulation (clamped bar,
# Gaussian-pulse initial condition) in which the middle subdomain starts as a
# full-order model (FOM) and is swapped mid-run to a linear OpInf reduced-order
# model (ROM), while the two outer subdomains run as ROMs throughout.
#
# Problem setup (see clamped-swap-middle.yaml):
#   domains: clamped-rom-1 (ROM), clamped-fom-2 (FOM), clamped-rom-3 (ROM)
#   Material: linear elastic, E = 1 GPa, ν = 0, ρ = 1000 kg/m³ (all subdomains)
#   Time integrator: Newmark (β = 0.25, γ = 0.5) for every subdomain
#   Initial condition (nsall, z-component): a Gaussian pulse
#       u_z(z, 0) = a * exp(-z² / (2 s²)),  a = 1.0e-3,  s = 0.02
#   All lateral faces (x±, y±) pinned in their normal direction.
#
# Swap: subdomain 2 (clamped-fom-2, a SolidMechanics FOM) is replaced by
# clamped-rom-2 (a LinearOpInfRom) once the controller's last completed step
# reaches t_swap = 3.0e-4.  This is the original single-swap scenario from
# issue #196 (FOM -> ROM, no return swap); the round-trip FOM -> ROM -> FOM
# chain from the issue is exercised separately by the multi-swap example in
# examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-multi-swap.
#
# Unlike clamped-swap-middle.yaml itself, this test only runs to 0.75 of the
# original final time (7.5e-4 instead of 1.0e-3) and disables CSV output
# (the YAML file on disk is left untouched; both overrides are applied
# in-memory to the loaded params before calling Norma.run). The swap at
# t_swap = 3.0e-4 still occurs well before the new, shorter final time, so
# the post-swap assertions below remain meaningful.
#
# 0.75 was chosen deliberately, after checking the exact solution's
# amplitude over the whole [0.5, 1.0] * T range, not just at one point:
# subdomain 2's z-window ([-1/3, 1/3]) sits in a near-cancellation valley
# for t roughly in [4.5e-4, 6.3e-4] — the two reflected wave packets pass
# through and largely cancel within that specific z-range over that
# specific time window (not just at the single point t = T/2 = 5.0e-4,
# where the cancellation is exact by symmetry; nearby times are merely
# small, not zero). Evaluating the relative-error check anywhere in that
# valley makes it pathologically sensitive: the true signal norm there is
# 3-4 orders of magnitude smaller than at the swap time itself, so even a
# small absolute discretization/coupling/ROM-projection difference — the
# same kind of difference that is negligible relative to the signal
# everywhere else — dominates the relative error. (An earlier version of
# this test used 0.6 * T = 6.0e-4, which is right at the edge of that valley
# and produced a relative error over 13, i.e. completely uninformative.)
# 0.75 * T = 7.5e-4 sits on the plateau well past the valley (the signal
# there is back to essentially the same amplitude as at the swap time),
# while still finishing comfortably before the original final time of
# 1.0e-3.
#
# The clamped bar's exact 1-D wave solution (reflections off both clamped
# ends superposed) is known in closed form; see schwarz-ahead-overlap-dynamic-
# clamped.jl for the same formula applied to a 2-subdomain, no-swap case.
# Here it is evaluated on subdomain 2 at the (shortened) final time, after
# the swap has replaced its model with a ROM, to confirm the swapped-in ROM
# is producing a physically correct field rather than e.g. stale or zeroed
# data.  Note T below is the *original* problem's final time (1.0e-3, the
# wave-reflection timescale baked into the closed-form solution), which is
# independent of how long this test actually runs the simulation for.

using LinearAlgebra
using YAML

@testset "AHeaD Overlap Dynamic Clamped 3SD FOM-ROM Swap" begin
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-single-swap"

    # ── Stage files ─────────────────────────────────────────────────────────
    # The subsim YAMLs reference their meshes as ../clamped-3sd-<k>.g
    # (relative to the test working directory), so the meshes are copied one
    # directory above, matching the convention used by the other ahead/
    # overlap/clamped tests.
    cp("../examples/ahead/overlap/clamped/clamped-3sd-1.g", "../clamped-3sd-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-2.g", "../clamped-3sd-2.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-3.g", "../clamped-3sd-3.g"; force=true)

    cp("$example_dir/clamped-fom-2.yaml", "clamped-fom-2.yaml"; force=true)
    cp("$example_dir/clamped-rom-1.yaml", "clamped-rom-1.yaml"; force=true)
    cp("$example_dir/clamped-rom-2.yaml", "clamped-rom-2.yaml"; force=true)
    cp("$example_dir/clamped-rom-3.yaml", "clamped-rom-3.yaml"; force=true)
    cp("$example_dir/linear-opinf-operator-M30-1.npz", "linear-opinf-operator-M30-1.npz"; force=true)
    cp("$example_dir/linear-opinf-operator-M30-2.npz", "linear-opinf-operator-M30-2.npz"; force=true)
    cp("$example_dir/linear-opinf-operator-M30-3.npz", "linear-opinf-operator-M30-3.npz"; force=true)

    # Load the parent YAML into a params dict (instead of running the file
    # directly) so the final time and CSV output can be overridden in-memory
    # without touching the original example file.
    input_file = "$example_dir/clamped-swap-middle.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    final_time = 0.75 * params["final time"]  # 7.5e-4; see header comment for why
    params["final time"] = final_time
    # CSV output interval = 0 disables write_stop's is_csv_step entirely, which
    # also gates the per-subsim "CSV write sidesets" branch (still set in the
    # individual subsim YAMLs) — so this one override suffices to suppress
    # all CSV output, sidesets included, without editing those files.
    params["CSV output interval"] = 0.0
    params["name"] = "clamped-swap-middle"

    sim = Norma.run(params)

    # ── Clean up ────────────────────────────────────────────────────────────
    rm("clamped-fom-2.yaml"; force=true)
    rm("clamped-rom-1.yaml"; force=true)
    rm("clamped-rom-2.yaml"; force=true)
    rm("clamped-rom-3.yaml"; force=true)
    rm("linear-opinf-operator-M30-1.npz"; force=true)
    rm("linear-opinf-operator-M30-2.npz"; force=true)
    rm("linear-opinf-operator-M30-3.npz"; force=true)
    rm("../clamped-3sd-1.g"; force=true)
    rm("../clamped-3sd-2.g"; force=true)
    rm("../clamped-3sd-3.g"; force=true)
    rm("clamped-3sd-1.e"; force=true)
    rm("clamped-3sd-fom-2.e"; force=true)
    rm("clamped-3sd-rom-2.e"; force=true)
    rm("clamped-3sd-3.e"; force=true)
    for f in filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ final_time rtol = 1.0e-9

    # ── No CSV output was written ────────────────────────────────────────────
    @test isempty(filter(f -> endswith(f, ".csv"), readdir()))

    # ── Swap fired exactly once, on slot 2 ───────────────────────────────────
    @test length(sim.swaps) == 1
    @test sim.swaps[1].applied == true
    @test sim.swaps[1].subsim_name == "clamped-fom-2"

    subsims = sim.subsims
    model1 = subsims[1].model
    model2 = subsims[2].model
    model3 = subsims[3].model

    # Slot 2 now holds the ROM replacement; slots 1 and 3 (never swapped)
    # keep their original ROM models throughout.
    @test subsims[2].name == "clamped-rom-2"
    @test model1 isa Norma.LinearOpInfRom
    @test model2 isa Norma.LinearOpInfRom
    @test model3 isa Norma.LinearOpInfRom

    # The slot ↔ name lookup is re-keyed to the post-swap occupant.
    @test sim.handle_by_name["clamped-rom-2"].id == 2
    @test subsims[2].handle.id == 2

    # ── Physical sanity check on the swapped-in subdomain ────────────────────
    # Subdomain 2's shadow FOM (model2.fom_model) holds the ROM's
    # reconstructed full-order kinematic fields after the final step.
    fom2 = model2.fom_model
    z2 = fom2.reference[3, :]
    disp2_x = fom2.displacement[1, :]
    disp2_y = fom2.displacement[2, :]
    disp2_z = fom2.displacement[3, :]

    # Exact 1-D wave solution for a Gaussian pulse on a bar clamped at both
    # ends (superposition of the two traveling-wave reflections), evaluated
    # at this test's (shortened) final time t = final_time = 7.5e-4.  T is the
    # *original* problem's final time (1.0e-3); it is a fixed parameter of the
    # closed-form solution (the reflection timescale) and does not change
    # just because this test stops the simulation early.  See
    # schwarz-ahead-overlap-dynamic-clamped.jl for the same closed-form
    # expression (disp/velo/acce) applied to a 2-subdomain, no-swap problem
    # with the same material properties and initial condition.
    c = sqrt(1.0e9 / 1.0e3)
    a = 1.0e-3
    b = 0.0
    s = 0.02
    T = 1.0e-3
    t = final_time

    n2 = length(z2)
    disp2_z_exact = zeros(Float64, n2)
    for i in 1:n2
        disp2_z_exact[i] =
            0.5 * a * (exp(-(z2[i] - c * t - b)^2 / 2 / s^2) + exp(-(z2[i] + c * t - b)^2 / 2 / s^2)) -
            0.5 * a * (exp(-(z2[i] - c * (T - t) - b)^2 / 2 / s^2) + exp(-(z2[i] + c * (T - t) - b)^2 / 2 / s^2))
    end

    disp2_z_relerr = norm(disp2_z_exact - disp2_z) / norm(disp2_z_exact)

    # Loose-ish tolerance: this checks that the post-swap ROM is tracking the
    # correct physical solution (catching e.g. a stale, zeroed, or
    # mis-projected state after the swap), not tight numerical accuracy.
    # The signal at t = 7.5e-4 is back to essentially the same amplitude as
    # at the swap time (t_swap = 3.0e-4; see the header comment for why this
    # particular time was chosen), so this tolerance is set with the same
    # order of margin as the no-swap, full-amplitude reference check in
    # schwarz-ahead-overlap-dynamic-clamped.jl (~2% there), loosened to
    # account for ROM truncation and the FOM -> ROM state-transfer projection
    # at the swap, neither of which that reference case has.
    @test disp2_z_relerr < 0.05
    @test norm(disp2_x) ≈ 0.0 atol = 1.0e-8
    @test norm(disp2_y) ≈ 0.0 atol = 1.0e-8
end

