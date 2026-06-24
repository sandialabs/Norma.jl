# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests a 3-subdomain overlapping-Schwarz dynamic simulation (clamped bar,
# Gaussian-pulse initial condition) in which the middle subdomain starts as a
# full-order model (FOM), is swapped mid-run to a linear OpInf reduced-order
# model (ROM), and is then swapped back to a FOM — the round-trip chain from
# issue #196 — while the two outer subdomains run as ROMs throughout.
#
# Problem setup (see clamped-swap-middle.yaml):
#   domains: clamped-rom-1 (ROM), clamped-fom-2 (FOM), clamped-rom-3 (ROM)
#   Material: linear elastic, E = 1 GPa, ν = 0, ρ = 1000 kg/m³ (all subdomains)
#   Time integrator: Newmark (β = 0.25, γ = 0.5) for every subdomain
#   Initial condition (nsall, z-component): a Gaussian pulse
#       u_z(z, 0) = a * exp(-z² / (2 s²)),  a = 1.0e-3,  s = 0.02
#   All lateral faces (x±, y±) pinned in their normal direction.
#
# Swaps on subdomain 2:
#   1. clamped-fom-2 (SolidMechanics FOM) -> clamped-rom-2 (LinearOpInfRom)
#      once the controller's last completed step reaches t_swap = 3.0e-4.
#   2. clamped-rom-2 -> clamped-fom-2 again, once t_swap = 6.0e-4 is reached.
# This is the scenario from issue #196: declaring both swaps up front in the
# parent file, where the second plan's `subsim: clamped-rom-2` only becomes
# resolvable after the first swap has fired.  It requires both the
# validate_swap_plans relaxation and the apply_swap! handle_by_name
# re-keying from the issue #196 fix in src/swap.jl.
#
# Unlike clamped-swap-middle.yaml itself, this test only runs to 0.8 of the
# original final time (8.0e-4 instead of 1.0e-3) and disables CSV output
# (the YAML file on disk is left untouched; both overrides are applied
# in-memory to the loaded params before calling Norma.run). 0.8 was chosen,
# rather than the originally-suggested 0.6 = t_swap2, deliberately: Norma's
# TimeSwapCriterion fires on the step *after* the controller's last completed
# time reaches t_swap (see should_swap's comment in src/swap.jl), so a run
# that stops at exactly t_swap2 = 6.0e-4 would end one step too early for the
# second swap to ever fire — the simulation would stay on the ROM and never
# observably swap back to a FOM.  Stopping at 8.0e-4 instead gives the second
# swap 200 steps of margin to fire and settle before the final-state checks
# below, while still finishing well before the original final time of 1.0e-3
# (and clear of t = T/2 = 5.0e-4, where the closed-form exact solution used
# for the physical check below vanishes identically by symmetry).
#
# The clamped bar's exact 1-D wave solution (reflections off both clamped
# ends superposed) is known in closed form; see schwarz-ahead-overlap-dynamic-
# clamped.jl for the same formula applied to a 2-subdomain, no-swap case, and
# schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap.jl for the same
# check applied to the single-swap (FOM -> ROM only) version of this problem.
# Here it is evaluated on subdomain 2 at the (shortened) final time, after
# both swaps have run, to confirm the final FOM is holding a physically
# correct field rather than e.g. stale, zeroed, or mis-projected data.

using LinearAlgebra
using YAML

@testset "AHeaD Overlap Dynamic Clamped 3SD FOM-ROM-FOM Multi-Swap" begin
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-multi-swap"

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
    final_time = 0.8 * params["final time"]  # 8.0e-4; see header comment for why
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
    # The second swap's replacement (clamped-fom-2.yaml) reuses slot 2's
    # original "output mesh file" (clamped-3sd-fom-2.e), already written by
    # the first FOM phase; uniquify_swap_output! detects the collision and
    # renames the output mesh file (and the subsim) with a "-phase2" suffix.
    rm("clamped-3sd-fom-2-phase2.e"; force=true)
    rm("clamped-3sd-3.e"; force=true)
    for f in filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ final_time rtol = 1.0e-9

    # ── No CSV output was written ────────────────────────────────────────────
    @test isempty(filter(f -> endswith(f, ".csv"), readdir()))

    # ── Both swaps fired, in order, on slot 2 ────────────────────────────────
    @test length(sim.swaps) == 2
    @test all(p.applied for p in sim.swaps)
    @test sim.swaps[1].subsim_name == "clamped-fom-2"
    @test sim.swaps[2].subsim_name == "clamped-rom-2"

    subsims = sim.subsims
    model1 = subsims[1].model
    model2 = subsims[2].model
    model3 = subsims[3].model

    # Slot 2 has round-tripped back to a FOM, renamed with a "-phase2" suffix
    # (see the uniquify_swap_output! note in the clean-up section above);
    # slots 1 and 3 (never swapped) keep their original ROM models throughout.
    @test subsims[2].name == "clamped-fom-2-phase2"
    @test model1 isa Norma.LinearOpInfRom
    @test model2 isa Norma.SolidMechanics
    @test model3 isa Norma.LinearOpInfRom

    # The slot ↔ name lookup is re-keyed to the post-swap occupant, and the
    # intermediate name from the first swap survives as an alias for the same
    # slot (see the apply_swap! comment in src/swap.jl for why old names are
    # kept rather than removed).
    @test sim.handle_by_name["clamped-fom-2-phase2"].id == 2
    @test sim.handle_by_name["clamped-rom-2"].id == 2
    @test sim.handle_by_name["clamped-fom-2"].id == 2
    @test subsims[2].handle.id == 2

    # ── Physical sanity check on the round-tripped subdomain ─────────────────
    # Slot 2 is a SolidMechanics FOM again after the second swap, so its
    # kinematic fields are read directly (no shadow FOM indirection needed,
    # unlike the single-swap test where slot 2 ends on a ROM).
    z2 = model2.reference[3, :]
    disp2_x = model2.displacement[1, :]
    disp2_y = model2.displacement[2, :]
    disp2_z = model2.displacement[3, :]

    # Exact 1-D wave solution for a Gaussian pulse on a bar clamped at both
    # ends (superposition of the two traveling-wave reflections), evaluated
    # at this test's (shortened) final time t = final_time = 8.0e-4.  T is the
    # *original* problem's final time (1.0e-3); it is a fixed parameter of the
    # closed-form solution (the reflection timescale) and does not change
    # just because this test stops the simulation early.
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
    # Loose tolerance: this checks that the round-tripped FOM is tracking the
    # correct physical solution (catching e.g. a stale, zeroed, or
    # mis-projected state after either swap), not tight numerical accuracy.
    @test disp2_z_relerr < 0.05
    @test norm(disp2_x) ≈ 0.0 atol = 1.0e-8
    @test norm(disp2_y) ≈ 0.0 atol = 1.0e-8
end
