# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests a 3-subdomain overlapping-Schwarz dynamic simulation (clamped bar,
# Gaussian-pulse initial condition) in which the two OUTER subdomains start
# as linear OpInf reduced-order models (ROM) and are swapped mid-run to
# full-order models (FOM), while the middle subdomain runs as a FOM
# throughout. This is the companion to schwarz-ahead-overlap-dynamic-clamped-
# 3sd-fom-rom-swap.jl (which swaps the *middle* subdomain FOM -> ROM); here
# it is the two *outer* subdomains that swap, ROM -> FOM, using
# clamped-swap-left-right.yaml. This is the exact reproduction scenario from
# issue #197 (artifacts in the acceleration field after a ROM-involved swap).
#
# Problem setup (see clamped-swap-left-right.yaml):
#   domains: clamped-rom-1 (ROM), clamped-fom-2 (FOM), clamped-rom-3 (ROM)
#   Material: linear elastic, E = 1 GPa, ν = 0, ρ = 1000 kg/m³ (all subdomains)
#   Time integrator: Newmark (β = 0.25, γ = 0.5) for every subdomain
#   Initial condition (nsall, z-component): a Gaussian pulse
#       u_z(z, 0) = a * exp(-z² / (2 s²)),  a = 1.0e-3,  s = 0.02
#   All lateral faces (x±, y±) pinned in their normal direction.
#
# Swaps: subdomains 1 and 3 (clamped-rom-1, clamped-rom-3 — both
# LinearOpInfRom) are replaced by clamped-fom-1 / clamped-fom-3 (both
# SolidMechanics FOMs) once the controller's last completed step reaches
# t_swap = 3.0e-4.
#
# Unlike clamped-swap-left-right.yaml itself, this test only runs to 0.75 of
# the original final time (7.5e-4 instead of 1.0e-3) and disables CSV output
# (the YAML file on disk is left untouched; both overrides are applied
# in-memory to the loaded params before calling Norma.run), matching the
# choice made in schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap.jl —
# see that test's header comment for why 0.75 * T sits on the signal plateau
# rather than in the near-cancellation valley around [4.5e-4, 6.3e-4].
#
# The clamped bar's exact 1-D wave solution (reflections off both clamped
# ends superposed) is known in closed form; see schwarz-ahead-overlap-dynamic-
# clamped.jl for the same formula (disp/velo/acce) applied to a 2-subdomain,
# no-swap case. Here it is evaluated on subdomain 1 (and, by symmetry,
# subdomain 3) at the (shortened) final time, after the swap has replaced
# their models with FOMs, to confirm the swapped-in FOM is producing a
# physically correct *acceleration* field — not just displacement — rather
# than carrying forward a stale or zeroed acceleration from the swap. T below
# is the *original* problem's final time (1.0e-3, the wave-reflection
# timescale baked into the closed-form solution), independent of how long
# this test actually runs the simulation for.

using LinearAlgebra
using YAML

@testset "AHeaD Overlap Dynamic Clamped 3SD ROM-FOM Swap" begin
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-single-swap"

    # ── Stage files ─────────────────────────────────────────────────────────
    # The subsim YAMLs reference their meshes as ../clamped-3sd-<k>.g
    # (relative to the test working directory), so the meshes are copied one
    # directory above, matching the convention used by the other ahead/
    # overlap/clamped tests.
    cp("../examples/ahead/overlap/clamped/clamped-3sd-1.g", "../clamped-3sd-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-2.g", "../clamped-3sd-2.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-3.g", "../clamped-3sd-3.g"; force=true)

    cp("$example_dir/clamped-fom-1.yaml", "clamped-fom-1.yaml"; force=true)
    cp("$example_dir/clamped-fom-2.yaml", "clamped-fom-2.yaml"; force=true)
    cp("$example_dir/clamped-fom-3.yaml", "clamped-fom-3.yaml"; force=true)
    cp("$example_dir/clamped-rom-1.yaml", "clamped-rom-1.yaml"; force=true)
    cp("$example_dir/clamped-rom-3.yaml", "clamped-rom-3.yaml"; force=true)
    cp("$example_dir/linear-opinf-operator-M30-1.npz", "linear-opinf-operator-M30-1.npz"; force=true)
    cp("$example_dir/linear-opinf-operator-M30-3.npz", "linear-opinf-operator-M30-3.npz"; force=true)

    # Load the parent YAML into a params dict (instead of running the file
    # directly) so the final time and CSV output can be overridden in-memory
    # without touching the original example file.
    input_file = "$example_dir/clamped-swap-left-right.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    final_time = 0.75 * params["final time"]  # 7.5e-4; see header comment for why
    params["final time"] = final_time
    # CSV output interval = 0 disables write_stop's is_csv_step entirely, which
    # also gates the per-subsim "CSV write sidesets" branch (still set in the
    # individual subsim YAMLs) — so this one override suffices to suppress
    # all CSV output, sidesets included, without editing those files.
    params["CSV output interval"] = 0.0
    params["name"] = "clamped-swap-left-right"

    sim = Norma.run(params)

    # ── Clean up ────────────────────────────────────────────────────────────
    rm("clamped-fom-1.yaml"; force=true)
    rm("clamped-fom-2.yaml"; force=true)
    rm("clamped-fom-3.yaml"; force=true)
    rm("clamped-rom-1.yaml"; force=true)
    rm("clamped-rom-3.yaml"; force=true)
    rm("linear-opinf-operator-M30-1.npz"; force=true)
    rm("linear-opinf-operator-M30-3.npz"; force=true)
    rm("../clamped-3sd-1.g"; force=true)
    rm("../clamped-3sd-2.g"; force=true)
    rm("../clamped-3sd-3.g"; force=true)
    rm("clamped-3sd-1.e"; force=true)
    rm("clamped-3sd-fom-2.e"; force=true)
    rm("clamped-3sd-3.e"; force=true)
    # clamped-fom-1.yaml / clamped-fom-3.yaml reuse their ROM phase's output
    # filename (clamped-3sd-1.e / clamped-3sd-3.e), already written before the
    # swap fires, so uniquify_swap_output! renames the post-swap output file
    # with a "-phase2" suffix (see the header comment above the name checks
    # below for details).
    rm("clamped-3sd-1-phase2.e"; force=true)
    rm("clamped-3sd-3-phase2.e"; force=true)
    for f in filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ final_time rtol = 1.0e-9

    # ── No CSV output was written ────────────────────────────────────────────
    @test isempty(filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir()))

    # ── Both swaps fired, on slots 1 and 3 ───────────────────────────────────
    @test length(sim.swaps) == 2
    @test all(p.applied for p in sim.swaps)
    @test sim.swaps[1].subsim_name == "clamped-rom-1"
    @test sim.swaps[2].subsim_name == "clamped-rom-3"

    subsims = sim.subsims
    model1 = subsims[1].model
    model2 = subsims[2].model
    model3 = subsims[3].model

    # Slots 1 and 3 now hold the FOM replacements; slot 2 (never swapped)
    # keeps its original FOM model throughout. Both replacement YAMLs
    # (clamped-fom-1.yaml, clamped-fom-3.yaml) reuse their respective ROM
    # phase's "output mesh file" name (clamped-3sd-1.e / clamped-3sd-3.e,
    # already written by the time the swap fires), so uniquify_swap_output!
    # renames the replacement subsim (and its output file) with a "-phase2"
    # suffix — see the apply_swap! / uniquify_swap_output! comments in
    # src/swap.jl and the analogous round-trip case in
    # schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-fom-multi-swap.jl.
    @test subsims[1].name == "clamped-fom-1-phase2"
    @test subsims[3].name == "clamped-fom-3-phase2"
    @test model1 isa Norma.SolidMechanics
    @test model2 isa Norma.SolidMechanics
    @test model3 isa Norma.SolidMechanics

    # The slot ↔ name lookup is re-keyed to the post-swap occupant, and the
    # original ROM name survives as an alias for the same slot (see the
    # apply_swap! comment in src/swap.jl for why old names are kept rather
    # than removed).
    @test sim.handle_by_name["clamped-fom-1-phase2"].id == 1
    @test sim.handle_by_name["clamped-fom-3-phase2"].id == 3
    @test sim.handle_by_name["clamped-rom-1"].id == 1
    @test sim.handle_by_name["clamped-rom-3"].id == 3
    @test subsims[1].handle.id == 1
    @test subsims[3].handle.id == 3

    # ── Physical sanity check on the swapped-in subdomains ───────────────────
    # Subdomains 1 and 3 are SolidMechanics FOMs again after their swaps, so
    # their kinematic fields — crucially including ACCELERATION, the field
    # reported as corrupted in issue #197 — are read directly.
    c = sqrt(1.0e9 / 1.0e3)
    a = 1.0e-3
    b = 0.0
    s = 0.02
    T = 1.0e-3
    t = final_time

    function exact_solution(z::Vector{Float64})
        n = length(z)
        disp_z = zeros(Float64, n)
        velo_z = zeros(Float64, n)
        acce_z = zeros(Float64, n)
        for i in 1:n
            disp_z[i] =
                0.5 * a * (exp(-(z[i] - c * t - b)^2 / 2 / s^2) + exp(-(z[i] + c * t - b)^2 / 2 / s^2)) -
                0.5 * a * (exp(-(z[i] - c * (T - t) - b)^2 / 2 / s^2) + exp(-(z[i] + c * (T - t) - b)^2 / 2 / s^2))
            velo_z[i] =
                0.5 * c * a / s^2 * (
                    (z[i] - c * t - b) * exp(-(z[i] - c * t - b)^2 / 2 / s^2) -
                    (z[i] + c * t - b) * exp(-(z[i] + c * t - b)^2 / 2 / s^2)
                ) +
                0.5 * c * a / s^2 * (
                    (z[i] - c * (T - t) - b) * exp(-(z[i] - c * (T - t) - b)^2 / 2 / s^2) -
                    (z[i] + c * (T - t) - b) * exp(-(z[i] + c * (T - t) - b)^2 / 2 / s^2)
                )
            acce_z[i] =
                0.5 * a * (
                    -c^2 / s^2 * exp(-0.5 * (z[i] - c * t - b)^2 / s^2) +
                    c^2 / s^4 * (z[i] - c * t - b)^2 * exp(-0.5 * (z[i] - c * t - b)^2 / s^2) -
                    c^2 / s^2 * exp(-0.5 * (z[i] + c * t - b)^2 / s^2) +
                    c^2 / s^4 * (z[i] + c * t - b)^2 * exp(-0.5 * (z[i] + c * t - b)^2 / s^2)
                ) -
                0.5 * a * (
                    -c^2 / s^2 * exp(-0.5 * (z[i] - c * (T - t) - b)^2 / s^2) +
                    c^2 / s^4 * (z[i] - c * (T - t) - b)^2 * exp(-0.5 * (z[i] - c * (T - t) - b)^2 / s^2) -
                    c^2 / s^2 * exp(-0.5 * (z[i] + c * (T - t) - b)^2 / s^2) +
                    c^2 / s^4 * (z[i] + c * (T - t) - b)^2 * exp(-0.5 * (z[i] + c * (T - t) - b)^2 / s^2)
                )
        end
        return disp_z, velo_z, acce_z
    end

    # Loose-ish tolerances throughout: these checks confirm the post-swap FOM
    # is tracking the correct physical solution (catching e.g. a stale,
    # zeroed, or mis-projected acceleration left over from the swap — the
    # symptom reported in issue #197), not tight numerical accuracy. The
    # tolerance is the same order of margin as the analogous FOM -> ROM swap
    # test (schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap.jl, ~10%),
    # loosened slightly to account for this being a ROM -> FOM transfer
    # (lifting a possibly-truncated reduced state into the full nodal space)
    # rather than the reverse.
    for (model, z) in ((model1, model1.reference[3, :]), (model3, model3.reference[3, :]))
        disp_z_exact, velo_z_exact, acce_z_exact = exact_solution(z)
        disp_z = model.displacement[3, :]
        velo_z = model.velocity[3, :]
        acce_z = model.acceleration[3, :]

        disp_z_relerr = norm(disp_z_exact - disp_z) / norm(disp_z_exact)
        velo_z_relerr = norm(velo_z_exact - velo_z) / norm(velo_z_exact)
        acce_z_relerr = norm(acce_z_exact - acce_z) / norm(acce_z_exact)

        @test disp_z_relerr < 0.05
        @test velo_z_relerr < 0.07
        @test acce_z_relerr < 0.1
        @test norm(model.displacement[1, :]) ≈ 0.0 atol = 1.0e-8
        @test norm(model.displacement[2, :]) ≈ 0.0 atol = 1.0e-8
    end
end
