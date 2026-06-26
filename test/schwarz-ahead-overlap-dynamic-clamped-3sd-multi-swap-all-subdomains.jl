# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests a 3-subdomain overlapping-Schwarz dynamic simulation (clamped bar,
# Gaussian-pulse initial condition) in which ALL THREE subdomains round-trip
# their model type mid-run: the two outer subdomains go ROM -> FOM -> ROM,
# while the middle subdomain goes FOM -> ROM -> FOM, via six swap plans
# declared up front in clamped-swap.yaml.
#
# Problem setup (see clamped-swap.yaml):
#   domains: clamped-rom-1 (ROM), clamped-fom-2 (FOM), clamped-rom-3 (ROM)
#   Material: linear elastic, E = 1 GPa, ν = 0, ρ = 1000 kg/m³ (all subdomains)
#   Time integrator: Newmark (β = 0.25, γ = 0.5) for every subdomain
#   Initial condition (nsall, z-component): a Gaussian pulse
#       u_z(z, 0) = a * exp(-z² / (2 s²)),  a = 1.0e-3,  s = 0.02
#   All lateral faces (x±, y±) pinned in their normal direction.
#
# Swaps (all six declared in clamped-swap.yaml):
#   1. clamped-rom-1 -> clamped-fom-1 (slot 1, ROM -> FOM) at t_swap = 3.0e-4
#   2. clamped-fom-2 -> clamped-rom-2 (slot 2, FOM -> ROM) at t_swap = 2.5e-4
#   3. clamped-rom-3 -> clamped-fom-3 (slot 3, ROM -> FOM) at t_swap = 3.0e-4
#   4. clamped-fom-1 -> clamped-rom-1 (slot 1, FOM -> ROM) at t_swap = 6.0e-4
#   5. clamped-rom-2 -> clamped-fom-2 (slot 2, ROM -> FOM) at t_swap = 5.5e-4
#   6. clamped-fom-3 -> clamped-rom-3 (slot 3, FOM -> ROM) at t_swap = 6.0e-4
# This is the same six-swap, all-subdomains round-trip chain exercised inline
# in schwarz-ahead-overlap-dynamic-clamped-3sd-round-trip-swap-chain.jl (see
# that test for the issue #196 silent-skip regression this chain originally
# targeted); here it is run from the materialized example directory itself
# rather than a YAML string built on the fly.
#
# Unlike clamped-swap.yaml itself, this test only runs to 0.75 of the
# original final time (7.5e-4 instead of 1.0e-3) and disables CSV output
# (the YAML file on disk is left untouched; both overrides are applied
# in-memory to the loaded params before calling Norma.run), matching the
# choice made in schwarz-ahead-overlap-dynamic-clamped-3sd-rom-fom-swap-sides.jl
# — see that test's header comment for why 0.75 * T sits on the signal
# plateau rather than in the near-cancellation valley around [4.5e-4, 6.3e-4].
# All six t_swap thresholds above are at or below 6.0e-4, so every swap has
# already fired by 7.5e-4, leaving slots 1 and 3 back on ROMs and slot 2 back
# on a FOM — the same final occupancy as the round-trip-chain test.
#
# The clamped bar's exact 1-D wave solution (reflections off both clamped
# ends superposed) is known in closed form; see schwarz-ahead-overlap-dynamic-
# clamped.jl for the same formula (disp/velo/acce) applied to a 2-subdomain,
# no-swap case. Here it is evaluated on all three subdomains at the
# (shortened) final time, after every swap has fired, to confirm each
# round-tripped subdomain holds physically correct displacement, velocity,
# AND acceleration fields rather than e.g. stale, zeroed, or mis-projected
# data — acceleration in particular is the field reported as corrupted in
# issue #197. Slots 1 and 3 end the run as ROMs, so their physical-space
# fields are read from the shadow FOM (model.fom_model); slot 2 ends as a
# FOM, so its fields are read directly.

using LinearAlgebra
using YAML

@testset "AHeaD Overlap Dynamic Clamped 3SD Multi-Swap All Subdomains" begin
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-multi-swap-all-subdomains"

    # ── Stage files ─────────────────────────────────────────────────────────
    # The subsim YAMLs reference their meshes as ../clamped-3sd-<k>.g
    # (relative to the test working directory), so the meshes are copied one
    # directory above, matching the convention used by the other ahead/
    # overlap/clamped tests.
    cp("../examples/ahead/overlap/clamped/clamped-3sd-1.g", "../clamped-3sd-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-2.g", "../clamped-3sd-2.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-3.g", "../clamped-3sd-3.g"; force=true)

    for f in (
        "clamped-fom-1.yaml",
        "clamped-fom-2.yaml",
        "clamped-fom-3.yaml",
        "clamped-rom-1.yaml",
        "clamped-rom-2.yaml",
        "clamped-rom-3.yaml",
        "linear-opinf-operator-M30-1.npz",
        "linear-opinf-operator-M30-2.npz",
        "linear-opinf-operator-M30-3.npz",
    )
        cp("$example_dir/$f", f; force=true)
    end

    # Load the parent YAML into a params dict (instead of running the file
    # directly) so the final time and CSV output can be overridden in-memory
    # without touching the original example file.
    input_file = "$example_dir/clamped-swap.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    final_time = 0.75 * params["final time"]  # 7.5e-4; see header comment for why
    params["final time"] = final_time
    # CSV output interval = 0 disables write_stop's is_csv_step entirely, which
    # also gates the per-subsim "CSV write sidesets" branch (still set in the
    # individual subsim YAMLs) — so this one override suffices to suppress
    # all CSV output, sidesets included, without editing those files.
    params["CSV output interval"] = 0.0
    params["name"] = "clamped-swap"

    sim = Norma.run(params)

    # ── Clean up ────────────────────────────────────────────────────────────
    for f in (
        "clamped-fom-1.yaml",
        "clamped-fom-2.yaml",
        "clamped-fom-3.yaml",
        "clamped-rom-1.yaml",
        "clamped-rom-2.yaml",
        "clamped-rom-3.yaml",
        "linear-opinf-operator-M30-1.npz",
        "linear-opinf-operator-M30-2.npz",
        "linear-opinf-operator-M30-3.npz",
        "../clamped-3sd-1.g",
        "../clamped-3sd-2.g",
        "../clamped-3sd-3.g",
        "clamped-3sd-1.e",
        "clamped-3sd-fom-2.e",
        "clamped-3sd-3.e",
        "clamped-3sd-1-phase2.e",
        "clamped-3sd-rom-2.e",
        "clamped-3sd-3-phase2.e",
        "clamped-3sd-1-phase3.e",
        "clamped-3sd-fom-2-phase2.e",
        "clamped-3sd-3-phase3.e",
    )
        rm(f; force=true)
    end
    for f in filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ final_time rtol = 1.0e-9

    # ── No CSV output was written ────────────────────────────────────────────
    # Scoped to this test's own "clamped-" prefix: other tests that ran
    # earlier in this shared directory may leave their own CSV files behind
    # (e.g. cantilever-*, cuboid-*), and this test has no business asserting
    # on those.
    @test isempty(filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir()))

    # ── All six swaps fired, in order, across all three slots ───────────────
    @test length(sim.swaps) == 6
    @test all(p.applied for p in sim.swaps)
    @test sim.swaps[1].subsim_name == "clamped-rom-1"
    @test sim.swaps[2].subsim_name == "clamped-fom-2"
    @test sim.swaps[3].subsim_name == "clamped-rom-3"
    @test sim.swaps[4].subsim_name == "clamped-fom-1"
    @test sim.swaps[5].subsim_name == "clamped-rom-2"
    @test sim.swaps[6].subsim_name == "clamped-fom-3"

    # ── Final occupant identity for each slot ────────────────────────────────
    # Slots 1 and 3: ROM -> FOM -> ROM, hitting the output-file-collision
    # rename on BOTH swaps (clamped-fom-k.yaml and clamped-rom-k.yaml both
    # write to clamped-3sd-k.e for k = 1, 3), so the final occupant name is
    # "-phase3" (uniquify_swap_output! increments past the already-used
    # "-phase2" from the first swap). Slot 2: FOM -> ROM -> FOM. The first
    # swap writes to a DIFFERENT output file (clamped-3sd-rom-2.e, distinct
    # from clamped-3sd-fom-2.e), so no rename happens there. The second swap
    # reuses clamped-3sd-fom-2.e (written by the very first, pre-swap phase),
    # so that one rename happens, landing on "-phase2".
    @test sim.subsims[1].name == "clamped-rom-1-phase3"
    @test sim.subsims[2].name == "clamped-fom-2-phase2"
    @test sim.subsims[3].name == "clamped-rom-3-phase3"
    @test sim.subsims[1].model isa Norma.LinearOpInfRom
    @test sim.subsims[2].model isa Norma.SolidMechanics
    @test sim.subsims[3].model isa Norma.LinearOpInfRom

    # The slot ↔ name lookup is re-keyed to the post-swap occupant, and every
    # intermediate name from earlier swaps survives as an alias for the same
    # slot (see the apply_swap! comment in src/swap.jl for why old names are
    # kept rather than removed).
    @test sim.handle_by_name["clamped-rom-1"].id == 1
    @test sim.handle_by_name["clamped-fom-1"].id == 1
    @test sim.handle_by_name["clamped-rom-1-phase3"].id == 1

    @test sim.handle_by_name["clamped-fom-2"].id == 2
    @test sim.handle_by_name["clamped-rom-2"].id == 2
    @test sim.handle_by_name["clamped-fom-2-phase2"].id == 2

    @test sim.handle_by_name["clamped-rom-3"].id == 3
    @test sim.handle_by_name["clamped-fom-3"].id == 3
    @test sim.handle_by_name["clamped-rom-3-phase3"].id == 3

    # ── Physical sanity check against the exact analytical solution ─────────
    # Exact 1-D wave solution for a Gaussian pulse on a bar clamped at both
    # ends (superposition of the two traveling-wave reflections), evaluated
    # at this test's (shortened) final time t = final_time = 7.5e-4. T is the
    # *original* problem's final time (1.0e-3); it is a fixed parameter of
    # the closed-form solution (the reflection timescale) and does not change
    # just because this test stops the simulation early.
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

    # Loose tolerance throughout: these checks confirm every subsim is still
    # tracking the correct physical solution after the full six-swap chain
    # (catching e.g. a stale, zeroed, or mis-projected state left over from
    # any of the swaps — the symptom reported in issue #197), not tight
    # numerical accuracy.
    for model in (sim.subsims[1].model, sim.subsims[2].model, sim.subsims[3].model)
        fom = model isa Norma.RomModel ? model.fom_model : model
        z = fom.reference[3, :]
        disp_z_exact, velo_z_exact, acce_z_exact = exact_solution(z)

        disp_z = fom.displacement[3, :]
        velo_z = fom.velocity[3, :]
        acce_z = fom.acceleration[3, :]

        disp_z_relerr = norm(disp_z_exact - disp_z) / norm(disp_z_exact)
        velo_z_relerr = norm(velo_z_exact - velo_z) / norm(velo_z_exact)
        acce_z_relerr = norm(acce_z_exact - acce_z) / norm(acce_z_exact)

        @test disp_z_relerr < 0.04
        @test velo_z_relerr < 0.07
        @test acce_z_relerr < 0.09
        @test norm(fom.displacement[1, :]) ≈ 0.0 atol = 1.0e-8
        @test norm(fom.displacement[2, :]) ≈ 0.0 atol = 1.0e-8
    end
end
