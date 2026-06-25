# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Regression test for a silent-skip bug in multi-link swap chains: a LATER
# plan in a `swaps:` chain, written to target a subsim by its ORIGINAL
# (un-renamed) name, must still fire even if an EARLIER plan's apply_swap!
# call triggered uniquify_swap_output! to rename that subsim's occupant.
#
# Root cause: build_replacement_subsim sets the replacement's intended name
# to stripped_name(plan.replacement_file) (e.g. "clamped-fom-1"), but if that
# replacement's `output mesh file` collides with one already written by an
# earlier phase of the same slot, uniquify_swap_output! renames the occupant
# (e.g. to "clamped-fom-1-phase2") before apply_swap! registers it in
# sim.handle_by_name. A later plan in the SAME `swaps:` chain that names the
# subsim by its original, un-renamed name (exactly what a person writing a
# round-trip chain like ROM -> FOM -> ROM would naturally write, and exactly
# what validate_swap_plans checks the name against) then has no matching key
# in handle_by_name. maybe_apply_swaps! silently skips any plan whose
# `subsim:` name isn't yet a registered key
# (`haskey(sim.handle_by_name, plan.subsim_name) || continue`), with no
# warning — so the plan is dropped forever, on every subsequent step, with
# no error and no indication anything went wrong.
#
# The fix (apply_swap! in src/swap.jl) registers the pre-rename intended name
# as a surviving alias for the slot alongside the actual (possibly renamed)
# occupant name, so later plans resolve correctly regardless of any rename.
#
# This test exercises a full round trip on all three subdomains of the
# clamped-bar overlap example (see the sibling
# dynamic-linear-elastic-opinf-3sd-fom-rom-single-swap example directory):
#   subdomain 1: ROM -> FOM (t=3e-4) -> ROM (t=6e-4)
#   subdomain 2: FOM -> ROM (t=2.5e-4) -> FOM (t=5.5e-4)
#   subdomain 3: ROM -> FOM (t=3e-4) -> ROM (t=6e-4)
# Subdomains 1 and 3 hit the rename (clamped-rom-*.yaml and clamped-fom-*.yaml
# both write to the same output mesh file, clamped-3sd-1.e / clamped-3sd-3.e),
# so their SECOND swap is exactly the case that silently never fired before
# the fix. Subdomain 2's two replacement files write to two DIFFERENT output
# files (clamped-3sd-rom-2.e vs clamped-3sd-fom-2.e), so its first swap never
# renames anything — included here as a same-chain control showing the
# no-rename case was never broken.
#
# Run to just past the second round of swaps (final time shortened to
# 7.0e-4, comfortably past the last t_swap of 6.0e-4). The primary check is
# that every plan in the chain actually FIRED (sim.swaps[i].applied); on
# top of that, each subsim's final kinematic state (displacement, velocity,
# acceleration) is also compared against the clamped bar's exact analytical
# solution, as a sanity check that the round trip didn't leave any subsim
# with stale, zeroed, or mis-projected state.

using LinearAlgebra
using YAML

@testset "AHeaD Overlap Dynamic Clamped 3SD Round-Trip Swap Chain" begin
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-single-swap"

    # ── Stage files ─────────────────────────────────────────────────────────
    cp("../examples/ahead/overlap/clamped/clamped-3sd-1.g", "../clamped-3sd-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-2.g", "../clamped-3sd-2.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-3sd-3.g", "../clamped-3sd-3.g"; force=true)

    for f in ("clamped-fom-1.yaml", "clamped-fom-2.yaml", "clamped-fom-3.yaml",
              "clamped-rom-1.yaml", "clamped-rom-2.yaml", "clamped-rom-3.yaml",
              "linear-opinf-operator-M30-1.npz", "linear-opinf-operator-M30-2.npz",
              "linear-opinf-operator-M30-3.npz")
        cp("$example_dir/$f", f; force=true)
    end

    t_swap_1 = 2.5e-4   # subdomain 2: FOM -> ROM
    t_swap_2 = 3.0e-4   # subdomains 1, 3: ROM -> FOM
    t_swap_3 = 5.5e-4   # subdomain 2: ROM -> FOM (back)
    t_swap_4 = 6.0e-4   # subdomains 1, 3: FOM -> ROM (back)
    final_time = 7.0e-4 # just past the last swap; see header comment
    time_step = 1.0e-6

    swap_yaml = """
    type: multi
    domains: ["clamped-rom-1.yaml", "clamped-fom-2.yaml", "clamped-rom-3.yaml"]
    swaps:
      - subsim: clamped-rom-1
        replacement: clamped-fom-1.yaml
        criterion:
          type: time
          t_swap: $t_swap_2
      - subsim: clamped-fom-2
        replacement: clamped-rom-2.yaml
        criterion:
          type: time
          t_swap: $t_swap_1
      - subsim: clamped-rom-3
        replacement: clamped-fom-3.yaml
        criterion:
          type: time
          t_swap: $t_swap_2
      - subsim: clamped-fom-1
        replacement: clamped-rom-1.yaml
        criterion:
          type: time
          t_swap: $t_swap_4
      - subsim: clamped-rom-2
        replacement: clamped-fom-2.yaml
        criterion:
          type: time
          t_swap: $t_swap_3
      - subsim: clamped-fom-3
        replacement: clamped-rom-3.yaml
        criterion:
          type: time
          t_swap: $t_swap_4
    Exodus output interval: 1.0e-5
    CSV output interval: 0.0
    initial time: 0.0
    final time: $final_time
    time step: $time_step
    minimum iterations: 1
    maximum iterations: 16
    relative tolerance: 1.0e-12
    absolute tolerance: 1.0e-08
    """
    write("clamped-swap-round-trip.yaml", swap_yaml)

    params = YAML.load_file("clamped-swap-round-trip.yaml"; dicttype=Norma.Parameters)
    params["name"] = "clamped-swap-round-trip"

    sim = Norma.run(params)

    # ── Clean up ────────────────────────────────────────────────────────────
    for f in ("clamped-fom-1.yaml", "clamped-fom-2.yaml", "clamped-fom-3.yaml",
              "clamped-rom-1.yaml", "clamped-rom-2.yaml", "clamped-rom-3.yaml",
              "clamped-swap-round-trip.yaml",
              "linear-opinf-operator-M30-1.npz", "linear-opinf-operator-M30-2.npz",
              "linear-opinf-operator-M30-3.npz",
              "../clamped-3sd-1.g", "../clamped-3sd-2.g", "../clamped-3sd-3.g",
              "clamped-3sd-1.e", "clamped-3sd-fom-2.e", "clamped-3sd-3.e",
              "clamped-3sd-1-phase2.e", "clamped-3sd-rom-2.e", "clamped-3sd-3-phase2.e",
              "clamped-3sd-1-phase3.e", "clamped-3sd-fom-2-phase2.e", "clamped-3sd-3-phase3.e")
        rm(f; force=true)
    end
    for f in filter(f -> startswith(f, "clamped-") && endswith(f, ".csv"), readdir())
        rm(f; force=true)
    end

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ final_time rtol = 1.0e-9

    # ── Every plan in the chain fired — this is the actual regression check.
    # Before the fix, the 4th and 6th plans here (subdomains 1 and 3's second
    # swap, back to ROM) were silently skipped forever because their
    # `subsim:` name ("clamped-fom-1" / "clamped-fom-3") never became a
    # resolvable key in handle_by_name once uniquify_swap_output! renamed
    # the occupant to "clamped-fom-1-phase2" / "clamped-fom-3-phase2" at the
    # first swap. ──────────────────────────────────────────────────────────
    @test length(sim.swaps) == 6
    @test all(p.applied for p in sim.swaps)

    # ── Final occupant identity for each slot ────────────────────────────────
    # Slots 1 and 3: ROM -> FOM -> ROM, hitting the output-file-collision
    # rename on BOTH swaps (clamped-3sd-1.e / clamped-3sd-3.e is reused by
    # every phase), so the final occupant is "-phase3" (uniquify_swap_output!
    # increments past the already-used "-phase2" from the first swap).
    @test sim.subsims[1].name == "clamped-rom-1-phase3"
    @test sim.subsims[3].name == "clamped-rom-3-phase3"
    @test sim.subsims[1].model isa Norma.LinearOpInfRom
    @test sim.subsims[3].model isa Norma.LinearOpInfRom

    # Slot 2: FOM -> ROM -> FOM. The first swap writes to a DIFFERENT output
    # file (clamped-3sd-rom-2.e, distinct from clamped-3sd-fom-2.e), so no
    # rename happens there ("clamped-rom-2", unrenamed). The second swap
    # reuses clamped-3sd-fom-2.e (written by the very first, pre-swap phase),
    # so THAT one rename happens, landing on "-phase2".
    @test sim.subsims[2].name == "clamped-fom-2-phase2"
    @test sim.subsims[2].model isa Norma.SolidMechanics

    # ── All the original AND intended-but-renamed names remain resolvable ──
    # (both kept as surviving aliases for their slot, per the apply_swap!
    # comments in src/swap.jl).
    @test sim.handle_by_name["clamped-rom-1"].id == 1
    @test sim.handle_by_name["clamped-fom-1"].id == 1
    @test sim.handle_by_name["clamped-fom-1-phase2"].id == 1
    @test sim.handle_by_name["clamped-rom-1-phase3"].id == 1

    @test sim.handle_by_name["clamped-fom-2"].id == 2
    @test sim.handle_by_name["clamped-rom-2"].id == 2
    @test sim.handle_by_name["clamped-fom-2-phase2"].id == 2

    @test sim.handle_by_name["clamped-rom-3"].id == 3
    @test sim.handle_by_name["clamped-fom-3"].id == 3
    @test sim.handle_by_name["clamped-fom-3-phase2"].id == 3
    @test sim.handle_by_name["clamped-rom-3-phase3"].id == 3

    # ── Physical sanity check against the exact analytical solution ─────────
    # Exact 1-D wave solution for a Gaussian pulse on a bar clamped at both
    # ends (superposition of the two traveling-wave reflections), evaluated
    # at this test's (shortened) final time t = final_time = 7.0e-4. T is
    # the *original* problem's final time (1.0e-3); it is a fixed parameter
    # of the closed-form solution (the reflection timescale) and does not
    # change just because this test stops the simulation early. See
    # schwarz-ahead-overlap-dynamic-clamped.jl for the same closed-form
    # expression (disp/velo/acce) applied to a 2-subdomain, no-swap problem
    # with the same material properties and initial condition.
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

    # Loose tolerance throughout: this checks that every subsim is still
    # tracking the correct physical solution after the full round trip
    # (catching e.g. a stale, zeroed, or mis-projected state left over from
    # any of the six swaps), not tight numerical accuracy. Slots 1 and 3 end
    # the run as ROMs, so their physical-space fields live on the shadow FOM
    # (model.fom_model); slot 2 ends as a FOM, so its fields are read
    # directly.
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

        @test disp_z_relerr < 0.045
        @test velo_z_relerr < 0.07
        @test acce_z_relerr < 0.09
    end
end
