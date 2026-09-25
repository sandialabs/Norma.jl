# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus

# Runs the AHeaD overlap dynamic ROM-FOM clamped-bar problem from
# examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-rom(-restart)
# (a single symmetric Gaussian displacement/velocity pulse, two overlapping
# Schwarz subdomains: "coarse" — a linear OpInf ROM (M20 basis) — and
# "fine" — a FOM) two ways and compares:
#
#   1. "dynamic"  — the plain (non-restart)
#                    dynamic-linear-elastic-single-gaussian-opinf-rom example,
#                    a single continuous coupled run from t = 0.0 to
#                    t = 1.0e-3.
#   2. "restart"  — the dynamic-linear-elastic-single-gaussian-opinf-rom-restart
#                    example, which resumes both subdomains from the
#                    checkpoint each subdomain wrote at t = 2.9e-4 (index 30
#                    in each subdomain's own Exodus output from run 1, since
#                    Exodus output interval = 1.0e-5 = 10 * time step), using
#                    a single shared `restart: index: 30` in the top-level
#                    multi-domain file, and continues the coupled iteration
#                    on to t = 1.0e-3.
#
# This is the mixed ROM-FOM counterpart of
# schwarz-ahead-overlap-dynamic-clamped-single-gaussian-fom-fom-restart.jl
# (see that test for the full history of the bug this exercises): it
# verifies that the restart-only Schwarz initial-acceleration refinement
# (see initialize(sim::MultiDomainSimulation) in src/simulation.jl and
# initialize(::Newmark, ...; trust_schwarz) / initialize(::RomNewmark, ...;
# trust_schwarz) in src/time_integrator.jl / src/opinf/opinf_time_integrator.jl)
# also produces a correct, non-placeholder-zero initial acceleration at the
# Schwarz-coupled boundary when one side of that boundary is a ROM. A ROM
# subdomain's restart is layered on top of the FOM restart machinery — its
# internal fom_model is restored from the checkpoint exactly like a
# standalone FOM model, and that restored state is then projected onto the
# reduced basis (see apply_ics(::Parameters, ::RomModel, ...) in
# opinf_ics_bcs.jl) — so is_restarted() (model.jl) reports true for the ROM
# subdomain via its fom_model, and the Schwarz refinement pass applies to it
# exactly as it does to a restarted FOM subdomain.
#
# As in the FOM-FOM version, the initial-acceleration-at-restart check below
# covers both subdomains directly (a targeted comparison at the Exodus
# level). The broader final-time kinematic comparison only covers subdomain
# 2 ("fine", the FOM, where the pulse actually lives at the point the two
# runs are compared) — subdomain 1 ("coarse", the ROM) is the far side of
# the bar the pulse has barely (or not yet) reached, so its true
# displacement/velocity/acceleration are themselves ~0 at any given instant;
# a relative-tolerance comparison there is dominated by floating-point-level
# noise accumulated independently over ~700 steps in two separately-run
# trajectories, and is not physically meaningful (same precedent as the
# FOM-FOM test and schwarz-ahead-overlap-dynamic-clamped-single-gaussian-rom-fom-multi-swap.jl).
@testset "Schwarz AHeaD Overlap Dynamic Clamped Single-Gaussian OpInf ROM Restart" begin
    # ── Phase 1: full, uninterrupted coupled run (t = 0.0 -> 1.0e-3) ───────
    example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-rom"
    restart_example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-rom-restart"

    cp("../examples/ahead/overlap/clamped/clamped-smaller-1.g", "../clamped-smaller-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-larger-2.g", "../clamped-larger-2.g"; force=true)
    cp("$example_dir/linear-opinf-operator-1-M20.npz", "linear-opinf-operator-1-M20.npz"; force=true)
    cp("$example_dir/clamped-1.yaml", "clamped-1.yaml"; force=true)
    cp("$example_dir/clamped-2.yaml", "clamped-2.yaml"; force=true)
    cp("$example_dir/clamped.yaml", "clamped.yaml"; force=true)
    # CSV output isn't used by this test (only the Exodus checkpoint and the
    # final in-memory model state are), so turn it off to save I/O.
    write("clamped.yaml", replace(read("clamped.yaml", String), "CSV output interval: 1.0e-5" => "CSV output interval: 0"))

    sim_dynamic = Norma.run("clamped.yaml")

    rm("clamped-1.yaml"; force=true)
    rm("clamped-2.yaml"; force=true)
    rm("clamped.yaml"; force=true)
    rm("../clamped-smaller-1.g"; force=true)
    rm("../clamped-larger-2.g"; force=true)

    # Each subdomain's own Exodus output becomes that subdomain's restart
    # checkpoint. With time step = 1.0e-6 and Exodus output interval =
    # 1.0e-5 (10 * time step), write index k (1-based) holds
    # t = (k - 1) * 1.0e-5, so index 30 is t = 2.9e-4. The ROM subdomain's
    # own Exodus output holds its internal fom_model's reconstructed
    # displacement/velocity/acceleration fields (refreshed from
    # reduced_state/reduced_velocity at every Exodus output step), so it is
    # a valid FOM-style restart checkpoint just like a standalone FOM
    # subdomain's output.
    mv("clamped-1.e", "clamped-1-in.e"; force=true)
    mv("clamped-2.e", "clamped-2-in.e"; force=true)

    # ── Phase 2: restart run resuming both subdomains from t = 2.9e-4 -> 1.0e-3 ──
    # These files already point at "clamped-1-in.e"/"clamped-2-in.e" (the
    # checkpoint names just produced above) and already omit `initial
    # conditions:` in favor of `restart: index: 30` -- no edits needed
    # beyond turning off CSV output.
    cp("$restart_example_dir/clamped-1.yaml", "clamped-1.yaml"; force=true)
    cp("$restart_example_dir/clamped-2.yaml", "clamped-2.yaml"; force=true)
    cp("$restart_example_dir/clamped.yaml", "clamped.yaml"; force=true)
    write("clamped.yaml", replace(read("clamped.yaml", String), "CSV output interval: 1.0e-5" => "CSV output interval: 0"))

    sim_restart = Norma.run("clamped.yaml")

    # ── Check the initial-acceleration-at-restart fix directly ─────────────
    # "clamped-1-in.e"/"clamped-2-in.e" (renamed copies of run 1's own
    # output) still hold the *true*, continuously-integrated acceleration at
    # t = 2.9e-4 (index 30) — ground truth. "clamped-1-out.e"/
    # "clamped-2-out.e" (this restart run's own fresh output, per "output
    # mesh file:" in the example) hold, at Exodus time_index 1, the
    # acceleration this run *solved for* right after restart, before any
    # time-stepping. These should agree, including on the Schwarz-coupled
    # (coarse/fine) overlap boundary specifically, which is exactly where
    # this problem previously showed a large spike -- for subdomain 1, that
    # means the FOM-space fields reconstructed from the ROM's reduced
    # state. A generous tolerance is used since a fresh one-shot solve need
    # not bit-match a value produced by ~30 real Newmark/Schwarz steps --
    # this is a smoke test against a full regression back to the pre-fix,
    # order-unity relative error at the Schwarz-coupled boundary, not a
    # tight accuracy bound.
    for domain in 1:2
        true_mesh = ExodusDatabase("clamped-$domain-in.e", "r")
        restart_mesh = ExodusDatabase("clamped-$domain-out.e", "r")
        local true_acce, fresh_acce
        try
            true_acce_x = Vector{Float64}(Exodus.read_values(true_mesh, NodalVariable, 30, "acce_x"))
            true_acce_y = Vector{Float64}(Exodus.read_values(true_mesh, NodalVariable, 30, "acce_y"))
            true_acce_z = Vector{Float64}(Exodus.read_values(true_mesh, NodalVariable, 30, "acce_z"))
            true_acce = vcat(true_acce_x, true_acce_y, true_acce_z)
            fresh_acce_x = Vector{Float64}(Exodus.read_values(restart_mesh, NodalVariable, 1, "acce_x"))
            fresh_acce_y = Vector{Float64}(Exodus.read_values(restart_mesh, NodalVariable, 1, "acce_y"))
            fresh_acce_z = Vector{Float64}(Exodus.read_values(restart_mesh, NodalVariable, 1, "acce_z"))
            fresh_acce = vcat(fresh_acce_x, fresh_acce_y, fresh_acce_z)
        finally
            Exodus.close(true_mesh)
            Exodus.close(restart_mesh)
        end
        @test fresh_acce ≈ true_acce rtol = 5.0e-02 atol = 1.0e-06 * maximum(abs.(true_acce))
    end

    rm("clamped-1.yaml"; force=true)
    rm("clamped-2.yaml"; force=true)
    rm("clamped.yaml"; force=true)
    rm("clamped-1-in.e"; force=true)
    rm("clamped-2-in.e"; force=true)
    rm("clamped-1-out.e"; force=true)
    rm("clamped-2-out.e"; force=true)
    rm("linear-opinf-operator-1-M20.npz"; force=true)

    # ── Both runs reached the same final time without failing ──────────────
    @test sim_dynamic.failed == false
    @test sim_restart.failed == false
    @test sim_dynamic.controller.time ≈ 1.0e-3 rtol = 1.0e-09
    @test sim_restart.controller.time ≈ 1.0e-3 rtol = 1.0e-09

    # ── Final-step kinematic fields agree between the two runs, subdomain 2 (FOM) only ──
    # See the comment at the top of this file for why subdomain 1 (the ROM)
    # is not checked quantitatively here. A looser tolerance than a typical
    # restart check is still used: Newmark here is non-dissipative (β=0.25,
    # γ=0.5), so any small mismatch between the restart run's fresh
    # one-shot initial acceleration and the continuous run's
    # continuously-integrated value at the same instant is not damped out
    # over the remaining ~700 steps, just carried along. The
    # initial-acceleration check above is the precise regression guard for
    # the actual bug that was fixed; this is a broader sanity check that the
    # two trajectories still end up close.
    model_dynamic = sim_dynamic.subsims[2].model
    model_restart = sim_restart.subsims[2].model

    disp_scale = maximum(abs.(model_dynamic.displacement))
    velo_scale = maximum(abs.(model_dynamic.velocity))
    acce_scale = maximum(abs.(model_dynamic.acceleration))

    @test model_restart.displacement ≈ model_dynamic.displacement rtol = 1.0e-02 atol = 1.0e-06 * disp_scale
    @test model_restart.velocity ≈ model_dynamic.velocity rtol = 1.0e-02 atol = 1.0e-06 * velo_scale
    @test model_restart.acceleration ≈ model_dynamic.acceleration rtol = 1.0e-02 atol = 1.0e-06 * acce_scale
end
