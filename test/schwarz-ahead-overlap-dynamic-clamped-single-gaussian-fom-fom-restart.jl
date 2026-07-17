# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus

# Runs the AHeaD overlap dynamic FOM-FOM clamped-bar problem (a single
# symmetric Gaussian displacement/velocity pulse, two overlapping Schwarz
# subdomains: "coarse" and "fine") two ways and compares:
#
#   1. "dynamic"  — a single continuous coupled run from t = 0.0 to t = 1.0e-3.
#   2. "restart"  — a restart run that resumes both subdomains from the
#                    checkpoint each subdomain wrote at t = 2.9e-4 (index 30
#                    in each subdomain's own Exodus output from run 1, since
#                    Exodus output interval = 1.0e-5 = 10 * time step), using
#                    a single shared `restart: index: 30` in the top-level
#                    multi-domain file, and continues the coupled iteration
#                    on to t = 1.0e-3.
#
# This is the FOM-FOM overlap configuration that originally exposed a bug in
# the Schwarz-restart initial-acceleration solve: apply_bc() for a
# Schwarz-coupled boundary interpolates from each subdomain's *history*
# (controller.*_hist), not its live integrator state, and that history's
# first entry is pushed before any subdomain has a real acceleration. Right
# after restart, the true acceleration at the coarse/fine coupling boundary
# is not small (the pulse is passing through it), so a stale placeholder
# there is very visibly wrong — first as a large spike at the shared
# boundary right at the restart instant, then propagating through the
# domain over subsequent steps (Newmark with β=0.25, γ=0.5 has no numerical
# damping to dissipate it). See initialize(sim::MultiDomainSimulation) and
# initialize(::Newmark, ...; trust_schwarz) in src/simulation.jl and
# src/time_integrator.jl for the fix.
#
# The initial-acceleration-at-restart check below covers both subdomains
# directly (it's a targeted comparison at the Exodus level, immune to the
# issue described next). The broader final-time kinematic comparison,
# though, covers subdomain 1 ("coarse") only implicitly via sim_dynamic
# failing to fail -- it is not checked quantitatively against subdomain 2,
# matching the precedent in
# schwarz-ahead-overlap-dynamic-clamped-single-gaussian-rom-fom-multi-swap.jl:
# subdomain 1 is the far side of the bar the pulse has barely (or not yet)
# reached, so its true displacement/velocity/acceleration are themselves
# ~0 at any given instant, not just small -- a relative-tolerance
# comparison there is dominated by floating-point-level noise accumulated
# independently over ~700 steps in two separately-run trajectories, and is
# not physically meaningful. Only subdomain 2 ("fine"), where the pulse
# actually lives, is checked quantitatively at final time.
@testset "Schwarz AHeaD Overlap Dynamic Clamped Single-Gaussian FOM-FOM Restart" begin
    # ── Phase 1: full, uninterrupted coupled run (t = 0.0 -> 1.0e-3) ───────
    cp("../examples/ahead/overlap/clamped/clamped-smaller-1.g", "../clamped-smaller-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-larger-2.g", "../clamped-larger-2.g"; force=true)
    cp(
        "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-fom/clamped.yaml",
        "clamped.yaml";
        force=true,
    )
    # CSV output isn't used by this test (only the Exodus checkpoint and the
    # final in-memory model state are), so turn it off to save I/O.
    write("clamped.yaml", replace(read("clamped.yaml", String), "CSV output interval: 1.0e-5" => "CSV output interval: 0"))
    cp(
        "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-fom/clamped-1.yaml",
        "clamped-1.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-fom/clamped-2.yaml",
        "clamped-2.yaml";
        force=true,
    )

    sim_dynamic = Norma.run("clamped.yaml")

    rm("clamped-1.yaml"; force=true)
    rm("clamped-2.yaml"; force=true)
    rm("clamped.yaml"; force=true)
    rm("../clamped-smaller-1.g"; force=true)
    rm("../clamped-larger-2.g"; force=true)

    # Each subdomain's own Exodus output becomes that subdomain's restart
    # checkpoint. With time step = 1.0e-6 and Exodus output interval =
    # 1.0e-5 (10 * time step), write index k (1-based) holds
    # t = (k - 1) * 1.0e-5, so index 30 is t = 2.9e-4.
    mv("clamped-1.e", "clamped-1-in.e"; force=true)
    mv("clamped-2.e", "clamped-2-in.e"; force=true)

    # ── Phase 2: restart run resuming both subdomains from t = 2.9e-4 -> 1.0e-3 ──
    # These files already point at "clamped-1-in.e"/"clamped-2-in.e" (the
    # checkpoint names just produced above), already omit `initial
    # conditions:` in favor of a `restart:` block (restart supplies the
    # initial displacement/velocity instead — the two can't coexist, see
    # process_restart!()), and already write to "clamped-1-out.e"/
    # "clamped-2-out.e" — no edits needed to the subdomain files. The
    # top-level file's own `restart: index: -11` is patched to `index: 30`
    # (t = 2.9e-4): the example's own choice of restart point isn't tied to
    # anything in particular, but this test deliberately restarts while the
    # pulse is actively crossing the Schwarz-coupled boundary — see the
    # comment at the top of this file — so that specific index is preserved
    # here rather than adopted from the example as-is.
    restart_example_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-single-gaussian-opinf-fom-restart"
    cp("$restart_example_dir/clamped-1.yaml", "clamped-1.yaml"; force=true)
    cp("$restart_example_dir/clamped-2.yaml", "clamped-2.yaml"; force=true)
    cp("$restart_example_dir/clamped.yaml", "clamped.yaml"; force=true)
    write(
        "clamped.yaml",
        replace(
            read("clamped.yaml", String),
            "index: -11" => "index: 30",
            "CSV output interval: 1.0e-5" => "CSV output interval: 0",
        ),
    )

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
    # this problem previously showed a large spike. A generous tolerance is
    # used since a fresh one-shot solve need not bit-match a value produced
    # by ~30 real Newmark/Schwarz steps -- this is a smoke test against a
    # full regression back to the pre-fix, order-unity relative error at the
    # Schwarz-coupled boundary, not a tight accuracy bound.
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

    # ── Both runs reached the same final time without failing ──────────────
    @test sim_dynamic.failed == false
    @test sim_restart.failed == false
    @test sim_dynamic.controller.time ≈ 1.0e-3 rtol = 1.0e-09
    @test sim_restart.controller.time ≈ 1.0e-3 rtol = 1.0e-09

    # ── Final-step kinematic fields agree between the two runs, subdomain 2 only ──
    # Subdomain 1 is deliberately not checked quantitatively here -- see the
    # comment at the top of this file. A looser tolerance than a typical
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
