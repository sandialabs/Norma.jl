# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# In-memory restart, check 1 of 3: repeated restarts that do not change the model
# must leave the answer at machine precision.
#
# sd1 is switched between two FOM variants of itself -- same mesh, same material,
# same integrator, only a different output file -- every 8 steps, 16 times over a
# 128-step march. FOM -> FOM `copy_model_state!` is an array copy, so none of
# those restarts may perturb the solution: the result must equal plain Norma.run
# to round-off, not merely closely.
#
# This exercises two things at once. The state transfer, clock alignment,
# integrator resync and rollback-buffer seeding run 16 times and must not drift.
# And `step!` in src/restart.jl is a hand-copy of the body of `evolve`'s loop
#
#     advance_control_time -> (hook) -> sync_control_time -> advance_control
#
# so if `evolve` ever gains a call, the copy silently diverges and every result
# produced through the restart API drifts from Norma's own. Nothing else in the
# suite covers that.
#
# 2-subdomain overlapping-Schwarz clamped bar, Gaussian pulse at z = 0 travelling
# left at c = sqrt(E/rho) = 1000 m/s. Meshes come from the shared meshes/ pool.

using LinearAlgebra

@testset "In-Memory Restart Fidelity" begin
    mesh_dir = "../examples/ahead/overlap/clamped/meshes"

    TFINAL, NSTEP = 4.0e-4, 128
    DT = TFINAL / NSTEP
    RESTART_EVERY = 8                    # 16 restarts over the march
    MESH = ("clamped-2sd-1.g", "clamped-2sd-2.g")
    BLOCK = ("coarse", "fine")
    ZCLAMP = ("nsz-", "nsz+")
    SIDESET = ("ssz+", "ssz-")

    ic = """
    initial conditions:
      displacement:
        - node set: nsall
          component: z
          function: "a=0.001; s=0.02; a*exp(-z*z/s/s/2)"
      velocity:
        - node set: nsall
          component: z
          function: "a=0.001; s=0.02; E=1.0e+09; rho=1000.0; c=sqrt(E/rho); -a*c*z/s/s*exp(-z*z/s/s/2)"
    """

    # `tag` distinguishes the two FOM variants of sd1. They differ only in the
    # output file name -- Norma deletes and rewrites its output .e by name, so
    # two live variants cannot share one.
    subdomain_yaml(i, tag) = """
    type: single
    input mesh file: $(MESH[i])
    output mesh file: sd$(i)$(tag).e
    model:
      type: solid mechanics
      material:
        blocks:
          $(BLOCK[i]): hyperelastic
        hyperelastic:
          model: linear elastic
          elastic modulus: 1.0e+09
          Poisson's ratio: 0.0
          density: 1000.0
    time integrator:
      type: Newmark
      β: 0.25
      γ: 0.5
      time step: $DT
    $(ic)boundary conditions:
      Dirichlet:
        - {node set: nsx-, component: x, function: "0.0"}
        - {node set: nsx+, component: x, function: "0.0"}
        - {node set: nsy-, component: y, function: "0.0"}
        - {node set: nsy+, component: y, function: "0.0"}
        - {node set: $(ZCLAMP[i]), component: z, function: "0.0"}
      Schwarz overlap:
        - side set: $(SIDESET[i])
          source: sd$(3 - i)fom
          source block: $(BLOCK[3 - i])
    solver:
      type: Hessian minimizer
      step: full Newton
      minimum iterations: 1
      maximum iterations: 16
      relative tolerance: 1.0e-10
      absolute tolerance: 1.0e-06
    """

    top_yaml = """
    type: multi
    domains: ["sd1fom.yaml", "sd2fom.yaml"]
    Exodus output interval: 1.0e+06
    CSV output interval: 0
    initial time: 0.0
    final time: $TFINAL
    time step: $DT
    minimum iterations: 1
    maximum iterations: 16
    relative tolerance: 1.0e-12
    absolute tolerance: 1.0e-08
    """

    # ── Stage ───────────────────────────────────────────────────────────────
    for m in MESH
        cp("$mesh_dir/$m", m; force=true)
    end
    write("sd1fom.yaml", subdomain_yaml(1, "fom"))
    write("sd2fom.yaml", subdomain_yaml(2, "fom"))
    write("sd1alt.yaml", subdomain_yaml(1, "alt"))     # second FOM variant of sd1
    write("top.yaml", top_yaml)

    final_displacement(sim) =
        vcat((vec(Norma.full_order_displacement(s)) for s in sim.subsims)...)

    # ── Run ─────────────────────────────────────────────────────────────────
    reference = final_displacement(Norma.run("top.yaml"))

    r = Norma.InMemoryRestart("top.yaml")
    restarted = try
        Norma.register_variant!(r, 1, "sd1alt", "sd1alt.yaml")
        Norma.march!(r; on_step = function (r, step)
            # Alternate between the two FOM variants. Both discretize the same
            # problem, so every one of these is a no-op physically -- and must
            # be a no-op numerically.
            if step % RESTART_EVERY == 0
                to = Norma.active_variant(r, 1) == "sd1fom" ? "sd1alt" : "sd1fom"
                Norma.switch!(r, 1, to)
            end
        end)
        Norma.displacement_vector(r)
    finally
        Norma.close_restart!(r)
    end

    nrestart = length(r.history)
    reached = r.sim.controller.time
    failed = r.sim.failed

    # ── Clean up ────────────────────────────────────────────────────────────
    for f in ("sd1fom.yaml", "sd2fom.yaml", "sd1alt.yaml", "top.yaml",
              "sd1fom.e", "sd2fom.e", "sd1alt.e", MESH...)
        rm(f; force=true)
    end

    # ── Assert ──────────────────────────────────────────────────────────────
    fidelity = norm(restarted - reference) / norm(reference)

    @test failed == false
    @test reached ≈ TFINAL rtol = 1.0e-9
    @test nrestart == NSTEP ÷ RESTART_EVERY       # steps 8, 16, ..., 128
    @test length(restarted) == length(reference)
    # Identical solves ⇒ identical answer, however many times the state was
    # handed over. Any nonzero value here means a restart is perturbing the
    # solution, or that the march loop no longer matches `evolve`.
    @test fidelity < 1.0e-13
end
