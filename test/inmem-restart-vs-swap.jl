# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# In-memory restart, check 2 of 3: restart and swap must agree.
#
# The same transition -- sd1 FOM -> 10-mode OpInf ROM at t = t_swap -- is driven
# two ways over identical inputs:
#
#   1. Norma's scheduled machinery: a YAML `swaps:` block with a time criterion,
#      run by plain Norma.run. `apply_swap!` rebuilds the replacement subsim from
#      disk, transfers state, and re-points the slot.
#   2. InMemoryRestart: the variant is built once up front, and `switch!` is
#      called from the march hook at the step whose end time is t_swap.
#
# They share `copy_model_state!`, `align_replacement_time!` and
# `_sync_integrator_from_model!`, and both construct the replacement through
# `build_replacement_subsim` -- so the physics after the transition should be the
# same. What differs is only when the variant is built and whether anything
# touches disk. If the two disagree, the in-memory path has dropped or reordered
# part of what `apply_swap!` does, which is exactly the failure this guards.
#
# t_swap is placed before the pulse reaches sd1 (z = -0.20 at t = 2.0e-4, step 64
# of 128), so the comparison is not dominated by M10 truncation over the arrival.

using LinearAlgebra

@testset "In-Memory Restart Matches Swap" begin
    mesh_dir = "../examples/ahead/overlap/clamped/meshes"
    rom_dir = "../examples/ahead/overlap/clamped/inmem_restart"

    TFINAL, NSTEP = 4.0e-4, 128
    DT = TFINAL / NSTEP
    # Deliberately BETWEEN step boundaries. `should_swap` tests
    # `prev_time >= t_swap`, and prev_time is accumulated by repeated addition
    # while t_swap makes a round trip through the YAML as a decimal string --
    # so a t_swap sitting exactly on a boundary (39*DT) differs from prev_time
    # by an ULP and the two routes fire one step apart, which shows up as a
    # ~1e-4 disagreement that has nothing to do with the mechanism.
    T_SWAP = 38.5 * DT             # both routes transition entering step 40
    MESH = ("clamped-2sd-1.g", "clamped-2sd-2.g")
    BLOCK = ("coarse", "fine")
    ZCLAMP = ("nsz-", "nsz+")
    SIDESET = ("ssz+", "ssz-")
    ROMFILE = "linear-opinf-operator-1-M10.npz"

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

    subdomain_yaml(i, kind=:fom) = """
    type: single
    input mesh file: $(MESH[i])
    output mesh file: sd$(i)$(kind).e
    model:
    $(kind === :rom ? "  type: linear opinf rom\n  model-file: $ROMFILE" : "  type: solid mechanics")
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

    # The two top files are identical apart from the `swaps:` block, so anything
    # that differs between the runs comes from the transition mechanism itself.
    top_yaml(swaps) = """
    type: multi
    domains: ["sd1fom.yaml", "sd2fom.yaml"]
    $(swaps)Exodus output interval: 1.0e+06
    CSV output interval: 0
    initial time: 0.0
    final time: $TFINAL
    time step: $DT
    minimum iterations: 1
    maximum iterations: 16
    relative tolerance: 1.0e-12
    absolute tolerance: 1.0e-08
    """

    swaps_block = """
    swaps:
      - subsim: sd1fom
        replacement: sd1rom.yaml
        criterion:
          type: time
          t_swap: $T_SWAP
    """

    # ── Stage ───────────────────────────────────────────────────────────────
    for m in MESH
        cp("$mesh_dir/$m", m; force=true)
    end
    cp("$rom_dir/$ROMFILE", ROMFILE; force=true)
    write("sd1fom.yaml", subdomain_yaml(1))
    write("sd2fom.yaml", subdomain_yaml(2))
    write("sd1rom.yaml", subdomain_yaml(1, :rom))
    write("top-swap.yaml", top_yaml(swaps_block))
    write("top-restart.yaml", top_yaml(""))

    final_displacement(sim) =
        vcat((vec(Norma.full_order_displacement(s)) for s in sim.subsims)...)

    # ── Run 1: scheduled swap ───────────────────────────────────────────────
    swap_sim = Norma.run("top-swap.yaml")
    swapped = final_displacement(swap_sim)
    swap_fired = length(swap_sim.swaps) == 1 && swap_sim.swaps[1].applied
    swap_model_is_rom = swap_sim.subsims[1].model isa Norma.LinearOpInfRom

    # ── Run 2: in-memory restart ────────────────────────────────────────────
    r = Norma.InMemoryRestart("top-restart.yaml")
    restarted = try
        Norma.register_variant!(r, 1, "sd1rom", "sd1rom.yaml")
        # Evaluate the SAME predicate the scheduled path does. Note
        # `maybe_apply_swaps!` calls `should_swap(plan.criterion, sim.subsims[slot])`
        # -- against the SUBSIM, not the MultiDomainSimulation -- and a subsim's
        # clock lags the outer controller by one step at this point in the loop,
        # because subsims are advanced inside `advance_control`, after the hook.
        # Using the outer `r.sim.controller` here instead fires one step early,
        # which shows up as a ~1e-4 disagreement that is a schedule difference,
        # not a mechanism difference.
        Norma.march!(r; on_step = (r, _) ->
            r.sim.subsims[1].controller.prev_time >= T_SWAP &&
                Norma.switch!(r, 1, "sd1rom"))
        Norma.displacement_vector(r)
    finally
        Norma.close_restart!(r)
    end

    nswitch = length(r.history)
    switch_step = nswitch == 1 ? r.history[1].step : -1
    restart_model_is_rom = r.sim.subsims[1].model isa Norma.LinearOpInfRom
    failed = swap_sim.failed || r.sim.failed

    # ── Clean up ────────────────────────────────────────────────────────────
    # The swap path renames the replacement's output mesh with a phase suffix
    # (sd1rom-phase2.e), so sweep by pattern rather than by a fixed list.
    for f in ("sd1fom.yaml", "sd2fom.yaml", "sd1rom.yaml",
              "top-swap.yaml", "top-restart.yaml",
              "sd1fom.e", "sd2fom.e", "sd1rom.e", ROMFILE, MESH...)
        rm(f; force=true)
    end
    for f in readdir()
        occursin(r"^sd\d+\w*-phase\d+\.e$", f) && rm(f; force=true)
    end

    # ── Assert ──────────────────────────────────────────────────────────────
    agreement = norm(restarted - swapped) / norm(swapped)

    @test failed == false
    # Both routes performed the transition exactly once, on the same slot.
    @test swap_fired == true
    @test nswitch == 1
    @test switch_step > 1 && switch_step < NSTEP
    # Both ended with sd1 reduced.
    @test swap_model_is_rom == true
    @test restart_model_is_rom == true
    # Same transition at the same step, same state transfer, same replacement
    # construction ⇒ the two marches must agree to round-off. A nonzero value
    # here means the in-memory path has dropped or reordered part of what
    # apply_swap! does.
    @test agreement < 1.0e-12
end
