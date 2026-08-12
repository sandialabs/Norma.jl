# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# In-memory restart, check 3 of 3: it works for any number of subdomains.
#
# The same FOM -> ROM -> FOM round trip is run on 2-, 3-, 4- and 5-subdomain
# clamped bars. Nothing in `switch!` is written for a particular slot count, but the
# things it re-points are shared across the whole simulation -- `subsims`,
# `handle_by_name`, `name_by_handle`, and the Schwarz BC wiring rebuilt by
# `pair_schwarz_bcs` -- so a subdomain count where a slot has TWO overlap
# partners rather than one exercises paths a 2-subdomain case never reaches.
# 3sd, 4sd and 5sd all have interior subdomains; 2sd does not.
#
# The switched subdomain is an interior one wherever there is a choice, since
# that is the case with two neighbours to re-couple.
#
# Assertions are deliberately structural rather than numerical: the march
# completes, the two switches land where asked, the model type flips and comes
# back, and the result stays finite and close to the untouched march. The
# quantitative checks live in inmem-restart-fidelity.jl (round-off under repeated
# restarts) and inmem-restart-vs-swap.jl (agreement with the scheduled path).
#
# Meshes come from the shared meshes/ pool, named clamped-<N>sd-<i>.g.

using LinearAlgebra

@testset "In-Memory Restart Multi-Subdomain" begin
    mesh_dir = "../examples/ahead/overlap/clamped/meshes"
    rom_dir = "../examples/ahead/overlap/clamped/inmem_restart"
    rom3_dir = "../examples/ahead/overlap/clamped/dynamic-linear-elastic-opinf-3sd-fom-rom-single-swap"

    TFINAL, NSTEP = 4.0e-4, 64
    DT = TFINAL / NSTEP
    ROM_WINDOW = (10, 30)

    # A ROM basis only spans the trajectory it was fitted on, so each case must
    # use the initial condition its operator was built for. The 2sd and 5sd
    # operators were fitted with a velocity consistent with a left-travelling
    # wave (v = c*u'); the 3sd operator that ships with the swap example was
    # fitted from rest, where the pulse instead splits into two half-amplitude
    # waves. Mixing them up makes the ROM meaningless -- measured: drift 0.9999,
    # i.e. the reduced solution is uncorrelated with the reference.
    ic_travelling = """
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
    ic_from_rest = """
    initial conditions:
      displacement:
        - node set: nsall
          component: z
          function: "a=0.001; s=0.02; a*exp(-z*z/s/s/2)"
    """

    # (label, meshes, blocks, ROM operator, switching slot, IC)
    cases = [
        ("2sd",
         ["clamped-2sd-1.g", "clamped-2sd-2.g"],
         ["coarse", "fine"],
         ["$rom_dir/linear-opinf-operator-1-M10.npz"],
         1, ic_travelling),
        ("3sd",
         ["clamped-3sd-1.g", "clamped-3sd-2.g", "clamped-3sd-3.g"],
         ["subdomain_1", "subdomain_2", "subdomain_3"],
         ["$rom3_dir/linear-opinf-operator-M30-2.npz"],
         2, ic_from_rest),
        ("4sd",
         ["clamped-4sd-$i.g" for i in 1:4],
         ["subdomain_$i" for i in 1:4],
         ["$rom_dir/linear-opinf-operator-4sd-2-M10.npz"],
         2, ic_travelling),
        ("5sd",
         ["clamped-5sd-$i.g" for i in 1:5],
         ["subdomain_$i" for i in 1:5],
         ["$rom_dir/linear-opinf-operator-5sd-3-M10.npz"],
         3, ic_travelling),
    ]

    """
    Subdomain YAML for a bar split into `n` pieces. The end subdomains carry the
    physical z clamp and one Schwarz overlap; interior ones carry no z clamp and
    two overlaps, which is the wiring a 2-subdomain case never produces.
    """
    function subdomain_yaml(i, n, meshes, blocks, romfile, ic)
        dir = [ "    - {node set: nsx-, component: x, function: \"0.0\"}",
                "    - {node set: nsx+, component: x, function: \"0.0\"}",
                "    - {node set: nsy-, component: y, function: \"0.0\"}",
                "    - {node set: nsy+, component: y, function: \"0.0\"}" ]
        i == 1 && push!(dir, "    - {node set: nsz-, component: z, function: \"0.0\"}")
        i == n && push!(dir, "    - {node set: nsz+, component: z, function: \"0.0\"}")
        ov = String[]
        i > 1 && push!(ov, "    - side set: ssz-\n      source: sd$(i-1)fom\n" *
                           "      source block: $(blocks[i-1])")
        i < n && push!(ov, "    - side set: ssz+\n      source: sd$(i+1)fom\n" *
                           "      source block: $(blocks[i+1])")
        model = romfile === nothing ? "  type: solid mechanics" :
                "  type: linear opinf rom\n  model-file: $romfile"
        tag = romfile === nothing ? "fom" : "rom"
        return """
        type: single
        input mesh file: $(meshes[i])
        output mesh file: sd$(i)$(tag).e
        model:
        $model
          material:
            blocks:
              $(blocks[i]): hyperelastic
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
        $(join(dir, "\n"))
          Schwarz overlap:
        $(join(ov, "\n"))
        solver:
          type: Hessian minimizer
          step: full Newton
          minimum iterations: 1
          maximum iterations: 16
          relative tolerance: 1.0e-10
          absolute tolerance: 1.0e-06
        """
    end

    top_yaml(n) = """
    type: multi
    domains: [$(join(["\"sd$(i)fom.yaml\"" for i in 1:n], ", "))]
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

    for (label, meshes, blocks, romfiles, slot, ic) in cases
        n = length(meshes)
        a, b = ROM_WINDOW
        romname = basename(romfiles[1])

        # ── Stage ───────────────────────────────────────────────────────────
        for m in meshes
            cp("$mesh_dir/$m", m; force=true)
        end
        cp(romfiles[1], romname; force=true)
        for i in 1:n
            write("sd$(i)fom.yaml", subdomain_yaml(i, n, meshes, blocks, nothing, ic))
        end
        write("sd$(slot)rom.yaml", subdomain_yaml(slot, n, meshes, blocks, romname, ic))
        write("top.yaml", top_yaml(n))

        # ── Run ─────────────────────────────────────────────────────────────
        base = Norma.run("top.yaml")
        baseline = vcat((vec(Norma.full_order_displacement(s)) for s in base.subsims)...)

        r = Norma.InMemoryRestart("top.yaml")
        switched = try
            Norma.register_variant!(r, slot, "sd$(slot)rom", "sd$(slot)rom.yaml")
            Norma.march!(r; on_step = (r, k) ->
                Norma.switch!(r, slot, a <= k <= b ? "sd$(slot)rom" : "sd$(slot)fom"))
            Norma.displacement_vector(r)
        finally
            Norma.close_restart!(r)
        end

        nswitch = length(r.history)
        steps = [h.step for h in r.history]
        slots = [h.slot for h in r.history]
        ended_fom = Norma.active_variant(r, slot) == "sd$(slot)fom"
        ended_model_is_fom = r.sim.subsims[slot].model isa Norma.SolidMechanics
        nslot = length(r.sim.subsims)
        reached = r.sim.controller.time
        failed = r.sim.failed
        drift = norm(switched - baseline) / norm(baseline)

        # ── Clean up ────────────────────────────────────────────────────────
        for i in 1:n
            rm("sd$(i)fom.yaml"; force=true)
            rm("sd$(i)fom.e"; force=true)
        end
        for f in ("sd$(slot)rom.yaml", "sd$(slot)rom.e", "top.yaml", romname, meshes...)
            rm(f; force=true)
        end

        # ── Assert ──────────────────────────────────────────────────────────
        @testset "$label" begin
            @test nslot == n
            @test failed == false
            @test reached ≈ TFINAL rtol = 1.0e-9
            # Into the ROM entering step a, back out entering step b+1.
            @test nswitch == 2
            @test steps == [a, b + 1]
            @test all(==(slot), slots)
            # The slot is a FOM again, by both the restart's bookkeeping and the
            # model actually sitting in it.
            @test ended_fom == true
            @test ended_model_is_fom == true
            # Finite, and in the neighbourhood of the untouched march -- the ROM
            # window introduces real error, so this is a sanity bound, not an
            # accuracy claim.
            @test isfinite(drift) && drift < 0.5
        end
    end
end
