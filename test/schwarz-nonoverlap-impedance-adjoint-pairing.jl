# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML
using LinearAlgebra

const cantilever_imp_nc_example = "../examples/nonoverlap/dynamic-same-step/cantilever-impedance-nonconforming"

# Run the nonconforming nonoverlap impedance cantilever for a few steps,
# optionally with adjoint pairing (the example enables it; pass false to
# strip the flag and get the legacy per-side transfer). With
# explicit_free=true the free subdomain runs central difference + explicit
# solver, exercising the IMEX interface treatment.
function run_cantilever_imp_nc(adjoint_pairing::Bool; num_steps=5, explicit_free=false, free_dt="5.0e-07", controller_dt=nothing, schwarz_tols=nothing, theta=nothing, aitken=false)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_imp_nc_example/$f", f; force=true)
    end
    if !adjoint_pairing
        # Pairing is the default; the legacy path needs an explicit opt-out.
        for f in ["cantilever-clamped.yaml", "cantilever-free.yaml"]
            doc = read(f, String)
            write(f, replace(doc, "adjoint pairing: true" => "adjoint pairing: false"))
        end
    end
    if explicit_free
        doc = read("cantilever-free.yaml", String)
        doc = replace(
            doc,
            "time integrator:\n  type: Newmark\n  β: 0.25\n  γ: 0.5\n  time step: 5.0e-07" =>
                "time integrator:\n  type: central difference\n  time step: $free_dt\n  CFL: 1.0\n  γ: 0.5",
        )
        doc = replace(
            doc,
            r"solver:\n  type: Hessian minimizer(\n  [^\n]+)+" =>
                "solver:\n  type: explicit solver\n  step: explicit",
        )
        write("cantilever-free.yaml", doc)
    end
    params = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cantilever-multi.yaml"
    params["final time"] = num_steps * 5.0e-07
    if controller_dt !== nothing
        params["time step"] = controller_dt
    end
    if schwarz_tols !== nothing
        params["relative tolerance"], params["absolute tolerance"] = schwarz_tols
    end
    if theta !== nothing
        params["relaxation parameter"] = theta
    end
    if aitken
        params["relaxation"] = "aitken recursive"
    end
    sim = Norma.run(params)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e",
              "cantilever-multi-energy.csv"]
        rm(f; force=true)
    end
    return sim
end

function impedance_bc_of(subsim)
    for bc in subsim.model.boundary_conditions
        if bc isa Norma.SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition
            return bc
        end
    end
    return nothing
end

@testset "Schwarz Nonoverlap Impedance Adjoint Pairing" begin
    sim = run_cantilever_imp_nc(true)
    @test sim.failed == false

    bc1 = impedance_bc_of(sim.subsims[1])
    bc2 = impedance_bc_of(sim.subsims[2])
    @test bc1 !== nothing && bc2 !== nothing
    @test bc1.adjoint_pairing && bc2.adjoint_pairing

    W1 = bc1.square_projector
    W2 = bc2.square_projector
    P1 = bc1.dirichlet_projector
    P2 = bc2.dirichlet_projector
    N1 = bc1.neumann_projector
    N2 = bc2.neumann_projector

    # The interface is genuinely nonconforming (different trace spaces).
    @test size(W1, 1) != size(W2, 1)

    # Discrete adjoint condition: W1*P1 = (W2*P2)' = B (single cross-mass).
    B1 = W1 * P1
    B2 = W2 * P2
    scale = max(norm(B1), norm(B2))
    @test norm(B1 - transpose(B2)) / scale < 1.0e-12

    # Force transfer is the adjoint of the partner's kinematic transfer.
    @test norm(N1 - transpose(P2)) / max(norm(N1), 1.0) < 1.0e-12
    @test norm(N2 - transpose(P1)) / max(norm(N2), 1.0) < 1.0e-12

    # Kinematic transfers reproduce constants (partition of unity) on this
    # flat, fully covered interface.
    @test maximum(abs.(P1 * ones(size(P1, 2)) .- 1.0)) < 1.0e-9
    @test maximum(abs.(P2 * ones(size(P2, 2)) .- 1.0)) < 1.0e-9

    # Shared pair impedance (equal materials: harmonic mean = one-sided value)
    # and shared Robin parameter on both sides.
    @test bc1.impedance ≈ bc2.impedance
    @test bc1.robin_parameter ≈ bc2.robin_parameter

    # The paired transfer preserves the total of a transferred nodal force
    # (conservation certificate of the shared cross-mass): a uniform unit
    # traction integrates to the same total force through either side.
    @test sum(W1 * ones(size(W1, 1))) ≈ sum(W2 * ones(size(W2, 1))) rtol = 1.0e-9

    # Explicit (central difference) free subdomain under pairing: the IMEX
    # interface treatment (implicit interface rows in the acceleration
    # update) keeps the paired consistent-traction exchange stable — it is
    # violently unstable when the interface rows are updated explicitly.
    sim_imex = run_cantilever_imp_nc(true; explicit_free=true, num_steps=20)
    @test sim_imex.failed == false
    u_imex = sim_imex.subsims[1].model.displacement
    @test all(isfinite, u_imex)
    @test maximum(abs.(u_imex)) < 0.1

    # Subcycled (multi-time-step) pairing: the explicit free subdomain runs
    # four substeps per Schwarz stop; the paired consistent traction is
    # evaluated on the partner's piecewise-linearly time-interpolated state
    # (the Gravouil-Combescure interpolated-traction structure).
    sim_sub = run_cantilever_imp_nc(true; explicit_free=true, free_dt="1.25e-07", num_steps=20)
    @test sim_sub.failed == false
    u_sub = sim_sub.subsims[1].model.displacement
    @test all(isfinite, u_sub)
    @test maximum(abs.(u_sub)) < 0.1

    # Windowed controller stop (both subdomains subcycle five substeps per
    # stop): the per-substep-time relaxation slots keep the windowed exchange
    # waveform-consistent, so the converged trajectory must match the
    # same-step run — the controller step is a cost knob, not a physics knob.
    # Before the slotted state, the per-pair relaxation blended iterates
    # across time (a causal low-pass on the exchanged traction) and the
    # windowed run drifted from the first stop.
    # The residual difference is Schwarz-truncation noise (it collapses as the
    # stopping tolerances tighten), so all runs use tightened tolerances to
    # keep the threshold crisp; the pre-fix per-pair state fails this by two
    # orders of magnitude regardless of tolerance (its windowed fixed point
    # genuinely differs). All three runs use relaxation parameter 1.0 so a
    # single θ semantics is compared: the converged paired-impedance solution
    # is measurably θ-sensitive (~4e-5 relative on this benchmark at 1e-13
    # tolerance, same-step and windowed alike — a pre-existing neutral mode
    # of the exchange, under separate investigation), and this block tests
    # windowing, not that sensitivity. θ = 1 is also the measured optimum for
    # windowed impedance stops (the dashpot supplies the contraction).
    tols = (1.0e-10, 1.0e-12)
    sim_ss = run_cantilever_imp_nc(true; num_steps=10, schwarz_tols=tols, theta=1.0)
    sim_win = run_cantilever_imp_nc(true; num_steps=10, controller_dt=2.5e-06, schwarz_tols=tols, theta=1.0)
    # With Aitken configured, windowed stops must fall back to the fixed
    # relaxation parameter (Aitken accelerates only single-slot stops; every
    # windowed Aitken policy measured on the 10 ms benchmark dies or loses to
    # fixed θ) and land on the same trajectory.
    sim_win_ait = run_cantilever_imp_nc(true; num_steps=10, controller_dt=2.5e-06, schwarz_tols=tols, theta=1.0, aitken=true)
    for i in 1:2
        u_ss = sim_ss.subsims[i].model.displacement
        u_win = sim_win.subsims[i].model.displacement
        @test norm(u_win - u_ss) / norm(u_ss) < 1.0e-6
        u_win_ait = sim_win_ait.subsims[i].model.displacement
        @test norm(u_win_ait - u_ss) / norm(u_ss) < 1.0e-6
    end

    # Temporal transfer pair of the subcycled pairing: the finer side receives
    # the piecewise-linear interpolant of the partner history (on the query's
    # own segment), and the coarser side receives the endpoint value of the
    # least-squares linear fit of the partner trajectory over its step window
    # — lag-free like endpoint sampling, filtering intra-window content
    # instead of aliasing it (the Schwarz-compatible form of the
    # Prakash-Hjelmstad linear-in-time interface traction).
    th2 = [0.0, 1.0]
    vh2 = [[0.0], [10.0]]
    @test Norma.interpolate(th2, vh2, 0.25) ≈ [2.5]        # not the endpoint
    th5 = collect(0.0:0.25:1.0)
    kinked = [[abs(t - 0.5)] for t in th5]
    @test Norma.interpolate(th5, kinked, 0.30) ≈ [0.2]     # own segment, no sign flip
    linear = [[t, 2.0 * t] for t in th5]
    @test Norma.time_average(th5, linear, 0.0, 1.0) ≈ [0.5, 1.0]
    @test Norma.time_average(th5, linear, 0.5, 0.5) ≈ [0.5, 1.0]      # degenerate window
    @test Norma.time_endpoint_fit(th5, linear, 0.0, 1.0) ≈ [1.0, 2.0]  # exact endpoint, no lag
    @test Norma.time_endpoint_fit(th5, linear, 0.1, 0.9) ≈ [0.9, 1.8]  # off-grid window
    @test Norma.time_endpoint_fit(th5, kinked, 0.0, 1.0) ≈ [0.25]      # symmetric kink filtered flat
    @test Norma.time_endpoint_fit(th5, linear, 0.5, 0.5) ≈ [0.5, 1.0]  # degenerate window

    # Mismatched Robin parameters abort under pairing (strict by design: the
    # Robin spring is conservative only when both sides' interface forces
    # pair through the same coefficient).
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_imp_nc_example/$f", f; force=true)
    end
    doc = read("cantilever-clamped.yaml", String)
    write("cantilever-clamped.yaml", replace(doc, "robin parameter: 2.0e+09" => "robin parameter: 4.0e+09"))
    params_bad = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
    params_bad["name"] = "cantilever-multi.yaml"
    params_bad["final time"] = 5.0e-07
    Norma.NORMA_TEST_MODE[] = true
    try
        @test_throws Norma.NormaAbortException Norma.run(params_bad)
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e",
              "cantilever-multi-energy.csv"]
        rm(f; force=true)
    end

    # `impedance scale` on the nonoverlap variant scales the dashpot; 0 is the
    # explicit pure-Robin opt-in (t + α W u = g, issue #176 warns) and must
    # produce a zero pair impedance without dividing by zero in the harmonic
    # mean. Both runs are short; pure Robin is not exercised at 1 ms here.
    for (scale, factor) in (("0.5", 0.5), ("0.0", 0.0))
        for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
                  "cantilever-clamped.g", "cantilever-free.g"]
            cp("$cantilever_imp_nc_example/$f", f; force=true)
        end
        for f in ["cantilever-clamped.yaml", "cantilever-free.yaml"]
            doc = read(f, String)
            write(f, replace(doc, "adjoint pairing: true" => "adjoint pairing: true\n      impedance scale: $scale"))
        end
        params_scaled = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
        params_scaled["name"] = "cantilever-multi.yaml"
        params_scaled["final time"] = 3 * 5.0e-07
        sim_scaled = Norma.run(params_scaled)
        @test sim_scaled.failed == false
        bc_scaled = impedance_bc_of(sim_scaled.subsims[1])
        @test bc_scaled.impedance ≈ factor * bc1.impedance
        for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
                  "cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e",
                  "cantilever-multi-energy.csv"]
            rm(f; force=true)
        end
    end

    # Legacy per-side transfer still runs on the same meshes (opt-out path).
    # On this 2:1 NESTED interface the legacy heuristics happen to select
    # adjoint operators too, so no non-adjointness assertion is made here;
    # the drift studies document the difference on non-nested interfaces.
    sim_legacy = run_cantilever_imp_nc(false)
    @test sim_legacy.failed == false
    bc1l = impedance_bc_of(sim_legacy.subsims[1])
    bc2l = impedance_bc_of(sim_legacy.subsims[2])
    @test !bc1l.adjoint_pairing && !bc2l.adjoint_pairing
end
