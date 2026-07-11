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
function run_cantilever_imp_nc(adjoint_pairing::Bool; num_steps=5, explicit_free=false, free_dt="5.0e-07")
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_imp_nc_example/$f", f; force=true)
    end
    if !adjoint_pairing
        for f in ["cantilever-clamped.yaml", "cantilever-free.yaml"]
            doc = read(f, String)
            write(f, replace(doc, r"\n *adjoint pairing: true" => ""))
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
