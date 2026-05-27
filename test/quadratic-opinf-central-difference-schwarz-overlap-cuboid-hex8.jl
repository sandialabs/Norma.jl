# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests the evaluate(::RomCentralDifference, ::RomExplicitSolver, ::QuadraticOpInfRom)
# specialization added to src/opinf/opinf_solver.jl.
#
# Setup: cuboid-1 (fine FOM, Newmark) and cuboid-2 (coarse QuadraticOpInfRom,
# central difference) are coupled via Schwarz overlap.  The applied loading is
# a z-displacement ramp at the top face of cuboid-2: z = 0.5 - 0.5*cos(π*t),
# which reaches exactly 1.0 at t = 1.0.
#
# Because the OpInf basis K has max eigenvalue ≈ 53.3 (ωmax ≈ 7.3 rad/s,
# critical Δt ≈ 0.274 s), the shared time step of 0.01 s is well within the
# CFL stability limit (CFL ≈ 0.037).
#
# Checked quantities (mirrors schwarz-ahead-nonoverlap-dynamic-cuboid-cubic-rom-rom.jl):
#   • min/max displacement components for both subdomains
#   • volume-averaged stress components for both subdomains
#   • Schwarz iteration counts per time step
#   • final reduced-state vector of the QuadraticOpInfRom
#
# Physics-based symmetry: modes 1 and 2 correspond to x- and y-aligned
# displacements on a geometry loaded purely in z; by problem symmetry they must
# be equal to machine precision.  Mode 3 is the dominant z-mode and must be
# strictly larger in magnitude.
#
# Regression note: the exact values below were recorded from the first passing
# run and must be updated whenever the operator file, mesh, or integrator
# parameters change.

@testset "Quadratic OpInf FOM-ROM Implicit-Explicit Schwarz Overlap Dynamic Cuboid Hex8" begin
    example_dir = "../examples/ahead/overlap/cuboid/dynamic-quadratic-opinf-rom-explicit"
    cp("$example_dir/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("$example_dir/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("$example_dir/cuboids.yaml", "cuboids.yaml"; force=true)
    cp(
        "$example_dir/quadratic-opinf-operator.npz",
        "quadratic-opinf-operator.npz";
        force=true,
    )
    cp("$example_dir/../cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("$example_dir/../cuboid-2.g", "../cuboid-2.g"; force=true)

    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine   = subsims[1].model   # SolidMechanics  (Newmark FOM)
    model_coarse = subsims[2].model   # QuadraticOpInfRom (central-difference ROM)

    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("quadratic-opinf-operator.npz"; force=true)

    # ── displacement statistics ───────────────────────────────────────────────
    # cuboid-1: FOM – displacement stored directly on model
    min_disp_x_cuboid1 = minimum(model_fine.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_fine.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_fine.displacement[3, :])

    # cuboid-2: ROM – physical displacement reconstructed in fom_model
    min_disp_x_cuboid2 = minimum(model_coarse.fom_model.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_coarse.fom_model.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_coarse.fom_model.displacement[3, :])

    # ── volume-averaged stress ────────────────────────────────────────────────
    avg_stress_cuboid1 = average_components(model_fine.stress)
    avg_stress_cuboid2 = average_components(model_coarse.fom_model.stress)

    # ── simulation must complete without failure ──────────────────────────────
    @test !model_coarse.failed

    # ── reduced-state sanity ──────────────────────────────────────────────────
    @test all(isfinite, model_coarse.reduced_state)
    @test norm(model_coarse.reduced_state) > 0.0

    # Symmetry: x-mode (1) and y-mode (2) are equal under pure z-loading
    @test model_coarse.reduced_state[1] ≈ model_coarse.reduced_state[2] rtol = 1.0e-06

    # The z-mode (3) must dominate the x/y modes
    @test abs(model_coarse.reduced_state[3]) > abs(model_coarse.reduced_state[1])

    # Regression: exact reduced-state values at t = 1.0
    @test model_coarse.reduced_state[1] ≈ 0.3282650287668367 rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ 0.3282650287657078 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ 2.0015111660534566 rtol = 1.0e-06

    # ── displacement regression ───────────────────────────────────────────────
    # cuboid-1 (FOM): Poisson contraction in x/y; z pushed up by Schwarz BC
    @test min_disp_x_cuboid1 ≈ -0.11926950567829919  atol = 1.0e-06
    @test min_disp_y_cuboid1 ≈ -0.11926950567794792 atol = 1.0e-06
    @test max_disp_z_cuboid1 ≈  0.6065838584133116 atol = 1.0e-06

    # cuboid-2 (ROM, fom_model): top face prescribed; z at top = 1.0 exactly
    @test min_disp_x_cuboid2 ≈ -0.11986552696734587 atol = 1.0e-06
    @test min_disp_y_cuboid2 ≈ -0.11986552696693364 atol = 1.0e-06
    @test max_disp_z_cuboid2 ≈  1.0                  atol = 1.0e-08

    # ── stress regression ─────────────────────────────────────────────────────
    @test avg_stress_cuboid1 ≈
        [-1.1863325261022714e6 -1.186332526035554e6 5.0812693386426896e8 -607310.0290748979 -607310.0291353256 -42857.96730489226] atol =
        1.0e1
    @test avg_stress_cuboid2 ≈
        [3.191982793806081e6 3.191982793961093e6 5.287770940238085e8 -661354.8217173232 -661354.821773954 -14387.794002634266] atol =
        1.0e1

    # ── Schwarz convergence history ───────────────────────────────────────────
    # Overlap coupling converges well within the 16-iteration cap.
    @test sim.controller.schwarz_iters ≈
        [2, 5, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 6, 8, 9, 9, 9, 9, 9] atol = 0
end
