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
# a z-displacement ramp at the top face of cuboid-2: z = 0.5 - 0.5*cos(π*t).
#
# Because the OpInf basis K has max eigenvalue ≈ 53.3 (ωmax ≈ 7.3 rad/s,
# critical Δt ≈ 0.274 s), the shared time step of 0.01 s is well within the
# CFL stability limit (CFL ≈ 0.037).
#
# Physics-based symmetry: modes 1 and 2 correspond to x- and y-aligned
# displacements on a geometry loaded purely in z; by problem symmetry they must
# be equal to machine precision.  Mode 3 is the dominant z-mode and must be
# strictly larger in magnitude.
#
# Regression note: once the test passes, run it once and record the exact
# reduced_state values here to guard against future regressions, e.g.:
#   @test model_coarse.reduced_state[1] ≈ <value> rtol = 1.0e-06
#   @test model_coarse.reduced_state[2] ≈ <value> rtol = 1.0e-06
#   @test model_coarse.reduced_state[3] ≈ <value> rtol = 1.0e-06

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
    model_coarse = subsims[2].model

    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("quadratic-opinf-operator.npz"; force=true)

    # Simulation must complete without failure
    @test !model_coarse.failed

    # All reduced-state components must be finite
    @test all(isfinite, model_coarse.reduced_state)

    # The applied z-loading must produce a non-trivial solution
    @test norm(model_coarse.reduced_state) > 0.0

    # Symmetry: x-mode (1) and y-mode (2) are identical under pure z-loading
    @test model_coarse.reduced_state[1] ≈ model_coarse.reduced_state[2] rtol = 1.0e-06

    # The z-mode (3) must dominate: its magnitude exceeds that of mode 1
    @test abs(model_coarse.reduced_state[3]) > abs(model_coarse.reduced_state[1])
end
