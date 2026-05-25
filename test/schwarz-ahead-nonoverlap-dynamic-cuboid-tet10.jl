# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Regression test for TET10-TET10 non-overlapping Schwarz coupling.
#
# This test exercises the fix for the SingularException that arose when
# compute_neumann_projector called get_square_projection_matrix with TRI6
# face elements (the faces of TET10 volumes).  The old default of 4
# quadrature points produced a rank-deficient (rank ≤ 4) 6×6 element mass
# matrix H = Σ_p w_p N_p N_p^T, making H \ I fail.  The fix raises the
# TRI6 default to the 7-point Dunavant degree-5 rule, which guarantees a
# full-rank, positive-definite H for any mesh topology.
#
# Geometry (z-axis is the coupling direction):
#   cuboid-1 (fine, TET10):   bottom fixed (nsz-), top coupled (ssz+)
#   cuboid-2 (coarse, TET10): top driven   (nsz+), bottom coupled (ssz-)
#
# Schwarz coupling type: Dirichlet-Neumann non-overlapping
#   cuboid-1 receives Dirichlet displacement from cuboid-2 at ssz+
#   cuboid-2 receives Neumann force       from cuboid-1 at ssz-
#
# Applied displacement BC on cuboid-2 top (nsz+):
#   u_z(t) = 0.25 - 0.25 * cos(π t)
#   At t = 1.0:  u_z = 0.25 - 0.25 * cos(π) = 0.5  (exact, Dirichlet)
#
# Subsim ordering follows the domains list in cuboids.yaml:
#   domains: ["cuboid-2.yaml", "cuboid-1.yaml"]
#   → subsims[1] = cuboid-2 (coarse, top, driven)
#   → subsims[2] = cuboid-1 (fine,   bottom, fixed at base)

const tet10_example_dir = "../examples/ahead/nonoverlap/cuboid"
const tet10_example_subdir = "$tet10_example_dir/dynamic-tet10"

@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid TET10-TET10 FOM-FOM" begin
    # --- set up: copy input files into the test working directory -----------
    cp("$tet10_example_subdir/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("$tet10_example_subdir/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("$tet10_example_subdir/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    # Mesh files are shared with other tet10 examples one level up.
    # The YAML files reference them as "../cuboid-tet10-{1,2}.g", i.e.
    # relative to the test/ directory → they go in the repo root.
    cp("$tet10_example_dir/cuboid-tet10-1.g", "../cuboid-tet10-1.g"; force=true)
    cp("$tet10_example_dir/cuboid-tet10-2.g", "../cuboid-tet10-2.g"; force=true)

    # --- run ----------------------------------------------------------------
    input_file = "cuboids.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["name"] = input_file
    sim = Norma.run(params)

    # --- clean up -----------------------------------------------------------
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-tet10-1.g"; force=true)
    rm("../cuboid-tet10-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)

    # --- extract results ----------------------------------------------------
    subsims = sim.subsims
    # domains list: ["cuboid-2.yaml", "cuboid-1.yaml"]
    u_coarse = subsims[1].model.displacement  # cuboid-2: coarse, top, driven
    u_fine   = subsims[2].model.displacement  # cuboid-1: fine, bottom, fixed

    # --- physics-based assertions -------------------------------------------

    # 1. Simulation completed without failure.
    @test sim.failed == false

    # 2. Dirichlet BC on top of coarse domain (cuboid-2, nsz+) is satisfied.
    #    u_z(t=1) = 0.25 - 0.25*cos(π) = 0.5
    @test maximum(u_coarse[3, :]) ≈ 0.5 atol = 1.0e-6

    # 3. Dirichlet BC on bottom of fine domain (cuboid-1, nsz-) is satisfied.
    #    u_z = 0 (fixed).
    @test minimum(u_fine[3, :]) ≈ 0.0 atol = 1.0e-6

    # 4. Interface displacement continuity in z.
    #    The top of the fine domain (max u_fine_z) receives u_z from the
    #    coarse domain as a Dirichlet BC, so it should match the bottom of
    #    the coarse domain (min u_coarse_z).  Loose tolerance: the Neumann
    #    side sees a projected force, not exact displacement continuity.
    @test maximum(u_fine[3, :]) ≈ minimum(u_coarse[3, :]) rtol = 1.0e-3

    # 5. Schwarz iterations are bounded (solver did not stall completely).
    @test all(sim.controller.schwarz_iters .≤ params["maximum iterations"])
end
