# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Cross-mesh single-domain swap: phase 1 is a hex8 cube and phase 2 is a tet4
# cube on the same physical [-0.5, 0.5]^3 domain.  Both phases are linear
# elastic with identical material parameters, so the post-swap solution
# should agree with the standard hex8 cube uniaxial-tension baseline.
#
# Neither YAML enables `stress recovery`, so the swap exercises the
# auto-build path that constructs a consistent recovery on demand to do the
# L2 transfer of the kinematic state.
@testset "Single Static Solid Cube Mid-Run Swap (cross-mesh)" begin
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube-hex8.g"; force=true)
    cp("../examples/element-types/tet4/cube/cube.g", "cube-tet4.g"; force=true)
    cp("../examples/single/static-solid/cube-swap/cross-mesh/cube-hex8.yaml", "cube-hex8.yaml"; force=true)
    cp("../examples/single/static-solid/cube-swap/cross-mesh/cube-tet4.yaml", "cube-tet4.yaml"; force=true)
    sim = Norma.run("cube-hex8.yaml")
    rm("cube-hex8.yaml"; force=true)
    rm("cube-tet4.yaml"; force=true)
    rm("cube-hex8.g"; force=true)
    rm("cube-tet4.g"; force=true)
    rm("cube-hex8.e"; force=true)
    rm("cube-tet4.e"; force=true)

    # Replacement was substituted in place; phase-2 has no swaps of its own.
    @test sim.name == "cube-tet4"
    @test isempty(sim.swaps)
    @test sim.failed == false
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-09

    # Both phases are linear elastic with the same parameters under uniaxial
    # tension; the final state is the exact linear solution
    #   u_z(Z) = Z + 0.5,  u_x(X) = -ν (X + 0.5),  u_y(Y) = -ν (Y + 0.5).
    # Averages depend on the mesh node distribution (tet4 nodes are not on
    # the same uniform grid as hex8), so we compute them analytically from
    # the reference coordinates rather than hard-coding hex8 values.
    avg_disp = average_components(sim.integrator.displacement)
    avg_stress = average_components(sim.model.stress)
    mean_x = mean(sim.model.reference[1, :])
    mean_y = mean(sim.model.reference[2, :])
    mean_z = mean(sim.model.reference[3, :])
    @test avg_disp[1] ≈ -0.25 * (mean_x + 0.5) rtol = 1.0e-3
    @test avg_disp[2] ≈ -0.25 * (mean_y + 0.5) rtol = 1.0e-3
    @test avg_disp[3] ≈ mean_z + 0.5 rtol = 1.0e-3
    @test avg_stress[3] ≈ 1.0e9 rtol = 1.0e-3
end
