# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Curved-interface regression for the nonoverlap DN Dirichlet transfer.
# A stiffer circular inclusion sits in a square matrix; the interface circle
# is discretized non-conformally (36 intervals on the inclusion side, 20 on
# the matrix side), so the two trace meshes are different chord polygons and
# the L2 projector does not map one side's reference geometry onto the
# other's. If the Dirichlet coupling transfers projected CURRENT POSITIONS
# (P(x + u) - x_dst), that geometric mismatch is injected as a spurious
# scalloped interface displacement at the coarse-facet frequency, with
# amplitude of the coarse-facet sagitta R(1 - cos(pi/20)) = 4.3e-3 -- as
# large as the applied loading (5.0e-3). Transferring projected
# DISPLACEMENTS (P u) leaves each side its own reference geometry; the two
# forms coincide on flat interfaces, which is why the cuboid tests are blind
# to this. Reported by Irina Tezaur on Giulia's 2D matrix-inclusion problem
# (2026-08-27): scalloped artifacts in the converged Schwarz solution that
# grow with interface non-conformality.
#
# The metric is the smoothness of the radial displacement along the
# inclusion's interface circle: the residual of a least-squares fit of low
# Fourier harmonics (k <= 4; the physics of this load is dominated by k = 0
# and 2, and any real mesh-induced content at k <= 4 is captured by the
# fit). Measured on this problem: position transfer gives a maximum residual
# of 3.4e-3 (the scallops at the 20-facet frequency), displacement transfer
# 2.2e-5. The threshold sits 9x above the fixed behavior and 17x below the
# defective one.
#
# The DN loop for a stiff inclusion fully embedded in a softer matrix has
# fixed-point gain > 1, so fixed relaxation diverges for every theta; the
# controller uses Aitken secant relaxation, which converges it.
@testset "Schwarz Nonoverlap Static Inclusion Curved Interface" begin
    src = "../examples/nonoverlap/static-same-step/matrix-inclusion/"
    cp(src * "inclusion.yaml", "inclusion.yaml"; force=true)
    cp(src * "circle.yaml", "circle.yaml"; force=true)
    cp(src * "square.yaml", "square.yaml"; force=true)
    cp(src * "circle.g", "circle.g"; force=true)
    cp(src * "square.g", "square.g"; force=true)
    sim = Norma.run("inclusion.yaml")
    rm("inclusion.yaml"; force=true)
    rm("circle.yaml"; force=true)
    rm("square.yaml"; force=true)
    rm("circle.g"; force=true)
    rm("square.g"; force=true)
    rm("circle.e"; force=true)
    rm("square.e"; force=true)
    model_square = sim.subsims[1].model
    model_circle = sim.subsims[2].model
    schwarz_bc_index = findfirst(
        b -> b isa Norma.SolidMechanicsNonOverlapSchwarzBoundaryCondition, model_circle.boundary_conditions
    )
    @test schwarz_bc_index !== nothing
    schwarz_bc = model_circle.boundary_conditions[schwarz_bc_index]
    interface_nodes = schwarz_bc.global_from_local_map
    angles = [atan(model_circle.reference[2, n], model_circle.reference[1, n]) for n in interface_nodes]
    radial_disp = [
        model_circle.displacement[1, n] * cos(angles[i]) + model_circle.displacement[2, n] * sin(angles[i]) for
        (i, n) in enumerate(interface_nodes)
    ]
    num_harmonics = 4
    fit_matrix = ones(length(angles), 2 * num_harmonics + 1)
    for k in 1:num_harmonics
        fit_matrix[:, 2 * k] = cos.(k .* angles)
        fit_matrix[:, 2 * k + 1] = sin.(k .* angles)
    end
    smoothness_residual = radial_disp - fit_matrix * (fit_matrix \ radial_disp)
    # Position transfer measures 3.4e-3 here; displacement transfer 2.2e-5.
    @test maximum(abs.(smoothness_residual)) < 2.0e-04
    # The geometric injection also inflates the interface motion itself:
    # 4.0e-3 with position transfer vs 1.2e-3 converged physics.
    @test maximum(abs.(radial_disp)) < 2.0e-03
    # Loading sanity: the matrix carries the prescribed 5.0e-3 boundary value.
    @test maximum(abs.(model_square.displacement)) ≈ 5.0e-03 rtol = 1.0e-06
end
