# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Quasistatic Dirichlet-Neumann cuboid pair, pulled to unit displacement on the
# far face. The reference values are shared with the overlapping decomposition
# in schwarz-ahead-overlap-quasistatic-cuboid.jl: the same physical problem is
# solved there with a different partition, so both must land on the same
# transverse contraction and the same axial stress.
@testset "Schwarz AHeaD Non-Overlap Quasistatic Cuboid HEX8-HEX8" begin
    src = "../examples/ahead/nonoverlap/cuboid"
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml"]
        cp("$src/quasistatic/$f", f; force=true)
    end
    cp("$src/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("$src/cuboid-2.g", "../cuboid-2.g"; force=true)
    input_file = "cuboids.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["name"] = input_file
    max_iters = params["maximum iterations"]
    sim = Norma.run(params)
    model_cuboid1 = sim.subsims[1].model
    model_cuboid2 = sim.subsims[2].model

    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)

    @test sim.failed == false

    # Every step must reach the tolerance rather than run out of iterations. The
    # cap was 16 while the hardest step needs 43, so every step used to end
    # unconverged and nothing said so.
    @test maximum(sim.controller.schwarz_iters) < max_iters
    @test sim.controller.schwarz_iters ≈ [43, 38, 32, 27, 23, 20, 19, 17, 16, 15] atol = 0

    avg_stress_cuboid1 = average_components(model_cuboid1.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    # Transverse contraction, identical in both subdomains and independent of
    # how the domain is partitioned.
    @test minimum(model_cuboid1.displacement[1, :]) ≈ -0.11889377666112658 rtol = 1.0e-6
    @test minimum(model_cuboid1.displacement[2, :]) ≈ -0.11889377666112656 rtol = 1.0e-6
    @test minimum(model_cuboid2.displacement[1, :]) ≈ -0.11889377701070683 rtol = 1.0e-6
    @test minimum(model_cuboid2.displacement[2, :]) ≈ -0.11889377701070669 rtol = 1.0e-6

    # Axial stretch: the interface sits at the half height, the far face at 1.
    @test maximum(model_cuboid1.displacement[3, :]) ≈ 0.5 atol = 1.0e-8
    @test maximum(model_cuboid2.displacement[3, :]) ≈ 1.0 atol = 1.0e-8

    # Uniaxial state: one axial component, everything else negligible beside it.
    @test avg_stress_cuboid1[3] ≈ 5.210645149300675e8 rtol = 1.0e-6
    @test avg_stress_cuboid2[3] ≈ 5.2106451505573684e8 rtol = 1.0e-6
    @test all(abs.(avg_stress_cuboid1[[1, 2, 4, 5, 6]]) .< 1.0)
    @test all(abs.(avg_stress_cuboid2[[1, 2, 4, 5, 6]]) .< 1.0)
end
