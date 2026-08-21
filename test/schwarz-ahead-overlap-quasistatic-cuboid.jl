# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Quasistatic overlapping cuboid pair, the same physical problem that
# schwarz-ahead-nonoverlap-quasistatic-cuboid.jl partitions without an overlap.
# The transverse contraction and the axial stress must agree between the two,
# since neither depends on how the domain was cut.
@testset "Schwarz AHeaD Overlap Quasistatic Cuboid HEX8-HEX8" begin
    src = "../examples/ahead/overlap/cuboid"
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
    # cap was 16 while every step needs 19 to 21, so every step used to end
    # unconverged and nothing said so.
    @test maximum(sim.controller.schwarz_iters) < max_iters
    @test sim.controller.schwarz_iters ≈ [19, 20, 21, 21, 21, 21, 21, 21, 21, 19] atol = 0

    avg_stress_cuboid1 = average_components(model_cuboid1.stress)
    avg_stress_cuboid2 = average_components(model_cuboid2.stress)

    # Transverse contraction, matching the non-overlapping decomposition.
    @test minimum(model_cuboid1.displacement[1, :]) ≈ -0.11889377676878005 rtol = 1.0e-6
    @test minimum(model_cuboid1.displacement[2, :]) ≈ -0.11889377676878016 rtol = 1.0e-6
    @test minimum(model_cuboid2.displacement[1, :]) ≈ -0.11889377673325015 rtol = 1.0e-6
    @test minimum(model_cuboid2.displacement[2, :]) ≈ -0.11889377673325018 rtol = 1.0e-6

    # Axial stretch: the subdomains overlap, so the first reaches past the half
    # height, and the second still ends at the pulled face.
    @test maximum(model_cuboid1.displacement[3, :]) ≈ 0.625 atol = 1.0e-8
    @test maximum(model_cuboid2.displacement[3, :]) ≈ 1.0 atol = 1.0e-8

    # Uniaxial state, and the same axial stress the non-overlapping run gives.
    @test avg_stress_cuboid1[3] ≈ 5.210645142916895e8 rtol = 1.0e-6
    @test avg_stress_cuboid2[3] ≈ 5.2106451526941085e8 rtol = 1.0e-6
    @test all(abs.(avg_stress_cuboid1[[1, 2, 4, 5, 6]]) .< 1.0)
    @test all(abs.(avg_stress_cuboid2[[1, 2, 4, 5, 6]]) .< 1.0)
end
