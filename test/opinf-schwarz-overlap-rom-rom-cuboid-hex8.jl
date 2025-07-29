# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Opinf Schwarz Overlap Dynamic Cuboid Hex8 Same Step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/opinf-operator-1.npz",
        "opinf-operator-1.npz";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/opinf-operator-2.npz",
        "opinf-operator-2.npz";
        force=true,
    )
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom-rom/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    rm("opinf-operator-1.npz")
    rm("opinf-operator-2.npz")

    @test model_coarse.reduced_state[1] ≈ -0.0008153342945627467 rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ -0.0008153342949887136 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ 0.19986541967787202 rtol = 1.0e-06
    @test model_fine.reduced_state[1] ≈ -0.0061291008333238 rtol = 1.0e-06
    @test model_fine.reduced_state[2] ≈ -0.006129100833376133 rtol = 1.0e-06
    @test model_fine.reduced_state[3] ≈ -0.34397149691049 rtol = 1.0e-06
end
