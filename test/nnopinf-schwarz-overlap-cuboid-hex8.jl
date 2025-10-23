# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "NN Opinf Schwarz Overlap Dynamic Cuboid Hex8 Same Step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/ml-models", "ml-models"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-2.g", "cuboid-2.g"; force=true)
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
    # Remove all files and subdirectories in the directory
    for file in readdir("ml-models")
        rm(joinpath("ml-models", file); force=true)
    end
    # Now remove the empty directory
    rm("ml-models")
    @test model_coarse.reduced_state[1] ≈ -0.024143687027586 rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ -0.014289943951844604 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ 0.06875121953836646 rtol = 1.0e-06
end
