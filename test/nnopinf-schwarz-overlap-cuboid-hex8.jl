# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Neural network Opinf Schwarz Overlap Static Cuboid Hex8 Same Step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-nn-rom/cuboid-1-for-test.yaml", "cuboid-1-for-test.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-nn-rom/cuboid-2-for-test.yaml", "cuboid-2-for-test.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-nn-rom/cuboids-for-test.yaml", "cuboids-for-test.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-nn-rom/neural-network-model", "neural-network-model"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids-for-test.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids-for-test.yaml")
    rm("cuboid-1-for-test.yaml")
    rm("cuboid-2-for-test.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    for file in readdir("neural-network-model")
      rm(joinpath("neural-network-model", file); force=true)  # Remove each file/directory
    end

    rm("neural-network-model"; force=true)
end
