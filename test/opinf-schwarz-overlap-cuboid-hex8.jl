@testset "schwarz-overlap-static-cuboid-hex8-same-step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-1.yaml", "cuboid-1.yaml", force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-2.yaml", "cuboid-2.yaml", force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboids.yaml", "cuboids.yaml", force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/opinf-operator.npz", "opinf-operator.npz", force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-1.g", "cuboid-1.g", force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-2.g", "cuboid-2.g", force=true)
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
    rm("opinf-operator.npz")
    print(model_coarse.reduced_state[1])
    @test model_coarse.reduced_state[1] ≈ 0.04374576440784916 rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ 0.04374576440786179 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ 0.19557108653526997 rtol = 1.0e-06
end
