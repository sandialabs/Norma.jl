# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Kernel ROM-ROM Schwarz Overlap Dynamic Cuboid Hex8 Same Step" begin
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/cuboid-1.yaml",
        "cuboid-1.yaml";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/cuboid-2.yaml",
        "cuboid-2.yaml";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/cuboids.yaml",
        "cuboids.yaml";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/linear-kernel-rom-1.npz",
        "linear-kernel-rom-1.npz";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/linear-matern-kernel-rom-2.npz",
        "linear-matern-kernel-rom-2.npz";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/cuboid-1.g",
        "cuboid-1.g";
        force = true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-kernel-rom-rom/cuboid-2.g",
        "cuboid-2.g";
        force = true,
    )
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml"; force = true)
    rm("cuboid-1.yaml"; force = true)
    rm("cuboid-2.yaml"; force = true)
    rm("cuboid-1.g"; force = true)
    rm("cuboid-2.g"; force = true)
    rm("cuboid-1.e"; force = true)
    rm("cuboid-2.e"; force = true)
    rm("linear-kernel-rom-1.npz"; force = true)
    rm("linear-matern-kernel-rom-2.npz"; force = true)

    @test model_coarse.reduced_state[1] ≈ -0.04408941424927044 rtol = 1.0e-05
    @test model_coarse.reduced_state[2] ≈ -0.0031640397058258734 rtol = 1.0e-05
    @test model_fine.reduced_state[1] ≈ -0.06992403693639167 rtol = 1.0e-05
    @test model_fine.reduced_state[2] ≈ -0.005910211088836505 rtol = 1.0e-05
end
