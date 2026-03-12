# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Opinf Schwarz Overlap Dynamic Cuboid Hex8 Same Step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/opinf-operator.npz", "opinf-operator.npz"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-opinf-rom/cuboid-2.g", "cuboid-2.g"; force=true)
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
    @test model_coarse.reduced_state[1] ≈ -1.235577703737587e-5 rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ -1.235577703737587e-5 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ -0.0008960511947529225 rtol = 1.0e-06
end

@testset "Opinf Schwarz Overlap L2 Error Dynamic Cuboid Hex8 Same Step" begin
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboid-1.yaml",
        "cuboid-1.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboid-2.yaml",
        "cuboid-2.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboids.yaml",
        "cuboids.yaml";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/opinf-operator.npz",
        "opinf-operator.npz";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboid-1.g",
        "cuboid-1.g";
        force=true,
    )
    cp(
        "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboid-2.g",
        "cuboid-2.g";
        force=true,
    )

    sim = Norma.create_simulation("cuboids.yaml")
    bc_fine = only([
        bc for bc in sim.subsims[1].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_coarse = only([
        bc for bc in sim.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test bc_fine.compute_overlap_l2_error == true
    @test bc_coarse.compute_overlap_l2_error == true
    @test !isempty(bc_fine.overlap_node_indices)
    @test !isempty(bc_coarse.overlap_node_indices)
    @test length(bc_fine.overlap_node_indices) > length(unique(bc_fine.side_set_node_indices))
    Exodus.close(sim.subsims[1].params["input_mesh"])
    Exodus.close(sim.subsims[1].params["output_mesh"])
    Exodus.close(sim.subsims[2].params["input_mesh"])
    Exodus.close(sim.subsims[2].params["output_mesh"])

    sim = Norma.run("cuboids.yaml")
    bc_fine = only([
        bc for bc in sim.subsims[1].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_coarse = only([
        bc for bc in sim.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test isfinite(bc_fine.overlap_l2_error)
    @test isfinite(bc_coarse.overlap_l2_error)
    @test bc_fine.overlap_l2_error >= 0.0
    @test bc_coarse.overlap_l2_error >= 0.0
    @test isfile("overlap-l2-errors-0001.csv")

    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    rm("opinf-operator.npz")
    for file in readdir()
        if startswith(file, "iterations-") || startswith(file, "overlap-l2-errors-")
            rm(file)
        end
    end
end
