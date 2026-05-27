# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

@testset "Opinf Schwarz Overlap L2 Error Dynamic Cuboid Hex8 Same Step" begin
    example_dir = "../examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error"
    input_files = (
        "cuboid-1.yaml",
        "cuboid-2.yaml",
        "cuboids.yaml",
        "opinf-operator.npz",
        "cuboid-1.g",
        "cuboid-2.g",
    )

    for file in input_files
        cp(joinpath(example_dir, file), file; force=true)
    end

    try
        sim = Norma.create_simulation("cuboids.yaml")
        bc_fine = only([
            bc for bc in sim.subsims[1].model.boundary_conditions if
            bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
        ])
        bc_coarse = only([
            bc for bc in sim.subsims[2].model.boundary_conditions if
            bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
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
            bc for bc in sim.subsims[1].model.boundary_conditions if
            bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
        ])
        bc_coarse = only([
            bc for bc in sim.subsims[2].model.boundary_conditions if
            bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
        ])
        @test isfinite(bc_fine.overlap_l2_error)
        @test isfinite(bc_coarse.overlap_l2_error)
        @test bc_fine.overlap_l2_error >= 0.0
        @test bc_coarse.overlap_l2_error >= 0.0
        @test isfile("overlap-l2-rel-errors-0001.csv")

        overlap_csv = read("overlap-l2-rel-errors-0001.csv", String)
        @test occursin("domain,side_set,overlap_l2_error", overlap_csv)
        @test occursin("cuboid-1,ssz+", overlap_csv)
        @test occursin("cuboid-2,ssz-", overlap_csv)
    finally
        for file in input_files
            rm(file; force=true)
        end
        rm("cuboid-1.e"; force=true)
        rm("cuboid-2.e"; force=true)
        for file in readdir()
            if startswith(file, "iterations-") || startswith(file, "overlap-l2-rel-errors-")
                rm(file; force=true)
            end
        end
    end
end
