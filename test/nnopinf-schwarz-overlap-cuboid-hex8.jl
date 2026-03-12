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
    @test model_coarse.reduced_state[1] ≈  -0.003507375025622216  rtol = 1.0e-06
    @test model_coarse.reduced_state[2] ≈ -0.002114866440831956 rtol = 1.0e-06
    @test model_coarse.reduced_state[3] ≈ 0.009165290330522333 rtol = 1.0e-06
end

@testset "NN Opinf Schwarz Overlap L2 Error Dynamic Cuboid Hex8 Same Step" begin
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/ml-models", "ml-models"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/ahead/overlap/cuboid/dynamic-nn-opinf-rom/cuboid-2.g", "cuboid-2.g"; force=true)

    sim_default = Norma.create_simulation("cuboids.yaml")
    bc_default = only([
        bc for bc in sim_default.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOpInfOverlapSchwarzBoundaryCondition
    ])
    @test bc_default isa Norma.SolidMechanicsOpInfOverlapSchwarzBoundaryCondition
    @test bc_default.fom_bc.compute_overlap_l2_error == false
    Exodus.close(sim_default.subsims[1].params["input_mesh"])
    Exodus.close(sim_default.subsims[1].params["output_mesh"])
    Exodus.close(sim_default.subsims[2].params["input_mesh"])
    Exodus.close(sim_default.subsims[2].params["output_mesh"])

    cuboids_text = read("cuboids.yaml", String)
    write("cuboids.yaml", replace(cuboids_text, "absolute tolerance: 1.0e-08" => "absolute tolerance: 1.0e-08\nCSV output interval: 1"))
    cuboid_2_text = read("cuboid-2.yaml", String)
    write(
        "cuboid-2.yaml",
        replace(
            cuboid_2_text,
            "source side set: ssz+" => "source side set: ssz+\n      compute overlap L2 error: true",
        ),
    )

    sim = Norma.run("cuboids.yaml")
    bc = only([
        bc for bc in sim.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOpInfOverlapSchwarzBoundaryCondition
    ])
    @test bc isa Norma.SolidMechanicsOpInfOverlapSchwarzBoundaryCondition
    @test bc.fom_bc.compute_overlap_l2_error == true
    @test length(bc.fom_bc.overlap_node_indices) > length(unique(bc.fom_bc.side_set_node_indices))
    @test isfinite(bc.fom_bc.overlap_l2_error)
    @test bc.fom_bc.overlap_l2_error >= 0.0
    @test Norma.get_overlap_l2_error(bc) == bc.fom_bc.overlap_l2_error
    @test isfile("overlap-l2-errors-0001.csv")

    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    for file in readdir()
        if startswith(file, "iterations-") || startswith(file, "overlap-l2-errors-")
            rm(file)
        end
    end
    for file in readdir("ml-models")
        rm(joinpath("ml-models", file); force=true)
    end
    rm("ml-models")
end
