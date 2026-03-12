# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Schwarz Overlap Static Cuboid Hex8 Same Step" begin
    cp("../examples/overlap/static-same-step/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
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
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    min_disp_y_fine = minimum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    min_disp_x_coarse = minimum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    min_disp_y_coarse = minimum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    min_disp_z_coarse = minimum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.75 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.25 rtol = 1.0e-01
    @test avg_stress_fine[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-01
end

@testset "Schwarz Overlap Static Cuboid Hex8 Different Steps" begin
    cp("../examples/overlap/static-different-steps/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
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
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    min_disp_y_fine = minimum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    min_disp_x_coarse = minimum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    min_disp_y_coarse = minimum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    min_disp_z_coarse = minimum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.75 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.25 rtol = 1.0e-01
    @test avg_stress_fine[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-01
end

@testset "Schwarz Overlap L2 Error Static Cuboid Hex8" begin
    cp("../examples/overlap/static-same-step/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)

    sim_default = Norma.create_simulation("cuboids.yaml")
    bc_default_1 = only([
        bc for bc in sim_default.subsims[1].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_default_2 = only([
        bc for bc in sim_default.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test bc_default_1.compute_overlap_l2_error == false
    @test bc_default_2.compute_overlap_l2_error == false
    Exodus.close(sim_default.subsims[1].params["input_mesh"])
    Exodus.close(sim_default.subsims[1].params["output_mesh"])
    Exodus.close(sim_default.subsims[2].params["input_mesh"])
    Exodus.close(sim_default.subsims[2].params["output_mesh"])

    sim_default = Norma.run("cuboids.yaml")
    @test !isfile("overlap-l2-errors-0001.csv")
    rm("cuboid-1.e")
    rm("cuboid-2.e")

    cuboids_text = read("cuboids.yaml", String)
    write("cuboids.yaml", replace(cuboids_text, "CSV output interval: 0" => "CSV output interval: 1"))
    cuboid_1_text = read("cuboid-1.yaml", String)
    write(
        "cuboid-1.yaml",
        replace(
            cuboid_1_text,
            "source side set: ssz-" => "source side set: ssz-\n      compute overlap L2 error: true",
        ),
    )
    cuboid_2_text = read("cuboid-2.yaml", String)
    write(
        "cuboid-2.yaml",
        replace(
            cuboid_2_text,
            "source side set: ssz+" => "source side set: ssz+\n      compute overlap L2 error: true",
        ),
    )

    sim_enabled = Norma.create_simulation("cuboids.yaml")
    bc_enabled_1 = only([
        bc for bc in sim_enabled.subsims[1].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_enabled_2 = only([
        bc for bc in sim_enabled.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test bc_enabled_1.compute_overlap_l2_error == true
    @test bc_enabled_2.compute_overlap_l2_error == true
    @test length(bc_enabled_1.overlap_node_indices) > length(unique(bc_enabled_1.side_set_node_indices))
    @test length(bc_enabled_2.overlap_node_indices) > length(unique(bc_enabled_2.side_set_node_indices))
    @test isnan(bc_enabled_1.overlap_l2_error)
    @test isnan(bc_enabled_2.overlap_l2_error)
    Exodus.close(sim_enabled.subsims[1].params["input_mesh"])
    Exodus.close(sim_enabled.subsims[1].params["output_mesh"])
    Exodus.close(sim_enabled.subsims[2].params["input_mesh"])
    Exodus.close(sim_enabled.subsims[2].params["output_mesh"])

    sim_enabled_ref = Ref{Any}(nothing)
    run_output = mktemp() do path, io
        redirect_stdout(io) do
            sim_enabled_ref[] = Norma.run("cuboids.yaml")
        end
        flush(io)
        seekstart(io)
        read(io, String)
    end
    sim_enabled = sim_enabled_ref[]
    bc_enabled_1 = only([
        bc for bc in sim_enabled.subsims[1].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_enabled_2 = only([
        bc for bc in sim_enabled.subsims[2].model.boundary_conditions if bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test isfinite(bc_enabled_1.overlap_l2_error)
    @test isfinite(bc_enabled_2.overlap_l2_error)
    @test bc_enabled_1.overlap_l2_error >= 0.0
    @test bc_enabled_2.overlap_l2_error >= 0.0
    @test isfile("overlap-l2-errors-0001.csv")
    overlap_csv = read("overlap-l2-errors-0001.csv", String)
    @test occursin("domain,side_set,overlap_l2_error", overlap_csv)
    @test occursin("cuboid-1,ssz+", overlap_csv)
    @test occursin("cuboid-2,ssz-", overlap_csv)
    @test occursin("Overlap L2 error [cuboid-1:ssz+]", run_output)
    @test occursin("Overlap L2 error [cuboid-2:ssz-]", run_output)

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
end
