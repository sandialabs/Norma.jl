# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using YAML

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
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
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
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
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

# A failed step is retried from the state it started from. Every Schwarz
# iteration restarts the subdomains at the start of the stop, so a failed first
# substep of an iteration must be retried from there too, with the internal
# variables of that state, and not from the end of the stop the previous
# iteration reached. The cuboids of Different Steps in J2 plasticity, with
# adaptive steps on cuboid-1, whose first substep fails in Schwarz iteration 2.
@testset "Schwarz Overlap Static Cuboid Hex8 Retry From Start Of Stop" begin
    folder = "../examples/overlap/static-different-steps/cuboids"
    for file in ("cuboid-1.g", "cuboid-2.g")
        cp(joinpath(folder, file), file; force=true)
    end
    plastic = Dict{String,Any}(
        "model" => "j2 plasticity",
        "elastic modulus" => 1.0e9,
        "Poisson's ratio" => 0.25,
        "density" => 1000.0,
        "yield stress" => 2.0e6,
        "hardening modulus" => 1.0e8,
    )
    for file in ("cuboid-1.yaml", "cuboid-2.yaml")
        params = YAML.load_file(joinpath(folder, file); dicttype=Norma.Parameters)
        params["model"]["material"]["elastic"] = plastic
        if file == "cuboid-1.yaml"
            merge!(
                params["time integrator"],
                Dict(
                    "minimum time step" => 0.0625,
                    "maximum time step" => 0.5,
                    "decrease factor" => 0.5,
                    "increase factor" => 1.0,
                ),
            )
        else
            params["boundary conditions"]["Dirichlet"][3]["function"] = "0.02 * t"
        end
        YAML.write_file(file, params)
    end
    params = YAML.load_file(joinpath(folder, "cuboids.yaml"); dicttype=Norma.Parameters)
    merge!(params, Dict("relative tolerance" => 1.0e-6, "absolute tolerance" => 1.0e-7, "maximum iterations" => 60))
    params["name"] = "cuboids.yaml"
    sim = Norma.create_simulation(params)
    fine = sim.subsims[1]
    controller = sim.controller
    failed = false
    retried = false
    offset = NaN
    same_state = false
    Norma.SOLVE_ITERATION_OBSERVER[] = function (model, iteration)
        (iteration == 0 && model === fine.model) || return nothing
        if failed && !retried
            # The first iterate of the retry: the free degrees of freedom keep
            # the displacement the step starts from.
            retried = true
            u = fine.integrator.displacement - controller.stop_disp[1]
            offset = maximum(abs, u[model.free_dofs])
            same_state = model.state_old == model.stop_state_old
        elseif !failed && controller.stop == 1 && controller.iteration_number == 2
            failed = true
            model.failed = true
        end
        return nothing
    end
    try
        Norma.run(sim)
    finally
        Norma.SOLVE_ITERATION_OBSERVER[] = nothing
        for file in ("cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e")
            rm(file; force=true)
        end
    end
    @test retried
    @test offset == 0.0
    @test same_state
    @test fine.integrator.time ≈ 1.0
    @test maximum(fine.model.displacement[3, :]) ≈ 0.015 rtol = 1.0e-3
end

@testset "Schwarz Overlap Static Cuboid Hex8 Weak DBC" begin
    cp("../examples/overlap/static-same-step/cuboids-weak/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids-weak/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids-weak/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids-weak/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids-weak/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
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

@testset "Schwarz Overlap Static Cuboid Hex8 Different Steps Weak DBC" begin
    cp("../examples/overlap/static-different-steps/cuboids-weak/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids-weak/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids-weak/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids-weak/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-different-steps/cuboids-weak/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
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
        bc for bc in sim_default.subsims[1].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_default_2 = only([
        bc for bc in sim_default.subsims[2].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test bc_default_1.compute_overlap_l2_error == ""
    @test bc_default_2.compute_overlap_l2_error == ""
    Exodus.close(sim_default.subsims[1].params["input_mesh"])
    Exodus.close(sim_default.subsims[1].params["output_mesh"])
    Exodus.close(sim_default.subsims[2].params["input_mesh"])
    Exodus.close(sim_default.subsims[2].params["output_mesh"])

    sim_default = Norma.run("cuboids.yaml")
    @test !isfile("overlap-l2-disp-rel-errors-0001.csv")
    rm("cuboid-1.e")
    rm("cuboid-2.e")

    cuboids_text = read("cuboids.yaml", String)
    write("cuboids.yaml", replace(cuboids_text, "CSV output interval: 0" => "CSV output interval: 1.0"))
    cuboid_1_text = read("cuboid-1.yaml", String)
    write(
        "cuboid-1.yaml",
        replace(
            cuboid_1_text,
            "source block: coarse" => "source block: coarse\n      compute overlap L2 relative error: disp",
        ),
    )
    cuboid_2_text = read("cuboid-2.yaml", String)
    write(
        "cuboid-2.yaml",
        replace(
            cuboid_2_text, "source block: fine" => "source block: fine\n      compute overlap L2 relative error: disp"
        ),
    )

    sim_enabled = Norma.create_simulation("cuboids.yaml")
    bc_enabled_1 = only([
        bc for bc in sim_enabled.subsims[1].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_enabled_2 = only([
        bc for bc in sim_enabled.subsims[2].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test bc_enabled_1.compute_overlap_l2_error == "disp"
    @test bc_enabled_2.compute_overlap_l2_error == "disp"
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
        bc for bc in sim_enabled.subsims[1].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    bc_enabled_2 = only([
        bc for bc in sim_enabled.subsims[2].model.boundary_conditions if
        bc isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition
    ])
    @test isfinite(bc_enabled_1.overlap_l2_error)
    @test isfinite(bc_enabled_2.overlap_l2_error)
    @test bc_enabled_1.overlap_l2_error >= 0.0
    @test bc_enabled_2.overlap_l2_error >= 0.0
    @test isfile("overlap-l2-disp-rel-errors-0001.csv")
    overlap_csv = read("overlap-l2-disp-rel-errors-0001.csv", String)
    @test occursin("domain,side_set,overlap_l2_error", overlap_csv)
    @test occursin("cuboid-1,ssz+", overlap_csv)
    @test occursin("cuboid-2,ssz-", overlap_csv)
    @test occursin("Overlap L2 displacement error [cuboid-1:ssz+]", run_output)
    @test occursin("Overlap L2 displacement error [cuboid-2:ssz-]", run_output)

    rm("cuboids.yaml")
    rm("cuboid-1.yaml")
    rm("cuboid-2.yaml")
    rm("cuboid-1.g")
    rm("cuboid-2.g")
    rm("cuboid-1.e")
    rm("cuboid-2.e")
    for file in readdir()
        if startswith(file, "iterations-") || startswith(file, "overlap-l2-")
            rm(file)
        end
    end
end

# The classical overlapping condition is iterated without relaxation, so
# relaxation requested for a simulation whose couplings are all classical
# overlapping conditions is rejected at setup; a relaxation parameter of one
# is the default and is accepted.
@testset "Schwarz Overlap Relaxation Rejected" begin
    src = "../examples/overlap/static-same-step/cuboids"
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"]
        cp("$src/$f", f; force=true)
    end
    function overlap_params(overrides)
        params = YAML.load_file("cuboids.yaml"; dicttype=Norma.Parameters)
        params["name"] = "cuboids.yaml"
        merge!(params, overrides)
        return params
    end
    Norma.NORMA_TEST_MODE[] = true
    try
        @test_throws Norma.NormaAbortException Norma.create_simulation(overlap_params(Dict("relaxation" => "anderson")))
        @test_throws Norma.NormaAbortException Norma.create_simulation(
            overlap_params(Dict("relaxation parameter" => 0.5))
        )
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
    sim = Norma.create_simulation(overlap_params(Dict("relaxation parameter" => 1.0)))
    @test sim isa Norma.MultiDomainSimulation
    for subsim in sim.subsims
        Norma.finalize_writing(subsim)
    end
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
end
