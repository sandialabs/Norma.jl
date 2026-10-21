# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using YAML

# The quasi-static `predictor` sets the iterate Newton starts a step from. It
# must leave the converged solution unchanged and start Newton closer to it.

# The linear elastic cube of examples/single/static-solid/cube/standard,
# pulled on rollers to u_z = t in steps of 0.1: its solution is the affine
# field of uniaxial stress, linear in time.
function predictor_cube(predictor; solver=nothing)
    cp("../examples/single/static-solid/cube/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    params = YAML.load_file("cube.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cube.yaml"
    predictor === nothing || (params["time integrator"]["predictor"] = predictor)
    solver === nothing || (params["solver"] = solver)
    try
        return Norma.run(params)
    finally
        rm("cube.yaml"; force=true)
        rm("cube.g"; force=true)
        rm("cube.e"; force=true)
    end
end

# The run of `simulate()` and the residual of the iterate every solve starts
# from, recorded by the solve observer.
function predicted_residuals(simulate)
    residuals = Float64[]
    Norma.SOLVE_ITERATION_OBSERVER[] = function (model, iteration)
        iteration == 0 || return nothing
        residual = model.internal_force - model.boundary_force
        push!(residuals, norm(residual[model.free_dofs]))
        return nothing
    end
    try
        return simulate(), residuals
    finally
        Norma.SOLVE_ITERATION_OBSERVER[] = nothing
    end
end

@testset "Quasi-Static Predictor Linear Elastic Cube" begin
    constant, constant_residuals = predicted_residuals(() -> predictor_cube(nothing))
    @test constant.integrator.predictor == :constant
    @test length(constant_residuals) == 10
    for predictor in ("linear", "tangent")
        sim, residuals = predicted_residuals(() -> predictor_cube(predictor))
        @test sim.integrator.predictor == Symbol(predictor)
        @test sim.model.displacement ≈ constant.model.displacement rtol = 1.0e-10
        # The constant predictor leaves the increment of the loaded face on
        # the top layer of elements. The other two predict the affine
        # solution, `tangent` from the elastic response and `linear` from the
        # increment of the previous step, which the first step does not have.
        predicted = predictor == "linear" ? (2:10) : (1:10)
        @test all(residuals[predicted] .< 1.0e-12 .* constant_residuals[predicted])
        if predictor == "linear"
            @test residuals[1] ≈ constant_residuals[1] rtol = 1.0e-12
        end
        # The relative tolerance, measured from the residual of the iterate
        # without the predictor, is met as with the constant predictor: in one
        # iteration per step.
        @test Norma.NEWTON_ITERATION_COUNT[] == 10
    end
end

@testset "Quasi-Static Predictor Input Errors" begin
    Norma.norma_log(0, :info, "Testing bad input. Mock code abort ...")
    Norma.NORMA_TEST_MODE[] = true
    try
        @test_throws Norma.NormaAbortException predictor_cube("quadratic")
        steepest_descent = Norma.Parameters(
            "type" => "steepest descent",
            "step" => "steepest descent",
            "step length" => 0.1,
            "minimum iterations" => 1,
            "maximum iterations" => 64,
            "relative tolerance" => 1.0e-15,
            "absolute tolerance" => 1.0e-08,
        )
        @test_throws Norma.NormaAbortException predictor_cube("tangent"; solver=steepest_descent)
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
end

# A matrix-free solver has no stiffness to estimate the residual of the
# iterate without the prediction, so it evaluates that residual, and the
# relative tolerance is measured from it as with the Hessian minimizer. Measured
# from the residual of the predicted iterate, which `linear` makes nearly the
# solution of the cube, the tolerance was below the floor that L-BFGS reaches,
# and every solve after the first stopped at `maximum iterations`.
@testset "Quasi-Static Predictor Matrix-Free Solver" begin
    lbfgs = Norma.Parameters(
        "type" => "steepest descent",
        "step" => "lbfgs",
        "memory" => 10,
        "step length" => 1.0,
        "minimum iterations" => 1,
        "maximum iterations" => 300,
        "relative tolerance" => 1.0e-8,
        "absolute tolerance" => 1.0e-10,
    )
    constant = predictor_cube("constant"; solver=lbfgs)
    @test Norma.UNCONVERGED_SOLVE_COUNT[] == 0
    constant_iterations = Norma.NEWTON_ITERATION_COUNT[]
    linear = predictor_cube("linear"; solver=lbfgs)
    @test Norma.UNCONVERGED_SOLVE_COUNT[] == 0
    @test Norma.NEWTON_ITERATION_COUNT[] < constant_iterations
    @test linear.model.displacement ≈ constant.model.displacement rtol = 1.0e-6
end

# LCM's RubiksCube with one element per grain (rubiks-cube-1.g), from the
# example input. At fixed time steps, a step that fails cannot be cut, and the
# run aborts; with `adaptive`, it is halved as in the example.
function predictor_rubiks_cube(predictor, time_step, final_time; adaptive=false)
    folder = "../examples/crystal-plasticity/rubiks-cube"
    cp(joinpath(folder, "rubiks-cube-1.g"), "rubiks-cube-1.g"; force=true)
    params = YAML.load_file(joinpath(folder, "rubiks-cube.yaml"); dicttype=Norma.Parameters)
    params["name"] = "rubiks-cube-predictor.yaml"
    params["input mesh file"] = "rubiks-cube-1.g"
    params["output mesh file"] = "rubiks-cube-predictor.e"
    params["Exodus output interval"] = 0.0
    integrator = params["time integrator"]
    if adaptive
        integrator["maximum time step"] = time_step
    else
        for key in ("minimum time step", "maximum time step", "decrease factor", "increase factor")
            delete!(integrator, key)
        end
    end
    integrator["time step"] = time_step
    integrator["final time"] = final_time
    integrator["predictor"] = predictor
    try
        sim = Norma.run(params)
        return sim, Norma.NEWTON_ITERATION_COUNT[]
    finally
        rm("rubiks-cube-1.g"; force=true)
        rm("rubiks-cube-predictor.e"; force=true)
    end
end

@testset "Quasi-Static Predictor Crystal Plasticity" begin
    # Steps of 2.5e-4 through the onset of flow, which the constant predictor
    # also takes: every predictor converges to the same solution.
    constant, constant_iterations = predictor_rubiks_cube("constant", 2.5e-4, 2.0e-3)
    eqps(sim) = [state[end] for block in sim.model.state_old for element in block for state in element]
    @test maximum(eqps(constant)) > 1.0e-4
    for predictor in ("linear", "tangent")
        sim, iterations = predictor_rubiks_cube(predictor, 2.5e-4, 2.0e-3)
        @test sim.integrator.time ≈ 2.0e-3
        @test sim.model.displacement ≈ constant.model.displacement rtol = 1.0e-8
        @test eqps(sim) ≈ eqps(constant) rtol = 1.0e-6
        @test iterations < constant_iterations
    end
    # LCM's steps of 5e-4: from the second step on, the constant predictor
    # strains the elements next to the loaded face beyond the range of the
    # flow rule, and the run aborts.
    Norma.norma_log(0, :info, "Testing a failed step at a fixed time step. Mock code abort ...")
    Norma.NORMA_TEST_MODE[] = true
    try
        @test_throws Norma.NormaAbortException predictor_rubiks_cube("constant", 5.0e-4, 1.0e-3)
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
    for predictor in ("linear", "tangent")
        sim, _ = predictor_rubiks_cube(predictor, 5.0e-4, 1.0e-3)
        @test sim.integrator.time ≈ 1.0e-3
    end
    # With adaptive time stepping, both predictors cut the step at the onset
    # of flow, [1e-3, 1.5e-3], and retry it with half the step from its start
    # and the history of the last accepted step: they take the same steps and
    # reach the same solution.
    adaptive = Dict{String,Any}()
    for predictor in ("linear", "tangent")
        times = Float64[]
        Norma.SOLVE_ITERATION_OBSERVER[] = function (model, iteration)
            iteration == 0 && push!(times, model.time)
            return nothing
        end
        sim = try
            first(predictor_rubiks_cube(predictor, 5.0e-4, 2.0e-3; adaptive=true))
        finally
            Norma.SOLVE_ITERATION_OBSERVER[] = nothing
        end
        @test sim.integrator.time ≈ 2.0e-3
        @test any(t -> isapprox(t, 1.25e-3; rtol=1.0e-9), times)
        adaptive[predictor] = sim
    end
    @test adaptive["linear"].model.displacement ≈ adaptive["tangent"].model.displacement rtol = 1.0e-8
    @test eqps(adaptive["linear"]) ≈ eqps(adaptive["tangent"]) rtol = 1.0e-6
end

# The overlapping Schwarz cuboids of examples/overlap/static-different-steps,
# whose subdomains take 2 and 3 steps per stop. Every Schwarz iteration
# restarts the subdomains from the start of the stop, which the predictors
# must start from too.
function predictor_cuboids(predictor)
    folder = "../examples/overlap/static-different-steps/cuboids"
    for file in ("cuboids.yaml", "cuboid-1.g", "cuboid-2.g")
        cp(joinpath(folder, file), file; force=true)
    end
    for file in ("cuboid-1.yaml", "cuboid-2.yaml")
        text = read(joinpath(folder, file), String)
        write(file, replace(text, "type: quasi static" => "type: quasi static\n  predictor: $predictor"))
    end
    try
        return Norma.run("cuboids.yaml")
    finally
        for file in ("cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g")
            rm(file; force=true)
        end
        rm("cuboid-1.e"; force=true)
        rm("cuboid-2.e"; force=true)
    end
end

@testset "Quasi-Static Predictor Schwarz" begin
    constant = predictor_cuboids("constant")
    for predictor in ("linear", "tangent")
        sim = predictor_cuboids(predictor)
        @test all(subsim.integrator.predictor == Symbol(predictor) for subsim in sim.subsims)
        for (subsim, reference) in zip(sim.subsims, constant.subsims)
            @test subsim.model.displacement ≈ reference.model.displacement rtol = 1.0e-6
        end
        fine, coarse = sim.subsims[1].model, sim.subsims[2].model
        @test maximum(fine.displacement[3, :]) ≈ 0.75 rtol = 1.0e-6
        @test minimum(coarse.displacement[1, :]) ≈ -0.125 rtol = 1.0e-6
        @test average_components(fine.stress)[3] ≈ 5.0e+08 rtol = 1.0e-6
    end
end

# An in-memory restart switches a subdomain between variants built once. A
# variant switched back in starts its predictor history again, as after a
# swap, instead of predicting from its last step before it was switched out.
# The cuboids with `linear`, in four stops: the fine cuboid switches to a copy
# of itself at the second stop and back at the fourth.
@testset "Quasi-Static Predictor In-Memory Restart" begin
    folder = "../examples/overlap/static-different-steps/cuboids"
    for file in ("cuboid-1.g", "cuboid-2.g")
        cp(joinpath(folder, file), file; force=true)
    end
    write(
        "cuboids.yaml", replace(read(joinpath(folder, "cuboids.yaml"), String), "time step: 1.0" => "time step: 0.25")
    )
    for file in ("cuboid-1.yaml", "cuboid-2.yaml")
        text = read(joinpath(folder, file), String)
        write(file, replace(text, "type: quasi static" => "type: quasi static\n  predictor: linear"))
    end
    write("cuboid-1-copy.yaml", replace(read("cuboid-1.yaml", String), "cuboid-1.e" => "cuboid-1-copy.e"))
    reference = Norma.run("cuboids.yaml")
    r = Norma.InMemoryRestart("cuboids.yaml")
    history = Float64[]
    restarted = try
        Norma.register_variant!(r, 1, "copy", "cuboid-1-copy.yaml")
        Norma.march!(r; on_step=function (r, step)
            if step == 2
                Norma.switch!(r, 1, "copy")
            elseif step == 4
                Norma.switch!(r, 1, "cuboid-1")
                push!(history, r.sim.subsims[1].integrator.last_time_step)
            end
            return nothing
        end)
        Norma.displacement_vector(r)
    finally
        Norma.close_restart!(r)
        for file in ("cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1-copy.yaml")
            rm(file; force=true)
        end
        for file in ("cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e", "cuboid-1-copy.e")
            rm(file; force=true)
        end
    end
    @test length(r.history) == 2
    @test history == [0.0]
    @test restarted ≈ vcat((vec(subsim.model.displacement) for subsim in reference.subsims)...) rtol = 1.0e-8
end
