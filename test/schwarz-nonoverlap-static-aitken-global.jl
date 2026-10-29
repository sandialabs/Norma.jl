# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Global Aitken factor (`aitken scope: global`): one factor for all relaxed
# interfaces, from their stacked residuals, against the default per-interface
# factors.

using YAML

# Run the multi-domain input `main` of example directory `src`, with `files`
# copied in, the controller keys of `overrides` set, and the outputs removed.
function run_aitken_scope(src, main, files, outputs; overrides=Dict{String,Any}())
    for file in vcat([main], files)
        cp(joinpath(src, file), file; force=true)
    end
    params = YAML.load_file(main; dicttype=Norma.Parameters)
    params["name"] = main
    for (k, v) in overrides
        params[k] = v
    end
    sim = Norma.run(params)
    for file in vcat([main], files, outputs)
        rm(file; force=true)
    end
    return sim
end

@testset "Schwarz Aitken Global Factor Formulas" begin
    ds = [[1.0, -2.0], [0.5, 0.25, -1.0]]
    δs = [[-0.5, 1.5], [2.0, -0.75, 0.5]]
    rs = [[0.2, -0.1], [1.0, 0.0, 0.5]]
    xs = [[3.0, 1.0], [1.0, 2.0, -1.0]]
    δ_sq = sum(δ -> dot(δ, δ), δs)

    # Secant form: the factor of the stacked vectors, not a combination of the
    # per-interface factors.
    θ = Norma.aitken_secant_stacked(ds, δs, xs, 0.5)
    @test θ ≈ -(dot(ds[1], δs[1]) + dot(ds[2], δs[2])) / δ_sq
    θ_1 = -dot(ds[1], δs[1]) / dot(δs[1], δs[1])
    θ_2 = -dot(ds[2], δs[2]) / dot(δs[2], δs[2])
    @test !(θ ≈ θ_1) && !(θ ≈ θ_2)
    # An unresolved difference falls back to the input factor.
    @test Norma.aitken_secant_stacked([zeros(2), zeros(3)], δs, xs, 0.5) == 0.5
    @test Norma.aitken_secant_stacked(ds, [zeros(2), zeros(3)], xs, 0.5) == 0.5

    # Recursive form; an unresolved residual difference keeps the previous factor.
    @test Norma.aitken_recursive_stacked(rs, δs, xs, 0.7) ≈ -0.7 * (dot(rs[1], δs[1]) + dot(rs[2], δs[2])) / δ_sq
    @test Norma.aitken_recursive_stacked(rs, [zeros(2), zeros(3)], xs, 0.7) == 0.7
end

@testset "Schwarz Aitken Global Factor Robin-Robin Single Interface" begin
    # One relaxed interface: the global factors are the per-interface ones, so
    # both scopes take the same Schwarz iterations to the same solution.
    # The Robin-Robin inputs use the meshes of the Dirichlet-Neumann example.
    function run_robin(overrides)
        meshes = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann"
        for mesh in ["cuboid-1.g", "cuboid-2.g"]
            cp(joinpath(meshes, mesh), mesh; force=true)
        end
        return run_aitken_scope(
            "../examples/nonoverlap/static-same-step/cuboids-robin-robin",
            "cuboids.yaml",
            ["cuboid-1.yaml", "cuboid-2.yaml"],
            ["cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e"];
            overrides=overrides,
        )
    end
    for relaxation in ["aitken secant", "aitken recursive"]
        per_interface = run_robin(Dict{String,Any}("relaxation" => relaxation, "aitken scope" => "per interface"))
        global_ = run_robin(Dict{String,Any}("relaxation" => relaxation, "aitken scope" => "global"))
        @test per_interface.controller.aitken_scope == :per_interface
        @test global_.controller.aitken_scope == :global
        @test per_interface.failed == false
        @test global_.failed == false
        @test global_.controller.schwarz_iters == per_interface.controller.schwarz_iters
        for (a, b) in zip(global_.subsims, per_interface.subsims)
            @test a.model.displacement ≈ b.model.displacement rtol = 1.0e-10
        end
    end

    # The per-interface secant factor of one interface, replayed on a key with
    # no history, is the stacked factor of that interface alone.
    controller = run_robin(Dict{String,Any}("relaxation" => "aitken secant")).controller
    key = (97, 96, 1)
    x_0, g_0 = [1.0, 2.0, 3.0], [1.5, 1.0, 3.5]
    x_1, g_1 = [1.2, 1.6, 3.2], [1.4, 1.5, 3.3]
    Norma.relaxation_aitken_secant_theta!(controller, key, 1, 1, g_0, x_0)
    θ = Norma.relaxation_aitken_secant_theta!(controller, key, 1, 2, g_1, x_1)
    @test θ == Norma.aitken_secant_stacked([x_1 .- x_0], [(g_1 .- x_1) .- (g_0 .- x_0)], [x_1], 0.0)
end

@testset "Schwarz Aitken Global Factor Dirichlet-Neumann Cuboids" begin
    # One unconstrained interface, whose global factor is formed from the
    # projected interface trace: the same solution as per-interface Aitken.
    src = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann"
    files = ["cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"]
    outputs = ["cuboid-1.e", "cuboid-2.e"]
    for relaxation in ["aitken secant", "aitken recursive"]
        sim = run_aitken_scope(
            src,
            "cuboids-aitken.yaml",
            files,
            outputs;
            overrides=Dict{String,Any}("relaxation" => relaxation, "aitken scope" => "global"),
        )
        @test sim.failed == false
        @test sim.controller.converged == true
        model_fine = sim.subsims[1].model
        @test minimum(model_fine.displacement[1, :]) ≈ -0.125 rtol = 1.0e-06
        @test maximum(model_fine.displacement[3, :]) ≈ 0.5 rtol = 1.0e-06
        @test average_components(model_fine.stress)[3] ≈ 5.0e+08 rtol = 1.0e-06
    end
end

@testset "Schwarz Aitken Global Factor Inclusion Split" begin
    # Four relaxed interfaces on the Dirichlet side (the inclusion split of
    # issue #233): every interface takes the same factor at every iteration.
    src = "../examples/nonoverlap/static-same-step/matrix-inclusion-split"
    files = ["disk.yaml", "matrix.yaml", "disk.g", "matrix.g"]
    outputs = ["disk.e", "matrix.e"]
    sim = run_aitken_scope(
        src,
        "inclusion.yaml",
        files,
        outputs;
        overrides=Dict{String,Any}(
            "relaxation" => "aitken secant", "aitken scope" => "global", "maximum iterations" => 256
        ),
    )
    @test sim.failed == false
    @test sim.controller.converged == true
    state = Norma.AITKEN_GLOBAL_STATE[sim.controller]
    @test length(state.theta) == 4
    @test length(unique(values(state.theta))) == 1
end

# Where the global factor does not apply, the run aborts instead of silently
# taking other factors.
@testset "Schwarz Aitken Global Factor Rejected" begin
    function expect_abort(f)
        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException f()
        finally
            Norma.NORMA_TEST_MODE[] = false
        end
    end
    function with_inputs(f, src, main, files, outputs, overrides)
        for file in vcat([main], files)
            cp(joinpath(src, file), file; force=true)
        end
        params = YAML.load_file(main; dicttype=Norma.Parameters)
        params["name"] = main
        merge!(params, overrides)
        try
            f(params)
        finally
            for file in vcat([main], files, outputs)
                rm(file; force=true)
            end
        end
    end

    # Controller keys: an unknown scope, and the global scope without an
    # Aitken relaxation.
    robin = "../examples/nonoverlap/static-same-step/cuboids-robin-robin"
    robin_files = ["cuboid-1.yaml", "cuboid-2.yaml"]
    for overrides in (
        Dict{String,Any}("relaxation" => "aitken secant", "aitken scope" => "globl"),
        Dict{String,Any}("aitken scope" => "global"),
        Dict{String,Any}("relaxation" => "anderson", "aitken scope" => "global"),
    )
        with_inputs(robin, "cuboids.yaml", robin_files, String[], overrides) do params
            expect_abort(() -> Norma.create_simulation(params))
        end
    end

    # Schwarz overlap matching: rejected at setup.
    matching = "../examples/overlap/dynamic-same-step/cantilever-matching"
    matching_files = ["cantilever-clamped.yaml", "cantilever-free.yaml", "cantilever-clamped.g", "cantilever-free.g"]
    with_inputs(
        matching,
        "cantilever.yaml",
        matching_files,
        ["cantilever-clamped.e", "cantilever-free.e"],
        Dict{String,Any}("relaxation" => "aitken secant", "aitken scope" => "global"),
    ) do params
        expect_abort(() -> Norma.create_simulation(params))
    end

    # A windowed stop (5 substeps of the relaxed side per stop): rejected when
    # it occurs, at the second substep of the first stop.
    dn = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"
    dn_files = ["cantilever-clamped.yaml", "cantilever-free.yaml", "cantilever-clamped.g", "cantilever-free.g"]
    with_inputs(
        dn,
        "cantilever-multi.yaml",
        dn_files,
        ["cantilever-clamped.e", "cantilever-free.e"],
        Dict{String,Any}(
            "relaxation" => "aitken secant",
            "aitken scope" => "global",
            "final time" => 5.0e-6,
            "time step" => 2.5e-6,
            "Exodus output interval" => 0,
            "CSV output interval" => 0,
        ),
    ) do params
        expect_abort(() -> Norma.run(params))
    end
end
