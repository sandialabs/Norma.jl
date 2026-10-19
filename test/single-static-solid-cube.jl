# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using YAML

@testset "Single Static Solid Cube" begin
    cp("../examples/single/static-solid/cube/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

@testset "Single Static Solid Cube Direct Solver" begin
    cp("../examples/single/static-solid/cube/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    params = YAML.load_file("cube.yaml"; dicttype=Norma.Parameters)
    params["solver"]["linear solver"] = "direct"
    params["name"] = "cube.yaml"
    simulation = Norma.run(params)
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
end

# The direct linear solver keeps sparse Cholesky for symmetric positive
# definite matrices, and switches to LU for a nonsymmetric or an indefinite
# matrix. A system with Inf or NaN, which a material can leave at a trial
# iterate, gives a NaN step, which fails the next evaluation, instead of the
# exception of LU that would end the run.
@testset "Direct Linear Solver Factorizations" begin
    n = 6
    T = SymTridiagonal(fill(4.0, n), fill(-1.0, n - 1))
    b = collect(1.0:n)
    cases = (
        ("symmetric positive definite", Norma.sparse(Matrix(T)), true),
        ("nonsymmetric", Norma.sparse(Matrix(T) + diagm(1 => fill(0.5, n - 1))), false),
        ("symmetric indefinite", Norma.sparse(Matrix(T) - 5.0 * I), false),
    )
    for (name, A, cholesky) in cases
        cache = Norma.DirectSolverCache(nothing, Int64[], Int64[])
        x = Norma.solve_linear_direct(cache, A, b)
        @test x ≈ Matrix(A) \ b rtol = 1.0e-12
        # Only a Cholesky factor is kept for reuse.
        @test (cache.factor !== nothing) == cholesky
    end
    for value in (NaN, Inf)
        A = Norma.sparse(Matrix(T))
        A[2, 3] = value
        cache = Norma.DirectSolverCache(nothing, Int64[], Int64[])
        Norma.solve_linear_direct(cache, Norma.sparse(Matrix(T)), b)
        @test all(isnan, Norma.solve_linear_direct(cache, A, b))
        @test cache.factor === nothing
        @test all(isnan, Norma.solve_linear_direct(cache, Norma.sparse(Matrix(T)), [NaN; b[2:end]]))
    end
end

@testset "Single Static Solid Cube Line Search" begin
    cp("../examples/single/static-solid/cube/standard/cube-line-search.yaml", "cube-line-search.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-line-search.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-line-search.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

@testset "Single Static Solid Cube Steepest Descent" begin
    cp(
        "../examples/single/static-solid/cube/standard/cube-steepest-descent.yaml",
        "cube-steepest-descent.yaml";
        force=true,
    )
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-steepest-descent.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-steepest-descent.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 3.0
    @test avg_stress[2] ≈ 0.0 atol = 3.0
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 3.0
    @test avg_stress[5] ≈ 0.0 atol = 3.0
    @test avg_stress[6] ≈ 0.0 atol = 3.0
end
