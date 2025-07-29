# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra: norm
using StaticArrays

@testset "Elastic Constants" begin
    @testset "Given E And Ν" begin
        params = Norma.Parameters("elastic modulus" => 210e9, "Poisson's ratio" => 0.3)
        E, ν, κ, λ, μ = Norma.elastic_constants(params)
        @test isapprox(E, 210e9)
        @test isapprox(ν, 0.3)
        @test isapprox(κ, E / 3(1 - 2ν))
        @test isapprox(λ, E * ν / ((1 + ν) * (1 - 2ν)))
        @test isapprox(μ, E / 2(1 + ν))
    end

    @testset "Given E And Κ" begin
        E = 210e9
        κ = 175e9
        params = Norma.Parameters("elastic modulus" => E, "bulk modulus" => κ)
        E2, ν, κ2, λ, μ = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(κ2, κ)
        @test isapprox(ν, (3κ - E) / (6κ))
        @test isapprox(λ, (3κ * (3κ - E)) / (9κ - E))
        @test isapprox(μ, 3κ * E / (9κ - E))
    end

    @testset "Given E And Λ" begin
        E = 210e9
        λ = 121e9
        params = Norma.Parameters("elastic modulus" => E, "Lamé's first constant" => λ)
        E2, ν, κ, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(λ2, λ)
        @test isapprox(μ, (E - 3λ + sqrt(E^2 + 9λ^2 + 2E * λ)) / 4)
    end

    @testset "Given E And Μ" begin
        E = 210e9
        μ = 80.77e9
        params = Norma.Parameters("elastic modulus" => E, "shear modulus" => μ)
        E2, ν, κ, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(μ2, μ)
        @test isapprox(ν, E / (2μ) - 1)
    end

    @testset "Given Ν And Κ" begin
        ν = 0.3
        κ = 150e9
        params = Norma.Parameters("Poisson's ratio" => ν, "bulk modulus" => κ)
        E, ν2, κ2, λ, μ = Norma.elastic_constants(params)
        @test isapprox(ν2, ν)
        @test isapprox(κ2, κ)
        @test isapprox(E, 3κ * (1 - 2ν))
    end

    @testset "Given Ν And Λ" begin
        ν = 0.3
        λ = 121e9
        params = Norma.Parameters("Poisson's ratio" => ν, "Lamé's first constant" => λ)
        E, ν2, κ, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(ν2, ν)
        @test isapprox(λ2, λ)
    end

    @testset "Given Ν And Μ" begin
        ν = 0.3
        μ = 80e9
        params = Norma.Parameters("Poisson's ratio" => ν, "shear modulus" => μ)
        E, ν2, κ, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(μ2, μ)
        @test isapprox(ν2, ν)
    end

    @testset "Given Κ And Λ" begin
        κ = 150e9
        λ = 90e9
        params = Norma.Parameters("bulk modulus" => κ, "Lamé's first constant" => λ)
        E, ν, κ2, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(κ2, κ)
        @test isapprox(λ2, λ)
    end

    @testset "Given Κ And Μ" begin
        κ = 150e9
        μ = 50e9
        params = Norma.Parameters("bulk modulus" => κ, "shear modulus" => μ)
        E, ν, κ2, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(κ2, κ)
        @test isapprox(μ2, μ)
    end

    @testset "Given Λ And Μ" begin
        λ = 100e9
        μ = 80e9
        params = Norma.Parameters("Lamé's first constant" => λ, "shear modulus" => μ)
        E, ν, κ, λ2, μ2 = Norma.elastic_constants(params)
        @test isapprox(λ2, λ)
        @test isapprox(μ2, μ)
    end

    @testset "Failure Cases" begin
        Norma.norma_log(0, :info, "Testing bad input. Mock code abort ...")
        Norma.NORMA_TEST_MODE[] = true
        @test_throws Norma.NormaAbortException Norma.elastic_constants(Norma.Parameters("elastic modulus" => 200e9))
        @test_throws Norma.NormaAbortException Norma.elastic_constants(Norma.Parameters("Poisson's ratio" => 0.25))
        @test_throws Norma.NormaAbortException Norma.elastic_constants(Norma.Parameters("bulk modulus" => 180e9))
        @test_throws Norma.NormaAbortException Norma.elastic_constants(
            Norma.Parameters("Lamé's first constant" => 90e9)
        )
        @test_throws Norma.NormaAbortException Norma.elastic_constants(Norma.Parameters("shear modulus" => 80e9))
        @test_throws Norma.NormaAbortException Norma.elastic_constants(Norma.Parameters())
        Norma.NORMA_TEST_MODE[] = false
    end
end

const I3 = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

@testset "Constitutive() For All Solids At Identity Deformation" begin
    F = I3

    @testset "Saint Venant Kirchhoff" begin
        params = Norma.Parameters("elastic modulus" => 100.0, "Poisson's ratio" => 0.3, "density" => 7800.0)
        mat = Norma.SaintVenant_Kirchhoff(params)
        W, P, AA = Norma.constitutive(mat, F)
        @test isapprox(W, 0.0; atol=1e-12)
        @test size(P) == (3, 3)
        @test size(AA) == (3, 3, 3, 3)
    end

    @testset "Linear Elastic" begin
        params = Norma.Parameters("elastic modulus" => 100.0, "Poisson's ratio" => 0.3, "density" => 7800.0)
        mat = Norma.Linear_Elastic(params)
        W, σ, CC = Norma.constitutive(mat, F)
        @test isapprox(W, 0.0; atol=1e-12)
        @test size(σ) == (3, 3)
        @test size(CC) == (3, 3, 3, 3)
    end

    @testset "Neohookean" begin
        params = Norma.Parameters("elastic modulus" => 100.0, "Poisson's ratio" => 0.3, "density" => 7800.0)
        mat = Norma.Neohookean(params)
        W, P, AA = Norma.constitutive(mat, F)
        @test isapprox(W, 0.0; atol=1e-12)
        @test size(P) == (3, 3)
        @test size(AA) == (3, 3, 3, 3)
    end

    @testset "Seth Hill" begin
        params = Norma.Parameters(
            "elastic modulus" => 100.0, "Poisson's ratio" => 0.3, "density" => 7800.0, "m" => 1, "n" => 1
        )
        mat = Norma.SethHill(params)
        W, P, AA = Norma.constitutive(mat, F)
        @test isapprox(W, 0.0; atol=1e-12)
        @test size(P) == (3, 3)
        @test size(AA) == (3, 3, 3, 3)
    end
end

@testset "J2 Stress" begin
    E = 200.0e+09
    nu = 0.25
    rho = 7800.0
    Y = 1.0e+09
    mu = E / 2 / (1 + nu)
    params = Norma.Parameters()
    params["elastic modulus"] = E
    params["Poisson's ratio"] = nu
    params["density"] = rho
    params["yield stress"] = Y
    material = Norma.J2(params)
    gamma = Y / mu / sqrt(3)
    F = [1.0 gamma 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Fp = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    eqps = 0.0
    dt = 1.0e-06
    Fe, Fp, eqps, sigma = Norma.stress_update(material, F, Fp, eqps, dt)
    @test norm(Fp - I(3)) ≈ 0.0 atol = 1.0e-12
    @test eqps ≈ 0.0 atol = 1.0e-12
    @test sqrt(1.5) * norm(Norma.dev(sigma)) / Y ≈ 1.0 rtol = 3.0e-06
    gamma = 2 * Y / mu / sqrt(3)
    F = [1.0 gamma 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Fp = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    eqps = 0.0
    dt = 1.0e-06
    Fe, Fp, eqps, sigma = Norma.stress_update(material, F, Fp, eqps, dt)
    @test sqrt(1.5) * norm(Norma.dev(sigma)) / Y ≈ 1.0 rtol = 2.0e-16
end
