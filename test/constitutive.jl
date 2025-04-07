# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra: norm

const Parameters = Dict{String, Any}

@testset "elastic_constants tests" begin

    @testset "Given E and ν" begin
        params = Parameters("elastic modulus" => 210e9, "Poisson's ratio" => 0.3)
        E, ν, κ, λ, μ = Norma.elastic_constants(params)
        @test isapprox(E, 210e9)
        @test isapprox(ν, 0.3)
        @test isapprox(κ, E / 3(1 - 2ν))
        @test isapprox(λ, E * ν / ((1 + ν) * (1 - 2ν)))
        @test isapprox(μ, E / 2(1 + ν))
    end

    @testset "Given E and κ" begin
        E = 210e9
        κ = 175e9
        params = Parameters("elastic modulus" => E, "bulk modulus" => κ)
        E2, ν, κ2, λ, μ = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(κ2, κ)
        @test isapprox(ν, (3κ - E) / (6κ))
        @test isapprox(λ, (3κ * (3κ - E)) / (9κ - E))
        @test isapprox(μ, 3κ * E / (9κ - E))
    end

    @testset "Given E and λ" begin
        E = 210e9
        λ = 121e9
        params = Parameters("elastic modulus" => E, "Lamé's first constant" => λ)
        E2, ν, κ, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(λ2, λ)
        @test isapprox(μ, (E - 3λ + sqrt(E^2 + 9λ^2 + 2E * λ)) / 4)
    end

    @testset "Given E and μ" begin
        E = 210e9
        μ = 80.77e9
        params = Parameters("elastic modulus" => E, "shear modulus" => μ)
        E2, ν, κ, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(E2, E)
        @test isapprox(μ2, μ)
        @test isapprox(ν, E / (2μ) - 1)
    end

    @testset "Given ν and κ" begin
        ν = 0.3
        κ = 150e9
        params = Parameters("Poisson's ratio" => ν, "bulk modulus" => κ)
        E, ν2, κ2, λ, μ = Norma.elastic_constants(params)
        @test isapprox(ν2, ν)
        @test isapprox(κ2, κ)
        @test isapprox(E, 3κ * (1 - 2ν))
    end

    @testset "Given ν and λ" begin
        ν = 0.3
        λ = 121e9
        params = Parameters("Poisson's ratio" => ν, "Lamé's first constant" => λ)
        E, ν2, κ, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(ν2, ν)
        @test isapprox(λ2, λ)
    end

    @testset "Given ν and μ" begin
        ν = 0.3
        μ = 80e9
        params = Parameters("Poisson's ratio" => ν, "shear modulus" => μ)
        E, ν2, κ, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(μ2, μ)
        @test isapprox(ν2, ν)
    end

    @testset "Given κ and λ" begin
        κ = 150e9
        λ = 90e9
        params = Parameters("bulk modulus" => κ, "Lamé's first constant" => λ)
        E, ν, κ2, λ2, μ = Norma.elastic_constants(params)
        @test isapprox(κ2, κ)
        @test isapprox(λ2, λ)
    end

    @testset "Given κ and μ" begin
        κ = 150e9
        μ = 50e9
        params = Parameters("bulk modulus" => κ, "shear modulus" => μ)
        E, ν, κ2, λ, μ2 = Norma.elastic_constants(params)
        @test isapprox(κ2, κ)
        @test isapprox(μ2, μ)
    end

    @testset "Given λ and μ" begin
        λ = 100e9
        μ = 80e9
        params = Parameters("Lamé's first constant" => λ, "shear modulus" => μ)
        E, ν, κ, λ2, μ2 = Norma.elastic_constants(params)
        @test isapprox(λ2, λ)
        @test isapprox(μ2, μ)
    end

    @testset "Failure cases" begin
        @test_throws ErrorException Norma.elastic_constants(Parameters("elastic modulus" => 200e9))
        @test_throws ErrorException Norma.elastic_constants(Parameters("Poisson's ratio" => 0.25))
        @test_throws ErrorException Norma.elastic_constants(Parameters("bulk modulus" => 180e9))
        @test_throws ErrorException Norma.elastic_constants(Parameters("Lamé's first constant" => 90e9))
        @test_throws ErrorException Norma.elastic_constants(Parameters("shear modulus" => 80e9))
        @test_throws ErrorException Norma.elastic_constants(Parameters())
    end
end

@testset "J2 stress" begin
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
