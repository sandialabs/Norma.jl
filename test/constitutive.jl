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

    @testset "Reciprocal_Neohookean" begin
        params = Norma.Parameters("elastic modulus" => 100.0, "Poisson's ratio" => 0.3, "density" => 7800.0)
        mat = Norma.Reciprocal_Neohookean(params)
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

@testset "Hyperelastic Objectivity Under Rigid Rotation" begin
    # For a hyperelastic material σ(F) = σ(R·U) should equal R·σ(U)·Rᵀ
    # (frame-indifference). This test would fail if F were stored as its
    # transpose anywhere in the constitutive → Cauchy pipeline.
    θ = π / 4
    R = SMatrix{3,3,Float64,9}([
        cos(θ) -sin(θ) 0.0
        sin(θ)  cos(θ) 0.0
        0.0     0.0    1.0
    ])
    U = SMatrix{3,3,Float64,9}([
        1.20 0.05 0.00
        0.05 0.95 0.00
        0.00 0.00 1.10
    ])
    F_nostretch_rotation_only = R
    F_stretch = U
    F_combined = R * U

    function cauchy(mat, F)
        W, P, _ = Norma.constitutive(mat, F)
        J = det(F)
        return F * P' / J
    end

    params = Norma.Parameters("elastic modulus" => 1.0e9, "Poisson's ratio" => 0.3, "density" => 1.0)
    for mat in (
        Norma.SaintVenant_Kirchhoff(params),
        Norma.Neohookean(params),
    )
        σ_U = cauchy(mat, F_stretch)
        σ_RU = cauchy(mat, F_combined)
        @test σ_RU ≈ R * σ_U * R' atol = 1.0e-08 * norm(σ_U)
        # Pure rotation → zero Cauchy
        σ_R = cauchy(mat, F_nostretch_rotation_only)
        @test norm(σ_R) < 1.0e-06
    end
end

@testset "J2Plasticity Constitutive At Identity" begin
    E = 200.0e9
    ν = 0.25
    σy = 1.0e9
    H = 20.0e9
    params = Norma.Parameters(
        "elastic modulus" => E, "Poisson's ratio" => ν,
        "density" => 7800.0, "yield stress" => σy, "hardening modulus" => H
    )
    mat = Norma.J2Plasticity(params)
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    state0 = Norma.initial_state(mat)
    W, P, AA, state1 = Norma.constitutive(mat, F, state0)
    @test isapprox(W, 0.0; atol=1e-12)
    @test isapprox(norm(P), 0.0; atol=1e-6)
    @test size(AA) == (3, 3, 3, 3)
    # No plastic flow at identity
    @test isapprox(state1[10], 0.0; atol=1e-12)
end

@testset "J2Plasticity Elastic Then Plastic Step" begin
    E = 200.0e9
    ν = 0.25
    σy = 1.0e9
    H = 0.0   # perfect plasticity
    μ = E / (2 * (1 + ν))
    params = Norma.Parameters(
        "elastic modulus" => E, "Poisson's ratio" => ν,
        "density" => 7800.0, "yield stress" => σy, "hardening modulus" => H
    )
    mat = Norma.J2Plasticity(params)
    state0 = Norma.initial_state(mat)

    # Pure shear at exactly the yield point (von Mises): γ = σy / (√3 μ)
    γ = σy / (sqrt(3.0) * μ)
    F_yield = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    W, P, AA, state1 = Norma.constitutive(mat, F_yield, state0)
    # Fᵖ should remain identity (just at yield, no plastic flow yet)
    @test isapprox(state1[10], 0.0; atol=1e-6)

    # Beyond yield: double the shear strain
    F_plastic = @SMatrix [1.0 2γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    W2, P2, AA2, state2 = Norma.constitutive(mat, F_plastic, state0)
    # Equivalent plastic strain should be positive
    @test state2[10] > 0.0
    # Cauchy stress von Mises should be ≈ σy (perfect plasticity)
    J = det(Matrix{Float64}(F_plastic))
    σ = Matrix{Float64}(F_plastic) * Matrix{Float64}(P2)' / J
    σdev = σ - tr(σ) / 3 * I(3)
    σvm = sqrt(1.5) * norm(σdev)
    @test isapprox(σvm, σy; rtol=1e-3)
end

@testset "J2Plasticity Consistent Tangent" begin
    # Verify the LIVE Simo-Hughes tangent (_sh_j2_tangent, reached through
    # constitutive) against a central-difference Jacobian of the stress update.
    #
    # Elastic branch: the tangent is the exact Jacobian.
    #
    # Plastic branch: it is exactly the MAJOR-SYMMETRIC PART of the exact
    # Jacobian, and not the Jacobian itself.  The Simo-Hughes return map
    # evaluates the effective shear modulus μ̄ = μ tr(b̄ᵉ_trial)/3 at the trial
    # state, so the discrete update is not exactly variational and its true
    # Jacobian carries a small antisymmetric part.  BOX 9.2 discards that part
    # by construction, which is what a symmetric solver wants.  Measured below:
    # the symmetric part matches to ~1e-8 while the discarded antisymmetric part
    # is O(1e-4) of the tangent norm.  Newton stays fast; asymptotic convergence
    # is not fully quadratic on steps that yield.
    params = Norma.Parameters(
        "elastic modulus" => 200.0e9, "Poisson's ratio" => 0.3,
        "density" => 7800.0, "yield stress" => 250.0e6, "hardening modulus" => 20.0e9,
    )
    mat = Norma.J2Plasticity(params)
    state0 = Norma.initial_state(mat)

    # Central-difference ∂P/∂F, holding the old state fixed, so this measures
    # the algorithmic (consistent) tangent and not the continuum one.
    function fd_tangent(mat, F, state_old; h=1.0e-7)
        AA = zeros(3, 3, 3, 3)
        for k in 1:3, l in 1:3
            Fp = MMatrix{3,3,Float64,9}(F)
            Fm = MMatrix{3,3,Float64,9}(F)
            Fp[k, l] += h
            Fm[k, l] -= h
            _, Pp, _, _ = Norma.constitutive(mat, SMatrix{3,3,Float64,9}(Fp), state_old; need_tangent=false)
            _, Pm, _, _ = Norma.constitutive(mat, SMatrix{3,3,Float64,9}(Fm), state_old; need_tangent=false)
            for i in 1:3, j in 1:3
                AA[i, j, k, l] = (Pp[i, j] - Pm[i, j]) / (2h)
            end
        end
        return AA
    end
    major_sym(A) = [0.5 * (A[i, j, k, l] + A[k, l, i, j]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    major_asym(A) = [0.5 * (A[i, j, k, l] - A[k, l, i, j]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

    # Elastic step: exact Jacobian, and no antisymmetric part to discard.
    F_el = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.001]
    _, _, AA_el, state_el = Norma.constitutive(mat, F_el, state0)
    @test state_el[10] == state0[10]          # did not yield
    AA_el_fd = fd_tangent(mat, F_el, state0)
    @test size(AA_el_fd) == (3, 3, 3, 3)
    @test norm(AA_el_fd - Array(AA_el)) / norm(AA_el_fd) < 1.0e-7
    @test norm(major_asym(AA_el_fd)) / norm(AA_el_fd) < 1.0e-7

    # Plastic steps: symmetric part exact, antisymmetric part discarded.
    F_pre = @SMatrix [0.98 0.02 0.0; 0.0 0.97 0.01; 0.0 0.0 1.05]
    _, _, _, s_pre = Norma.constitutive(mat, F_pre, state0; need_tangent=false)
    plastic_cases = (
        ("isochoric shear", @SMatrix([1.0 0.032 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]), state0),
        ("triaxial + shear", @SMatrix([1.10 0.05 0.0; 0.02 0.95 0.01; 0.0 0.01 0.97]), state0),
        ("preloaded state", @SMatrix([0.97 0.03 0.0; 0.0 0.96 0.02; 0.0 0.0 1.08]), s_pre),
    )
    for (name, F, s_old) in plastic_cases
        @testset "$name" begin
            _, _, AA, s_new = Norma.constitutive(mat, F, s_old)
            @test s_new[10] > s_old[10]           # this step yielded
            AA_fd = fd_tangent(mat, F, s_old)
            scale = norm(AA_fd)
            # The coded tangent IS the symmetric part, to FD accuracy.
            @test norm(major_sym(AA_fd) - Array(AA)) / scale < 1.0e-5
            # ... and it is NOT the full Jacobian: the discarded antisymmetric
            # part is real and measurable, which is why the bound above is on
            # major_sym rather than on the difference itself.
            asym = norm(major_asym(AA_fd)) / scale
            @test 1.0e-6 < asym < 1.0e-2
            @test isapprox(norm(AA_fd - Array(AA)) / scale, asym; rtol=0.05)
        end
    end
end
