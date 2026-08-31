# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Three-field Hu-Washizu element with a projected volumetric strain.
#
# The element is exercised directly rather than through a simulation, because
# every property that defines it is visible at the element level and none of
# them need a solver: that an affine field is reproduced exactly, that the
# projection engages on a non-affine one, and that the tangent is the true
# Hessian of the energy the residual comes from.

@testset "Three-Field Element" begin
    hencky() = Norma.Hencky(
        Dict{String,Any}("elastic modulus" => 1.0e9, "Poisson's ratio" => 0.3, "density" => 1000.0)
    )

    "A straight-edged reference TETRA10, mid-side nodes at edge midpoints."
    function tet10_reference()
        v = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
        X = zeros(3, 10)
        X[:, 1:4] = v
        for (k, (a, b)) in enumerate(((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)))
            X[:, 4 + k] = 0.5 * (v[:, a] + v[:, b])
        end
        return X
    end

    "Energy of the same element under the ordinary displacement formulation."
    function displacement_energy(material, dN, w, X, u)
        x = X + reshape(u, 3, size(X, 2))
        E = 0.0
        for q in 1:length(w)
            dNdξ = dN[:, :, q]
            dXdξ = Norma.SMatrix{3,3,Float64,9}(dNdξ * X')
            F = Norma.SMatrix{3,3,Float64,9}(x * (dXdξ \ dNdξ)')
            E += Norma.strain_energy(material, F) * det(dXdξ) * w[q]
        end
        return E
    end

    element_type = Norma.element_type_from_string("TETRA10")
    X = tet10_reference()
    material = hencky()

    # A field with genuinely varying J, so the projection has something to do.
    nonaffine = 1.0e-3 * [
        sin(3.1 * X[1, i] + 1.7 * X[2, i]) * (j == 1) +
        cos(2.3 * X[3, i]) * (j == 2) +
        X[1, i] * X[2, i] * (j == 3) for i in 1:10 for j in 1:3
    ]

    @testset "affine fields are reproduced exactly" begin
        # Under an affine motion F is constant, so J is constant, so the
        # projection of log J onto ANY pressure space returns log J itself and
        # the element must coincide with the displacement element to roundoff.
        # This is the property that makes the formulation consistent: it can
        # only change the answer where the volumetric strain actually varies.
        A = [0.03 0.01 -0.02; -0.015 0.02 0.01; 0.005 -0.01 0.025]
        b = [0.001, -0.002, 0.003]
        u = vec(A * X .+ b)
        for (order, G) in ((0, 4), (1, 5))
            _, dN, w, ξ = Norma.isoparametric(element_type, G)
            E_disp = displacement_energy(material, dN, w, X, u)
            E_three, _, _ = Norma.three_field_element(material, order, dN, w, ξ, X, u)
            @test E_three ≈ E_disp rtol = 1.0e-12
        end
    end

    @testset "the projection engages on a non-affine field" begin
        # Non-vacuity.  A projection that never changes anything would satisfy
        # every other test in this file.
        _, dN, w, ξ = Norma.isoparametric(element_type, 4)
        E_disp = displacement_energy(material, dN, w, X, nonaffine)
        E_three, _, _ = Norma.three_field_element(material, 0, dN, w, ξ, X, nonaffine)
        @test abs(E_three - E_disp) / abs(E_disp) > 0.1
    end

    @testset "quadrature must exceed the pressure dimension" begin
        # With as many quadrature points as pressure degrees of freedom the L²
        # projection is exact interpolation at those points, so θ̄ == θ wherever
        # the energy is sampled and the element degenerates to a displacement
        # element -- silently.  Measured: relative energy difference exactly
        # 0.0 at order 1 with the four-point rule.
        _, dN, w, ξ = Norma.isoparametric(element_type, 4)
        E_disp = displacement_energy(material, dN, w, X, nonaffine)
        E_deg = Norma.three_field_energy(material, 1, dN, w, ξ, X, nonaffine)
        @test E_deg ≈ E_disp rtol = 1.0e-12      # the degeneracy is real...
        # ...so the entry point refuses it rather than running it.
        @test_throws ErrorException Norma.three_field_element(material, 1, dN, w, ξ, X, nonaffine)
        @test_throws ErrorException Norma.check_three_field_quadrature(1, 4)
        @test Norma.check_three_field_quadrature(1, 5) === nothing
        @test Norma.check_three_field_quadrature(0, 4) === nothing
    end

    @testset "residual and tangent are consistent with the energy" begin
        # Residual and tangent both come from the same energy by automatic
        # differentiation, so this checks the energy against finite differences
        # rather than checking AD against itself.
        _, dN, w, ξ = Norma.isoparametric(element_type, 5)
        for order in (0, 1)
            E0, R, K = Norma.three_field_element(material, order, dN, w, ξ, X, nonaffine)
            @test norm(K - K') / norm(K) < 1.0e-12

            h = 1.0e-7
            for i in (1, 7, 19, 30)
                up = copy(nonaffine); up[i] += h
                um = copy(nonaffine); um[i] -= h
                Ep = Norma.three_field_energy(material, order, dN, w, ξ, X, up)
                Em = Norma.three_field_energy(material, order, dN, w, ξ, X, um)
                @test (Ep - Em) / (2h) ≈ R[i] rtol = 1.0e-5
            end
        end
    end

    @testset "rigid translation costs nothing" begin
        _, dN, w, ξ = Norma.isoparametric(element_type, 5)
        t = repeat([0.01, -0.02, 0.03], 10)
        for order in (0, 1)
            E, R, _ = Norma.three_field_element(material, order, dN, w, ξ, X, t)
            @test E ≈ 0.0 atol = 1.0e-16
            @test norm(R) ≈ 0.0 atol = 1.0e-6
        end
    end

    @testset "unsupported materials are refused" begin
        # The subtraction of the pointwise volumetric energy is exact only for a
        # material whose energy splits additively in the logarithmic measure.
        # For anything else it would leave a volumetric remainder inside the
        # deviatoric term -- a quiet inaccuracy, not a visible failure.
        neo = Norma.Neohookean(
            Dict{String,Any}("elastic modulus" => 1.0e9, "Poisson's ratio" => 0.3, "density" => 1000.0)
        )
        @test Norma.three_field_supported(hencky()) == true
        @test Norma.three_field_supported(neo) == false
        _, dN, w, ξ = Norma.isoparametric(element_type, 5)
        @test_throws ErrorException Norma.three_field_element(neo, 0, dN, w, ξ, X, nonaffine)
    end

    @testset "pressure order is validated" begin
        @test_throws ErrorException Norma.pressure_shape(2, [0.25, 0.25, 0.25])
        @test length(Norma.pressure_shape(0, [0.25, 0.25, 0.25])) == 1
        @test length(Norma.pressure_shape(1, [0.25, 0.25, 0.25])) == 4
    end
end
