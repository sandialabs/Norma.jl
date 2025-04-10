# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

@testset "Barycentric Shape Functions" begin
    @testset "D=2, N=3 (tri3)" begin
        ξ = @SVector [1 / 3, 1 / 3]
        N, dN, ddN = Norma.barycentric_shape_functions(Val(2), Val(3), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (2, 3)
        @test size(ddN) == (2, 2, 3)
    end

    @testset "D=3, N=4 (tetra4)" begin
        ξ = @SVector [1 / 4, 1 / 4, 1 / 4]
        N, dN, ddN = Norma.barycentric_shape_functions(Val(3), Val(4), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (3, 4)
        @test size(ddN) == (3, 3, 4)
    end

    @testset "D=3, N=10 (tetra10)" begin
        ξ = @SVector [0.1, 0.2, 0.3]
        N, dN, ddN = Norma.barycentric_shape_functions(Val(3), Val(10), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (3, 10)
        @test size(ddN) == (3, 3, 10)
    end
end

@testset "Lagrangian Shape Functions" begin
    @testset "D=1, N=2 (bar2)" begin
        ξ = @SVector [0.0]
        N, dN, ddN = Norma.lagrangian_shape_functions(Val(1), Val(2), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (1, 2)
        @test size(ddN) == (1, 1, 2)
    end

    @testset "D=2, N=4 (quad4)" begin
        ξ = @SVector [0.0, 0.0]
        N, dN, ddN = Norma.lagrangian_shape_functions(Val(2), Val(4), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (2, 4)
        @test size(ddN) == (2, 2, 4)
    end

    @testset "D=3, N=8 (hex8)" begin
        ξ = @SVector [0.0, 0.0, 0.0]
        N, dN, ddN = Norma.lagrangian_shape_functions(Val(3), Val(8), ξ)
        @test isapprox(sum(N), 1.0; atol=1e-12)
        @test size(dN) == (3, 8)
        @test size(ddN) == (3, 3, 8)
    end
end

@testset "Barycentric Quadrature And Interpolation" begin
    @testset "Tri3 1 Point" begin
        N, dN, w, ξ = Norma.barycentric(Val(2), Val(3), Val(1))
        @test size(N) == (3, 1)
        @test size(dN) == (2, 3, 1)
        @test length(w) == 1
        @test size(ξ) == (2, 1)
    end

    @testset "Tri3 3 Point" begin
        N, dN, w, ξ = Norma.barycentric(Val(2), Val(3), Val(3))
        @test size(N) == (3, 3)
        @test size(dN) == (2, 3, 3)
        @test length(w) == 3
        @test size(ξ) == (2, 3)
    end

    @testset "Tetra10 5 Point" begin
        N, dN, w, ξ = Norma.barycentric(Val(3), Val(10), Val(5))
        @test size(N) == (10, 5)
        @test size(dN) == (3, 10, 5)
    end
end

@testset "Lagrangian Quadrature And Interpolation" begin
    @testset "Bar2 2 Point" begin
        N, dN, w, ξ = Norma.lagrangian(Val(1), Val(2), Val(2))
        @test size(N) == (2, 2)
        @test size(dN) == (1, 2, 2)
        @test length(w) == 2
    end

    @testset "Quad4 4 Point" begin
        N, dN, w, ξ = Norma.lagrangian(Val(2), Val(4), Val(4))
        @test size(N) == (4, 4)
        @test size(dN) == (2, 4, 4)
        @test length(w) == 4
    end

    @testset "Quad4 9 Point" begin
        N, dN, w, ξ = Norma.lagrangian(Val(2), Val(4), Val(9))
        @test size(N) == (4, 9)
        @test size(dN) == (2, 4, 9)
        @test length(w) == 9
    end

    @testset "Hex8 8 Point" begin
        N, dN, w, ξ = Norma.lagrangian(Val(3), Val(8), Val(8))
        @test size(N) == (8, 8)
        @test size(dN) == (3, 8, 8)
        @test length(w) == 8
    end
end

@testset "Tetrahedron Interpolation" begin
    vertices = [
        0 1 0 0
        0 0 1 0
        0 0 0 1
    ] * 1.0
    element_type = Norma.TETRA4
    x1 = zeros(3)
    ξ1 = Norma.map_to_parametric(element_type, vertices, x1)
    x2 = ones(3) / 3.0
    ξ2 = Norma.map_to_parametric(element_type, vertices, x2)
    x3 = [1.0, 0.0, 0.0]
    ξ3 = Norma.map_to_parametric(element_type, vertices, x3)
    ξ4 = [0.4, 0.3, 0.2]
    ξ5 = [0.0, 0.5, 0.5]
    ξ6 = [1.0, 1.0, 1.0001]
    x4 = [0.5, 0.5, 0.5]
    x5 = [0.0, 0.0, 0.5]
    x6 = x2 * (1.0 + 1.0e-07)
    @test norm(ξ1) ≈ 0.0 atol = 1.0e-08
    @test norm(ξ2) ≈ sqrt(3.0) / 3.0 atol = 1.0e-08
    @test norm(ξ3) ≈ 1.0 atol = 1.0e-08
    @test Norma.is_inside_parametric(element_type, ξ4) == true
    @test Norma.is_inside_parametric(element_type, ξ5) == true
    @test Norma.is_inside_parametric(element_type, ξ6) == false
    @test Norma.is_inside(element_type, vertices, x4)[2] == false
    @test Norma.is_inside(element_type, vertices, x5)[2] == true
    @test Norma.is_inside(element_type, vertices, x6)[2] == true
end

@testset "Hexahedron Interpolation" begin
    vertices = [
        -1 1 1 -1 -1 1 1 -1
        -1 -1 1 1 -1 -1 1 1
        -1 -1 -1 -1 1 1 1 1
    ] * 0.5
    element_type = Norma.HEX8
    x1 = zeros(3)
    ξ1 = Norma.map_to_parametric(element_type, vertices, x1)
    x2 = 0.5 * ones(3)
    ξ2 = Norma.map_to_parametric(element_type, vertices, x2)
    x3 = -0.5 * ones(3)
    ξ3 = Norma.map_to_parametric(element_type, vertices, x3)
    ξ4 = [0.1, 0.5, 0.9]
    ξ5 = [-1.0, 1.0, 1.0]
    ξ6 = [1.0, 1.0, 1.0001]
    x4 = [0.4, 0.5, 0.6]
    x5 = [-0.5, 0.0, 0.5]
    x6 = [0.5, 0.5, 0.5] * (1.0 + 1.0e-07)
    @test norm(ξ1) ≈ 0.0 atol = 1.0e-08
    @test norm(ξ2 - ones(3)) ≈ 0.0 atol = 1.0e-08
    @test norm(ξ3 + ones(3)) ≈ 0.0 atol = 1.0e-08
    @test Norma.is_inside_parametric(element_type, ξ4) == true
    @test Norma.is_inside_parametric(element_type, ξ5) == true
    @test Norma.is_inside_parametric(element_type, ξ6) == false
    @test Norma.is_inside(element_type, vertices, x4)[2] == false
    @test Norma.is_inside(element_type, vertices, x5)[2] == true
    @test Norma.is_inside(element_type, vertices, x6)[2] == true
end
