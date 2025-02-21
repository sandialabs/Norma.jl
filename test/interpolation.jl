# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "interpolation tet" begin
    vertices = [
        0 1 0 0
        0 0 1 0
        0 0 0 1
    ] * 1.0
    element_type = "TETRA4"
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

@testset "interpolation hex" begin
    vertices = [
        -1 1 1 -1 -1 1 1 -1
        -1 -1 1 1 -1 -1 1 1
        -1 -1 -1 -1 1 1 1 1
    ] * 0.5
    element_type = "HEX8"
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
