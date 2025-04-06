# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra: norm

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
