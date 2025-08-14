# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using Random
using StaticArrays

function finite_difference(material::Norma.Solid, F::SMatrix{3,3,Float64}, dF::SMatrix{3,3,Float64}, h::Float64)
    return (
        Norma.constitutive(material, F - 2 * h * dF)[1] - 8 * Norma.constitutive(material, F - h * dF)[1] +
        8 * Norma.constitutive(material, F + h * dF)[1] - Norma.constitutive(material, F + 2 * h * dF)[1]
    ) / (12 * h)
end

@testset "Constitutive Model Energy Gradient" begin
    Random.seed!(0)

    base_params = Norma.Parameters("elastic modulus" => 1.0, "Poisson's ratio" => 0.3, "density" => 1.0)
    sh_params = merge(base_params, Norma.Parameters("m" => 2, "n" => 2))
    models = [Norma.Neohookean(base_params), Norma.SaintVenant_Kirchhoff(base_params), Norma.SethHill(sh_params)]
    F_n = 10
    dF_n = 10

    for model in models
        for _ in 1:F_n
            F = SMatrix{3,3,Float64}(randn(3,3) * 0.1 + I)
            dWdF = Norma.constitutive(model, F)[2]
            for _ in 1:dF_n
                dF = SMatrix{3,3,Float64}(randn(9) * 0.1)
                dWdF_fd = finite_difference(model, F, dF, 1e-6)
                @test isapprox(sum(dWdF .* dF), dWdF_fd, atol=1e-6)
            end
        end
    end
end
