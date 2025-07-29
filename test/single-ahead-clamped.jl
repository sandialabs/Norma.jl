# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Clamped HEX8" begin
    cp("../examples/ahead/single/clamped/dynamic/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/ahead/single/clamped/clamped.g", "../clamped.g"; force=true)
    input_file = "clamped.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    t = 1.0e-4
    params["time integrator"]["time step"] = 1.0e-6
    params["time integrator"]["final time"] = t

    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("clamped.yaml")
    rm("../clamped.g")
    rm("clamped.e")

    z = model.reference[3, :]
    disp_x = model.current[1, :] - model.reference[1, :]
    disp_y = model.current[2, :] - model.reference[2, :]
    disp_z = model.current[3, :] - model.reference[3, :]
    velo_x = model.velocity[1, :]
    velo_y = model.velocity[2, :]
    velo_z = model.velocity[3, :]
    acce_x = model.acceleration[1, :]
    acce_y = model.acceleration[2, :]
    acce_z = model.acceleration[3, :]

    c = sqrt(1e9 / 1e3)
    a = 0.001
    b = 0.0
    s = 0.02
    T = 1.0e-3

    #Create and populate exact solution vectors 
    n = size(z)[1]
    disp_z_exact = zeros(Float64, n)
    velo_z_exact = zeros(Float64, n)
    acce_z_exact = zeros(Float64, n)
    for i in 1:n
        disp_z_exact[i] =
            1 / 2 * a * (exp(-(z[i] - c * t - b)^2 / 2 / s^2) + exp(-(z[i] + c * t - b)^2 / 2 / s^2)) -
            1 / 2 * a * (exp(-(z[i] - c * (T - t) - b)^2 / 2 / s^2) + exp(-(z[i] + c * (T - t) - b)^2 / 2 / s^2))
        velo_z_exact[i] =
            c / 2 * a / s^2 * (
                (z[i] - c * t - b) * exp(-(z[i] - c * t - b)^2 / 2 / s^2) -
                (z[i] + c * t - b) * exp(-(z[i] + c * t - b)^2 / 2 / s^2)
            ) +
            c / 2 * a / s^2 * (
                (z[i] - c * (T - t) - b) * exp(-(z[i] - c * (T - t) - b)^2 / 2 / s^2) -
                (z[i] + c * (T - t) - b) * exp(-(z[i] + c * (T - t) - b)^2 / 2 / s^2)
            )
        acce_z_exact[i] =
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z[i] - c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z[i] - c * t - b)^2) * exp(-1 / 2 * (z[i] - c * t - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z[i] + c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z[i] + c * t - b)^2) * exp(-1 / 2 * (z[i] + c * t - b)^2 / s^2)
            ) -
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z[i] - c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z[i] - c * (T - t) - b)^2) * exp(-1 / 2 * (z[i] - c * (T - t) - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z[i] + c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z[i] + c * (T - t) - b)^2) * exp(-1 / 2 * (z[i] + c * (T - t) - b)^2 / s^2)
            )
    end
    disp_z_err = norm(disp_z_exact - disp_z)
    disp_z_norm = norm(disp_z_exact)
    disp_z_relerr = disp_z_err / disp_z_norm
    velo_z_err = norm(velo_z_exact - velo_z)
    velo_z_norm = norm(velo_z_exact)
    velo_z_relerr = velo_z_err / velo_z_norm
    acce_z_err = norm(acce_z_exact - acce_z)
    acce_z_norm = norm(acce_z_exact)
    acce_z_relerr = acce_z_err / acce_z_norm

    @test disp_z_relerr ≈ 0.02203179834481467 atol = 1e-12
    @test velo_z_relerr ≈ 0.048747445420986746 atol = 1e-12
    @test acce_z_relerr ≈ 0.11433498445981101 atol = 1e-6
    @test norm(disp_x) ≈ 0.0 atol = 0.0
    @test norm(disp_y) ≈ 0.0 atol = 0.0
    @test norm(velo_x) ≈ 0.0 atol = 0.0
    @test norm(velo_y) ≈ 0.0 atol = 0.0
    @test norm(acce_x) ≈ 0.0 atol = 0.0
    @test norm(acce_y) ≈ 0.0 atol = 0.0
end
