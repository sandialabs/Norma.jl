# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Overlap Dynamic Clamped HEX8" begin
    cp("../examples/ahead/overlap/clamped/dynamic/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/ahead/overlap/clamped/dynamic/clamped-1.yaml", "clamped-1.yaml"; force=true)
    cp("../examples/ahead/overlap/clamped/dynamic/clamped-2.yaml", "clamped-2.yaml"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-1.g", "../clamped-1.g"; force=true)
    cp("../examples/ahead/overlap/clamped/clamped-2.g", "../clamped-2.g"; force=true)
    input_file = "clamped.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    t = 1.0e-4
    params["time step"] = 1.0e-6
    params["final time"] = t

    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_clamped0 = subsims[1].model
    model_clamped1 = subsims[2].model

    rm("clamped.yaml")
    rm("clamped-1.yaml")
    rm("clamped-2.yaml")
    rm("../clamped-1.g")
    rm("../clamped-2.g")
    rm("clamped-1.e")
    rm("clamped-2.e")

    z0 = model_clamped0.reference[3, :]
    disp0_x = model_clamped0.current[1, :] - model_clamped0.reference[1, :]
    disp0_y = model_clamped0.current[2, :] - model_clamped0.reference[2, :]
    disp0_z = model_clamped0.current[3, :] - model_clamped0.reference[3, :]
    velo0_x = model_clamped0.velocity[1, :]
    velo0_y = model_clamped0.velocity[2, :]
    velo0_z = model_clamped0.velocity[3, :]
    acce0_x = model_clamped0.acceleration[1, :]
    acce0_y = model_clamped0.acceleration[2, :]
    acce0_z = model_clamped0.acceleration[3, :]

    z1 = model_clamped1.reference[3, :]
    disp1_x = model_clamped1.current[1, :] - model_clamped1.reference[1, :]
    disp1_y = model_clamped1.current[2, :] - model_clamped1.reference[2, :]
    disp1_z = model_clamped1.current[3, :] - model_clamped1.reference[3, :]
    velo1_x = model_clamped1.velocity[1, :]
    velo1_y = model_clamped1.velocity[2, :]
    velo1_z = model_clamped1.velocity[3, :]
    acce1_x = model_clamped1.acceleration[1, :]
    acce1_y = model_clamped1.acceleration[2, :]
    acce1_z = model_clamped1.acceleration[3, :]

    c = sqrt(1e9 / 1e3)
    a = 0.001
    b = 0.0
    s = 0.02
    T = 1.0e-3

    #Create and populate exact solution vectors  
    n0 = size(z0)[1]
    disp0_z_exact = zeros(Float64, n0)
    velo0_z_exact = zeros(Float64, n0)
    acce0_z_exact = zeros(Float64, n0)
    for i in 1:n0
        disp0_z_exact[i] =
            1 / 2 * a * (exp(-(z0[i] - c * t - b)^2 / 2 / s^2) + exp(-(z0[i] + c * t - b)^2 / 2 / s^2)) -
            1 / 2 * a * (exp(-(z0[i] - c * (T - t) - b)^2 / 2 / s^2) + exp(-(z0[i] + c * (T - t) - b)^2 / 2 / s^2))
        velo0_z_exact[i] =
            c / 2 * a / s^2 * (
                (z0[i] - c * t - b) * exp(-(z0[i] - c * t - b)^2 / 2 / s^2) -
                (z0[i] + c * t - b) * exp(-(z0[i] + c * t - b)^2 / 2 / s^2)
            ) +
            c / 2 * a / s^2 * (
                (z0[i] - c * (T - t) - b) * exp(-(z0[i] - c * (T - t) - b)^2 / 2 / s^2) -
                (z0[i] + c * (T - t) - b) * exp(-(z0[i] + c * (T - t) - b)^2 / 2 / s^2)
            )
        acce0_z_exact[i] =
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z0[i] - c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z0[i] - c * t - b)^2) * exp(-1 / 2 * (z0[i] - c * t - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z0[i] + c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z0[i] + c * t - b)^2) * exp(-1 / 2 * (z0[i] + c * t - b)^2 / s^2)
            ) -
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z0[i] - c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z0[i] - c * (T - t) - b)^2) * exp(-1 / 2 * (z0[i] - c * (T - t) - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z0[i] + c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z0[i] + c * (T - t) - b)^2) * exp(-1 / 2 * (z0[i] + c * (T - t) - b)^2 / s^2)
            )
    end
    disp0_z_err = norm(disp0_z_exact - disp0_z)
    disp0_z_norm = norm(disp0_z_exact)
    disp0_z_relerr = disp0_z_err / disp0_z_norm
    velo0_z_err = norm(velo0_z_exact - velo0_z)
    velo0_z_norm = norm(velo0_z_exact)
    velo0_z_relerr = velo0_z_err / velo0_z_norm
    acce0_z_err = norm(acce0_z_exact - acce0_z)
    acce0_z_norm = norm(acce0_z_exact)
    acce0_z_relerr = acce0_z_err / acce0_z_norm

    @test disp0_z_relerr ≈ 0.022031767756555413 atol = 1e-12
    @test velo0_z_relerr ≈ 0.04874736761742479 atol = 1e-12
    @test acce0_z_relerr ≈ 0.11433470809918797 atol = 1e-6
    @test norm(disp0_x) ≈ 0.0 atol = 0.0
    @test norm(disp0_y) ≈ 0.0 atol = 0.0
    @test norm(velo0_x) ≈ 0.0 atol = 0.0
    @test norm(velo0_y) ≈ 0.0 atol = 0.0
    @test norm(acce0_x) ≈ 0.0 atol = 0.0
    @test norm(acce0_y) ≈ 0.0 atol = 0.0

    #Create and populate exact solution vectors  
    n1 = size(z1)[1]
    disp1_z_exact = zeros(Float64, n1)
    velo1_z_exact = zeros(Float64, n1)
    acce1_z_exact = zeros(Float64, n1)
    for i in 1:n1
        disp1_z_exact[i] =
            1 / 2 * a * (exp(-(z1[i] - c * t - b)^2 / 2 / s^2) + exp(-(z1[i] + c * t - b)^2 / 2 / s^2)) -
            1 / 2 * a * (exp(-(z1[i] - c * (T - t) - b)^2 / 2 / s^2) + exp(-(z1[i] + c * (T - t) - b)^2 / 2 / s^2))
        velo1_z_exact[i] =
            c / 2 * a / s^2 * (
                (z1[i] - c * t - b) * exp(-(z1[i] - c * t - b)^2 / 2 / s^2) -
                (z1[i] + c * t - b) * exp(-(z1[i] + c * t - b)^2 / 2 / s^2)
            ) +
            c / 2 * a / s^2 * (
                (z1[i] - c * (T - t) - b) * exp(-(z1[i] - c * (T - t) - b)^2 / 2 / s^2) -
                (z1[i] + c * (T - t) - b) * exp(-(z1[i] + c * (T - t) - b)^2 / 2 / s^2)
            )
        acce1_z_exact[i] =
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z1[i] - c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z1[i] - c * t - b)^2) * exp(-1 / 2 * (z1[i] - c * t - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z1[i] + c * t - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z1[i] + c * t - b)^2) * exp(-1 / 2 * (z1[i] + c * t - b)^2 / s^2)
            ) -
            1 / 2 *
            a *
            (
                -c^2 / s^2 * exp(-1 / 2 * (z1[i] - c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z1[i] - c * (T - t) - b)^2) * exp(-1 / 2 * (z1[i] - c * (T - t) - b)^2 / s^2) -
                c^2 / s^2 * exp(-1 / 2 * (z1[i] + c * (T - t) - b)^2 / s^2) +
                1 / s^4 * c^2 * ((z1[i] + c * (T - t) - b)^2) * exp(-1 / 2 * (z1[i] + c * (T - t) - b)^2 / s^2)
            )
    end
    disp1_z_err = norm(disp1_z_exact - disp1_z)
    disp1_z_norm = norm(disp1_z_exact)
    disp1_z_relerr = disp1_z_err / disp1_z_norm
    velo1_z_err = norm(velo1_z_exact - velo1_z)
    velo1_z_norm = norm(velo1_z_exact)
    velo1_z_relerr = velo1_z_err / velo1_z_norm
    acce1_z_err = norm(acce1_z_exact - acce1_z)
    acce1_z_norm = norm(acce1_z_exact)
    acce1_z_relerr = acce1_z_err / acce1_z_norm

    @test disp1_z_relerr ≈ 0.02203176775652283 atol = 1e-12
    @test velo1_z_relerr ≈ 0.0487473676174047 atol = 1e-12
    @test acce1_z_relerr ≈ 0.11433470809930901 atol = 1e-6
    @test norm(disp1_x) ≈ 0.0 atol = 0.0
    @test norm(disp1_y) ≈ 0.0 atol = 0.0
    @test norm(velo1_x) ≈ 0.0 atol = 0.0
    @test norm(velo1_y) ≈ 0.0 atol = 0.0
    @test norm(acce1_x) ≈ 0.0 atol = 0.0
    @test norm(acce1_y) ≈ 0.0 atol = 0.0

    @test sim.controller.schwarz_iters ≈ [
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
    ] atol = 0
end
