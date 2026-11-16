# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
@testset "Schwarz Nonoverlap Static Cuboid Hex8 Aitken Secant" begin
    src = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/"
    cp(src * "cuboids-aitken.yaml", "cuboids-aitken-secant.yaml"; force=true)
    cp(src * "cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp(src * "cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp(src * "cuboid-1.g", "cuboid-1.g"; force=true)
    cp(src * "cuboid-2.g", "cuboid-2.g"; force=true)

    # Select the Aitken secant (non-recursive, paper eq. 9) instead of the
    # Aitken-recursive Irons-Tuck form.
    write(
        "cuboids-aitken-secant.yaml",
        replace(
            read("cuboids-aitken-secant.yaml", String), "relaxation: aitken recursive" => "relaxation: aitken secant"
        ),
    )

    sim = Norma.run("cuboids-aitken-secant.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model

    rm("cuboids-aitken-secant.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)

    @test sim.controller.relaxation_method == :aitken_secant
    @test sim.failed == false

    # The first secant pair on the DN and overlap path is degenerate: a fresh
    # slot seeds the stored iterate with the incoming datum, so the next Schwarz iteration
    # finds g^(n) == g^(n-1), d = 0, and the raw factor -dot(d, δ)/‖δ‖² is
    # exactly zero. A zero factor freezes the interface, and in a chain of three
    # or more subdomains that freezes the whole Schwarz iteration and is read as convergence
    # after one iteration (issue #219). The degenerate pair must fall back to
    # the input theta instead. Replayed here on a key with no history.
    theta0 = sim.controller.relaxation_parameter
    degenerate_key = (97, 96, 1)
    v = [1.0, 2.0, 3.0]
    w = [1.5, 2.5, 3.5]
    Norma.relaxation_aitken_secant_theta!(sim.controller, degenerate_key, 1, 0, v, v)
    theta_second = Norma.relaxation_aitken_secant_theta!(sim.controller, degenerate_key, 1, 1, w, v)
    @test theta_second == theta0
    @test theta_second != 0.0

    # The smallest |θ| applied in a Schwarz iteration is recorded, so that the
    # convergence criterion can check the update divided by it.
    c = sim.controller
    saved = (c.converged, c.relative_error, c.absolute_error, c.relaxation_theta_min)
    c.relaxation_theta_min = Inf
    Norma.record_relaxation_factor!(c, 0.6)
    Norma.record_relaxation_factor!(c, -6.2e-8)
    Norma.record_relaxation_factor!(c, 1.3)
    @test c.relaxation_theta_min == 6.2e-8

    # The guard against a converged update produced by a collapsed factor
    # (refuse_relaxed_convergence!). The cases are final iterations of the
    # static 3D 5-subdomain beam study (November 2026), with the update
    # expressed as a fraction of the relative tolerance of that study (1e-8),
    # here applied to this run's tolerances. Each case enters with the update
    # criterion met.
    rtol = c.relative_tolerance
    atol = c.absolute_tolerance
    function guard(θ, relative, absolute)
        c.relaxation_theta_min = θ
        c.relative_error = relative
        c.absolute_error = absolute
        c.converged = true
        return Norma.refuse_relaxed_convergence!(sim)
    end
    # False convergences, whose error against the monolithic solution was 0.875
    # (global factor 6.2e-8) and 0.023 (global factor 1.27e-5, the mildest of
    # the 28), are refused, also when the absolute update meets its tolerance.
    @test guard(6.2e-8, 0.279 * rtol, 0.5 * atol) == true
    @test c.converged == false
    @test guard(1.27e-5, 0.907 * rtol, 0.5 * atol) == true
    @test c.converged == false
    # A factor of exactly zero, the case of the guard of August 2026.
    @test guard(0.0, 0.0, 0.0) == true
    # Genuine convergences are accepted: the global factor 0.187 at CR = 2, and
    # the smallest factor of any genuinely converged run, 7.06e-4.
    @test guard(0.187, 0.266 * rtol, 10.0 * atol) == false
    @test c.converged == true
    @test guard(7.06e-4, 0.582 * rtol, 10.0 * atol) == false
    @test c.converged == true
    # Fixed relaxation and Anderson acceleration (no factor recorded) are not affected.
    @test guard(0.2, rtol, 10.0 * atol) == false
    @test guard(Inf, rtol, 10.0 * atol) == false
    # An iteration that has not converged is left as it is.
    c.relaxation_theta_min = 1.0e-9
    c.converged = false
    @test Norma.refuse_relaxed_convergence!(sim) == false
    @test c.converged == false
    c.converged, c.relative_error, c.absolute_error, c.relaxation_theta_min = saved

    # Aitken secant must reach the same physical solution as classical
    # relaxation and the Aitken-recursive Irons-Tuck form.
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.5 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.5 rtol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
end
