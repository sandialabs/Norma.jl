# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Clamped bar (z ∈ [-L/2, L/2]) initialised with the Gaussian pulse
# u(z, 0) = α·exp(-z²/(2s²)), v(z, 0) = 0.  The d'Alembert split solution
#   u(z, t) = (α/2)·[g(z-ct) + g(z+ct)],   g(ξ) = exp(-ξ²/(2s²))
# is exact for t < (L/2)/c (before the front first hits a clamped end).
# Derivatives in t give v and a.  Sampling on the mesh node z-coordinates
# yields the reference extrema this test compares the simulation against.
@testset "Single Implicit Dynamic Solid Clamped" begin
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    simulation = Norma.run("clamped.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("clamped.yaml"; force=true)
    rm("clamped.g"; force=true)
    rm("clamped.e"; force=true)

    α = 0.01
    s = 0.02
    c = sqrt(1.0e9 / 1000.0)
    t = 1.0e-5
    ct = c * t

    g(ξ) = exp(-ξ^2 / (2 * s^2))
    u_ref(z) = (α / 2) * (g(z - ct) + g(z + ct))
    v_ref(z) = (α * c / (2 * s^2)) * ((z - ct) * g(z - ct) - (z + ct) * g(z + ct))
    a_ref(z) = (α * c^2 / (2 * s^4)) *
               (((z - ct)^2 - s^2) * g(z - ct) + ((z + ct)^2 - s^2) * g(z + ct))

    z_nodes = model.reference[3, :]
    u_field = u_ref.(z_nodes)
    v_field = v_ref.(z_nodes)
    a_field = a_ref.(z_nodes)

    max_disp = maximum_components(integrator.displacement)
    max_velo = maximum_components(integrator.velocity)
    max_acce = maximum_components(integrator.acceleration)
    min_disp = minimum_components(integrator.displacement)
    min_velo = minimum_components(integrator.velocity)
    min_acce = minimum_components(integrator.acceleration)

    # Lateral motion is identically zero.
    @test max_disp[1] ≈ 0.0 atol = 1.0e-06
    @test max_disp[2] ≈ 0.0 atol = 1.0e-06
    @test min_disp[1] ≈ 0.0 atol = 1.0e-06
    @test min_disp[2] ≈ 0.0 atol = 1.0e-06
    @test max_velo[1] ≈ 0.0 atol = 1.0e-06
    @test max_velo[2] ≈ 0.0 atol = 1.0e-06
    @test min_velo[1] ≈ 0.0 atol = 1.0e-06
    @test min_velo[2] ≈ 0.0 atol = 1.0e-06
    @test max_acce[1] ≈ 0.0 atol = 1.0e-06
    @test max_acce[2] ≈ 0.0 atol = 1.0e-06
    @test min_acce[1] ≈ 0.0 atol = 1.0e-06
    @test min_acce[2] ≈ 0.0 atol = 1.0e-06

    # Axial extrema match the d'Alembert reference sampled at the mesh nodes.
    @test max_disp[3] ≈ maximum(u_field) rtol = 8.0e-05
    @test min_disp[3] ≈ minimum(u_field) atol = 1.0e-06
    @test max_velo[3] ≈ maximum(v_field) rtol = 5.0e-05
    @test min_velo[3] ≈ minimum(v_field) rtol = 2.0e-04
    @test max_acce[3] ≈ maximum(a_field) rtol = 5.0e-04
    @test min_acce[3] ≈ minimum(a_field) rtol = 9.0e-04
end
