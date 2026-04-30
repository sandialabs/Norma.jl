# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Single Static Solid Robin Bc" begin
    cp("../examples/single/static-solid/cube/robin-bc/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/neumann-bc/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    # Analytical solution for Robin BC: alpha * T_z + beta * u_z = g on z=1 face
    # With alpha=1, beta=E=1e9, g=1e6 at t=1:
    #   eps_z = g / (alpha * E + beta) = 1e6 / (1e9 + 1e9) = 5e-4
    #   u_z = eps_z * z,  avg u_z = eps_z * 0.5 = 2.5e-4
    #   eps_x = eps_y = -nu * eps_z = -0.25 * 5e-4 = -1.25e-4
    #   avg u_x = eps_x * 0.5 = -6.25e-5,  avg u_y = -6.25e-5
    #   sigma_zz = E * eps_z = 1e9 * 5e-4 = 5e5
    @test avg_disp[1] ≈ -6.25e-05 rtol = 1.0e-06
    @test avg_disp[2] ≈ -6.25e-05 rtol = 1.0e-06
    @test avg_disp[3] ≈  2.50e-04 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress[3] ≈ 5.0e+05 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-01
end
