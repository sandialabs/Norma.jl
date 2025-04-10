# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Tet10 Static Solid Cube" begin
    cp("../examples/element-types/tet10/cube/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/element-types/tet10/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    min_disp_x = minimum(integrator.displacement[1:3:end])
    min_disp_y = minimum(integrator.displacement[2:3:end])
    max_disp_z = maximum(integrator.displacement[3:3:end])
    avg_stress = average_components(model.stress)
    @test min_disp_x ≈ -0.25 rtol = 1.0e-06
    @test min_disp_y ≈ -0.25 rtol = 1.0e-06
    @test max_disp_z ≈ 1.00 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end
