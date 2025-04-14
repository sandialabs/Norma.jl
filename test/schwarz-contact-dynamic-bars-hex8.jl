# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using YAML

@testset "Schwarz Contact Dynamic Bars Hex8" begin
    cp("../examples/contact/implicit-dynamic/2-bars/bars.yaml", "bars.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/2-bars/bar-1.yaml", "bar-1.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/2-bars/bar-2.yaml", "bar-2.yaml"; force=true)
    cp("../examples/contact/implicit-dynamic/2-bars/bar-1.g", "bar-1.g"; force=true)
    cp("../examples/contact/implicit-dynamic/2-bars/bar-2.g", "bar-2.g"; force=true)
    input_file = "bars.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = -1.0e-06
    params["final time"] = 1.0e-06
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_fine = subsims[1].model
    integrator_fine = subsims[1].integrator
    model_coarse = subsims[2].model
    integrator_coarse = subsims[2].integrator
    rm("bars.yaml")
    rm("bar-1.yaml")
    rm("bar-2.yaml")
    rm("bar-1.g")
    rm("bar-2.g")
    rm("bar-1.e")
    rm("bar-2.e")
    min_disp_x_fine = minimum(model_fine.current[1, :] - model_fine.reference[1, :])
    max_disp_y_fine = maximum(model_fine.current[2, :] - model_fine.reference[2, :])
    max_disp_z_fine = maximum(model_fine.current[3, :] - model_fine.reference[3, :])
    max_disp_x_coarse = maximum(model_coarse.current[1, :] - model_coarse.reference[1, :])
    max_disp_y_coarse = maximum(model_coarse.current[2, :] - model_coarse.reference[2, :])
    max_disp_z_coarse = maximum(model_coarse.current[3, :] - model_coarse.reference[3, :])
    potential_energy_fine = integrator_fine.stored_energy
    potential_energy_coarse = integrator_coarse.stored_energy
    kinetic_energy_fine = integrator_fine.kinetic_energy
    kinetic_energy_coarse = integrator_coarse.kinetic_energy
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ 0.0 atol = 2.0e-05
    @test max_disp_y_fine ≈ 1.0e-6 atol = 1.0e-06
    @test max_disp_z_fine ≈ 1.0e-6 atol = 1.0e-06
    @test max_disp_x_coarse ≈ 0.0 atol = 2.0e-05
    @test max_disp_y_coarse ≈ 1.0e-6 atol = 1.0e-06
    @test max_disp_z_coarse ≈ 1.0e-6 atol = 1.0e-06
    @test potential_energy_fine ≈ 5.0e-5 atol = 3.0e-06
    @test potential_energy_coarse ≈ 5.0e-5 atol = 4.0e-06
    @test kinetic_energy_fine ≈ 0.0 atol = 3.0e-06
    @test kinetic_energy_coarse ≈ 0.0 atol = 3.0e-06
    @test avg_stress_fine[1] ≈ -1.0e08 rtol = 1.0e-01
    @test avg_stress_coarse[1] ≈ -1.0e08 rtol = 1.0e-01
end
