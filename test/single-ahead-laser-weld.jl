# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Laser Weld HEX8 Clamped BCs" begin
    cp("../examples/ahead/single/laser-weld/dynamic/clamped/laser-weld.yaml", "laser-weld.yaml"; force=true)
    cp("../examples/ahead/single/laser-weld/laser-weld.g", "../../laser-weld.g"; force=true)
    input_file = "laser-weld.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 0.01
    params["time integrator"]["final time"] = 0.05
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("laser-weld.yaml")
    rm("../../laser-weld.g")
    rm("laser-weld.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -6.297730171953009e-5 atol = 1e-12
    @test min_disp_y ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z ≈ 0.0007221807237683814 atol = 1e-12
    @test avg_stress ≈
        [15742.885247375629 1.6479602309168996e6 20705.484327184666 -14.635248791516172 1.5557909225756466e-8 3.578215426191114e-8] atol =
        1e-6
end

@testset "AHeaD Single Quasistatic Laser Weld HEX8 Symmetry BCs" begin
    cp("../examples/ahead/single/laser-weld/quasistatic/symmetry/laser-weld.yaml", "laser-weld.yaml"; force=true)
    cp("../examples/ahead/single/laser-weld/laser-weld.g", "../../laser-weld.g"; force=true)
    input_file = "laser-weld.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 0.1
    params["time integrator"]["final time"] = 0.5
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("laser-weld.yaml")
    rm("../../laser-weld.g")
    rm("laser-weld.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.006324585400834512 atol = 1e-12
    @test min_disp_y ≈ -0.04999999999999999 atol = 1e-12
    @test max_disp_z ≈ 0.0014201639297437146 atol = 1e-12
    @test avg_stress ≈
        [-124951.56298187797 1.3542079941106066e8 1.5670276966246015e6 1132.7823840653277 920346.061777385 2765.970872725516] atol =
        1.0e1
end
