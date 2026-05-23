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

    rm("laser-weld.yaml"; force=true)
    rm("../../laser-weld.g"; force=true)
    rm("laser-weld.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -6.297730171953009e-5 atol = 1e-12
    @test min_disp_y ≈ -0.0006155829702431115 atol = 1e-12
    @test max_disp_z ≈ 0.0007221807237683814 atol = 1e-12
    @test avg_stress ≈
        [15753.432582384918 1.6487160372148945e6 19939.130694208845 -14.023455872528954 1.7834705077557913e-8 3.44378979382038e-8] atol =
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

    rm("laser-weld.yaml"; force=true)
    rm("../../laser-weld.g"; force=true)
    rm("laser-weld.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.006324585400834512 atol = 1e-12
    @test min_disp_y ≈ -0.04999999999999999 atol = 1e-12
    @test max_disp_z ≈ 0.0014201639297437146 atol = 1e-12
    @test avg_stress ≈
        [-69354.06156634973 1.3762557588580975e8 -693346.2795400913 -1336.0754968765782 940957.2800525569 4167.462436451759] atol =
        1.0e1
end
