# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Bracket HEX8" begin
    cp("../examples/ahead/single/bracket/dynamic/bracket.yaml", "bracket.yaml"; force=true)
    cp("../examples/ahead/single/bracket/bracket.g", "../bracket.g"; force=true)
    input_file = "bracket.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 1.0e-5
    params["time integrator"]["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("bracket.yaml")
    rm("../bracket.g")
    rm("bracket.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    println(avg_stress)

    @test min_disp_x ≈ -0.0008933290654759424 atol = 1e-12
    @test min_disp_y ≈ -0.0001655247963407644 atol = 1e-12
    @test max_disp_z ≈ 0.0101486664023163 atol = 1e-12
    @test avg_stress ≈
        [-8.971263149258394e6 6.383924802446201e6 -476409.7156016053 304437.4306066531 3.0513844880209036e7 -4.0955010966046983e6] atol =
        1.0e1
end
