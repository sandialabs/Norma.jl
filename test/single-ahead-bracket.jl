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

    rm("bracket.yaml"; force=true)
    rm("../bracket.g"; force=true)
    rm("bracket.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    println(avg_stress)

    @test min_disp_x ≈ -0.0008933290654759424 atol = 1e-12
    @test min_disp_y ≈ -0.0001655247963407644 atol = 1e-12
    @test max_disp_z ≈ 0.0101486664023163 atol = 1e-12
    @test avg_stress ≈
        [-1.3069569134033617e7 6.19466287887055e6 3.811158192718337e6 -74052.105767898 2.9738669739486277e7 -4.072007814346489e6] atol =
        1.0e1
end
