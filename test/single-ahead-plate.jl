# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Plate HEX8" begin
    cp("../examples/ahead/single/plate/dynamic/plate.yaml", "plate.yaml"; force=true)
    cp("../examples/ahead/single/plate/plate.g", "../plate.g"; force=true)
    input_file = "plate.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 1.0e-5
    params["time integrator"]["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("plate.yaml"; force=true)
    rm("../plate.g"; force=true)
    rm("plate.e"; force=true)

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)
    println("avg_stress = ", avg_stress)

    @test min_disp_x ≈ -0.00018834575405897889 atol = 1e-12
    @test min_disp_y ≈ -3.0412993194159355e-5 atol = 1e-12
    @test max_disp_z ≈ 0.002069427539933659 atol = 1e-12
    #@test avg_stress ≈
    #    [1.1020826350682063e6 -35067.15388944201 … 2.3344589081990052e7 5.46817767220394e-] atol =
    #    1.0e1
end
