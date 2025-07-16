# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "AHeaD Single Dynamic Torsion HEX8" begin
    cp("../examples/ahead/single/torsion/dynamic/torsion.yaml", "torsion.yaml"; force=true)
    cp("../examples/ahead/single/torsion/torsion.g", "../torsion.g"; force=true)
    input_file = "torsion.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["time integrator"]["initial time"] = 0.0
    params["time integrator"]["time step"] = 1.0e-6
    params["time integrator"]["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    model = sim.model

    rm("torsion.yaml")
    rm("../torsion.g")
    rm("torsion.e")

    min_disp_x = minimum(model.current[1, :] - model.reference[1, :])
    min_disp_y = minimum(model.current[2, :] - model.reference[2, :])
    max_disp_z = maximum(model.current[3, :] - model.reference[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.044380857564717574 atol = 1e-12
    @test min_disp_y ≈ -0.04438085756471812 atol = 1e-12
    @test max_disp_z ≈ 8.033797118034425e-5 atol = 1e-12
    @test avg_stress ≈
        [-1.3396757511504798e6 -1.33967575115108e6 457701.97380569484 1.0842104529729113e-7 3.0323481041705237e-7 1.7932685807409144e-7] atol =
        1.0e1
end
