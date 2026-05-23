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

    rm("torsion.yaml"; force=true)
    rm("../torsion.g"; force=true)
    rm("torsion.e"; force=true)

    min_disp_x = minimum(model.displacement[1, :])
    min_disp_y = minimum(model.displacement[2, :])
    max_disp_z = maximum(model.displacement[3, :])
    avg_stress = average_components(model.stress)

    @test min_disp_x ≈ -0.044380857564717574 atol = 1e-12
    @test min_disp_y ≈ -0.04438085756471812 atol = 1e-12
    @test max_disp_z ≈ 8.033797118034425e-5 atol = 1e-12
    @test avg_stress ≈
        [-569805.695585395 -569805.6955840682 -1.0820381373333016e6 -5.802448868053034e-8 -5.752426659455523e-8 -8.791887495362972e-7] atol =
        1.0e1
end
