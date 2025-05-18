# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Overlap Dynamic Torsion" begin
    cp("../examples/ahead/overlap/torsion/dynamic/torsion.yaml", "torsion.yaml"; force=true)
    cp("../examples/ahead/overlap/torsion/dynamic/torsion-1.yaml", "torsion-1.yaml"; force=true)
    cp("../examples/ahead/overlap/torsion/dynamic/torsion-2.yaml", "torsion-2.yaml"; force=true)
    cp("../examples/ahead/overlap/torsion/torsion-1.g", "../torsion-1.g"; force=true)
    cp("../examples/ahead/overlap/torsion/torsion-2.g", "../torsion-2.g"; force=true)
    input_file = "torsion.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-6
    params["final time"] = 7.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_torsion1 = subsims[1].model
    model_torsion2 = subsims[2].model

    rm("torsion.yaml")
    rm("torsion-1.yaml")
    rm("torsion-2.yaml")
    rm("../torsion-1.g")
    rm("../torsion-2.g")
    rm("torsion-1.e")
    rm("torsion-2.e")

    min_disp_x_torsion1 = minimum(model_torsion1.current[1, :] - model_torsion1.reference[1, :])
    min_disp_y_torsion1 = minimum(model_torsion1.current[2, :] - model_torsion1.reference[2, :])
    max_disp_z_torsion1 = maximum(model_torsion1.current[3, :] - model_torsion1.reference[3, :])
    min_disp_x_torsion2 = minimum(model_torsion2.current[1, :] - model_torsion2.reference[1, :])
    min_disp_y_torsion2 = minimum(model_torsion2.current[2, :] - model_torsion2.reference[2, :])
    max_disp_z_torsion2 = maximum(model_torsion2.current[3, :] - model_torsion2.reference[3, :])
    avg_stress_torsion1 = average_components(model_torsion1.stress)
    avg_stress_torsion2 = average_components(model_torsion2.stress)

    println("min_disp_x_torsion1 = ", min_disp_x_torsion1, "\n")
    println("min_disp_y_torsion1 = ", min_disp_y_torsion1, "\n")
    println("max_disp_z_torsion1 = ", max_disp_z_torsion1, "\n")
    println("min_disp_x_torsion2 = ", min_disp_x_torsion2, "\n")
    println("min_disp_y_torsion2 = ", min_disp_y_torsion2, "\n")
    println("max_disp_z_torsion2 = ", max_disp_z_torsion2, "\n")
    println("avg_stress_torsion1 = ", avg_stress_torsion1, "\n")
    println("avg_stress_torsion2 = ", avg_stress_torsion2, "\n")

    @test min_disp_x_torsion1 ≈ 0.0 atol = 1e-8
    @test min_disp_y_torsion1 ≈ 0.0 atol = 1e-8
    @test max_disp_z_torsion1 ≈ 0.0 atol = 1e-8
    @test min_disp_x_torsion2 ≈ 0.0 atol = 1e-8
    @test min_disp_y_torsion2 ≈ 0.0 atol = 1e-8
    @test max_disp_z_torsion2 ≈ 0.0 atol = 1e-8
    @test avg_stress_torsion1 ≈ [1.4889709835472307e-8 6.412694449530423e-9 3.5272710678153343e-10 -2.9469083079595386e-8 1.4680773495310078e-8 -1.0532666845550373e-9] atol = 1e-6
    @test avg_stress_torsion2 ≈ [-5.669538912417606e-7 -5.865215720508659e-7 -4.190166732106221e-7 -6.759932128000284e-8 1.5515847522107731e-7 -1.8285069934175328e-7] atol = 1e-6

end

