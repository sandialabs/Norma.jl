# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

@testset "Schwarz AHeaD Non-Overlap Dynamic Torsion HEX8-HEX8" begin
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion.yaml", "torsion.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion-1.yaml", "torsion-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/dynamic/torsion-2.yaml", "torsion-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/torsion-1.g", "../torsion-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/torsion/torsion-2.g", "../torsion-2.g"; force=true)
    input_file = "torsion.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 1.0e-6
    # Long enough for several transits of the torsional wave across the
    # non-conforming interface; the short 5.0e-5 horizon could not see the
    # interface effects investigated in issue #176.
    params["final time"] = 5.0e-4
    params["name"] = input_file
    sim = Norma.run(params)
    subsims = sim.subsims
    model_torsion1 = subsims[1].model
    model_torsion2 = subsims[2].model

    rm("torsion.yaml"; force=true)
    rm("torsion-1.yaml"; force=true)
    rm("torsion-2.yaml"; force=true)
    rm("../torsion-1.g"; force=true)
    rm("../torsion-2.g"; force=true)
    rm("torsion-1.e"; force=true)
    rm("torsion-2.e"; force=true)

    min_disp_x_torsion1 = minimum(model_torsion1.displacement[1, :])
    min_disp_y_torsion1 = minimum(model_torsion1.displacement[2, :])
    max_disp_z_torsion1 = maximum(model_torsion1.displacement[3, :])
    min_disp_x_torsion2 = minimum(model_torsion2.displacement[1, :])
    min_disp_y_torsion2 = minimum(model_torsion2.displacement[2, :])
    max_disp_z_torsion2 = maximum(model_torsion2.displacement[3, :])
    avg_stress_torsion1 = average_components(model_torsion1.stress)
    avg_stress_torsion2 = average_components(model_torsion2.stress)

    @test min_disp_x_torsion1 ≈ -0.04545670147801213 atol = 1e-10
    @test min_disp_y_torsion1 ≈ -0.04545670147800828 atol = 1e-10
    @test max_disp_z_torsion1 ≈ 0.00040594143734673033 atol = 1e-10
    @test min_disp_x_torsion2 ≈ -0.04438065776439695 atol = 1e-10
    @test min_disp_y_torsion2 ≈ -0.04438065776439266 atol = 1e-10
    @test max_disp_z_torsion2 ≈ -6.92733021142234e-7 atol = 1e-10
    @test avg_stress_torsion1 ≈
        [637892.225556833 637892.2255569943 -350476.29504338105 1.5970434787959676e-7 1.354154392174678e-7 -4.1442580425155027e-7] atol =
        1.0e1
    @test avg_stress_torsion2 ≈
        [-556091.8441462617 -556091.8441454791 -962694.439218634 3.046938218176365e-7 -6.516056600958109e-7 3.259494769736193e-6] atol =
        1.0e1
    @test sum(sim.controller.schwarz_iters) == 979
    @test maximum(sim.controller.schwarz_iters) == 2

    # Pure torsion about z keeps the bar axis on the z axis: any lateral
    # centerline displacement is symmetry error injected by the interface
    # transfer (issue #176 — the nearest-node facet search used to seed a
    # growing lateral mode here).
    for model in (model_torsion1, model_torsion2)
        onaxis = [
            i for i in 1:size(model.reference, 2) if
            abs(model.reference[1, i]) < 1.0e-12 && abs(model.reference[2, i]) < 1.0e-12
        ]
        @test maximum(abs.(model.displacement[1:2, onaxis])) < 1.0e-10
    end
end
