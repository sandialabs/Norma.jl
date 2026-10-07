# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using YAML

@testset "Schwarz Nonoverlap Static Cuboid Hex8 Robin-Robin Same Step" begin
    cp("../examples/nonoverlap/static-same-step/cuboids-robin-robin/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/nonoverlap/static-same-step/cuboids-robin-robin/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/nonoverlap/static-same-step/cuboids-robin-robin/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.run("cuboids.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    min_disp_x_fine = minimum(model_fine.displacement[1, :])
    min_disp_y_fine = minimum(model_fine.displacement[2, :])
    max_disp_z_fine = maximum(model_fine.displacement[3, :])
    min_disp_x_coarse = minimum(model_coarse.displacement[1, :])
    min_disp_y_coarse = minimum(model_coarse.displacement[2, :])
    min_disp_z_coarse = minimum(model_coarse.displacement[3, :])
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test min_disp_x_fine ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_fine ≈ -0.125 rtol = 1.0e-06
    @test max_disp_z_fine ≈ 0.5 rtol = 1.0e-06
    @test min_disp_x_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_y_coarse ≈ -0.125 rtol = 1.0e-06
    @test min_disp_z_coarse ≈ 0.5 rtol = 1.0e-01
    @test avg_stress_fine[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_fine[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_fine[6] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[1] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[2] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[4] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[5] ≈ 0.0 atol = 1.0e-01
    @test avg_stress_coarse[6] ≈ 0.0 atol = 1.0e-01
end

# The relaxed side of a Robin-Robin pair relaxes the Robin datum only. A
# Neumann load on its interface nodes is applied unrelaxed, so the converged
# solution does not depend on the relaxation factor; relaxing the whole
# boundary force scaled such a load by 1/θ at the fixed point.
@testset "Schwarz Nonoverlap Robin-Robin Relaxation Leaves Interface Loads Unscaled" begin
    example = "../examples/nonoverlap/static-same-step/cuboids-robin-robin"
    mesh_dir = "../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann"
    function run_with(theta)
        cp("$mesh_dir/cuboid-1.g", "cuboid-1.g"; force=true)
        cp("$mesh_dir/cuboid-2.g", "cuboid-2.g"; force=true)
        cp("$example/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
        sub = YAML.load_file("$example/cuboid-2.yaml"; dicttype=Norma.Parameters)
        sub["boundary conditions"]["Neumann"] = [
            Norma.Parameters("side set" => "ssz-", "component" => "x", "function" => "1.0e+06 * t")
        ]
        YAML.write_file("cuboid-2.yaml", sub)
        top = YAML.load_file("$example/cuboids.yaml"; dicttype=Norma.Parameters)
        top["relaxation parameter"] = theta
        top["maximum iterations"] = 256
        YAML.write_file("cuboids.yaml", top)
        sim = Norma.run("cuboids.yaml")
        displacements = [copy(subsim.model.displacement) for subsim in sim.subsims]
        for f in ("cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g", "cuboid-1.e",
                  "cuboid-2.e")
            rm(f; force=true)
        end
        return displacements
    end
    full = run_with(1.0)
    half = run_with(0.5)
    for k in 1:2
        @test norm(half[k] - full[k]) ≤ 1.0e-6 * norm(full[k])
    end
    # The load moves the interface in x, so the test is not vacuous.
    @test maximum(abs.(full[2][1, :])) > 1.0e-4
end
