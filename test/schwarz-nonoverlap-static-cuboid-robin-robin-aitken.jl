# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Run the Robin-Robin DD cuboid with the requested relaxation method.
function run_robin_robin(relaxation)
    src = "../examples/nonoverlap/static-same-step/cuboids-robin-robin"
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml"]
        cp("$src/$f", f; force=true)
    end
    cp("../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/nonoverlap/static-same-step/cuboids-dirichlet-neumann/cuboid-2.g", "cuboid-2.g"; force=true)
    params = YAML.load_file("cuboids.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cuboids.yaml"
    relaxation !== nothing && (params["relaxation"] = relaxation)
    sim = Norma.run(params)
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml",
              "cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
    return sim
end

# Both Aitken variants must reach the same physical solution as the fixed
# relaxation and converge in no more Schwarz iterations.
function check_robin_physics(sim)
    model_fine = sim.subsims[1].model
    model_coarse = sim.subsims[2].model
    @test sim.failed == false
    @test minimum(model_fine.displacement[1, :]) ≈ -0.125 rtol = 1.0e-06
    @test minimum(model_fine.displacement[2, :]) ≈ -0.125 rtol = 1.0e-06
    @test maximum(model_fine.displacement[3, :]) ≈ 0.5 rtol = 1.0e-06
    @test minimum(model_coarse.displacement[3, :]) ≈ 0.5 rtol = 1.0e-01
    avg_stress_fine = average_components(model_fine.stress)
    avg_stress_coarse = average_components(model_coarse.stress)
    @test avg_stress_fine[3] ≈ 5.0e+08 rtol = 1.0e-06
    @test avg_stress_coarse[3] ≈ 5.0e+08 rtol = 1.0e-06
end

@testset "Schwarz Nonoverlap Static Cuboid Hex8 Robin-Robin Aitken" begin
    sim_fixed = run_robin_robin(nothing)
    sim_aitken = run_robin_robin("aitken")
    sim_secant = run_robin_robin("aitken secant")

    @test sim_fixed.controller.relaxation_method == :fixed
    @test sim_aitken.controller.relaxation_method == :aitken
    @test sim_secant.controller.relaxation_method == :aitken_secant

    # Dynamic Aitken relaxation of the Robin RHS must not change the converged
    # physical solution.
    check_robin_physics(sim_aitken)
    check_robin_physics(sim_secant)

    # Aitken should not converge slower than the fixed relaxation on this
    # problem (it accelerates the Robin-Robin fixed-point iteration).
    iters_fixed = sum(sim_fixed.controller.schwarz_iters)
    @test sum(sim_aitken.controller.schwarz_iters) ≤ iters_fixed
    @test sum(sim_secant.controller.schwarz_iters) ≤ iters_fixed
end
