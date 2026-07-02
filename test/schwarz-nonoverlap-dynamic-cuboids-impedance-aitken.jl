# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Run the non-overlapping impedance (Robin-type) DD cuboid for a few steps with
# the requested relaxation method.
function run_cuboids_impedance(relaxation; num_steps=5)
    src = "../examples/nonoverlap/dynamic-same-step/cuboids-impedance"
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"]
        cp("$src/$f", f; force=true)
    end
    params = YAML.load_file("cuboids.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cuboids.yaml"
    params["final time"] = num_steps * 0.01
    relaxation !== nothing && (params["relaxation"] = relaxation)
    sim = Norma.run(params)
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml",
              "cuboid-1.g", "cuboid-2.g", "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cuboids Impedance Aitken" begin
    sim_fixed = run_cuboids_impedance(nothing)
    sim_aitken = run_cuboids_impedance("aitken recursive")
    sim_secant = run_cuboids_impedance("aitken secant")

    @test sim_fixed.controller.relaxation_method == :fixed
    @test sim_aitken.controller.relaxation_method == :aitken_recursive
    @test sim_secant.controller.relaxation_method == :aitken_secant
    @test sim_aitken.failed == false
    @test sim_secant.failed == false

    # The impedance BC relaxes a force RHS (not a kinematic field), so every
    # relaxation method drives the same fixed point to the same Schwarz
    # tolerance: the converged solution must be identical.
    u_fixed = sim_fixed.subsims[1].model.displacement
    @test sim_aitken.subsims[1].model.displacement ≈ u_fixed rtol = 1.0e-06
    @test sim_secant.subsims[1].model.displacement ≈ u_fixed rtol = 1.0e-06

    # Aitken should not converge slower than the fixed relaxation.
    iters_fixed = sum(sim_fixed.controller.schwarz_iters)
    @test sum(sim_aitken.controller.schwarz_iters) ≤ iters_fixed
    @test sum(sim_secant.controller.schwarz_iters) ≤ iters_fixed
end
