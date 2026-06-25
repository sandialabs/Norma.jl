# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

const cantilever_imp_example = "../examples/nonoverlap/dynamic-same-step/cantilever-impedance"

# Run the non-overlapping impedance DD cantilever (Newmark) for a few steps with
# the requested relaxation method.
function run_cantilever_imp_aitken(relaxation; num_steps=5)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_imp_example/$f", f; force=true)
    end
    params = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cantilever-multi.yaml"
    params["final time"] = num_steps * 5.0e-07
    relaxation !== nothing && (params["relaxation"] = relaxation)
    sim = Norma.run(params)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cantilever Impedance Aitken" begin
    sim_fixed = run_cantilever_imp_aitken(nothing)
    sim_aitken = run_cantilever_imp_aitken("aitken")
    sim_secant = run_cantilever_imp_aitken("aitken secant")

    @test sim_fixed.controller.relaxation_method == :fixed
    @test sim_aitken.controller.relaxation_method == :aitken
    @test sim_secant.controller.relaxation_method == :aitken_secant
    @test sim_aitken.failed == false
    @test sim_secant.failed == false

    # Impedance relaxes a single force RHS, so every method drives the same fixed
    # point to the same Schwarz tolerance: the solution must agree.
    u_fixed = sim_fixed.subsims[1].model.displacement
    @test sim_aitken.subsims[1].model.displacement ≈ u_fixed rtol = 1.0e-04
    @test sim_secant.subsims[1].model.displacement ≈ u_fixed rtol = 1.0e-04

    # Aitken accelerates vs the (under-relaxed) fixed default on this problem.
    iters_fixed = sum(sim_fixed.controller.schwarz_iters)
    @test sum(sim_aitken.controller.schwarz_iters) ≤ iters_fixed
    @test sum(sim_secant.controller.schwarz_iters) ≤ iters_fixed
end
