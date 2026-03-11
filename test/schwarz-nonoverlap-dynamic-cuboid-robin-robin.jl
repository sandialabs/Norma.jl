# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

const rr_dynamic_example = "../examples/nonoverlap/dynamic-same-step/cuboids-robin-robin"

# Run 5 steps (t = 0.05) of the given input file and return the simulation.
function run_rr_dynamic(input_filename, extra_params=Dict())
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"]
        cp("$rr_dynamic_example/$(f == "cuboids.yaml" ? input_filename : f)", f; force=true)
    end
    params = YAML.load_file("cuboids.yaml"; dicttype=Norma.Parameters)
    params["final time"] = 0.05
    params["name"] = "cuboids.yaml"
    merge!(params, extra_params)
    sim = Norma.run(params)
    for f in ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g",
              "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cuboid Hex8 Robin-Robin Baseline" begin
    sim = run_rr_dynamic("cuboids.yaml")
    m1 = sim.subsims[1].model
    m2 = sim.subsims[2].model
    u1z = m1.current[3, :] .- m1.reference[3, :]
    u2z = m2.current[3, :] .- m2.reference[3, :]
    iters = sim.controller.schwarz_iters

    # Interface displacement continuity: top of domain 1 ≈ bottom of domain 2
    @test maximum(u1z) ≈ minimum(u2z) rtol = 1.0e-03

    # Applied BC at t=0.05: u_z(z=1) = 0.5 - 0.5*cos(π*0.05) ≈ 6.156e-3
    @test maximum(u2z) ≈ 6.156e-03 rtol = 1.0e-02

    # Schwarz iteration counts for first 5 steps are bounded
    @test all(iters[1:5] .≤ 20)
    @test sum(iters[1:5]) ≤ 80
end

@testset "Schwarz Nonoverlap Dynamic Cuboid Hex8 Robin-Robin Predictor" begin
    sim_base = run_rr_dynamic("cuboids.yaml")
    sim_pred = run_rr_dynamic("cuboids-predictor.yaml")
    m1 = sim_pred.subsims[1].model
    m2 = sim_pred.subsims[2].model
    u1z = m1.current[3, :] .- m1.reference[3, :]
    u2z = m2.current[3, :] .- m2.reference[3, :]

    # Same physical result
    @test maximum(u1z) ≈ minimum(u2z) rtol = 1.0e-03

    # Predictor reduces total Schwarz iterations over 5 steps
    iters_base = sum(sim_base.controller.schwarz_iters[1:5])
    iters_pred = sum(sim_pred.controller.schwarz_iters[1:5])
    @test iters_pred < iters_base
end

