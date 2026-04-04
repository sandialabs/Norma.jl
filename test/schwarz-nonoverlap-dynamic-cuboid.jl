# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

const dn_example = "../examples/nonoverlap/dynamic-same-step/cuboids-dirichlet-neumann"
const rr_example = "../examples/nonoverlap/dynamic-same-step/cuboids-robin-robin"

# Run the single-domain reference solution for 5 steps.
function run_single_domain()
    for f in ["cuboid.yaml", "cuboid.g"]
        cp("$dn_example/$f", f; force=true)
    end
    params = YAML.load_file("cuboid.yaml"; dicttype=Norma.Parameters)
    params["time integrator"]["final time"] = 0.05
    params["name"] = "cuboid.yaml"
    sim = Norma.run(params)
    for f in ["cuboid.yaml", "cuboid.g", "cuboid.e"]
        rm(f; force=true)
    end
    return sim
end

# Run a Schwarz simulation for 5 steps from the given example directory.
function run_schwarz(example_dir, input_filename, extra_params=Dict())
    for f in [input_filename, "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"]
        cp("$example_dir/$f", f; force=true)
    end
    params = YAML.load_file(input_filename; dicttype=Norma.Parameters)
    params["final time"] = 0.05
    params["name"] = input_filename
    merge!(params, extra_params)
    sim = Norma.run(params)
    for f in [input_filename, "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g",
              "cuboid-1.e", "cuboid-2.e"]
        rm(f; force=true)
    end
    return sim
end

# Extract interface displacements (z-component) from two subdomains.
function interface_displacements(sim)
    m1 = sim.subsims[1].model
    m2 = sim.subsims[2].model
    u1z = m1.current[3, :] .- m1.reference[3, :]
    u2z = m2.current[3, :] .- m2.reference[3, :]
    return u1z, u2z
end

@testset "Single Domain Reference" begin
    sim = run_single_domain()
    m = sim.model
    uz = m.current[3, :] .- m.reference[3, :]
    # Applied BC at t=0.05: u_z(z=0.75) = 0.5 - 0.5*cos(π*0.05) ≈ 6.156e-3
    @test maximum(uz) ≈ 6.156e-03 rtol = 1.0e-02
end

@testset "Schwarz Nonoverlap Dynamic Cuboid DN" begin
    sim = run_schwarz(dn_example, "cuboids.yaml")
    u1z, u2z = interface_displacements(sim)
    iters = sim.controller.schwarz_iters

    # Interface displacement continuity
    @test maximum(u1z) ≈ minimum(u2z) rtol = 1.0e-03

    # Applied BC matches single domain
    @test maximum(u2z) ≈ 6.156e-03 rtol = 1.0e-02

    # Schwarz convergence bounded
    @test all(iters[1:5] .≤ 30)
end

@testset "Schwarz Nonoverlap Dynamic Cuboid RR" begin
    sim = run_schwarz(rr_example, "cuboids.yaml")
    u1z, u2z = interface_displacements(sim)
    iters = sim.controller.schwarz_iters

    # Interface displacement continuity
    @test maximum(u1z) ≈ minimum(u2z) rtol = 1.0e-03

    # Applied BC matches single domain
    @test maximum(u2z) ≈ 6.156e-03 rtol = 1.0e-02

    # Schwarz convergence bounded
    @test all(iters[1:5] .≤ 20)
end
