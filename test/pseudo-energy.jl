# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Pseudo-energy Ẽ, the energy functional applied to the velocity and the acceleration (subdomain_pseudo_energy in
# overlap_energy.jl), Ẽ = 1/2 a^T A a + 1/2 v^T K v with
# A = M + (dt^2/2)(2 beta - gamma) K, on the undecomposed linear elastic
# cantilever written as a one-subdomain multidomain run. The Newmark scheme
# with beta = 1/4, gamma = 1/2 and the central difference scheme both conserve
# Ẽ of an undecomposed linear problem without loads, so its relative change
# over 50 steps measures roundoff. For Newmark the value written is also
# compared with the quadratic forms of the assembled mass and stiffness.

using YAML
using LinearAlgebra

const pseudo_example = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"

function run_monolithic_e2(explicit::Bool; num_steps=50, dt=5.0e-7)
    cp("$pseudo_example/cantilever.g", "cantilever.g"; force=true)
    sub = YAML.load_file("$pseudo_example/cantilever.yaml"; dicttype=Norma.Parameters)
    delete!(sub["time integrator"], "initial time")
    delete!(sub["time integrator"], "final time")
    if explicit
        sub["time integrator"] = Norma.Parameters(
            "type" => "central difference", "time step" => dt, "CFL" => 1.0, "γ" => 0.5
        )
        sub["solver"] = Norma.Parameters("type" => "explicit solver", "step" => "explicit")
    else
        sub["time integrator"]["time step"] = dt
        sub["solver"]["linear solver"] = "direct"
        sub["solver"]["linear solver relative tolerance"] = 1.0e-14
        sub["solver"]["absolute tolerance"] = 1.0e-6
    end
    YAML.write_file("cantilever.yaml", sub)
    params = Norma.Parameters(
        "type" => "multi",
        "name" => "mono-pseudo",
        "domains" => ["cantilever.yaml"],
        "Exodus output interval" => 1.0,
        "initial time" => 0.0,
        "final time" => num_steps * dt,
        "time step" => dt,
        "minimum iterations" => 1,
        "maximum iterations" => 4,
        "relative tolerance" => 1.0e-12,
        "absolute tolerance" => 1.0e-14,
        "blended energy output" => true,
    )
    sim = Norma.run(params)
    rows = readlines("mono-pseudo-energy.csv")
    header = split(rows[1], ",")
    data = [parse.(Float64, split(row, ",")) for row in rows[2:end]]
    for f in ["cantilever.g", "cantilever.e", "cantilever.yaml", "mono-pseudo-energy.csv"]
        rm(f; force=true)
    end
    return sim, header, data
end

@testset "Differentiated-System Energy: Newmark" begin
    sim, header, data = run_monolithic_e2(false)
    @test sim.failed == false
    @test length(data) == 51
    i_pseudo = findfirst(==("pseudo_energy_total"), header)
    i_e1 = findfirst(==("total_energy"), header)
    pseudo = [row[i_pseudo] for row in data]
    e1 = [row[i_e1] for row in data]
    @test maximum(abs.(pseudo ./ pseudo[1] .- 1.0)) < 1.0e-11
    @test maximum(abs.(e1 ./ e1[1] .- 1.0)) < 1.0e-11
    # The quadrature form of Ẽ equals the quadratic forms of the assembled
    # matrices (A = M for beta = 1/4, gamma = 1/2).
    model = sim.subsims[1].model
    a = vec(model.acceleration)
    v = vec(model.velocity)
    pseudo_matrix = 0.5 * dot(a, model.mass * a) + 0.5 * dot(v, model.stiffness * v)
    pseudo_written, _, _ = Norma.subdomain_pseudo_energy(sim.subsims[1])
    @test pseudo_written ≈ pseudo_matrix rtol = 1.0e-12
    @test pseudo_written ≈ pseudo[end] rtol = 1.0e-14
end

@testset "Differentiated-System Energy: Central Difference" begin
    sim, header, data = run_monolithic_e2(true)
    @test sim.failed == false
    @test length(data) == 51
    i_pseudo = findfirst(==("pseudo_energy_total"), header)
    pseudo = [row[i_pseudo] for row in data]
    @test maximum(abs.(pseudo ./ pseudo[1] .- 1.0)) < 1.0e-12
    # Without Schwarz side sets every row of the staggered kinetic energy is
    # an interior row.
    i_interface = findfirst(==("staggered_kinetic_interface"), header)
    @test all(row[i_interface] == 0.0 for row in data)
end
