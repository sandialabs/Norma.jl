# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The two estimates of the critical step of central difference (issue #243,
# src/stable_time_step.jl) on the undistorted hexahedral mesh of the cube
# example: they agree there, the trilinear hexahedron with lumped mass at
# Poisson's ratio zero has the critical step h/c_p, flattening the mesh scales
# the length estimate with the altitude, the estimate is recomputed every
# `stable time step interval` evaluations, and an unknown method aborts.

using YAML

function stable_step_simulation(overrides)
    input = YAML.load_file("../examples/single/explicit-dynamic-solid/cube/cube.yaml"; dicttype=Norma.Parameters)
    merge!(input["time integrator"], overrides)
    input["name"] = "cube-stable.yaml"
    return Norma.create_simulation(input)
end

@testset "Explicit Stable Time Step" begin
    cp("../examples/single/explicit-dynamic-solid/cube/cube.g", "cube.g"; force=true)

    sim = stable_step_simulation(Dict{String,Any}())
    integrator, model = sim.integrator, sim.model
    integrator.stable_method = :length
    length_step = Norma.critical_time_step(integrator, model)
    integrator.stable_method = :eigenvalue
    eigenvalue_step = Norma.critical_time_step(integrator, model)
    @test length_step ≈ eigenvalue_step rtol = 1.0e-10

    # Flattening the mesh by one half in z halves the altitude of every
    # element, and the length estimate with it; the eigenvalue estimate
    # cannot rise.
    model.reference[3, :] .*= 0.5
    integrator.stable_method = :length
    @test Norma.critical_time_step(integrator, model) ≈ 0.5 * length_step rtol = 1.0e-10
    integrator.stable_method = :eigenvalue
    @test Norma.critical_time_step(integrator, model) ≤ eigenvalue_step
    model.reference[3, :] .*= 2.0

    # Recomputed every third evaluation: a change of the mesh between
    # recomputations does not change the step.
    integrator.stable_method = :length
    integrator.stable_interval = 3
    integrator.stable_count = 0
    integrator.stable_critical = NaN
    first = Norma.stable_time_step(integrator, model)
    model.reference[3, :] .*= 0.5
    @test Norma.stable_time_step(integrator, model) == first
    @test Norma.stable_time_step(integrator, model) == first
    @test Norma.stable_time_step(integrator, model) ≈ 0.5 * first rtol = 1.0e-10
    model.reference[3, :] .*= 2.0
    Norma.finalize_writing(sim)

    # Trilinear hexahedron, lumped mass, Poisson's ratio zero: 2/ω_max = h/c_p.
    input = YAML.load_file("../examples/single/explicit-dynamic-solid/cube/cube.yaml"; dicttype=Norma.Parameters)
    input["model"]["material"]["elastic"]["Poisson's ratio"] = 0.0
    input["name"] = "cube-stable.yaml"
    sim = Norma.create_simulation(input)
    block = sim.model.blocks[1]
    material = sim.model.materials[1]
    c_p = sqrt(Norma.get_p_wave_modulus(material) / material.ρ)
    h = Norma.minimum_block_altitude(sim.model, block, block.N, block.dN, block.weights)
    sim.integrator.stable_method = :eigenvalue
    @test Norma.critical_time_step(sim.integrator, sim.model) ≈ h / c_p rtol = 1.0e-10
    Norma.finalize_writing(sim)

    Norma.NORMA_TEST_MODE[] = true
    try
        @test_throws Norma.NormaAbortException stable_step_simulation(
            Dict{String,Any}("stable time step" => "diameter")
        )
        @test_throws Norma.NormaAbortException stable_step_simulation(
            Dict{String,Any}("stable time step interval" => 0)
        )
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
    rm("cube.g"; force=true)
    rm("cube-stable.e"; force=true)
    rm("cube.e"; force=true)
end
