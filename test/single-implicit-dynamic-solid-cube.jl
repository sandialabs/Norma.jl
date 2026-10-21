# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using YAML

@testset "Single Implicit Dynamic Solid Cube" begin
    cp("../examples/single/implicit-dynamic-solid/cube/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/cube/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ 0.0 atol = 1.0e-06
    @test avg_disp[2] ≈ 0.0 atol = 1.0e-06
    @test avg_disp[3] ≈ 1.0 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end

# The consistent mass of the standard TETRA10, integrated with its default
# 5-point rule, whose centroid weight is negative, is indefinite. With
# `linear solver: direct`, the initial acceleration factors it with LU, where
# Cholesky threw a PosDefException. The TETRA10 cube of
# examples/element-types/tet10/cube, at rest, with its upper face accelerated
# along z.
@testset "Single Implicit Dynamic Solid Tet10 Cube Direct Solver" begin
    cp("../examples/element-types/tet10/cube/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/element-types/tet10/cube/cube.g", "cube.g"; force=true)
    params = YAML.load_file("cube.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cube.yaml"
    params["time integrator"] = Dict{String,Any}(
        "type" => "Newmark",
        "initial time" => 0.0,
        "final time" => 2.0e-4,
        "time step" => 1.0e-4,
        "β" => 0.25,
        "γ" => 0.5,
    )
    upper = only(bc for bc in params["boundary conditions"]["Dirichlet"] if bc["node set"] == "nsz+")
    upper["function"] = "1.0e3 * t^2"
    # The tolerances of the implicit dynamic cube: the residual of a Newmark
    # step stops at about 1e-11 relative.
    merge!(
        params["solver"],
        Dict{String,Any}("linear solver" => "direct", "relative tolerance" => 1.0e-10, "absolute tolerance" => 1.0e-06),
    )
    sim = Norma.create_simulation(params)
    try
        Norma.sync_control_time(sim)
        Norma.initialize(sim)
        model = sim.model
        free = model.free_dofs
        @test isposdef(Symmetric(Matrix(model.mass[free, free]))) == false
        # At rest, M_ff a_f = -M_fc a_c.
        a = sim.integrator.acceleration
        @test maximum(abs, a[.!free]) ≈ 2.0e3
        residual = model.mass[free, free] * a[free] + model.mass[free, .!free] * a[.!free]
        @test norm(residual) < 1.0e-10 * norm(model.mass[free, .!free] * a[.!free])
    finally
        Norma.finalize_writing(sim)
    end
    simulation = Norma.run(params)
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    @test simulation.failed == false
    @test simulation.integrator.time ≈ 2.0e-4
end
