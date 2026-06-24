# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Dirichlet BCs applied on side sets (ssx-, ssy-, ssz-, ssz+) instead of node
# sets.  Those side sets cover the same cube faces as the node sets in
# cube.yaml, so the solution must be identical to the Single Static Solid Cube
# reference.  This exercises SolidMechanicsSideSetDirichletBoundaryCondition.
@testset "Single Static Solid Cube Side Set DBC" begin
    cp("../examples/single/static-solid/cube/standard/cube-sideset-dbc.yaml", "cube-sideset-dbc.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube-sideset-dbc.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube-sideset-dbc.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube-sideset-dbc.e"; force=true)

    # The four BCs must all be the side-set variant, resolving their nodes from
    # side sets.
    ss_dbcs = filter(bc -> bc isa Norma.SolidMechanicsSideSetDirichletBoundaryCondition, model.boundary_conditions)
    @test length(ss_dbcs) == 4

    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)
    @test avg_disp[1] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[2] ≈ -0.125 rtol = 1.0e-06
    @test avg_disp[3] ≈ 0.500 rtol = 1.0e-06
    @test avg_stress[1] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[2] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test avg_stress[4] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[5] ≈ 0.0 atol = 1.0e-06
    @test avg_stress[6] ≈ 0.0 atol = 1.0e-06
end
