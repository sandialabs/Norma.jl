# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

@testset "Smoothing Energy Stagnation" begin
    # Integration test for the energy stagnation exit (issue #220) on a mesh
    # that genuinely stagnates: the awful cube smoothed with L-BFGS over two
    # pseudo-time steps.  The first solve starts from the badly distorted mesh
    # and converges on the relative residual tolerance, which is reachable only
    # because the initial residual is enormous.  The second solve starts from
    # the smoothed mesh, where the energy reaches the floor its topology
    # permits within a few hundred iterations while the residual stays many
    # orders of magnitude above both tolerances, so without the stagnation exit
    # it rides the iteration cap doing nothing.  Iteration counts are
    # deliberately not pinned; they are one-ulp sensitive (see the history of
    # test 71).
    mesh = joinpath(@__DIR__, "..", "examples", "ems", "awful-cube", "awful-cube.g")
    out = "awful_cube_stagnation.e"

    params = Dict{String,Any}(
        "name" => "awful_cube_stagnation",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => out,
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing",
            "smooth reference" => "max",
            "material" => Dict{String,Any}(
                "blocks" => Dict{String,Any}("awful" => "elastic"),
                "elastic" => Dict{String,Any}(
                    "model" => "seth-hill", "m" => 2, "n" => 2,
                    "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                ),
            ),
        ),
        "time integrator" => Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 2.0, "time step" => 1.0
        ),
        "boundary conditions" => Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => "nsx-", "component" => "x", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsy-", "component" => "y", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsz-", "component" => "z", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsx+", "component" => "x", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsy+", "component" => "y", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsz+", "component" => "z", "function" => "0.0"),
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent", "step" => "lbfgs", "memory" => 10,
            "minimum iterations" => 1, "maximum iterations" => 1024,
            "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
            "step length" => 1.0e-3, "use line search" => true,
            "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
            "line search maximum iterations" => 16,
            "energy stagnation window" => 10,
            "energy stagnation tolerance" => 1.0e-6,
        ),
    )

    try
        rm(out; force=true)
        sim = Norma.run(params)

        # The second solve exited by declaring stagnation, not by the residual
        # tolerances and not by the iteration cap (a cap exit leaves the flag
        # down), and the run is accepted as successful.
        @test sim.solver.stagnated == true
        @test sim.solver.converged == true
        @test sim.model.failed == false
        @test sim.failed == false

        # The residual is far above both tolerances at exit, which is the point:
        # the gradient cannot know the energy floor has been reached.
        @test sim.solver.absolute_error > 1.0e2
        @test sim.solver.relative_error > 1.0e-4

        # The energy sits at the topology-limited floor: reduced by orders of
        # magnitude from the awful mesh, and strictly positive because the
        # reference is an ideal regular tet and regular tets do not tile space.
        @test sim.model.strain_energy < 5.0e5
        @test sim.model.strain_energy > 4.0e5
    finally
        rm(out; force=true)
    end
end
