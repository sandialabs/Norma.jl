# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Anderson acceleration over several interfaces (issue #233). A stiffer disk
# in a softer matrix, two subdomains with a nonconforming interface, coupled
# through four side sets (one per quadrant of the cut at x = 0.1, y = -0.08),
# so that the Dirichlet side relaxes four interfaces. Anderson acceleration
# solves one least-squares problem for the stacked traces of the four and
# converges in 23 Schwarz iterations (relaxation parameter 0.2, depth 20).
# With one least-squares problem per interface, the residual of each
# interface depends on the data of the other three, the per-interface secant
# models are wrong, and the iteration inverts an element (depths 20 and 50)
# or does not converge in 500 iterations (depth 10). The same interface as one
# side set needs 21 iterations.
@testset "Schwarz Nonoverlap Static Inclusion Split Anderson" begin
    src = "../examples/nonoverlap/static-same-step/matrix-inclusion-split/"
    files = ["inclusion.yaml", "disk.yaml", "matrix.yaml", "disk.g", "matrix.g"]
    for file in files
        cp(src * file, file; force=true)
    end
    sim = Norma.run("inclusion.yaml")
    for file in vcat(files, ["disk.e", "matrix.e"])
        rm(file; force=true)
    end
    @test sim.failed == false
    @test sim.controller.converged == true
    @test 0 < sim.controller.schwarz_iters[1] <= 30
    num_relaxed = count(
        bc -> bc isa Norma.SolidMechanicsNonOverlapSchwarzBoundaryCondition && bc.is_dirichlet,
        sim.subsims[1].model.boundary_conditions,
    )
    @test num_relaxed == 4
end
