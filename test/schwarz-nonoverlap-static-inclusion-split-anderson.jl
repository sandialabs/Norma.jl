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

# Order of the subdomains under Anderson acceleration (issue #233): a Neumann
# side solved after the first subdomain with a relaxed Dirichlet side and before
# its own Dirichlet partner is flagged. Pairs are (Dirichlet index, Neumann
# index) in the order of `domains`. The cases are the four orders of the
# five-subdomain inclusion of the issue (circle, pp, pn, nn, np): only the
# order with every Dirichlet side first converged.
@testset "Schwarz Nonoverlap Anderson Neumann Before Dirichlet" begin
    flagged(pairs) = Set(Norma.neumann_before_dirichlet(pairs))
    # Circle first, each rectangle the Dirichlet side toward the rectangles after it.
    @test isempty(flagged([(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]))
    # Circle first, cyclic roles pp -> np -> nn -> pn -> pp.
    @test flagged([(1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (5, 4), (4, 3), (3, 2)]) == Set([(5, 4), (4, 3), (3, 2)])
    # Circle last: every rectangle is a Neumann side solved before the circle.
    @test flagged([(5, 1), (5, 2), (5, 3), (5, 4), (1, 2), (1, 4), (2, 3), (3, 4)]) ==
        Set([(5, 1), (5, 2), (5, 3), (5, 4)])
    # Two subdomains with the Neumann side first are not affected.
    @test isempty(flagged([(2, 1)]))
    @test isempty(flagged(Tuple{Int,Int}[]))
end
