# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Robin-Robin Schwarz with Anderson acceleration over five subdomains (issue
# #245): the circle-in-square inclusion with eight relaxed interfaces. The
# Robin datum is relaxed on the side listed first in `domains`, which receives
# the partner's datum of the previous Schwarz iteration in the Gauss-Seidel
# pass; Anderson acceleration then converges in 59 iterations. With the side
# listed later relaxed, the accelerated data were computed in the same pass,
# the iteration was not a function of them alone, and Anderson acceleration
# stalled with |ΔU| between 1e-4 and 1e-3 until the iteration limit.
@testset "Schwarz Nonoverlap Static Inclusion 5SD Robin-Robin Anderson" begin
    src = "../examples/nonoverlap/static-same-step/inclusion-5sd-robin-robin/"
    names = ["circle", "square_pp", "square_pn", "square_nn", "square_np"]
    files = vcat(["circle_in_square.yaml"], ["$n.yaml" for n in names], ["$(n)_tet4.g" for n in names])
    for file in files
        cp(src * file, file; force=true)
    end
    sim = Norma.run("circle_in_square.yaml")
    for file in vcat(files, ["$n.e" for n in names])
        rm(file; force=true)
    end
    @test sim.failed == false
    @test sim.controller.converged == true
    @test 0 < sim.controller.schwarz_iters[1] <= 70
end
