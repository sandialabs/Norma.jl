# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Dirichlet-Neumann Schwarz across a stiffness ratio of 100 (linear elastic
# cuboids). The gains of the iteration are the spectrum of -S_N^-1 S_D, which
# a ratio r scales by r: with the soft side as the Dirichlet side every gain
# is below 0.034 and the unrelaxed iteration converges in 8 iterations; with
# the stiff side as the Dirichlet side the gains reach -340 and Aitken
# recursive relaxation converges in 73. The example records the measurements.
@testset "Schwarz Nonoverlap Static Cuboids Stiffness Ratio" begin
    src = "../examples/nonoverlap/static-same-step/cuboids-stiffness-ratio/"
    for (input, max_iters) in [("soft-dirichlet", 10), ("stiff-dirichlet", 80)]
        domains = ["$input-soft", "$input-stiff"]
        files = vcat(["$input.yaml"], ["$d.yaml" for d in domains], ["cuboid-1.g", "cuboid-2.g"])
        for file in files
            cp(src * file, file; force=true)
        end
        sim = Norma.run("$input.yaml")
        for file in vcat(files, ["$d.e" for d in domains])
            rm(file; force=true)
        end
        @test sim.failed == false
        @test sim.controller.converged == true
        @test 0 < sim.controller.schwarz_iters[1] <= max_iters
    end
end
