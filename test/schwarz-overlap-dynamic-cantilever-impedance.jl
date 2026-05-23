# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

const cantilever_impedance_example = "../examples/overlap/dynamic-same-step/cantilever"

@testset "Schwarz Overlap Dynamic Cantilever Impedance" begin
    files = [
        "cantilever-impedance.yaml",
        "cantilever-free-impedance.yaml",
        "cantilever-clamped-impedance.yaml",
        "cantilever-free.g",
        "cantilever-clamped.g",
    ]
    for f in files
        cp(joinpath(cantilever_impedance_example, f), f; force=true)
    end
    sim = Norma.run("cantilever-impedance.yaml")
    model_free = sim.subsims[1].model
    model_clamped = sim.subsims[2].model
    for f in files
        rm(f; force=true)
    end
    rm("cantilever-free.e"; force=true)
    rm("cantilever-clamped.e"; force=true)

    @test sim.failed == false

    # Initial displacement comes from IC u_y = 0.393701 * x^2. At t = 3e-6 the
    # solution is essentially the initial parabola — just verify it's finite and nonzero.
    max_uy_free = maximum(model_free.displacement[2, :])
    max_uy_clamped = maximum(model_clamped.displacement[2, :])
    @test isfinite(max_uy_free)
    @test isfinite(max_uy_clamped)
    @test max_uy_free > 0.0
    @test max_uy_clamped > 0.0

    # Schwarz convergence bounded
    iters = sim.controller.schwarz_iters
    @test all(iters[1:3] .≤ 128)
end
