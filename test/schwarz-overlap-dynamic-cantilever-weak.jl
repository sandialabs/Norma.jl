# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

const cantilever_weak_example = "../examples/overlap/dynamic-same-step/cantilever-weak"

@testset "Schwarz Overlap Dynamic Cantilever Weak DBC" begin
    for f in readdir(cantilever_weak_example)
        cp(joinpath(cantilever_weak_example, f), f; force=true)
    end
    sim = Norma.run("cantilever.yaml")
    subsims = sim.subsims
    model_clamped = subsims[2].model
    model_free = subsims[1].model
    for f in readdir(cantilever_weak_example)
        rm(f; force=true)
    end
    rm("cantilever-clamped.e"; force=true)
    rm("cantilever-free.e"; force=true)
    max_disp_y_clamped = maximum(model_clamped.displacement[2, :])
    max_disp_y_free = maximum(model_free.displacement[2, :])
    @test max_disp_y_clamped ≈ 9.147e-03 rtol = 1.0e-02
    @test max_disp_y_free ≈ 2.509e-02 rtol = 1.0e-02
end
