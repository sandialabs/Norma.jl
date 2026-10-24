# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Bytes allocated by one explicit step. On 110 592 hexahedra a step allocated
# 549 MB, and the garbage collection of those bytes took half of the wall time
# on 12 threads: the stress was stored as a vector per quadrature point, the
# stable time step read the untyped shape function tables inside its element
# loop, and the predictor, corrector, step, and residual norm indexed with the
# mask of free degrees of freedom, which copies the vector. The bounds here are
# about five times the allocation measured after those were removed, and a
# tenth of the allocation before, so that a return of any of them fails.
using Test

@testset "Allocations Of An Explicit Step" begin
    cp("../examples/single/explicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/single/explicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    sim = Norma.create_simulation("clamped.yaml")
    Norma.initialize(sim)
    model, integrator, solver = sim.model, sim.integrator, sim.solver
    num_dofs = 3 * size(model.reference, 2)
    @test num_dofs == 12012
    # Warm up: compile every path of a step before measuring.
    for _ in 1:2
        Norma.advance_time(sim)
        Norma.advance_one_step(sim)
    end
    # The pieces that must not allocate at all.
    @test (@allocated Norma.predict(integrator, solver, model)) == 0
    @test (@allocated Norma.correct(integrator, solver, model)) == 0
    @test (@allocated Norma.compute_step(integrator, model, solver, solver.step)) == 0
    @test (@allocated Norma.masked_norm(solver.gradient, model.free_dofs)) == 0
    @test (@allocated Norma.stable_time_step(integrator, model)) < 4096
    # The whole step: the per-evaluation force vectors and the thread-local
    # element arrays remain, about 0.4 MB on this mesh with one thread.
    Norma.advance_time(sim)
    step_bytes = @allocated Norma.advance_one_step(sim)
    @test step_bytes < 2_000_000
    Norma.finalize_writing(sim)
    rm("clamped.yaml"; force=true)
    rm("clamped.g"; force=true)
    rm("clamped.e"; force=true)
end
