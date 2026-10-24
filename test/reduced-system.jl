# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The system on the free degrees of freedom is gathered through a cached
# structure (ReducedSystem) instead of indexing the assembled matrix with the
# mask; it must equal that indexing entry for entry, reuse its structure while
# the pattern and the mask are unchanged, and be rebuilt when the mask changes.
using SparseArrays
using Test

@testset "Reduced System" begin
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    sim = Norma.create_simulation("clamped.yaml")
    Norma.initialize(sim)
    model, integrator, solver = sim.model, sim.integrator, sim.solver
    Norma.predict(integrator, solver, model)
    Norma.evaluate(integrator, solver, model)
    free = model.free_dofs
    @test count(free) < length(free)
    A = Norma.reduced_matrix!(solver, solver.hessian, free)
    @test A == solver.hessian[free, free]
    @test A.colptr == solver.hessian[free, free].colptr
    @test A.rowval == solver.hessian[free, free].rowval
    b = Norma.reduced_vector!(solver.reduced, solver.gradient)
    @test b == solver.gradient[free]
    # Same pattern and mask: the structure arrays are the same objects.
    A2 = Norma.reduced_matrix!(solver, solver.hessian, free)
    @test A2.colptr === A.colptr && A2.rowval === A.rowval
    # The mass shares the pattern of the stiffness, so it is reduced through
    # the same structure.
    M = Norma.reduced_matrix!(solver, model.mass, free)
    @test M == model.mass[free, free]
    @test M.colptr === A.colptr
    # A changed mask rebuilds the structure.
    changed = copy(free)
    changed[findfirst(free)] = false
    A3 = Norma.reduced_matrix!(solver, solver.hessian, changed)
    @test A3 == solver.hessian[changed, changed]
    @test A3.colptr !== A.colptr
    @test Norma.symmetric_tangent(model) === true
    Norma.finalize_writing(sim)
    rm("clamped.yaml"; force=true)
    rm("clamped.g"; force=true)
    rm("clamped.e"; force=true)
end
