# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Arlequin (distance-ratio) partition-of-unity energy for overlapping Schwarz.
# The weights remove the double count of the physical overlap region.  These
# tests pin the two exact identities the construction must satisfy (with unit
# weights the blended energy reduces to the monodomain strain and kinetic
# energy), the outright rejection of a degenerate overlap, and the end-to-end
# CSV output showing the double count removed.

using LinearAlgebra

# Overwrite every cached weight of a subdomain with 1 (bypassing the overlap
# correction) so the blended energy must reproduce the monodomain energy.
function force_unit_arlequin_weights!(subsim)
    aw = Norma.compute_arlequin_weights(subsim)   # exercises the real classifier
    for block_weights in aw.weights, element_weights in block_weights
        fill!(element_weights, 1.0)
    end
    Norma.ARLEQUIN_WEIGHT_CACHE[subsim.model] = aw
    return nothing
end

@testset "Overlap Blended Energy: Unit-Weight Stored Identity" begin
    cp("../examples/single/static-solid/cube/standard/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    sim = Norma.create_simulation("cube.yaml")
    model = sim.model
    # A subdomain with no overlap coupling has unit weights everywhere, so the
    # blended stored energy must equal model.strain_energy exactly.
    model.displacement .= 0.003 .* model.reference
    Norma.evaluate(model, sim.integrator, sim.solver)
    stored, kinetic = Norma.blended_subdomain_energy(sim)
    @test stored ≈ model.strain_energy rtol = 1.0e-10
    @test kinetic == 0.0
    Exodus.close(sim.params["input_mesh"])
    Exodus.close(sim.params["output_mesh"])
    empty!(Norma.ARLEQUIN_WEIGHT_CACHE)
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
end

@testset "Overlap Blended Energy: Unit-Weight Kinetic Identity" begin
    cp("../examples/overlap/dynamic-same-step/cantilever/cantilever.yaml", "cantilever.yaml"; force=true)
    cp("../examples/overlap/dynamic-same-step/cantilever/cantilever-free.yaml", "cantilever-free.yaml"; force=true)
    cp("../examples/overlap/dynamic-same-step/cantilever/cantilever-clamped.yaml", "cantilever-clamped.yaml"; force=true)
    cp("../examples/overlap/dynamic-same-step/cantilever/cantilever-free.g", "cantilever-free.g"; force=true)
    cp("../examples/overlap/dynamic-same-step/cantilever/cantilever-clamped.g", "cantilever-clamped.g"; force=true)
    sim = Norma.create_simulation("cantilever.yaml")
    subsim = sim.subsims[1]
    model = subsim.model
    # Impose a nonuniform velocity and force unit weights: the per-quadrature
    # re-integration of 0.5 rho |v|^2 must reproduce the consistent-mass kinetic
    # energy 0.5 vᵀ M v exactly (Newmark uses the consistent mass).
    model.displacement .= 0.001 .* model.reference
    model.velocity .= 0.5 .* model.reference .+ 0.2
    model.compute_mass = true
    Norma.evaluate(model, subsim.integrator, subsim.solver)
    force_unit_arlequin_weights!(subsim)
    stored, kinetic = Norma.blended_subdomain_energy(subsim)
    v = vec(model.velocity)
    reference_kinetic = 0.5 * dot(v, model.mass, v)
    @test kinetic ≈ reference_kinetic rtol = 1.0e-10
    @test stored ≈ model.strain_energy rtol = 1.0e-10
    for s in sim.subsims
        Exodus.close(s.params["input_mesh"])
        Exodus.close(s.params["output_mesh"])
    end
    empty!(Norma.ARLEQUIN_WEIGHT_CACHE)
    rm("cantilever.yaml"; force=true)
    rm("cantilever-free.yaml"; force=true)
    rm("cantilever-clamped.yaml"; force=true)
    rm("cantilever-free.g"; force=true)
    rm("cantilever-clamped.g"; force=true)
    rm("cantilever-free.e"; force=true)
    rm("cantilever-clamped.e"; force=true)
end

@testset "Overlap Blended Energy: Degenerate Overlap Rejected" begin
    cp("../examples/overlap/static-same-step/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
    sim = Norma.create_simulation("cuboids.yaml")
    subsim = sim.subsims[1]
    model = subsim.model
    partners = Norma.overlap_partners(subsim)
    @test !isempty(partners)

    # Locate an integration point that lies in the physical overlap.
    blocks = Exodus.read_sets(model.mesh, Exodus.Block)
    overlap_point = nothing
    for (block_index, block) in enumerate(blocks)
        element_type = Norma.element_type_from_string(Exodus.read_block_parameters(model.mesh, block.id)[1])
        num_points = model.num_int_pts[block_index]
        N, _, _ = Norma.isoparametric(element_type, num_points)
        conn = Norma.get_block_connectivity(model.mesh, block.id)
        num_elements, num_nodes = size(conn)
        for element_index in 1:num_elements
            node_indices = conn[((element_index - 1) * num_nodes + 1):(element_index * num_nodes)]
            element_reference = model.reference[:, node_indices]
            for point in 1:num_points
                x = Vector{Float64}(element_reference * N[:, point])
                if any(Norma.point_in_model_reference(x, p.partner_model, p.tol) for p in partners)
                    overlap_point = x
                    break
                end
            end
            overlap_point === nothing || break
        end
        overlap_point === nothing || break
    end
    @test overlap_point !== nothing

    # A legitimate overlap yields a valid partition-of-unity weight.
    @test 0.0 <= Norma.arlequin_weight(overlap_point, model, partners, 1.0, subsim.name) <= 1.0
    # Forcing the characteristic length huge makes the band width fall below the
    # degeneracy tolerance, which must be rejected outright, not regularized.
    Norma.NORMA_TEST_MODE[] = true
    @test_throws Norma.NormaAbortException Norma.arlequin_weight(
        overlap_point, model, partners, 1.0e12, subsim.name
    )
    Norma.NORMA_TEST_MODE[] = false

    for s in sim.subsims
        Exodus.close(s.params["input_mesh"])
        Exodus.close(s.params["output_mesh"])
    end
    empty!(Norma.ARLEQUIN_WEIGHT_CACHE)
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
end

@testset "Overlap Blended Energy: Double Count Removed" begin
    cp("../examples/overlap/static-same-step/cuboids/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-1.g", "cuboid-1.g"; force=true)
    cp("../examples/overlap/static-same-step/cuboids/cuboid-2.g", "cuboid-2.g"; force=true)
    write("cuboids.yaml", read("cuboids.yaml", String) * "\nblended energy output: true\n")
    sim = Norma.run("cuboids.yaml")
    # Naive sum double-counts the shared overlap region.
    naive_stored = sim.subsims[1].model.strain_energy + sim.subsims[2].model.strain_energy

    @test isfile("cuboids-energy.csv")
    rows = readlines("cuboids-energy.csv")
    @test rows[1] == "time,stored_energy,kinetic_energy,total_energy"
    fields = split(rows[end], ",")
    stored = parse(Float64, fields[2])
    kinetic = parse(Float64, fields[3])
    total = parse(Float64, fields[4])
    @test stored > 0.0
    @test kinetic == 0.0                       # static: no kinetic contribution
    @test total ≈ stored rtol = 1.0e-12
    @test stored < naive_stored                # overlap double count removed
    @test stored ≈ 2.496168897869401e8 rtol = 1.0e-3

    empty!(Norma.ARLEQUIN_WEIGHT_CACHE)
    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("cuboid-1.g"; force=true)
    rm("cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("cuboids-energy.csv"; force=true)
end
