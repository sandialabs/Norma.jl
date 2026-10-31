# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Transfer operators on curved interfaces (issue #258): the three material
# boundaries of the quarter of the nested cylinders, each meshed on its two
# sides with different chords of the same circle. The rim facets of the flat
# caps are covered only in part by the projected facets of the other side,
# the sliver between the chords of the two sides. Every facet takes the
# segmented rule, with the remainder integrated on its own, the two
# integration directions of the cross-mass matrix agree to a few percent
# (0.11 on sscore and 0.32 on ssshell when a rim facet took the point rule as
# a whole), and each side's projector reproduces constants exactly.
@testset "Transfer Operators Curved Interface" begin
    S2 = Norma.SVector{2,Float64}
    square = [S2(0.0, 0.0), S2(1.0, 0.0), S2(1.0, 1.0), S2(0.0, 1.0)]
    piece = [S2(0.25, 0.25), S2(0.75, 0.25), S2(0.75, 0.75), S2(0.25, 0.75)]
    rest = Norma.subtract_convex_polygon([square], piece, 1.0e-12, sqrt(2.0))
    @test sum(Norma.polygon_area.(rest)) ≈ 0.75 atol = 1.0e-14
    # A piece along an edge of the polygon, with a repeated vertex as the
    # clipping leaves it: the zero-length edge must clip nothing.
    piece_edge = [S2(0.0, 0.0), S2(0.5, 0.0), S2(0.5, 0.0), S2(0.5, 1.0), S2(0.0, 1.0)]
    rest_edge = Norma.subtract_convex_polygon([square], piece_edge, 1.0e-12, sqrt(2.0))
    @test sum(Norma.polygon_area.(rest_edge)) ≈ 0.5 atol = 1.0e-14
    # The inverse map of a distorted planar quadrilateral near its corners,
    # where the closest-point projection in space converged to a spurious
    # stationary point.
    P = Norma.SMatrix{2,4,Float64,8}([0.0 1.0 1.2 -0.3; 0.0 0.1 1.0 0.8])
    for ξ_true in (S2(0.95, -0.9), S2(-0.98, 0.97), S2(0.3, 0.2))
        p = P * Norma.interpolate(Norma.QUAD4, ξ_true)[1]
        ξ, found = Norma.plane_inverse_map(Norma.QUAD4, P, p, 1.5)
        @test found
        @test ξ ≈ ξ_true atol = 1.0e-10
    end
    T = Norma.SMatrix{2,3,Float64,6}([0.0 2.0 0.5; 0.0 0.3 1.5])
    p = T * Norma.interpolate(Norma.TRI3, S2(0.2, 0.7))[1]
    ξ, found = Norma.plane_inverse_map(Norma.TRI3, T, p, 2.0)
    @test found
    @test ξ ≈ S2(0.2, 0.7) atol = 1.0e-12

    src = "../examples/jmp/concentric-cylinders/quasistatic/quarter/nonoverlap/"
    names = ["core-4", "filler_in-4", "filler_out-4", "shell-4"]
    files = vcat("cylinders-4.yaml", names .* ".yaml", ["core.g", "filler_in.g", "filler_out.g", "shell.g"])
    for file in files
        cp(src * file, file; force=true)
    end
    sim = Norma.create_simulation("cylinders-4.yaml")
    # Largest admissible relative difference between the two integrations of
    # B (its largest entry) and between the nodal areas of one side integrated
    # over the facets of either side (the largest nodal area); the remaining
    # difference is that of the areas of the two chord surfaces.
    bounds = Dict("sscore" => 0.05, "ssint" => 1.0e-4, "ssshell" => 0.02)
    seen = Set{String}()
    for subsim in sim.subsims, bc in subsim.model.boundary_conditions
        bc isa Norma.SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition || continue
        bc.name in seen && continue
        push!(seen, bc.name)
        psim = Norma.coupled_subsim_of(bc)
        pbc = psim.model.boundary_conditions[bc.coupled_bc_index]
        m1, m2 = subsim.model, psim.model
        all_segmented = true
        for (m, b, pm, pb) in ((m1, bc, m2, pbc), (m2, pbc, m1, bc))
            search = Norma.build_side_set_search(pm, pb.side_set_id)
            nd, ns = length(b.local_from_global_map), length(pb.local_from_global_map)
            k = 1
            for nn in b.num_nodes_sides
                nodes = b.side_set_node_indices[k:(k + nn - 1)]
                k += nn
                segmented, _ = Norma.add_facet_segmented_quadrature!(
                    zeros(nd, ns),
                    zeros(nd, nd),
                    get.(Ref(b.local_from_global_map), nodes, 0),
                    m.reference[:, nodes],
                    search,
                    pb.local_from_global_map,
                )
                all_segmented &= segmented
            end
        end
        @test all_segmented
        B1, W1 = Norma.get_rectangular_projection_matrix(m1, bc, m2, pbc)
        B2t, W2 = Norma.get_rectangular_projection_matrix(m2, pbc, m1, bc)
        bound = bounds[bc.name]
        @test maximum(abs.(B1 - transpose(B2t))) / maximum(abs.(B1)) < bound
        @test (W1 \ B1) * ones(size(B1, 2)) ≈ ones(size(B1, 1)) atol = 1.0e-12
        @test (W2 \ B2t) * ones(size(B2t, 2)) ≈ ones(size(B2t, 1)) atol = 1.0e-12
        areas_1 = vec(sum(W1; dims=2))
        areas_2 = vec(sum(W2; dims=2))
        @test maximum(abs.(vec(sum(B1; dims=1)) - areas_2)) / maximum(areas_2) < bound
        @test maximum(abs.(vec(sum(B2t; dims=1)) - areas_1)) / maximum(areas_1) < bound
    end
    @test length(seen) == 3
    for subsim in sim.subsims
        Exodus.close(subsim.params["input_mesh"])
        Exodus.close(subsim.params["output_mesh"])
    end
    for file in vcat(files, names .* ".e")
        rm(file; force=true)
    end
end
