# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML
using Exodus

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

# The weak overlap transfer, W⁻¹L with L the cross-mass matrix of the Schwarz
# boundary against the partner's volume shape functions, is integrated by the
# segmented rule: each facet is cut by the sections of the partner's elements.
# On the beam overlaps every facet is covered, Π reproduces constants and
# linear fields to roundoff, and on the 1:0.75 overlap Π equals an
# independent exact integration (each facet cut at the partner's node planes
# in y and z, 2 × 2 Gauss points per piece, exact for the piecewise bilinear
# integrand) to roundoff; the facet's own rule differed from it by 8.8e-2.
@testset "Schwarz Overlap Weak DBC: Segmented Transfer" begin
    overlap_bc_of(subsim) =
        only(filter(b -> b isa Norma.SolidMechanicsOverlapSchwarzBoundaryCondition, subsim.model.boundary_conditions))
    function cut_reference(model, bc, partner)
        block_id = Norma.block_id_from_name(bc.coupled_block_name, partner.mesh)
        element_type = only(partner.blocks).element_type
        L = zeros(length(bc.global_from_local_map), size(partner.reference, 2))
        g = [-1.0, 1.0] ./ sqrt(3.0)
        ycuts = unique(round.(partner.reference[2, :]; digits=12))
        zcuts = unique(round.(partner.reference[3, :]; digits=12))
        k = 1
        for ns in bc.num_nodes_sides
            facet_nodes = bc.side_set_node_indices[k:(k + ns - 1)]
            k += ns
            rows = get.(Ref(bc.local_from_global_map), facet_nodes, 0)
            X = model.reference[:, facet_nodes]
            y0, y1 = extrema(X[2, :])
            z0, z1 = extrema(X[3, :])
            ys = sort(unique(vcat(y0, y1, filter(y -> y0 + 1.0e-9 < y < y1 - 1.0e-9, ycuts))))
            zs = sort(unique(vcat(z0, z1, filter(z -> z0 + 1.0e-9 < z < z1 - 1.0e-9, zcuts))))
            for a in 1:(length(ys) - 1), b in 1:(length(zs) - 1), i in 1:2, j in 1:2
                y = (ys[a] + ys[a + 1]) / 2 + (ys[a + 1] - ys[a]) / 2 * g[i]
                z = (zs[b] + zs[b + 1]) / 2 + (zs[b + 1] - zs[b]) / 2 * g[j]
                weight = (ys[a + 1] - ys[a]) * (zs[b + 1] - zs[b]) / 4
                x = [X[1, 1], y, z]
                _, ξ, _, _ = Norma.closest_point_projection(X, x; strict=false)
                N = Vector(Norma.interpolate(Norma.QUAD4, ξ)[1])
                cols, η, found = Norma.find_point_in_mesh(x, partner, block_id, bc.search_tolerance)
                @test found
                Ns = Vector(Norma.interpolate(element_type, η)[1])
                L[rows, cols] .+= N * Ns' .* weight
            end
        end
        return Norma.get_square_projection_matrix(model, bc) \ L
    end
    for example in ("cantilever-conforming", "cantilever-nonconforming")
        dir = "../examples/overlap/dynamic-same-step/$example"
        for name in ("cantilever-free", "cantilever-clamped")
            cp("$dir/$name.g", "$name.g"; force=true)
            sub = YAML.load_file("$dir/$name.yaml"; dicttype=Norma.Parameters)
            sub["boundary conditions"]["Schwarz overlap"][1]["weak"] = true
            YAML.write_file("$name.yaml", sub)
        end
        params = YAML.load_file("$dir/cantilever.yaml"; dicttype=Norma.Parameters)
        params["name"] = "cantilever"
        sim, out = mktemp() do path, io
            sim = redirect_stdout(io) do
                sim = Norma.create_simulation(params)
                Norma.sync_control_time(sim)
                Norma.initialize(sim)
                sim
            end
            flush(io)
            seekstart(io)
            sim, replace(read(io, String), r"\s+" => " ")
        end
        @test !occursin("Overlap transfer", out)
        for subsim in sim.subsims
            bc = overlap_bc_of(subsim)
            partner = Norma.coupled_subsim_of(bc).model
            P = bc.dirichlet_projector
            X = subsim.model.reference[:, bc.global_from_local_map]
            @test maximum(abs.(P * ones(size(P, 2)) .- 1.0)) ≤ 1.0e-12
            @test maximum(abs.(X' - P * partner.reference')) ≤ 1.0e-12 * maximum(abs.(X))
            if example == "cantilever-nonconforming"
                @test maximum(abs.(P - cut_reference(subsim.model, bc, partner))) ≤ 1.0e-12
            end
        end
        for subsim in sim.subsims
            Exodus.close(subsim.params["input_mesh"])
            Exodus.close(subsim.params["output_mesh"])
        end
        foreach(
            f -> rm(f; force=true),
            [
                "cantilever-free.g",
                "cantilever-clamped.g",
                "cantilever-free.yaml",
                "cantilever-clamped.yaml",
                "cantilever-free.e",
                "cantilever-clamped.e",
            ],
        )
    end
end
