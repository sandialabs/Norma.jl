# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Variational (L2-projection) transfer for the impedance-overlap Schwarz BC on a
# NON-CONFORMING mesh pair (h = 8.47 vs 6.35 mm, 4:3, non-nested — the
# instability-study benchmark). The variational projector P = W \ L (mortar,
# in the domain-decomposition literature) projects the
# partner fields onto this side's trace space; it is non-expansive in L2 for
# any quadrature, reproduces constants exactly (partner partition of unity),
# and reproduces linear fields exactly on the flat interface facets. Those
# two reproduction properties are asserted here to machine precision, along
# with the stage-2a consistent traction: on a non-aligned interface with
# variational transfer, the partner traction is the one-sided weak flux on
# the offset partner facet surface, carried over by the surface-to-surface
# variational operator.

@testset "Schwarz Overlap Impedance Variational Transfer" begin
    example = "../examples/overlap/dynamic-same-step/cantilever-nonconforming"
    for f in ("cantilever-clamped.g", "cantilever-free.g", "cantilever-impedance.yaml")
        cp(joinpath(example, f), f; force=true)
    end
    for f in ("cantilever-clamped-impedance.yaml", "cantilever-free-impedance.yaml")
        y = read(joinpath(example, f), String)
        # Subdivided facet quadrature exercises the piecewise-integrand rule;
        # the reproduction invariants below hold for any rule.
        y = replace(y, r"(source block: [\w+-]+)" => s"\1
      transfer: variational
      transfer quadrature subdivisions: 3
      content aware absorption: true")
        write(f, y)
    end
    multi = read("cantilever-impedance.yaml", String)
    multi = replace(multi, "final time: 3.0e-4" => "final time: 3.0e-6")
    write("cantilever-impedance.yaml", multi)

    sim = Norma.run("cantilever-impedance.yaml")

    for f in (
        "cantilever-impedance.yaml",
        "cantilever-clamped-impedance.yaml",
        "cantilever-free-impedance.yaml",
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
    )
        rm(f; force=true)
    end

    @test sim.failed == false

    for subsim in sim.subsims
        model = subsim.model
        bc = first(
            b for b in model.boundary_conditions if
            b isa Norma.SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
        )
        # Non-conforming + variational transfer: the consistent traction is
        # active through the offset-surface patch (stage 2a), with a
        # surface-to-surface transfer onto this side's boundary nodes.
        patch = bc.traction_patch
        @test patch isa Norma.ConsistentTractionPatch
        @test patch.num_targets > 0
        @test size(patch.transfer) == (length(bc.global_from_local_map), patch.num_targets)
        # Variational projector present with the right shape.
        P = bc.variational_projector
        num_boundary_nodes = length(bc.global_from_local_map)
        coupled = Norma.get_fom_model(Norma.coupled_subsim_of(bc))
        @test size(P) == (num_boundary_nodes, size(coupled.reference, 2))
        # Constants transfer exactly (partner partition of unity).
        @test maximum(abs.(P * ones(size(P, 2)) .- 1.0)) < 1.0e-12
        # Linear fields transfer exactly on the flat interface: the projected
        # partner reference coordinates equal this side's boundary node
        # coordinates, even though the meshes do not align.
        own = model.reference[:, unique(bc.side_set_node_indices)]
        @test maximum(abs.(coupled.reference * P' - own)) < 1.0e-12
        # Content-aware absorption filter: annihilates constants (rigid
        # motions are transferable) and is dissipative in the W inner
        # product per scalar channel.
        F = bc.content_filter
        W = bc.square_projector
        @test size(F) == (num_boundary_nodes, num_boundary_nodes)
        @test maximum(abs.(F * ones(num_boundary_nodes))) < 1.0e-10
        for _ in 1:20
            v = randn(num_boundary_nodes)
            @test dot(v, W * (F * v)) > -1.0e-10
        end
    end
end
