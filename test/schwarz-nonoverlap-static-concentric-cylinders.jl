# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Quasi-static Dirichlet-Neumann Schwarz with Anderson acceleration on the
# quarter of the three nested cylinders (issue #246), in two subdomains (inner,
# outer) and in a chain of three (inner, filler_out, shell), against the
# monolithic run. Each unsupported subdomain is the Dirichlet side toward the
# supported shell. The chain of four that also separates the core from the
# filler diverges with Dirichlet-Neumann coupling (cylinders-3.yaml records
# why) and converges with Robin-Robin coupling on all interfaces
# (cylinders-4.yaml); neither is tested.
function mean_core_displacement_z(model)
    nodes = unique(vec(model.blocks[1].connectivity))
    return sum(model.displacement[3, nodes]) / length(nodes)
end

@testset "Schwarz Nonoverlap Static Concentric Cylinders" begin
    src = "../examples/nonoverlap/static-same-step/concentric-cylinders/quarter/"
    mono_files = ["cylinders.yaml", "cylinders.g"]
    for file in mono_files
        cp(src * "monolithic/" * file, file; force=true)
    end
    mono = Norma.run("cylinders.yaml")
    for file in vcat(mono_files, ["cylinders.e"])
        rm(file; force=true)
    end
    @test mono.failed == false
    reference = mean_core_displacement_z(mono.model)
    @test reference ≈ -5.797e-4 rtol = 1.0e-3

    cases = [
        ("cylinders-2.yaml", ["inner-2", "outer-2"], ["inner", "outer"], 24),
        ("cylinders-3.yaml", ["inner-3", "filler_out-3", "shell-3"], ["inner", "filler_out", "shell"], 27),
    ]
    for (input, domains, meshes, max_iters) in cases
        files = vcat([input], ["$d.yaml" for d in domains], ["$m.g" for m in meshes])
        for file in files
            cp(src * "nonoverlap/" * file, file; force=true)
        end
        sim = Norma.run(input)
        for file in vcat(files, ["$d.e" for d in domains])
            rm(file; force=true)
        end
        @test sim.failed == false
        @test sim.controller.converged == true
        @test all(0 .< sim.controller.schwarz_iters .<= max_iters)
        # The meshes are nonconforming and coarser than the monolithic one away
        # from the region: 0.08% (two) and 0.33% (three) at refinement 0.5.
        @test mean_core_displacement_z(sim.subsims[1].model) ≈ reference rtol = 5.0e-3
    end
end
