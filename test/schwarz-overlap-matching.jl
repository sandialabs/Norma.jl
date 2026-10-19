# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML
using Exodus
using LinearAlgebra

# The overlap matching coupling (`Schwarz overlap matching`): each side is
# constrained on its own Schwarz boundary to the L2 projection of the partner's
# velocity and receives the partner's d'Alembert reaction on all of its nodes
# through the transposed projector, with the overlap weighted by one half in
# both assemblies. With equal steps the window balance of the rate invariant
# cancels, so Ẽ is conserved to the Schwarz tolerance for the implicit and the
# explicit pair, and the energy with it for identical integrators; the blended
# energy output counts the overlap once with the same weights. Measured on the
# conforming and the 1:0.75 beam overlaps over 50 steps of 1 µs.
const overlap_matching_example = "../examples/overlap/dynamic-same-step/cantilever-matching"
const overlap_conforming_example = "../examples/overlap/dynamic-same-step/cantilever-conforming"

function run_overlap_matching(explicit::Bool, conforming::Bool; num_steps::Int=50)
    meshes = conforming ? overlap_conforming_example : overlap_matching_example
    for name in ("cantilever-free", "cantilever-clamped")
        cp("$meshes/$name.g", "$name.g"; force=true)
        sub = YAML.load_file("$overlap_matching_example/$name.yaml"; dicttype=Norma.Parameters)
        if explicit
            sub["time integrator"] = Norma.Parameters(
                "type" => "central difference", "time step" => 1.0e-6, "CFL" => 1.0, "γ" => 0.5
            )
            sub["solver"] = Norma.Parameters("type" => "explicit solver", "step" => "explicit")
        end
        YAML.write_file("$name.yaml", sub)
    end
    params = YAML.load_file("$overlap_matching_example/cantilever.yaml"; dicttype=Norma.Parameters)
    params["name"] = "overlap-matching"
    params["final time"] = num_steps * 1.0e-6
    sim, out = mktemp() do path, io
        sim = redirect_stdout(io) do
            return Norma.run(params)
        end
        flush(io)
        seekstart(io)
        return sim, read(io, String)
    end
    rows = readlines("overlap-matching-energy.csv")
    header = split(rows[1], ",")
    data = [parse.(Float64, split(row, ",")) for row in rows[2:end]]
    for f in (
        "cantilever-free.g",
        "cantilever-clamped.g",
        "cantilever-free.yaml",
        "cantilever-clamped.yaml",
        "cantilever-free.e",
        "cantilever-clamped.e",
        "overlap-matching-energy.csv",
    )
        rm(f; force=true)
    end
    return sim, out, header, data
end

function overlap_matching_bc_of(subsim)
    return only(
        filter(b -> b isa Norma.SolidMechanicsOverlapMatchingSchwarzBoundaryCondition, subsim.model.boundary_conditions)
    )
end

@testset "Schwarz Overlap Matching: Conservation (II and EE, conforming and 1:0.75)" begin
    for explicit in (false, true), conforming in (true, false)
        sim, out, header, data = run_overlap_matching(explicit, conforming)
        @test !sim.failed
        @test !occursin("Overlap transfer", out)
        @test occursin("Overlap matching interface", out)
        # Transfer operators: both projections reproduce constants and linear
        # fields, and each force transfer is the transpose of the partner's
        # projector.
        for subsim in sim.subsims
            bc = overlap_matching_bc_of(subsim)
            partner_sim = Norma.coupled_subsim_of(bc)
            P = bc.dirichlet_projector
            X = subsim.model.reference[:, bc.global_from_local_map]
            @test maximum(abs.(P * ones(size(P, 2)) .- 1.0)) ≤ 1.0e-12
            @test maximum(abs.(X' - P * partner_sim.model.reference')) ≤ 1.0e-12 * maximum(abs.(X))
            @test bc.neumann_projector == transpose(overlap_matching_bc_of(partner_sim).dirichlet_projector)
        end
        # Conservation over the run, relative to the initial values.
        i_E = findfirst(==("total_energy"), header)
        i_R = findfirst(==("rate_invariant_total"), header)
        E0 = data[1][i_E]
        R0 = data[1][i_R]
        δE = maximum(abs(row[i_E] / E0 - 1.0) for row in data)
        δR = maximum(abs(row[i_R] / R0 - 1.0) for row in data)
        @test δE ≤ 1.0e-11
        @test δR ≤ 1.0e-11
        # The iteration converges in a moderate number of iterations per window.
        iterations = [parse(Int, m.captures[1]) for m in eachmatch(r"Performed (\d+) Schwarz Iteration", out)]
        @test length(iterations) == 50
        @test maximum(iterations) ≤ 40
    end
end

# The weights of the overlap: the blended energy at t = 0 counts the overlap
# once, so the total stored energy equals the strain energy of the undecomposed
# initial deflection of the beam, which the nonoverlap matching exchange gives
# as the sum of its two parts.
@testset "Schwarz Overlap Matching: Blended Energy at t = 0" begin
    _, _, header, data = run_overlap_matching(false, true; num_steps=1)
    i_S = findfirst(==("stored_energy"), header)
    # Strain energy of the beam's initial deflection u_y = 0.393701 x², linear
    # elastic, E = 6.895e9, ν = 0.25, 254 × 25.4 × 12.7 mm (the undecomposed
    # mesh of the conforming example gives 1506.19 J).
    @test data[1][i_S] ≈ 1506.19 rtol = 2.0e-3
end
