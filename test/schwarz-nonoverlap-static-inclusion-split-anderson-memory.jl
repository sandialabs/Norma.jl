# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Memory of the Anderson histories against the depth m, on the inclusion with
# its interface split into four side sets (as in the Inclusion Split Anderson
# test). Each Anderson step stores one column per interface: the iterate and
# the residual of every relaxed field on the footprint of the interface, the
# partner degrees of freedom that its transfer reads, and the residual trace
# once. Both subdomains are quasi-static, so the displacement is the only
# relaxed field, and a column of the four interfaces holds
#   8 × Σ (2 × footprint + trace) bytes = 8 × (2 × 114 + 258) = 3888 bytes,
# against 401.5 KB when the histories held the whole partner field for all
# three kinematic fields. The history keeps at most m + 1 columns.
# The peak of a stop is therefore
#   peak = bytes per column × min(m + 1, columns stored),
# and the columns stored are the Schwarz iterations of the stop, one step per
# iteration, unless the condition limit drops some. With m = 5, 10, 20 the
# stop needs more iterations than m + 1 (27, 25, 23 when this test was
# written), the window fills, and the peak equals the full-window estimate
# bytes per column × (m + 1): it grows linearly in m.
using Random
using Test

@testset "Schwarz Nonoverlap Static Inclusion Split Anderson Memory" begin
    src = "../examples/nonoverlap/static-same-step/matrix-inclusion-split/"
    files = ["inclusion.yaml", "disk.yaml", "matrix.yaml", "disk.g", "matrix.g"]
    depths = [5, 10, 20]
    peaks = Dict{Int,Int}()
    iterations = Dict{Int,Int}()
    expected_bytes_per_column = 0
    for m in depths
        for file in files
            cp(src * file, file; force=true)
        end
        input = read("inclusion.yaml", String)
        @test occursin(r"anderson depth: *\d+", input)
        write("inclusion.yaml", replace(input, r"anderson depth: *\d+" => "anderson depth: $m"))
        sim = Norma.run("inclusion.yaml")
        for file in vcat(files, ["disk.e", "matrix.e"])
            rm(file; force=true)
        end
        @test sim.failed == false
        @test sim.controller.converged == true
        @test sim.controller.anderson_depth == m
        iterations[m] = sim.controller.schwarz_iters[1]
        peaks[m] = sim.controller.anderson_peak_bytes[1]
        @test 0 < iterations[m] <= 30
        @test peaks[m] > 0
        # The column of the relaxed interfaces from their footprints and traces.
        footprint_entries, trace_entries = 0, 0
        for subsim in sim.subsims, bc in subsim.model.boundary_conditions
            bc isa Norma.SolidMechanicsNonOverlapSchwarzBoundaryCondition && bc.is_dirichlet || continue
            @test Norma.relaxes_displacement_only(bc)
            footprint_entries += length(Norma.anderson_footprint(bc))
            trace_entries += 3 * size(bc.dirichlet_projector, 1)
        end
        @test (footprint_entries, trace_entries) == (114, 258)
        expected_bytes_per_column = 8 * (2 * footprint_entries + trace_entries)
    end

    # Bytes per column, from the smallest depth, whose window fills first,
    # against the count from the footprints and traces.
    m0 = first(depths)
    columns0 = min(m0 + 1, iterations[m0])
    @test peaks[m0] % columns0 == 0
    bytes_per_column = peaks[m0] ÷ columns0
    @test bytes_per_column == expected_bytes_per_column == 3888

    for m in depths
        full_window = bytes_per_column * (m + 1)
        # The window filled: the stop took at least m + 1 Anderson steps.
        @test iterations[m] >= m + 1
        # The measured peak never exceeds the full-window estimate, and here
        # equals it exactly (no column dropped by the condition limit).
        @test peaks[m] <= full_window
        @test peaks[m] == bytes_per_column * min(m + 1, iterations[m])
        @test peaks[m] == full_window
    end

    # Linear in m: equal increments per added column across the depths.
    for (a, b) in zip(depths[1:(end - 1)], depths[2:end])
        @test peaks[b] - peaks[a] == bytes_per_column * (b - a)
    end
end

# The step on the footprint gives, on the footprint, the iterates of the step
# on the whole field, when the data of the footprint depend on the iterates of
# the footprint alone, as they do for the exchange of an interface: the
# coefficients come from the traces, and every other operation is entrywise.
@testset "Anderson Step On The Footprint" begin
    rng = Random.MersenneTwister(0)
    n, footprint = 40, [3, 4, 5, 17, 18, 19, 31, 32, 33]
    B = 0.3 * randn(rng, length(footprint), length(footprint))
    c = randn(rng, length(footprint))
    trace_map = randn(rng, 6, length(footprint))
    # A map of the whole field whose footprint entries depend on the footprint
    # entries alone and whose other entries are arbitrary.
    function partner_data(x)
        g = randn(rng, n)
        g[footprint] = c .+ B * x[footprint]
        return g
    end
    full, restricted = Norma.AndersonHistory(), Norma.AndersonHistory()
    x = zeros(n)
    x_r = x[footprint]
    for _ in 1:12
        g = partner_data(x)
        t = trace_map * (g[footprint] .- x[footprint])
        x = Norma.anderson_step!(full, x, g, t, 0.5, 4)
        x_r = Norma.anderson_step!(restricted, x_r, g[footprint], t, 0.5, 4)
        @test x_r == x[footprint]
        @test length(full.T) <= 5
    end
    @test all(length(columns) <= 5 for columns in vcat(restricted.X, restricted.F, [restricted.T]))
end
