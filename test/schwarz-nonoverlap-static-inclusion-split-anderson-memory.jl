# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Memory of the Anderson histories against the depth m, on the inclusion with
# its interface split into four side sets (as in the Inclusion Split Anderson
# test). Each Anderson step stores one column (iterate, residual, and trace of
# every relaxed interface and field), all of the same size, and the history
# keeps at most m + 1 columns. The peak of a stop is therefore
#   peak = bytes per column × min(m + 1, columns stored),
# and the columns stored are the Schwarz iterations of the stop, one step per
# iteration, unless the condition limit drops some. With m = 5, 10, 20 the
# stop needs more iterations than m + 1 (27, 25, 23 when this test was
# written), the window fills, and the peak equals the full-window estimate
# bytes per column × (m + 1): it grows linearly in m.
@testset "Schwarz Nonoverlap Static Inclusion Split Anderson Memory" begin
    src = "../examples/nonoverlap/static-same-step/matrix-inclusion-split/"
    files = ["inclusion.yaml", "disk.yaml", "matrix.yaml", "disk.g", "matrix.g"]
    depths = [5, 10, 20]
    peaks = Dict{Int,Int}()
    iterations = Dict{Int,Int}()
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
    end

    # Bytes per column, from the smallest depth, whose window fills first.
    m0 = first(depths)
    columns0 = min(m0 + 1, iterations[m0])
    @test peaks[m0] % columns0 == 0
    bytes_per_column = peaks[m0] ÷ columns0
    @test bytes_per_column % 8 == 0

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
