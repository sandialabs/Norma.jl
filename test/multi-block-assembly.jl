# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Assembly on meshes of many small blocks (issue #251): one inverse map for all
# blocks, whose gather adds the sources of every target in block and element
# order and so equals the scatter bit for bit, and element chunks small enough
# that a block of a few hundred elements is divided among the threads.
using Random
using Test

@testset "Multi-Block Assembly" begin
    rng = MersenneTwister(251)
    num_targets = 1000
    targets = [rand(rng, 1:num_targets, n) for n in (37, 590, 1, 0, 212)]
    values = [randn(rng, length(t)) for t in targets]
    ptr, block, source = Norma.inverse_map(targets, num_targets)
    # One map for all blocks: a single pointer array over the targets.
    @test length(ptr) == num_targets + 1
    @test length(block) == length(source) == sum(length, targets)
    scattered = zeros(num_targets)
    for (b, block_targets) in enumerate(targets), (k, t) in enumerate(block_targets)
        scattered[t] += values[b][k]
    end
    gathered = zeros(num_targets)
    for t in 1:num_targets, p in ptr[t]:(ptr[t + 1] - 1)
        gathered[t] += values[block[p]][source[p]]
    end
    @test gathered == scattered
    # Element chunks: at least 32 per block down to 16 elements, at most 256,
    # from the number of elements alone.
    @test Norma.element_chunk(590) == 19
    @test cld(590, Norma.element_chunk(590)) >= 32
    @test Norma.element_chunk(100) == 16
    @test Norma.element_chunk(110_592) == Norma.CHUNK_ELEMENTS
end
