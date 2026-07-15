# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# SingleDomainSimulation() does `rm(output_mesh_file; force=true)` and then
# recreates it as a fresh copy of `input mesh file`, before the input mesh is
# ever read. If `input mesh file` and `output mesh file` resolve to the same
# path, that rm deletes the input mesh before the copy that's supposed to
# preserve it -- most dangerous on restart, where the input mesh is often
# the only surviving copy of the checkpoint. This should abort instead of
# silently running (and destroying data) in place.
@testset "In-Place Mesh Abort" begin
    @testset "Non-restart run, same path" begin
        cp("../examples/ahead/single/cuboid/cuboid-hex.g", "cuboid-hex-inplace.g"; force=true)
        original_bytes = read("cuboid-hex-inplace.g")

        doc = YAML.load_file(
            "../examples/ahead/single/cuboid/dynamic/cuboid.yaml"; dicttype=Norma.Parameters
        )
        doc["name"] = "restart-inplace-non-restart"
        doc["input mesh file"] = "cuboid-hex-inplace.g"
        doc["output mesh file"] = "cuboid-hex-inplace.g"

        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException Norma.run(doc)
            # The would-be checkpoint must survive the aborted attempt untouched.
            @test isfile("cuboid-hex-inplace.g")
            @test read("cuboid-hex-inplace.g") == original_bytes
        finally
            Norma.NORMA_TEST_MODE[] = false
            rm("cuboid-hex-inplace.g"; force=true)
        end
    end

    @testset "Restart run, same path (the motivating case)" begin
        cp(
            "../examples/ahead/single/cuboid/dynamic-restart/cuboid-restart-in.e",
            "cuboid-restart-inplace.e";
            force=true,
        )
        original_bytes = read("cuboid-restart-inplace.e")

        doc = YAML.load_file(
            "../examples/ahead/single/cuboid/dynamic-restart/cuboid.yaml"; dicttype=Norma.Parameters
        )
        doc["name"] = "restart-inplace-restart"
        doc["input mesh file"] = "cuboid-restart-inplace.e"
        doc["output mesh file"] = "cuboid-restart-inplace.e"

        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException Norma.run(doc)
            # The one and only copy of the checkpoint must still be there,
            # byte-for-byte, after the aborted attempt.
            @test isfile("cuboid-restart-inplace.e")
            @test read("cuboid-restart-inplace.e") == original_bytes
        finally
            Norma.NORMA_TEST_MODE[] = false
            rm("cuboid-restart-inplace.e"; force=true)
        end
    end

    @testset "Different-looking but equal paths still caught" begin
        cp("../examples/ahead/single/cuboid/cuboid-hex.g", "cuboid-hex-inplace2.g"; force=true)
        original_bytes = read("cuboid-hex-inplace2.g")

        doc = YAML.load_file(
            "../examples/ahead/single/cuboid/dynamic/cuboid.yaml"; dicttype=Norma.Parameters
        )
        doc["name"] = "restart-inplace-relative"
        doc["input mesh file"] = "cuboid-hex-inplace2.g"
        doc["output mesh file"] = "./cuboid-hex-inplace2.g"

        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException Norma.run(doc)
            @test isfile("cuboid-hex-inplace2.g")
            @test read("cuboid-hex-inplace2.g") == original_bytes
        finally
            Norma.NORMA_TEST_MODE[] = false
            rm("cuboid-hex-inplace2.g"; force=true)
        end
    end
end
