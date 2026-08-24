# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Unit tests for the input keyword validation (input_validation.jl): every
# loaded input file is walked against the known-key sets, and unknown keys
# warn with a Levenshtein "did you mean" suggestion instead of being
# silently ignored.
@testset "Input Validation" begin
    @testset "Levenshtein distance" begin
        @test Norma.levenshtein_distance("", "") == 0
        @test Norma.levenshtein_distance("abc", "abc") == 0
        @test Norma.levenshtein_distance("abc", "") == 3
        @test Norma.levenshtein_distance("kitten", "sitting") == 3
        @test Norma.levenshtein_distance("adjoint paring", "adjoint pairing") == 1
        # Non-ASCII keys must be measured per character, not per byte.
        @test Norma.levenshtein_distance("β", "γ") == 1
        @test Norma.levenshtein_distance("beta", "β") == 4
        @test Norma.levenshtein_distance("HHT alpha", "HHT α") == 5
    end

    @testset "Suggestions" begin
        known = Set(["adjoint pairing", "robin parameter", "impedance scale"])
        @test Norma.suggest_key("adjoint paring", known) == "adjoint pairing"
        @test Norma.suggest_key("Robin Parameter", known) == "robin parameter"
        @test Norma.suggest_key("completely unrelated key name", known) == ""
    end

    @testset "Misspelled BC entry key" begin
        params = Norma.Parameters(
            "type" => "single",
            "boundary conditions" => Norma.Parameters(
                "Schwarz impedance nonoverlap" => [
                    Norma.Parameters(
                        "side set" => "ssz+",
                        "source" => "cuboid-2",
                        "source side set" => "ssz-",
                        "adjoint paring" => false,
                    ),
                ],
            ),
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("adjoint paring", messages[1])
        @test occursin("Did you mean \"adjoint pairing\"?", messages[1])
    end

    @testset "Misspelled BC type key" begin
        params = Norma.Parameters(
            "type" => "single", "boundary conditions" => Norma.Parameters("Schwarz nonoverlap" => [])
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("Did you mean \"Schwarz DN nonoverlap\"?", messages[1])
    end

    @testset "Misspelled controller and integrator keys" begin
        params = Norma.Parameters("type" => "multi", "domains" => ["cuboid-1.yaml"], "relaxation type" => "classical")
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("Did you mean \"relaxation\"?", messages[1])

        params = Norma.Parameters(
            "type" => "single",
            "time integrator" => Norma.Parameters("type" => "Newmark", "beta" => 0.25, "β" => 0.25, "γ" => 0.5),
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("\"beta\"", messages[1])
    end

    @testset "Newmark keys valid only for Newmark" begin
        params = Norma.Parameters(
            "type" => "single",
            "time integrator" =>
                Norma.Parameters("type" => "central difference", "CFL" => 0.5, "γ" => 0.5, "β" => 0.25),
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("\"β\"", messages[1])
    end

    @testset "Unreferenced material property dict" begin
        params = Norma.Parameters(
            "type" => "single",
            "model" => Norma.Parameters(
                "type" => "solid mechanics",
                "material" => Norma.Parameters(
                    "blocks" => Norma.Parameters("fine" => "hyperelastic"),
                    "hyperelastic" => Norma.Parameters("model" => "neohookean", "elastic modulus" => 1.0e9),
                    "stale" => Norma.Parameters("model" => "neohookean"),
                ),
            ),
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("\"stale\"", messages[1])
    end

    @testset "Per-model material keys" begin
        props = Norma.Parameters("model" => "neohookean", "elastic modulus" => 1.0e9, "yield stress" => 1.0e9)
        params = Norma.Parameters(
            "type" => "single",
            "model" => Norma.Parameters(
                "type" => "solid mechanics",
                "material" => Norma.Parameters("blocks" => Norma.Parameters("fine" => "mat"), "mat" => props),
            ),
        )
        messages = Norma.validate_input_parameters(params, "test.yaml")
        @test length(messages) == 1
        @test occursin("\"yield stress\"", messages[1])
        props["model"] = "j2 plasticity"
        @test isempty(Norma.validate_input_parameters(params, "test.yaml"))
    end

    @testset "Shipped examples are clean" begin
        for input_file in (
            "../examples/contact/implicit-dynamic/2-bars/bars.yaml",
            "../examples/contact/implicit-dynamic/2-bars/bar-1.yaml",
            "../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/cuboids.yaml",
            "../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/cuboid-1.yaml",
            "../examples/single/static-solid/cube/plastic/cube.yaml",
        )
            params = YAML.load_file(input_file; dicttype=Norma.Parameters)
            @test isempty(Norma.validate_input_parameters(params, input_file))
        end
    end
end
