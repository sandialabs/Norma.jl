# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using Test
using Logging

@testset "Utils Tests" begin
    @testset "Format Time" begin
        @test Norma.format_time(0.0) == "0.00s"
        @test Norma.format_time(59.9) == "59.90s"
        @test Norma.format_time(60.0) == "1m 0.00s"
        @test Norma.format_time(3600.0) == "1h 0m 0.00s"
        @test Norma.format_time(86400.0) == "1d 0h 0m 0.00s"
        @test Norma.format_time(90061.5) == "1d 1h 1m 1.50s"
    end

    @testset "Configure Logger" begin
        # Save current ENV
        original_debug = get(ENV, "JULIA_DEBUG", nothing)

        try
            ENV["JULIA_DEBUG"] = "Norma"
            Norma.configure_logger()
            @test true  # Just ensure no error thrown

            ENV["JULIA_DEBUG"] = ""
            Norma.configure_logger()
            @test true
        finally
            if original_debug === nothing
                pop!(ENV, "JULIA_DEBUG", nothing)
            else
                ENV["JULIA_DEBUG"] = original_debug
            end
        end
    end

    @testset "norma_abort logging" begin
        Norma.norma_log(0, :info, "Mock code abort ...")
        @test Norma._norma_abort_message("Testing aborting") === nothing
        @test Norma._norma_abort_messagef("Testing failure code = %d", 42) === nothing
    end

    @testset "Parse Args" begin
        let saved_args = copy(ARGS)
            try
                # Success case
                empty!(ARGS)
                push!(ARGS, "input.yaml")
                @test Norma.parse_args() == "input.yaml"

                # Failure case
                empty!(ARGS)
                Norma.NORMA_TEST_MODE[] = true
                @test_throws Norma.NormaAbortException Norma.parse_args()
            finally
                empty!(ARGS)
                append!(ARGS, saved_args)
                Norma.NORMA_TEST_MODE[] = false
            end
        end
    end

    @testset "Thread Configuration" begin
        using LinearAlgebra: BLAS
        saved_blas_threads = BLAS.get_num_threads()
        saved_env = [(var, get(ENV, var, nothing)) for var in Norma.NORMA_BLAS_THREAD_VARS]
        try
            for (var, _) in saved_env
                pop!(ENV, var, nothing)
            end

            # Julia's `-t` sets the task pool; BLAS must follow it, not the
            # hardware thread count.
            BLAS.set_num_threads(saved_blas_threads + 1)
            Norma.configure_threads()
            @test BLAS.get_num_threads() == Base.Threads.threadpoolsize(:default)

            # An explicit request in the environment wins.
            for (var, _) in saved_env
                ENV[var] = "1"
                BLAS.set_num_threads(2)
                Norma.configure_threads()
                @test BLAS.get_num_threads() == 2
                pop!(ENV, var, nothing)
            end

            @test occursin("BLAS", Norma.thread_report())
        finally
            for (var, value) in saved_env
                value === nothing ? pop!(ENV, var, nothing) : (ENV[var] = value)
            end
            BLAS.set_num_threads(saved_blas_threads)
        end
    end

end
