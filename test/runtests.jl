# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using LinearAlgebra
using Printf
using Statistics
using Test

include("../src/Norma.jl")
include("helpers.jl")

# List of all test files (ordered)
const all_test_files = [
    "minitensor.jl", #1
    "interpolation.jl", #2
    "constitutive.jl", #3
    "single-static-solid-cube.jl", #4
    "single-static-solid-neumann-bc.jl", #5
    "single-implicit-dynamic-solid-cube.jl", #6
    "single-implicit-dynamic-solid-sho.jl", #7
    "single-implicit-dynamic-solid-clamped.jl", #8
    "single-explicit-dynamic-solid-cube.jl", #9
    "single-explicit-dynamic-solid-sho.jl", #10
    "single-explicit-dynamic-solid-clamped.jl", #11
    "tet4-static-solid-cube.jl", #12
    "tet10-static-solid-cube.jl", #13
    "schwarz-overlap-static-cuboid-hex8.jl", #14
    "schwarz-nonoverlap-static-cuboid-hex8.jl", #15
    "transfer-operators.jl", #16
    "schwarz-contact-static-cubes.jl", #17
    "schwarz-contact-dynamic-cubes.jl", #18
    "solid-inclined-displacement.jl", #19
    "opinf-schwarz-overlap-cuboid-hex8.jl", #20
    "quadratic-opinf-schwarz-overlap-cuboid-hex8.jl", #21
    "cubic-opinf-schwarz-overlap-cuboid-hex8.jl", #22
    "adaptive-time-stepping.jl", #23
    "schwarz-ahead-overlap-dynamic-notched-cylinder.jl", #24
    "schwarz-ahead-overlap-dynamic-laser-weld.jl", #25
    "schwarz-ahead-overlap-dynamic-torsion.jl", #26
    "schwarz-ahead-overlap-dynamic-bracket.jl", #27
    "schwarz-ahead-overlap-dynamic-plate.jl", #28
    "ahead-single-notched-cylinder.jl", #29
    "ahead-single-laser-weld.jl", #30
    "ahead-single-torsion.jl", #31
    "ahead-single-bracket.jl", #32
    "ahead-single-plate.jl", #33
    "schwarz-ahead-nonoverlap-dynamic-laser-weld.jl", #34
    "schwarz-ahead-nonoverlap-dynamic-torsion.jl", #35
    # Must go last due to FPE traps
    "utils.jl", #36
]

# Enumerated test files
const indexed_test_files = collect(enumerate(all_test_files))

# Optional test indices
const optional_test_indices = Int[]

# Default test indices: all except optional ones
const default_test_indices = [i for (i, _) in indexed_test_files if i ∉ optional_test_indices]

function print_available_tests()
    Norma.norma_log(0, :info, "Available tests:")
    for (i, file) in indexed_test_files
        Norma.norma_log(0, :info, rpad("[$i]", 5) * file)
    end
end

function parse_args(args)
    if "--list" in args
        print_available_tests()
        exit(0)
    end

    # Extract optional filter string
    filter_idx = findfirst(isequal("--filter"), args)
    name_filter = filter_idx !== nothing && filter_idx < length(args) ? lowercase(args[filter_idx + 1]) : ""

    run_all = "--all" in args

    # Extract indices
    selected_indices = try
        parse.(Int, filter(x -> occursin(r"^\d+$", x), args))
    catch
        Norma.norma_log(0, :error, "Invalid test index provided.")
        exit(1)
    end

    # Determine candidate set
    candidate_tests = if run_all
        indexed_test_files
    elseif !isempty(selected_indices)
        valid_indices = Set(1:length(all_test_files))
        for i in selected_indices
            if i ∉ valid_indices
                Norma.norma_log(0, :error, "Invalid test index: $i")
                exit(1)
            end
        end
        filter(t -> t[1] in selected_indices, indexed_test_files)
    else
        Norma.norma_log(0, :info, "No tests specified. Running default set.")
        filter(t -> t[1] in default_test_indices, indexed_test_files)
    end

    # Apply filter if present
    if !isempty(name_filter)
        candidate_tests = filter(t -> occursin(name_filter, lowercase(t[2])), candidate_tests)
        if isempty(candidate_tests)
            Norma.norma_log(0, :warn, "No tests match filter \"$name_filter\".")
            exit(1)
        end
    end

    return candidate_tests
end

# Get test list to run
test_files_to_run = parse_args(ARGS)

# Run test suite
start_time = time()
Norma.norma_log(0, :norma, "BEGIN TESTS")

@testset verbose = true "Norma.jl Test Suite" begin
    for (i, file) in test_files_to_run
        Norma.norma_log(0, :test, "[$i] Running $file...")
        include(file)
    end
end

elapsed_time = time() - start_time
Norma.norma_log(0, :done, "Tests Complete")
Norma.norma_log(0, :time, "Tests Run Time = " * Norma.format_time(elapsed_time))
Norma.norma_log(0, :norma, "END TESTS")

# Cleanup artifacts
for ext in ["yaml", "e", "g", "csv"]
    for file in filter(f -> endswith(f, ".$ext"), readdir())
        rm(file; force=true)
    end
end
