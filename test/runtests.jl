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
const indexed_test_files = [
    (1, "minitensor.jl"),
    (2, "interpolation.jl"),
    (3, "constitutive.jl"),
    (4, "single-static-solid-cube.jl"),
    (5, "single-static-solid-neumann-bc.jl"),
    (6, "single-implicit-dynamic-solid-cube.jl"),
    (7, "single-implicit-dynamic-solid-sho.jl"),
    (8, "single-implicit-dynamic-solid-clamped.jl"),
    (9, "single-explicit-dynamic-solid-cube.jl"),
    (10, "single-explicit-dynamic-solid-sho.jl"),
    (11, "single-explicit-dynamic-solid-clamped.jl"),
    (12, "tet4-static-solid-cube.jl"),
    (13, "tet10-static-solid-cube.jl"),
    (14, "schwarz-overlap-static-cuboid-hex8.jl"),
    (15, "schwarz-nonoverlap-static-cuboid-hex8.jl"),
    (16, "transfer-operators.jl"),
    (17, "schwarz-contact-static-cubes.jl"),
    (18, "schwarz-contact-dynamic-cubes.jl"),
    (19, "solid-inclined-displacement.jl"),
    (20, "opinf-schwarz-overlap-cuboid-hex8.jl"),
    (21, "opinf-schwarz-overlap-rom-rom-cuboid-hex8.jl"),
    (22, "quadratic-opinf-schwarz-overlap-cuboid-hex8.jl"),
    (23, "cubic-opinf-schwarz-overlap-cuboid-hex8.jl"),
    (24, "adaptive-time-stepping.jl"),
    (25, "schwarz-ahead-overlap-dynamic-clamped.jl"),
    (26, "schwarz-ahead-overlap-dynamic-notched-cylinder.jl"),
    (26, "schwarz-ahead-overlap-dynamic-laser-weld.jl"),
    (28, "schwarz-ahead-overlap-dynamic-torsion.jl"),
    (29, "schwarz-ahead-overlap-dynamic-bracket.jl"),
    (30, "schwarz-ahead-overlap-dynamic-plate.jl"),
    (31, "single-ahead-clamped.jl"),
    (32, "single-ahead-notched-cylinder.jl"),
    (33, "single-ahead-laser-weld.jl"),
    (34, "single-ahead-torsion.jl"),
    (35, "single-ahead-bracket.jl"),
    (36, "single-ahead-plate.jl"),
    (37, "schwarz-ahead-nonoverlap-dynamic-clamped.jl"),
    (38, "schwarz-ahead-nonoverlap-dynamic-laser-weld.jl"),
    (39, "schwarz-ahead-nonoverlap-dynamic-torsion.jl"),
    (40, "schwarz-ahead-nonoverlap-dynamic-plate.jl"),
    (41, "schwarz-ahead-nonoverlap-dynamic-bracket.jl"),
    (42, "single-static-solid-pressure-bc.jl"), 
    (43, "single-implicit-dynamic-solid-cube-pressure-nbc-stretch.jl"), 
    (44, "single-implicit-dynamic-solid-cube-pressure-nbc-expand.jl"), 
    (45, "single-implicit-dynamic-solid-can-pressure-nbc.jl"), 
    (46, "single-static-solid-cube-sd-dbc.jl"), 
    (47, "constitutive-model-energy-gradient.jl"),
    (48, "smoothing.jl"),
    (49, "utils.jl"), # Must go last due to FPE traps
]

# Extract test file names
const all_test_files = [file for (_, file) in indexed_test_files]

# Optional test indices (excluded from quick runs)
const optional_test_indices = Int[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]

# Quick test set (subset of all tests)
const quick_test_indices = [i for (i, _) in indexed_test_files if i ∉ optional_test_indices]

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

    # Optional name filter
    filter_idx = findfirst(isequal("--filter"), args)
    name_filter = filter_idx !== nothing && filter_idx < length(args) ? lowercase(args[filter_idx + 1]) : ""

    quick_only = "--quick" in args

    # Parse integer indices
    selected_indices = try
        parse.(Int, filter(x -> occursin(r"^\d+$", x), args))
    catch
        Norma.norma_log(0, :error, "Invalid test index provided.")
        exit(1)
    end

    candidate_tests = if !isempty(selected_indices)
        valid_indices = Set(i for (i, _) in indexed_test_files)
        for i in selected_indices
            if i ∉ valid_indices
                Norma.norma_log(0, :error, "Invalid test index: $i")
                exit(1)
            end
        end
        filter(t -> t[1] in selected_indices, indexed_test_files)
    elseif quick_only
        Norma.norma_log(0, :info, "Running quick test set.")
        filter(t -> t[1] in quick_test_indices, indexed_test_files)
    else
        Norma.norma_log(0, :info, "Running full test suite.")
        indexed_test_files
    end

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
