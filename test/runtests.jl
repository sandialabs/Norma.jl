# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using LinearAlgebra
using Printf
using Statistics
using Test
using Norma

include("helpers.jl")

# Capture the whole suite in a single log file (test/runtests.log) using Norma's
# own logging system, instead of the many per-input .log files.  Every norma_log
# message is written to it as sanitized ASCII (color codes stripped), matching
# what an individual run writes to its own .log file.
const runtests_log_path = joinpath(@__DIR__, "runtests.log")
Norma.open_session_log_file(runtests_log_path)

# Recursively count passes/fails/errors in a testset, descending into nested
# testsets that a test file may define.  Used to record a per-file result line
# in the log, since the Test stdlib's own summary is not routed through norma_log.
function count_test_results(ts::Test.AbstractTestSet)
    passes = fails = errors = 0
    if isdefined(ts, :n_passed)
        passes += ts.n_passed
    end
    for r in ts.results
        if r isa Test.Fail
            fails += 1
        elseif r isa Test.Error
            errors += 1
        elseif r isa Test.Pass
            passes += 1
        elseif r isa Test.AbstractTestSet
            p, f, e = count_test_results(r)
            passes += p
            fails += f
            errors += e
        end
    end
    return passes, fails, errors
end

# List of all test files (ordered)
const indexed_test_files = [
    (1,  "minitensor.jl"),
    (2,  "interpolation.jl"),
    (3,  "constitutive.jl"),
    (4,  "single-static-solid-cube.jl"),
    (5,  "single-static-solid-neumann-bc.jl"),
    (6,  "single-static-solid-robin-bc.jl"),
    (7,  "single-implicit-dynamic-solid-cube.jl"),
    (8,  "single-implicit-dynamic-solid-sho.jl"),
    (9,  "single-implicit-dynamic-solid-clamped.jl"),
    (10, "single-explicit-dynamic-solid-cube.jl"),
    (11, "single-explicit-dynamic-solid-sho.jl"),
    (12, "single-explicit-dynamic-solid-clamped.jl"),
    (13, "tet4-static-solid-cube.jl"),
    (14, "tet10-static-solid-cube.jl"),
    (15, "schwarz-overlap-static-cuboid-hex8.jl"),
    (16, "schwarz-nonoverlap-static-cuboid-hex8.jl"),
    (17, "schwarz-nonoverlap-static-cuboid-robin-robin.jl"),
    (18, "schwarz-nonoverlap-dynamic-cuboid.jl"),
    (19, "transfer-operators.jl"),
    (20, "schwarz-contact-static-cubes.jl"),
    (21, "schwarz-contact-dynamic-cubes.jl"),
    (22, "opinf-schwarz-overlap-l2-error-cuboid-hex8.jl"),
    (23, "opinf-schwarz-overlap-cuboid-hex8.jl"),
    (24, "opinf-schwarz-overlap-rom-rom-cuboid-hex8.jl"),
    (25, "quadratic-opinf-schwarz-overlap-cuboid-hex8.jl"),
    (26, "cubic-opinf-schwarz-overlap-cuboid-hex8.jl"),
    (27, "adaptive-time-stepping.jl"),
    (28, "schwarz-ahead-overlap-dynamic-clamped.jl"),
    (29, "schwarz-ahead-overlap-dynamic-notched-cylinder.jl"),
    (29, "schwarz-ahead-overlap-dynamic-laser-weld.jl"),
    (30, "schwarz-ahead-overlap-dynamic-torsion.jl"),
    (31, "schwarz-ahead-overlap-dynamic-bracket.jl"),
    (32, "schwarz-ahead-overlap-dynamic-plate.jl"),
    (33, "single-ahead-clamped.jl"),
    (34, "single-ahead-notched-cylinder.jl"),
    (35, "single-ahead-laser-weld.jl"),
    (36, "single-ahead-torsion.jl"),
    (37, "single-ahead-bracket.jl"),
    (38, "single-ahead-plate.jl"),
    (39, "schwarz-ahead-nonoverlap-dynamic-clamped.jl"),
    (40, "schwarz-ahead-nonoverlap-dynamic-laser-weld.jl"),
    (41, "schwarz-ahead-nonoverlap-dynamic-torsion.jl"),
    (42, "schwarz-ahead-nonoverlap-dynamic-plate.jl"),
    (43, "schwarz-ahead-nonoverlap-dynamic-bracket.jl"),
    (44, "single-static-solid-pressure-bc.jl"),
    (45, "single-implicit-dynamic-solid-cube-pressure-nbc-stretch.jl"),
    (46, "single-implicit-dynamic-solid-cube-pressure-nbc-expand.jl"),
    (47, "single-implicit-dynamic-solid-can-pressure-nbc.jl"),
    (48, "single-static-solid-cube-sd-dbc.jl"),
    (49, "constitutive-model-energy-gradient.jl"),
    (50, "smoothing.jl"),
    (51, "nnopinf-schwarz-overlap-cuboid-hex8.jl"),
    (52, "single-static-solid-j2-plasticity.jl"),
    (53, "utils.jl"),
    (54, "schwarz-nonoverlap-dynamic-cantilever-dn.jl"),
    (55, "schwarz-overlap-dynamic-cantilever-weak.jl"),
    (56, "schwarz-overlap-static-cuboid-hex8-swap.jl"),
    (57, "schwarz-nonoverlap-static-cuboid-hex8-aitken.jl"),
    (58, "schwarz-overlap-dynamic-cantilever-impedance.jl"),
    (59, "single-static-solid-cube-time-swap.jl"),
    (60, "linear-krom-schwarz-overlap-cuboid-hex8.jl"),
    (61, "rbf-krom-schwarz-overlap-cuboid-hex8.jl"),
    (62, "krom-schwarz-overlap-rom-rom-cuboid-hex8.jl"),
    (63, "recovery-cube.jl"),
    (64, "transfer-volumetric.jl"),
    (65, "single-static-solid-cube-time-swap-mesh-change.jl"),
    (66, "schwarz-ahead-nonoverlap-dynamic-cuboid.jl"),
    (67, "schwarz-ahead-nonoverlap-dynamic-cuboid-fom-rom.jl"),
    (68, "schwarz-ahead-nonoverlap-dynamic-cuboid-rom-fom.jl"),
    (69, "schwarz-ahead-nonoverlap-dynamic-cuboid-rom-rom.jl"),
    (70, "schwarz-ahead-nonoverlap-dynamic-cuboid-quadratic-rom-fom.jl"),
    (71, "schwarz-ahead-nonoverlap-dynamic-cuboid-cubic-rom-fom.jl"),
    (72, "schwarz-ahead-nonoverlap-dynamic-cuboid-cubic-rom-rom.jl"),
    (73, "single-static-solid-cube-stress-recovery-swap.jl"),
    (74, "schwarz-overlap-static-cuboid-hex8-overlap-l2-swap.jl"),
    (75, "quadratic-opinf-central-difference-schwarz-overlap-cuboid-hex8.jl"),
    (76, "cubic-opinf-central-difference-schwarz-overlap-cuboid-hex8.jl"),
    (77, "single-static-solid-notched-cylinder-j2-elastic-to-plastic-swap.jl"),
    (78, "single-dynamic-tension-specimen-j2-elastic-to-plastic-swap.jl"),
    (79, "single-dynamic-tension-specimen-j2-elastic-to-plastic-swap-long.jl"),
    (80, "schwarz-ahead-overlap-dynamic-notched-cylinder-j2-swap-fom.jl"),
    (81, "single-implicit-dynamic-solid-clamped-bc.jl"),
    (82, "single-explicit-dynamic-solid-clamped-bc.jl"),
    (83, "single-dynamic-opinf-rom-to-fom-time-swap.jl"),
    (84, "schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap-middle.jl"),
    (85, "schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-fom-multi-swap.jl"),
    (86, "schwarz-ahead-overlap-dynamic-notched-cylinder-opinf-swap-rom-fom.jl"),
    (87, "single-static-solid-cube-sideset-dbc.jl"),
    (88, "schwarz-ahead-overlap-dynamic-clamped-3sd-rom-fom-swap-sides.jl"),
    (89, "schwarz-ahead-overlap-dynamic-clamped-3sd-round-trip-swap-chain.jl"),
    (90, "schwarz-nonoverlap-static-cuboid-hex8-aitken-secant.jl"),
    (91, "schwarz-nonoverlap-dynamic-cantilever-dn-aitken-secant.jl"),
    (92, "schwarz-nonoverlap-static-cuboid-robin-robin-aitken.jl"),
    (93, "schwarz-nonoverlap-dynamic-cuboids-impedance-aitken.jl"),
    (94, "schwarz-nonoverlap-dynamic-cantilever-rr-aitken.jl"),
    (95, "schwarz-nonoverlap-dynamic-cantilever-impedance-aitken.jl"),
    (96, "schwarz-ahead-overlap-dynamic-clamped-3sd-multi-swap-all-subdomains.jl"),
    (97, "schwarz-ahead-overlap-dynamic-clamped-single-gaussian-rom-fom-multi-swap.jl"),
    (98, "single-ahead-cuboid-dynamic-restart.jl"),
    (99, "single-ahead-clamped-opinf-fom-restart.jl"),
    (100, "single-ahead-clamped-opinf-rom-restart.jl"),
    (101, "schwarz-ahead-overlap-dynamic-clamped-single-gaussian-fom-fom-restart.jl"),
    (102, "schwarz-ahead-overlap-dynamic-clamped-single-gaussian-opinf-rom-restart.jl"),
    (103, "schwarz-overlap-dynamic-cantilever-impedance-energy.jl"),
    (104, "schwarz-overlap-impedance-tensor-crossscale.jl"),
    (105, "schwarz-overlap-impedance-variational-transfer.jl"),
    (106, "newmark-hht-alpha.jl"),
    (107, "schwarz-ahead-nonoverlap-dynamic-notched-cylinder.jl"),
    (108, "schwarz-overlap-blended-energy.jl"),
    (109, "schwarz-nonoverlap-impedance-adjoint-pairing.jl"),
    (110, "restart-past-final-time.jl"),
    (111, "restart-inplace-checkpoint.jl"),
    (112, "schwarz-overlap-static-cuboids-restart.jl"),
]

# Extract test file names
const nnopinf_test_indices = Int[51]

const all_test_files = [file for (_, file) in indexed_test_files]
const standard_test_indices = [i for (i, _) in indexed_test_files if i ∉ nnopinf_test_indices]
# Optional test indices (excluded from quick runs):
# - 27..32 schwarz-ahead-overlap-dynamic-*
# - 33..38 single-ahead-*
# - 39..43 schwarz-ahead-nonoverlap-dynamic-*
# - 51     nnopinf-schwarz-overlap-cuboid-hex8
# - 54     schwarz-nonoverlap-dynamic-cantilever-dn
# - 84     schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap
# - 85     schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-fom-multi-swap
# - 86     schwarz-ahead-overlap-dynamic-notched-cylinder-opinf-swap-rom-fom
# - 107    schwarz-ahead-nonoverlap-dynamic-notched-cylinder
const optional_test_indices = Int[27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 51, 54, 77, 79, 80, 83, 84, 85, 86, 88, 89, 96, 97, 99, 101, 102, 103, 107]

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
    nnopinf = "--with-nnopinf" in args

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
    elseif nnopinf
        indexed_test_files
    else
        Norma.norma_log(0, :info, "Running standard test suite.")
        filter(t -> t[1] in standard_test_indices, indexed_test_files)
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
        ts = @testset "[$i] $file" begin
            include(file)
        end
        passes, fails, errors = count_test_results(ts)
        if fails == 0 && errors == 0
            Norma.norma_log(0, :done, "[$i] $file: $passes passed")
        else
            Norma.norma_log(0, :error, "[$i] $file: $passes passed, $fails failed, $errors errored")
        end
    end
end

elapsed_time = time() - start_time
Norma.norma_log(0, :done, "Tests Complete")
Norma.norma_log(0, :time, "Tests Run Time = " * Norma.format_time(elapsed_time))
Norma.norma_log(0, :norma, "END TESTS")

# Finish writing the suite log and report where it was written.
Norma.close_session_log_file()
Norma.norma_log(0, :info, "Test log written to $runtests_log_path")

# Cleanup artifacts
for ext in ["yaml", "e", "g", "csv"]
    for file in filter(f -> endswith(f, ".$ext"), readdir())
        rm(file; force=true)
    end
end
