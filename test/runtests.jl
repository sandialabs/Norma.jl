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
    (113, "schwarz-ahead-overlap-dynamic-notched-cylinder.jl"),
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

# A top-level @testset prints its own summary to stdout (bypassing norma_log,
# so it never reaches runtests.log) and, if ANY test anywhere in the tree
# failed or errored, throws a Test.TestSetException right at the closing
# `end` of the block below. That exception would otherwise abort the script
# before any of the summary/cleanup code runs, which is why the log was
# missing its final summary (and "END TESTS") on any run with at least one
# failure -- the common case, not the exception.
#
# We catch that exception so the log always gets finalized. The aggregate
# pass/fail/error totals are accumulated as we go, from the same per-file
# count_test_results(ts) call already used for each file's log line.
#
# An earlier version of this function tried to hand Test.jl's own outer
# testset object to Test.print_test_results a second time, so the log could
# get a byte-for-byte copy of the table @testset prints to stdout. That
# meant reconstructing the outer testset by hand with Test's private
# push_testset/pop_testset/finish functions instead of the `@testset`
# macro, to keep a live reference to it even when finishing it throws. That
# turned out to be too fragile in practice -- it crashed partway through
# rendering the table on a real run, after every test file had already
# completed, with nothing left standing to catch the error or log it,
# silently truncating runtests.log right after the last test file with no
# summary and no "END TESTS" at all. Rather than lean further on that
# undocumented internal API, we go back to the plain `@testset` macro here
# (proven to reliably reach the end of a run, failures included) and build
# our own summary table below from data we already collect per file, which
# has no such dependency.
#
# This all lives in a function rather than directly at the top level of the
# script: assigning to `overall_passes` etc. inside the nested for/try below
# would otherwise hit Julia's soft-scope-vs-hard-scope ambiguity for
# non-interactive scripts (confirmed: this raised
# `UndefVarError: overall_passes not defined in local scope`). Inside a
# function body, scoping is unambiguous and no `global` declarations are
# needed.
function run_test_suite(test_files_to_run)
    overall_passes = 0
    overall_fails = 0
    overall_errors = 0
    suite_failed = false
    # (label, passes, fails, errors, elapsed) per file, in run order -- used
    # to render a Test-Summary-style table after the run, mirroring what
    # Julia's own Test stdlib prints to the console (see below).
    file_summaries = Tuple{String,Int,Int,Int,Float64}[]

    try
        @testset verbose = true "Norma.jl Test Suite" begin
            for (i, file) in test_files_to_run
                Norma.norma_log(0, :test, "[$i] Running $file...")
                file_start = time()
                ts = @testset "[$i] $file" begin
                    include(file)
                end
                file_elapsed = time() - file_start
                passes, fails, errors = count_test_results(ts)
                overall_passes += passes
                overall_fails += fails
                overall_errors += errors
                push!(file_summaries, ("[$i] $file", passes, fails, errors, file_elapsed))
                if fails == 0 && errors == 0
                    Norma.norma_log(0, :done, "[$i] $file: $passes passed")
                else
                    Norma.norma_log(0, :error, "[$i] $file: $passes passed, $fails failed, $errors errored")
                end
            end
        end
    catch e
        if e isa Test.TestSetException
            # The per-file totals accumulated above are already complete and
            # correct at this point -- the outer @testset only throws once
            # its whole body (the entire for loop) has finished running, so
            # every file has already been counted.
            suite_failed = true
        else
            # Not a test failure -- a real error (e.g. a bug in runtests.jl
            # itself). Still try to close the log below, then rethrow so it
            # isn't silently swallowed.
            Norma.norma_log(0, :error, "Unexpected error while running the test suite: $e")
            Norma.close_session_log_file()
            rethrow()
        end
    end

    return overall_passes, overall_fails, overall_errors, suite_failed, file_summaries
end

# Format a duration the way the Test stdlib's own summary table does: plain
# minutes (no rollover to hours) plus one decimal place of seconds, e.g.
# "45.2s" under a minute, "9m10.0s" or "1m03.7s" over one, with the seconds
# field zero-padded to two digits once there's a minutes part.
function format_duration_like_testset(seconds::Real)
    if seconds < 60
        return @sprintf("%.1fs", seconds)
    end
    minutes = floor(Int, seconds / 60)
    secs = seconds - minutes * 60
    return @sprintf("%dm%04.1fs", minutes, secs)
end

# Render a table in the same visual style as the Test stdlib's own
# "Test Summary:" table -- Pass/Fail/[Error]/Total/Time columns, with the
# Error column dropped entirely when nothing in the run actually errored
# (as opposed to failed), same as the real one -- and write it straight to
# the log file as plain text, not through norma_log, so it isn't prefixed
# with a [TAG] the way every other log line is; that's what makes it read
# as one clean table instead of another wall of tagged lines. Each row is
# (label, passes, fails, errors, elapsed_seconds); the first row is treated
# as the grand total and is not indented.
#
# Earlier this used Norma.format_time, which formats "1h 7m 45.24s" -- a
# string longer than the column width the header line reserved for it. The
# lpad() used to right-align the Time column silently does nothing when the
# content is already wider than the requested width, so that string landed
# glued directly onto the end of the previous (Total) column with no space
# between them, e.g. "20491h 7m 45.24s" for a Total of 2049. Sizing every
# column to the actual widest value it will hold (as done below), rather
# than a width chosen in advance, avoids that class of bug entirely.
function log_test_summary_table(io, rows)
    show_errors = any(errors > 0 for (_, _, _, errors, _) in rows)
    # Blank, rather than "0", for a zero-count cell -- Total and Time are
    # the only columns always populated -- matching how the real Test
    # stdlib table reads (a row of all-Pass has nothing at all under Fail).
    zero_blank(n) = n == 0 ? "" : string(n)

    labels = [i == 1 ? label : "  " * label for (i, (label, _, _, _, _)) in enumerate(rows)]
    pass_strs = [zero_blank(p) for (_, p, _, _, _) in rows]
    fail_strs = [zero_blank(f) for (_, _, f, _, _) in rows]
    error_strs = [zero_blank(e) for (_, _, _, e, _) in rows]
    total_strs = [string(p + f + e) for (_, p, f, e, _) in rows]
    time_strs = [format_duration_like_testset(t) for (_, _, _, _, t) in rows]

    label_w = maximum(length, labels)
    pass_w = max(length("Pass"), maximum(length, pass_strs))
    fail_w = max(length("Fail"), maximum(length, fail_strs))
    error_w = show_errors ? max(length("Error"), maximum(length, error_strs)) : 0
    total_w = max(length("Total"), maximum(length, total_strs))
    time_w = max(length("Time"), maximum(length, time_strs))

    function render(label, pass_s, fail_s, error_s, total_s, time_s)
        cols = [lpad(pass_s, pass_w), lpad(fail_s, fail_w)]
        show_errors && push!(cols, lpad(error_s, error_w))
        push!(cols, lpad(total_s, total_w), lpad(time_s, time_w))
        return rpad(label, label_w) * " | " * join(cols, "  ")
    end

    println(io, render("Test Summary:", "Pass", "Fail", "Error", "Total", "Time"))
    for (i, _) in enumerate(rows)
        println(io, render(labels[i], pass_strs[i], fail_strs[i], error_strs[i], total_strs[i], time_strs[i]))
    end
    return nothing
end

overall_passes, overall_fails, overall_errors, suite_failed, file_summaries =
    run_test_suite(test_files_to_run)

elapsed_time = time() - start_time

summary_rows = vcat([("Norma.jl Test Suite", overall_passes, overall_fails, overall_errors, elapsed_time)], file_summaries)
let io = Norma.NORMA_LOG_FILE[]
    if io !== nothing
        log_test_summary_table(io, summary_rows)
        flush(io)
    end
end

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

# Preserve the previous behavior of exiting non-zero when the suite failed,
# now that the exception is caught rather than propagating on its own.
suite_failed && exit(1)
