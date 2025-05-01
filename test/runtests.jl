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

all_test_files = [
    "minitensor.jl",
    "interpolation.jl",
    "constitutive.jl",
    "single-static-solid-cube.jl",
    "single-static-solid-neumann-bc.jl",
    "single-implicit-dynamic-solid-cube.jl",
    "single-implicit-dynamic-solid-sho.jl",
    "single-implicit-dynamic-solid-clamped.jl",
    "single-explicit-dynamic-solid-cube.jl",
    "single-explicit-dynamic-solid-sho.jl",
    "single-explicit-dynamic-solid-clamped.jl",
    "tet4-static-solid-cube.jl",
    "tet10-static-solid-cube.jl",
    "schwarz-overlap-static-cuboid-hex8.jl",
    "schwarz-nonoverlap-static-cuboid-hex8.jl",
    "transfer-operators.jl",
    "schwarz-contact-static-cubes.jl",
    "schwarz-contact-dynamic-cubes.jl",
    "solid-inclined-displacement.jl",
    "opinf-schwarz-overlap-cuboid-hex8.jl",
    "quadratic-opinf-schwarz-overlap-cuboid-hex8.jl",
    "adaptive-time-stepping.jl",
]

nn_opinf_test_files = [
    "nnopinf-schwarz-overlap-cuboid-hex8.jl",
]

if "run-nn-opinf-tests" in ARGS
   append!(all_test_files,nn_opinf_test_files)
end

# Must go last for now due to FPE trapping
"utils.jl",
push!(all_test_files,"utils.jl")

const indexed_test_files = collect(enumerate(all_test_files))

# Display test options
Norma.norma_log(0, :info, "Available tests (e.g. ./runtests.jl 2 4 7):")
for (i, file) in indexed_test_files
    Norma.norma_log(0, :info, rpad("[$i]", 5) * file)
end

# Parse command-line arguments as indices (if any)

selected_test_indices = parse.(Int, filter(arg-> arg!= "run-nn-opinf-tests",ARGS))
       
# Determine which tests to run
test_files_to_run = isempty(selected_test_indices) ? indexed_test_files :
    filter(t -> t[1] in selected_test_indices, indexed_test_files)

# Start test run
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

# Cleanup: WARNING — Do not leave files here!
for ext in ["yaml", "e", "g", "csv"]
    for file in filter(f -> endswith(f, ".$ext"), readdir())
        rm(file; force=true)
    end
end
