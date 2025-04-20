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

const test_files = [
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
    "cubic-opinf-schwarz-overlap-cuboid-hex8.jl",
    "adaptive-time-stepping.jl",
    # Must go last for now due to FPE trapping
    "utils.jl",
]

@testset verbose = true "Norma.jl Test Suite" begin
    for file in test_files
        println("▶ \e[1;32mRunning $file...\e[0m")
        include(file)
    end
end

# WARNING: Do not leave output, meshes or inputs here.
# They will be removed.
for ext in ["yaml", "e", "g", "csv"]
    for file in filter(f -> endswith(f, ".$ext"), readdir())
        rm(file; force=true)
    end
end
