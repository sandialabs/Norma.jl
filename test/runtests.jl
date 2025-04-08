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

@testset verbose = true "Norma.jl Test Suite" begin
    include("minitensor.jl")
    include("interpolation.jl")
    include("constitutive.jl")
    include("single-static-solid-cube.jl")
    include("single-static-solid-neumann-bc.jl")
    include("single-implicit-dynamic-solid-cube.jl")
    include("single-implicit-dynamic-solid-sho.jl")
    include("single-implicit-dynamic-solid-clamped.jl")
    include("single-explicit-dynamic-solid-cube.jl")
    include("single-explicit-dynamic-solid-sho.jl")
    include("single-explicit-dynamic-solid-clamped.jl")
    include("tet4-static-solid-cube.jl")
    include("tet10-static-solid-cube.jl")
    include("schwarz-overlap-static-cuboid-hex8.jl")
    include("schwarz-nonoverlap-static-cuboid-hex8.jl")
    include("transfer-operators.jl")
    include("schwarz-contact-static-cubes.jl")
    include("schwarz-contact-dynamic-cubes.jl")
    include("solid-inclined-displacement.jl")
    include("opinf-schwarz-overlap-cuboid-hex8.jl")
    include("quadratic-opinf-schwarz-overlap-cuboid-hex8.jl")
    include("adaptive-time-stepping.jl")
    # Must go last for now due to FPE trapping
    include("utils.jl")
end
# WARNING: Do not leave output, meshes or inputs here.
# They will be removed.
for ext in ["yaml", "e", "g", "csv"]
    for file in filter(f -> endswith(f, ".$ext"), readdir())
        rm(file; force=true)
    end
end
