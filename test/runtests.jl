# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using Statistics
using Test

include("../src/Norma.jl")
include("helpers.jl")

include("interpolation.jl")
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
#include("schwarz-contact-dynamic-cubes.jl")
include("solid-cube-inclined-support.jl")
include("opinf-schwarz-overlap-cuboid-hex8.jl")
include("adaptive-time-stepping.jl")

