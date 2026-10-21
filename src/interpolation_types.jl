# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

@enum ElementType begin
    BAR2
    TRI3
    TRI6
    QUAD4
    TETRA4
    TETRA10
    HEX8
    # The 10-node composite tetrahedron of LCM (`Use Composite Tet 10`): a
    # TETRA10 mesh block whose interpolation is piecewise linear on 12
    # sub-tetrahedra. It is never read from a mesh; a block is switched to it
    # with the `composite tet 10` model key.
    COMPOSITE_TETRA10
end
