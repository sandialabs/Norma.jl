# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Options of the topological phase of the adaptivity loop
# (docs/notes/ems-adaptivity), from the `adaptivity` block of the input.
struct AdaptivityOptions
    desired_density::Float64     # elements above this energy density are candidates
    allowed_density::Float64     # no accepted operation may create an element above this
    minimum_decrease::Float64    # relative decrease of the cavity energy an operation must achieve
    adjacency_layers::Int        # rings of adjacent elements added to the candidate set
    maximum_passes::Int          # passes of a topology phase
    outer_iterations::Int        # alternations of smoothing and topology
    swaps::Bool
    collapses::Bool
end

# Result of one topological operation: the cavity replaced, the energies
# before and after, and for a collapse the node removed and the node it was
# moved onto (zero otherwise).
struct CavityProposal
    old_elements::Vector{Int}
    new_connectivity::Matrix{Int}
    block::Int
    energy_before::Float64
    energy_after::Float64
    removed_node::Int
    surviving_node::Int
end
