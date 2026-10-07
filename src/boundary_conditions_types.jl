# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

abstract type BoundaryCondition end
abstract type InitialCondition end
abstract type SolidMechanicsBoundaryCondition <: BoundaryCondition end
abstract type SolidMechanicsRegularBoundaryCondition <: SolidMechanicsBoundaryCondition end
abstract type SolidMechanicsNeumannRobinBoundaryCondition <: SolidMechanicsRegularBoundaryCondition end
abstract type SolidMechanicsSchwarzBoundaryCondition <: SolidMechanicsBoundaryCondition end
abstract type SolidMechanicsCouplingSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition end
# Coupling Schwarz BCs split along two orthogonal axes: geometry (overlap vs.
# non-overlap) and transmission condition (Dirichlet vs. Robin).  The
# geometry axis is captured here so "is this an overlapping coupling?" is a
# single `isa` on the abstract type rather than an enumeration of concretes that
# every call site must keep in sync (the omission of one concrete from such an
# enumeration is exactly what double-counted the overlap energy).  An overlapping
# coupling meshes a shared volume twice; a non-overlapping one partitions the
# domain across a measure-zero interface.
abstract type SolidMechanicsOverlapCouplingSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition end
abstract type SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition <: SolidMechanicsCouplingSchwarzBoundaryCondition end

using Symbolics

mutable struct SolidMechanicsDirichletBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_fun::Function
    velo_fun::Function
    acce_fun::Function
    # An expression without t is evaluated once per node on first application;
    # its velocity and acceleration are zero by construction.
    time_dependent::Bool
    cached_displacement::Vector{Float64}
end

# Dirichlet BC whose constrained nodes are the nodes of a side set rather than
# a node set.  Behaves exactly like SolidMechanicsDirichletBoundaryCondition
# (prescribes displacement/velocity/acceleration and constrains the DOFs), but
# resolves its node list from a side set via read_side_set_node_list.
mutable struct SolidMechanicsSideSetDirichletBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    node_indices::Vector{Int64}  # unique nodes belonging to the side set
    disp_fun::Function
    velo_fun::Function
    acce_fun::Function
    # An expression without t is evaluated once per node on first application;
    # its velocity and acceleration are zero by construction.
    time_dependent::Bool
    cached_displacement::Vector{Float64}
end

mutable struct SolidMechanicsNeumannBoundaryCondition <: SolidMechanicsNeumannRobinBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_fun::Function
end

# Analytic level-set surface constraint for energetic mesh smoothing.  The nodes
# of a side set are held to a smooth boundary surface g(x) = 0 while sliding
# freely within it — an inclined roller.  The exact gradient ∇g is obtained
# automatically from the symbolic g via Symbolics.gradient, so the user supplies
# only g.  A node shared by two such side sets accumulates both constraints,
# realizing the edge (curve) case; three, a pinned vertex — the dimensional
# hierarchy falls out of side-set membership with no special case.
#
# Two enforcement modes share this type and its YAML (selected by `enforcement`):
#   :exact   — local-frame roller (the production mechanism).  The residual is
#              projected onto the tangent subspace (its normal component is the
#              constraint reaction) and, after each step, the node is
#              closest-point-projected back onto the surface.  Surface-exact up
#              to the return-to-surface tolerance; needs a matrix-free solver.
#   :penalty — a quadratic penalty P = ½ κ g² whose force κ g ∇g is added to the
#              internal force and whose Gauss–Newton Hessian κ ∇g ∇gᵀ is added to
#              the stiffness.  Approximate (residual ~ 1/κ) but works with the
#              Newton solver too; a derisking/validation tool.
mutable struct SolidMechanicsSurfaceBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    side_set_id::Int64
    node_indices::Vector{Int64}  # unique nodes belonging to the side set
    level_set_fun::Function      # g(t, x, y, z)
    level_set_grad::Function     # ∇g(t, x, y, z) -> 3-vector
    enforcement::Symbol          # :exact or :penalty
    penalty::Float64             # used only when enforcement == :penalty
end

mutable struct SolidMechanicsRobinBoundaryCondition <: SolidMechanicsNeumannRobinBoundaryCondition
    name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_fun::Function
    robin_parameter::Float64
end

mutable struct SolidMechanicsNeumannPressureBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    pressure_fun::Function
end

mutable struct SolidMechanicsContactSchwarzBoundaryCondition <: SolidMechanicsSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    active_contact::Bool
    friction_type::Int64
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

mutable struct SolidMechanicsOverlapSchwarzBoundaryCondition <: SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    compute_overlap_l2_error::String
    overlap_node_indices::Vector{Int64}
    overlap_coupled_nodes_indices::Vector{Vector{Int64}}
    overlap_interpolation_function_values::Vector{Vector{Float64}}
    overlap_l2_error::Float64
    coupled_block_name::String
    search_tolerance::Float64
    dirichlet_projector::Matrix{Float64}
    use_weak::Bool
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

mutable struct SolidMechanicsNonOverlapSchwarzBoundaryCondition <: SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    is_dirichlet::Bool
    swap_bcs::Bool
    # Constrained exchange (input key `constrained`): the Dirichlet side imposes
    # one interface quantity, the projected velocity or displacement of the
    # partner (`constraint`), and derives the other two kinematic fields from
    # its own Newmark relations; the Neumann side receives the d'Alembert
    # reaction of the Dirichlet side through the transpose of the Dirichlet
    # projector, and both transfer operators are built from one cross mass
    # matrix. The fixed point of the Schwarz iteration is then the constrained
    # problem of Gravouil and Combescure (2001). `constraint` is :velocity or
    # :displacement once the pair is resolved, :unset while neither side of the
    # pair has named it.
    constrained::Bool
    constraint::Symbol
    # Direct interface solve (input key `interface solve: direct`): for a
    # constrained pair of two central difference subdomains at equal steps the
    # interface reaction is computed in one solve per stop instead of by the
    # Schwarz iteration (direct_interface_stop! in schwarz.jl).
    direct_solve::Bool
    # Interface force, in interleaved components of the interface nodes, that
    # this side received at its last application as the Neumann side. Read by
    # the interface force residual of the Schwarz loop.
    transferred_force::Vector{Float64}
    # The same force at every substep of the current stop, with its time, for
    # the per-substep force residual of a subcycled constrained pair.
    transferred_force_history::Vector{Tuple{Float64,Vector{Float64}}}
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# Robin-Robin nonoverlap Schwarz (input key `Schwarz RR nonoverlap`): the
# classical Robin transmission condition t + α W u = g on each side of the
# interface, with t the interface traction, u the interface displacement, W the
# boundary mass matrix of the side, α the Robin parameter, and g the datum
# assembled from the traction and displacement of the partner. Each side
# transfers the partner fields with its own Dirichlet and Neumann projectors,
# so the two sides of an interface may use different values of α.
mutable struct SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition <: SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_bc_name::String
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    robin_parameter::Float64
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end


include("opinf/opinf_ics_bcs_types.jl")
