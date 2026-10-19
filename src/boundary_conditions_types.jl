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

using SparseArrays

# Periodic condition (input key `Periodic`) on a box cell whose opposite faces
# are pairs of node sets, the second the image of the first under a
# translation along a coordinate axis. Every node of the second set follows its
# image on the first, u(X⁺) = u(X⁻) + (F̄ - I)(X⁺ - X⁻), with F̄ the macroscopic
# deformation gradient. The condition eliminates the displacements of these
# dependent nodes: each is that of its independent image plus the relative
# displacements of the corners it is away from it, Σ_k n_k (u(C_k) - u(O)). The
# origin O, the node common to the first sets of all pairs, is held fixed,
# and its image C_k across pair k carries the column of F̄ - I along the axis
# of the pair: u(C_k) - u(O) = (F̄ - I) L_k e_k. A component of F̄ is either
# prescribed, F̄_ij(t), or free with the matching component of the volume
# average of the first Piola-Kirchhoff stress prescribed, P̄_ij(t), through a
# force on the corner. The transformation A maps the independent degrees of
# freedom to all of them, u = A u, and condenses the residual and the tangent
# to Aᵀ r and Aᵀ K A.
mutable struct SolidMechanicsPeriodicBoundaryCondition <: SolidMechanicsRegularBoundaryCondition
    name::String
    axes::Vector{Int64}             # coordinate axis of each pair of node sets
    pair_names::Vector{Vector{String}}  # names of the node sets of each pair
    lengths::Vector{Float64}        # period L_k along that axis
    volume::Float64                 # volume of the cell
    origin::Int64
    corners::Vector{Int64}          # C_k, the image of the origin across pair k
    dependent_nodes::Vector{Int64}
    independent_nodes::Vector{Int64}  # image on the first sets of each dependent node
    periods::Matrix{Int64}          # n_k of each dependent node, pairs × dependent nodes
    components::Vector{CartesianIndex{2}}  # components (i, j) of F̄ along the periodic axes
    component_pairs::Vector{Int64}  # pair k whose axis is j, for each component
    stress_controlled::Vector{Bool}
    component_funs::Vector{Function}  # F̄_ij(t), or P̄_ij(t) if stress controlled
    displacement_gradient::Matrix{Float64}  # F̄ - I at the last application
    # The increment of the prescribed components of F̄ since the last solve,
    # which the solve's prediction applies to the whole cell (see
    # periodic_affine_increment!).
    pending_increment::Matrix{Float64}
    # The translation of the origin since the last solve, which the
    # prediction applies to the whole cell too, and the displacement of the
    # origin before the last application, from which it is measured.
    pending_translation::Vector{Float64}
    origin_displacement::Vector{Float64}
    transformation::SparseMatrixCSC{Float64,Int64}  # A
    transformation_transpose::SparseMatrixCSC{Float64,Int64}
    checked::Bool  # whether the prescribed conditions were checked against it
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

mutable struct SolidMechanicsNonOverlapSchwarzBoundaryCondition <:
               SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition
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
    # Cholesky factorization of square_projector, the boundary mass matrix W of
    # this side, which is symmetric positive definite and constant during a run.
    # It is computed once with the projectors (set_square_projector! in
    # schwarz.jl) and used for every solve with W, among them the W⁻¹ norm of
    # the interface force residual. `nothing` until the projectors are built.
    square_factor::Union{Nothing,Cholesky{Float64,Matrix{Float64}}}
    is_dirichlet::Bool
    swap_bcs::Bool
    # Matching exchange (input key `matching`): the Dirichlet side imposes one
    # interface quantity, the projected velocity or displacement of the
    # partner, and derives the other two kinematic fields from its own Newmark
    # relations; the Neumann side receives the d'Alembert reaction of the
    # Dirichlet side through the transpose of the Dirichlet projector, and both
    # transfer operators are built from one cross mass matrix. The fixed point
    # of the Schwarz iteration is then the constrained problem of Gravouil and
    # Combescure (2001), whose constraint is the matching. `constrained` is true
    # for a matching pair; `constraint` is the matched quantity, :velocity or
    # :displacement, and :unset for the standard exchange.
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
    # Input key `matching at: window end`: a subcycled Dirichlet side imposes
    # the matched quantity at the end of the window only, and at its interior
    # substeps its interface rows carry the linear ramp in time between its
    # reactions at the two ends of the window (window_start_reaction, the
    # converged reaction of the previous stop, and window_end_reaction, the
    # reaction at the end of its last pass). The force increments of its
    # substeps then add up to the change of its reaction and the window balance
    # of the rate invariant cancels, and the impulse of the ramp equals the
    # impulse the coarse side receives: a fine-step Dirichlet side conserves
    # both, as the coarse-step Dirichlet side does. With `every substep` (the
    # default) the matched quantity is imposed at every substep from the
    # partner's interpolated history.
    window_end_constraint::Bool
    window_start_reaction::Vector{Float64}
    # The window-end reaction that forms the ramp of the current iteration
    # (relaxed or Anderson-accelerated) and the one measured at the end of
    # the last pass.
    window_end_reaction::Vector{Float64}
    window_end_reaction_raw::Vector{Float64}
    # The end value of the ramp that the interior substeps of the current
    # iteration applied; the stopping criterion compares it with the
    # reaction measured at the end of the window.
    window_ramp_end::Vector{Float64}
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# Energy-conserving two-sided overlap coupling (input key `Schwarz overlap
# matching`): each subdomain is constrained on its own Schwarz boundary Γ_k to
# the L2 projection W_k⁻¹ L_k of the partner's matched quantity, with L_k the
# cross-mass matrix of Γ_k against the partner's volume shape functions, and
# receives on all of its nodes the partner's d'Alembert reaction on Γ_j through
# the transpose of the partner's projector, -Π_jᵀ ρ_j. Mass, internal force,
# stiffness, and energy of every quadrature point inside the partner are
# weighted by one half (ASSEMBLY_WEIGHT_CACHE in model.jl), so that the fixed
# point of the Schwarz iteration is the blended constrained problem and the
# window balance of the rate invariant cancels as for a nonoverlap matching
# pair. Both sides are constrained. With different time steps the fine-step
# side sets `matching at: window end`: it imposes the matched quantity at the
# end of the window only and carries the ramp of its own reaction at its
# interior substeps, while the coarse-step side takes one step per window.
# `interface solve: direct` (both sides explicit) applies at equal steps only.
# The fields mirror those of the nonoverlap condition that the matching
# machinery reads: `is_dirichlet`, `constrained`, and `swap_bcs` are fixed,
# and `transferred_force` holds the force received on all nodes.
mutable struct SolidMechanicsOverlapMatchingSchwarzBoundaryCondition <:
               SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    coupled_block_name::String
    search_tolerance::Float64
    coupled_bc_index::Int64
    dirichlet_projector::Matrix{Float64}
    neumann_projector::Matrix{Float64}
    square_projector::Matrix{Float64}
    square_factor::Union{Nothing,Cholesky{Float64,Matrix{Float64}}}
    is_dirichlet::Bool
    swap_bcs::Bool
    constrained::Bool
    constraint::Symbol
    direct_solve::Bool
    transferred_force::Vector{Float64}
    transferred_force_history::Vector{Tuple{Float64,Vector{Float64}}}
    # As for the nonoverlap condition: the fine-step side of an overlap matching
    # pair with different steps must set `matching at: window end`.
    window_end_constraint::Bool
    window_start_reaction::Vector{Float64}
    window_end_reaction::Vector{Float64}
    window_end_reaction_raw::Vector{Float64}
    # The end value of the ramp that the interior substeps of the current
    # iteration applied; the stopping criterion compares it with the
    # reaction measured at the end of the window.
    window_ramp_end::Vector{Float64}
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# The two conditions whose pairs are constrained problems: the nonoverlap
# matching pair and the two-sided overlap matching coupling.
const SolidMechanicsConstrainedSchwarzBoundaryCondition = Union{
    SolidMechanicsNonOverlapSchwarzBoundaryCondition,SolidMechanicsOverlapMatchingSchwarzBoundaryCondition
}

# Robin-Robin nonoverlap Schwarz (input key `Schwarz RR nonoverlap`): the
# classical Robin transmission condition t + α W u = g on each side of the
# interface, with t the interface traction, u the interface displacement, W the
# boundary mass matrix of the side, α the Robin parameter, and g the datum
# assembled from the traction and displacement of the partner. Each side
# transfers the partner fields with its own Dirichlet and Neumann projectors,
# so the two sides of an interface may use different values of α.
mutable struct SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition <:
               SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition
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
