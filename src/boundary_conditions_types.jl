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
# non-overlap) and transmission condition (Dirichlet vs. impedance/Robin).  The
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
    rotation_matrix::Matrix{Float64}
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

# Impedance-matching overlap Schwarz: replaces DBC-DBC with absorbing
# conditions on the overlap boundaries. Same strong interpolation as
# regular overlap, but applies t + Z u̇ = g as a force (not a constraint).
# Setup data for variationally consistent partner-traction extraction on a
# node-aligned (conforming) impedance-overlap interface. The weak partner
# traction ∫_Γ t_p·φ_i is obtained by partial assembly of the partner's
# discrete momentum residual (internal force + inertia) over the single layer
# of partner elements on the exterior side of the interface — the exact
# discrete traction, including mesh-scale content that nodal stress recovery
# cannot represent. Built in compute_impedance_overlap_schwarz_projectors!.
struct ConsistentTractionPatch
    element_nodes::Vector{Vector{Int64}}              # global partner node indices per patch element
    element_block::Vector{Int64}                      # partner block index per patch element
    element_index::Vector{Int64}                      # block-local element index (for material state)
    element_rows::Vector{Vector{Tuple{Int64,Int64}}}  # (element-local node, target-index) pairs
    block_element_type::Vector{ElementType}           # indexed by partner block index
    lumped_mass::Bool                                 # partner integrator uses lumped mass (explicit)
    # Target space of the accumulated weak flux. On a node-aligned interface
    # the targets are the BC-local boundary nodes and transfer is empty (the
    # weak flux lands on this side's trace space directly). On a
    # non-conforming interface the targets are the nodes of the offset
    # partner facet surface Γ̃ (the exterior/interior element boundary within
    # one partner element of Γ), and transfer = R W̃⁻¹ carries the weak flux
    # from Γ̃ to this side's trace space (surface-to-surface variational
    # transfer; stage 2a of the variational-transfer design note).
    num_targets::Int64
    transfer::Matrix{Float64}                         # (BC-local nodes) x (Γ̃ nodes); empty if identity
    # Stage 2c (characteristic transfer) data, empty on node-aligned
    # interfaces: the partner's full Robin datum g̃ = t̃ + Z u̇ + α u is
    # assembled on Γ̃ from the partner's own nodal trace values and carried
    # over by the single operator `transfer`, so the mesh-scale cancellation
    # between traction and velocity happens on the partner side, before
    # transfer.
    tilde_mass::Matrix{Float64}                       # W̃, surface mass on Γ̃
    tilde_nodes::Vector{Int64}                        # Γ̃ global partner node indices (target order)
    tilde_normals::Matrix{Float64}                    # 3 x (Γ̃ nodes), receiving-side normals
end

mutable struct SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition <: SolidMechanicsOverlapCouplingSchwarzBoundaryCondition
    name::String
    side_set_id::Int64
    side_set_node_indices::Vector{Int64}
    num_nodes_sides::Vector{Int64}
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
    local_from_global_map::Dict{Int64,Int64}
    global_from_local_map::Vector{Int64}
    square_projector::Matrix{Float64}
    # P/S-split tensor impedance Z = Z_p n⊗n + Z_s (I - n⊗n), with n the
    # interface normal. Both impedances are the NEIGHBOR subdomain's
    # characteristic values (optimized-Schwarz cross-scaling: the optimal
    # transmission operator approximates the neighbor's DtN map).
    impedance::Float64           # Z_p = ρ c_p = √(ρ(λ + 2μ)) of the neighbor
    impedance_shear::Float64     # Z_s = ρ c_s = √(ρμ) of the neighbor
    robin_parameter::Float64     # α for displacement penalty (0 = pure impedance)
    impedance_scale::Vector{Float64}  # multiplier on Z per step (default [1.0])
    # Partner-traction evaluation: "auto" (consistent traction when the interface
    # is node-aligned, recovered stress otherwise), "consistent traction", or
    # "recovered stress". traction_patch is nothing when recovered stress is used.
    partner_traction_mode::String
    traction_patch::Union{ConsistentTractionPatch,Nothing}
    # Transfer of the partner fields to this boundary: "pointwise"
    # interpolation (default) or "variational" L2 projection onto this side's
    # trace space (mortar, in the domain-decomposition literature), which is
    # non-expansive in L2 (the contractivity that pointwise interpolation
    # lacks on non-conforming interfaces).
    # variational_projector is (num boundary nodes) x (num partner nodes),
    # empty when pointwise.
    transfer_mode::String
    transfer_subdivisions::Int64  # facet-quadrature subdivisions for the L assembly
    variational_projector::Matrix{Float64}
    # Content-aware absorption: LK-dissipate the component of this side's
    # boundary velocity that the partner's trace space cannot represent,
    # (I - Π) u̇ with Π the W-orthogonal projection onto the span of the
    # variational projector's columns. Constants are transferable (P·1 = 1),
    # so rigid motions are untouched; on node-aligned interfaces the filter
    # vanishes identically. content_filter = I - Π, empty when disabled.
    content_absorption::Bool
    content_filter::Matrix{Float64}
    # Representable dashpot: restrict the impedance term to the component of
    # the velocity jump that BOTH trace spaces can represent,
    # Z Π (u̇_p − u̇), with the same W-orthogonal Π as above. On a
    # nonconforming interface the unfiltered jump cannot vanish at Schwarz
    # convergence (the coarse trace space cannot represent the fine side's
    # content), so the full dashpot persistently absorbs that content —
    # measured as the dominant spurious interface dissipation, and of
    # indefinite sign under non-adjoint (pointwise) transfer. Filtered, the
    # dashpot does zero work at convergence; the unrepresentable content sees
    # traction-only transmission and reflects (conserving) instead of being
    # silently absorbed. Absorption of that content remains available,
    # explicitly and sign-controlled, via content_absorption. On node-aligned
    # interfaces Π = I and this is a no-op. representable_projector = Π,
    # empty when disabled.
    representable_dashpot::Bool
    representable_projector::Matrix{Float64}
    coupled_block_name::String
    search_tolerance::Float64
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
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end

# Impedance-matching Robin-Robin Schwarz: t + Z u̇ + α W u = g
# Z = ρ c_p (characteristic impedance) absorbs outgoing waves at the interface,
# preventing reflections that cause energy growth with mixed integrators.
# The displacement penalty α W u provides quasi-static stability.
mutable struct SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition <: SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition
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
    impedance::Float64           # Z = ρ c_p = √(ρ(λ + 2μ))
    robin_parameter::Float64     # α for displacement penalty (0 = pure impedance)
    # Adjoint (variationally paired) transfer: both sides of the interface
    # derive their transfer operators from ONE shared cross-mass matrix
    # B_mn = ∫_Γ φ¹_m φ²_n dS, so that W₁ Π₁ = (W₂ Π₂)ᵀ = B and each side's
    # force transfer is the adjoint of the partner's kinematic transfer
    # (N₁ = Π₂ᵀ, N₂ = Π₁ᵀ). With a shared impedance Z and Robin α this makes
    # the dashpot's interface power telescope to -Z ∫_Γ [[u̇ʰ]]² dS ≤ 0 and
    # the Robin term a conservative interface spring — the discrete
    # energy-stability condition of the DG/mortar literature (see
    # doc/notes/schwarz-interface-energy, Eq. (adjoint)), available on this
    # nonoverlap variant precisely because the two sides share a single
    # interface Γ. Both sides of a pair must set it.
    adjoint_pairing::Bool
    parent::Simulation
    self_handle::DomainHandle
    coupled_handle::DomainHandle
end


include("opinf/opinf_ics_bcs_types.jl")
