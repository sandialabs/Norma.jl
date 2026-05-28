# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Criterion for triggering a mid-run sub-simulation swap.
abstract type SwapCriterion end

# Trigger when the multi-domain controller's time first exceeds t_swap.
struct TimeSwapCriterion <: SwapCriterion
    t_swap::Float64
end

# Trigger a swap from a stress-recovery error indicator: the relative
# Frobenius-norm difference between the lumped and consistent L2-projected
# nodal stress fields.  The two projections agree on a well-resolved field and
# diverge on an under-resolved one, so their difference is a cheap surrogate
# for the discretization error.  `direction` selects which way the swap goes:
#
#   :refine   fire when the difference EXCEEDS `tolerance` — the error
#             indicator is large, so a coarse model should be replaced by a
#             finer one.
#   :coarsen  fire when the difference is BELOW `tolerance` — the error
#             indicator is small, so a fine model may be replaced by a
#             coarser one.
#
# The lumped inverse-mass vector and the Cholesky factorization of the
# consistent mass matrix depend only on the mesh, so they are built once and
# cached.  The cache is keyed by model identity: a multi-domain criterion is
# evaluated against every subsim, so each model keeps its own operators rather
# than sharing one slot.  Per-step cost is then one shared RHS assembly plus
# two cheap solves (one diagonal scaling, one triangular back-substitution).
# When a model already carries a `BothRecovery`, its operators are reused.
#
# A zero stress field is reported as zero difference: under :refine that does
# not fire (no error), under :coarsen it does (trivially well-resolved).
mutable struct StressRecoverySwapCriterion <: SwapCriterion
    tolerance::Float64
    direction::Symbol   # :refine (swap above tolerance) or :coarsen (below)
    # Per-model recovery operators, built lazily on first evaluation.
    _recovery::IdDict{SolidMechanics,Tuple{LumpedRecovery,ConsistentRecovery}}
end

# Convenience constructor: the recovery-operator cache starts empty.
StressRecoverySwapCriterion(tolerance::Float64, direction::Symbol) =
    StressRecoverySwapCriterion(
        tolerance, direction, IdDict{SolidMechanics,Tuple{LumpedRecovery,ConsistentRecovery}}()
    )

# Trigger a swap based on the Schwarz overlap L2 displacement error: the L2
# norm of the displacement mismatch between subdomains at the overlapping
# interface nodes.
#
# RESTRICTIONS — validation aborts with a descriptive message if violated:
#   • Only valid for MultiDomainSimulation with overlapping Schwarz coupling.
#     Using it with non-overlapping Schwarz, Schwarz contact, or a monolithic
#     (single-domain) formulation is a configuration error.
#   • At least one SolidMechanicsOverlapSchwarzBoundaryCondition must carry
#     `compute overlap L2 relative error: disp` (or `velo` or `acce`).  Without that
#     setting the error is never computed and the criterion cannot function.
#
# `direction` selects which way the swap fires:
#   :refine   fire when the maximum overlap L2 error EXCEEDS `tolerance` —
#             the coupling mismatch is too large, so swap to a finer mesh.
#   :coarsen  fire when the maximum overlap L2 error falls BELOW `tolerance`
#             — the coupling is sufficiently tight, so swap to a coarser mesh.
#
# For a MultiDomainSimulation the criterion evaluates the maximum error over
# all overlap BCs (with the flag enabled) across every subsim.
struct SchwarzOverlapL2SwapCriterion <: SwapCriterion
    tolerance::Float64
    direction::Symbol   # :refine (swap above tolerance) or :coarsen (below)
end

# A scheduled swap: replace a simulation with the one described by
# `replacement_file` once `criterion` fires.  `subsim_name` is required for
# multi-domain swaps (it selects which subsim to replace) and `nothing` for
# single-domain swaps (where there is no choice).  One-shot.
mutable struct SwapPlan
    subsim_name::Union{Nothing,String}
    criterion::SwapCriterion
    replacement_file::String
    applied::Bool
end
