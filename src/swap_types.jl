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

# Trigger when the relative Frobenius-norm difference between the lumped and
# consistent L2-projected nodal stress fields falls below `tolerance`.
#
# On the first evaluation the lumped inverse-mass vector and the Cholesky
# factorization of the consistent mass matrix are built from the model's
# current mesh and cached for all subsequent time steps, so the mesh-traversal
# cost is paid only once.  Per-step cost is one shared RHS assembly plus two
# cheap solves (one diagonal scaling, one triangular back-substitution).
#
# A small relative difference indicates that both projection methods agree,
# which is characteristic of a well-resolved stress field.  When the stress
# is identically zero (both norms vanish) the method returns `true`
# (converged) by convention.
#
# The default tolerance is 1 % (1.0e-2).
mutable struct StressRecoverySwapCriterion <: SwapCriterion
    tolerance::Float64
    # Populated on the first should_swap evaluation; nothing until then.
    _lumped::Union{Nothing,LumpedRecovery}
    _consistent::Union{Nothing,ConsistentRecovery}
end

# Convenience constructor: specify only the tolerance (or rely on the default).
StressRecoverySwapCriterion(tol::Float64 = 1.0e-2) =
    StressRecoverySwapCriterion(tol, nothing, nothing)

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
