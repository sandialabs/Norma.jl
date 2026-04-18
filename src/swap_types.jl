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

# A scheduled swap: replace the subsim named `subsim_name` with the sim
# described by `replacement_file` once `criterion` fires. One-shot.
mutable struct SwapPlan
    subsim_name::String
    criterion::SwapCriterion
    replacement_file::String
    applied::Bool
end
