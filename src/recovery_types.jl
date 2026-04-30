# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using SparseArrays

abstract type AbstractRecoveryData end

struct NoRecovery <: AbstractRecoveryData end

# Geometry-only (density-free), scalar (1 DOF/node) lumped projection mass.
# inv_m[i] = 1 / Σ_e Σ_q N_i(ξ_q) w_q |J|, with inv_m[i] = 0 for nodes whose
# geometric tributary is non-positive (orphan nodes, or sign cancellation in
# higher-order quadrature rules).
struct LumpedRecovery <: AbstractRecoveryData
    inv_m::Vector{Float64}
end

# Geometry-only (density-free), scalar (1 DOF/node) consistent projection
# mass: M_ij = Σ_e Σ_q N_i(ξ_q) N_j(ξ_q) w_q |J|.  The matrix is symmetric
# positive definite, so we cache a Cholesky factorization once at setup and
# back-substitute per component at output time.
struct ConsistentRecovery{F} <: AbstractRecoveryData
    M::SparseMatrixCSC{Float64,Int64}
    factor::F
end
