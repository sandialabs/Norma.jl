# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Benchmark: finite-difference vs analytical consistent tangent for J2Plasticity.
#
# For each of two representative states (elastic step, plastic step):
#   • Verifies that the analytical and FD tangents agree to within FD truncation error.
#   • Times both implementations over N_REPS repetitions and reports the speedup.
#
# Run standalone:
#   cd test && julia --project=.. benchmark-j2-tangent.jl

using LinearAlgebra
using Printf
using StaticArrays
using Test

include("../src/Norma.jl")

const N_REPS = 500

# ── helpers ────────────────────────────────────────────────────────────────────

function time_reps(f, args...; n=N_REPS)
    f(args...)            # warm-up / JIT
    t = @elapsed for _ in 1:n
        f(args...)
    end
    return t / n          # seconds per call
end

function max_rel_err(A::SArray{Tuple{3,3,3,3}}, B::SArray{Tuple{3,3,3,3}})
    return norm(A - B) / max(norm(A), norm(B), 1.0)
end

# ── material ───────────────────────────────────────────────────────────────────

params = Norma.Parameters(
    "model"             => "j2 plasticity",
    "elastic modulus"   => 200.0e9,
    "Poisson's ratio"   => 0.3,
    "density"           => 7800.0,
    "yield stress"      => 250.0e6,
    "hardening modulus" => 20.0e9,
)
mat = Norma.J2Plasticity(params)

# ── deformation states ─────────────────────────────────────────────────────────

# State 1 – elastic: small uniaxial stretch (well below yield, σy = 250 MPa)
F_el = SMatrix{3,3,Float64,9}(
    1.0,  0.0,  0.0,
    0.0,  1.0,  0.0,
    0.0,  0.0,  1.001,
)
s_el = Norma.initial_state(mat)   # Fp = I, eqps = 0

# State 2 – plastic: larger strain, pre-loaded plastic state
F_pre = SMatrix{3,3,Float64,9}(
    0.98,  0.02,  0.0,
    0.0,   0.97,  0.01,
    0.0,   0.0,   1.05,
)
_, _, s_pre = Norma._j2_stress(mat, F_pre, s_el)

F_pl = SMatrix{3,3,Float64,9}(
    0.97,  0.03,  0.0,
    0.0,   0.96,  0.02,
    0.0,   0.0,   1.08,
)

# ── benchmark ──────────────────────────────────────────────────────────────────

println("\n── J2 Plasticity Tangent Benchmark  (N_REPS = $N_REPS) ──\n")
@printf("  %-16s  %12s  %12s  %8s  %12s\n",
        "case", "FD (μs)", "Analytic (μs)", "Speedup", "Rel error")
println("  ", "─"^68)

@testset "J2 Tangent: analytical vs FD" begin
    for (label, F, s_old) in [
            ("elastic step", F_el, s_el),
            ("plastic step", F_pl, s_pre),
        ]

        W, P, s_new = Norma._j2_stress(mat, F, s_old)

        AA_fd = Norma._j2_tangent_fd(mat, F, s_old, P)
        AA_an = Norma._j2_tangent_analytical(mat, F, s_old, P, s_new)

        t_fd = time_reps(Norma._j2_tangent_fd, mat, F, s_old, P)
        t_an = time_reps(Norma._j2_tangent_analytical, mat, F, s_old, P, s_new)

        speedup = t_fd / t_an
        err     = max_rel_err(AA_fd, AA_an)

        @printf("  %-16s  %12.2f  %12.2f  %7.1fx  %12.2e\n",
                label, t_fd*1e6, t_an*1e6, speedup, err)

        # The analytical tangent uses the frozen-Fᵖ approximation (Fᵖ_new treated
        # as constant when differentiating P).  For elastic steps this is exact
        # (no plastic flow, so frozen-Fᵖ = true tangent); for plastic steps the
        # error is O(Δεᵖ) relative to the full consistent tangent.
        tol = label == "elastic step" ? 1e-4 : 1e-2
        @test err < tol
    end
end
println()
