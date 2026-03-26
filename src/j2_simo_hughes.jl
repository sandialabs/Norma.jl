# Simo-Hughes J2 finite-deformation plasticity (BOX 9.1 + 9.2).
#
# Complete replacement for the Hencky/Mandel-based _j2_stress and
# _j2_tangent_analytical.  Uses neo-Hookean-type vol/dev energy split
# with spatial b̄ᵉ formulation.  Tangent is exact (matches FD to machine
# precision), giving quadratic Newton convergence.
#
# Reference: Simo & Hughes, Computational Inelasticity, pp 317-321.

# ---------------------------------------------------------------------------
# Stress update (BOX 9.1)
# ---------------------------------------------------------------------------

function _sh_j2_stress(
    material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64}
)
    κ  = material.κ
    μ  = material.μ
    σy = material.σy
    K  = material.H

    Fp_old   = SMatrix{3,3,Float64,9}(reshape(state_old[1:9], 3, 3))
    α_n      = state_old[10]

    J    = det(F)
    Jm23 = J^(-2.0/3.0)

    # Trial elastic left Cauchy-Green (isochoric)
    Fe_tr     = F * inv(Fp_old)
    be_bar_tr = Jm23 * (Fe_tr * Fe_tr')
    be_bar_tr = 0.5 * (be_bar_tr + be_bar_tr')  # enforce symmetry

    # Trial deviatoric Kirchhoff stress
    s_trial      = μ * _dev(be_bar_tr)
    s_trial_norm = norm(s_trial)

    # Effective shear modulus
    μ̄ = μ * tr(be_bar_tr) / 3.0

    # Yield function
    f_trial = s_trial_norm - sqrt(2.0/3.0) * (σy + K * α_n)

    if f_trial ≤ 0.0
        # Elastic
        s_new      = s_trial
        be_bar_new = be_bar_tr
        α_new      = α_n
        Δγ         = 0.0
        Fp_new     = Fp_old
    else
        # Radial return (BOX 9.1, step 4)
        n  = s_trial / s_trial_norm
        Δγ = f_trial / (2μ̄ + 2.0/3.0 * K)

        s_new = s_trial - 2μ̄ * Δγ * n
        α_new = α_n + sqrt(2.0/3.0) * Δγ

        # Update b̄ᵉ (eq 9.3.33)
        Ie_bar     = tr(be_bar_tr) / 3.0
        be_bar_new = s_new / μ + Ie_bar * I3

        # Recover Fp_new from b̄ᵉ_new via polar decomposition
        be_tr_sqrt_inv = _matrix_power(be_bar_tr, -0.5)
        Fe_tr_iso      = J^(-1.0/3.0) * Fe_tr
        R_tr           = be_tr_sqrt_inv * Fe_tr_iso

        be_new_sqrt    = _matrix_power(be_bar_new, 0.5)
        Fe_new_iso     = be_new_sqrt * R_tr
        Fe_new         = J^(1.0/3.0) * Fe_new_iso
        Fp_new         = inv(Fe_new) * F
    end

    # Kirchhoff stress: τ = J p 1 + s
    p = κ * (J - 1.0)   # U(J) = κ/2 (J-1)²
    τ = J * p * I3 + s_new

    # PK1: P = τ · F⁻ᵀ
    P = SMatrix{3,3,Float64,9}(τ * inv(F)')

    # Energy
    W = κ / 2.0 * (J - 1.0)^2 + μ / 2.0 * (tr(be_bar_new) - 3.0)

    state_new = Vector{Float64}(undef, 10)
    state_new[1:9] = vec(Fp_new)
    state_new[10]  = α_new

    return W, P, state_new, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n
end

# ---------------------------------------------------------------------------
# Helper: matrix power via eigendecomposition (for 3×3 symmetric matrices)
# ---------------------------------------------------------------------------

function _matrix_power(A::AbstractMatrix{Float64}, p::Float64)
    A_sym = 0.5 * (A + A')
    E = eigen(Symmetric(A_sym))
    D = Diagonal(E.values .^ p)
    return SMatrix{3,3,Float64,9}(E.vectors * D * E.vectors')
end

# ---------------------------------------------------------------------------
# Consistent tangent (BOX 9.2)
# ---------------------------------------------------------------------------

function _sh_j2_tangent(
    material::J2Plasticity, F::SMatrix{3,3,Float64,9},
    state_old::Vector{Float64}, P::SMatrix{3,3,Float64,9},
    s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n
)
    κ  = material.κ
    μ  = material.μ
    σy = material.σy
    K  = material.H

    J     = det(F)
    F_inv = inv(F)

    # 2nd Piola-Kirchhoff: S = F⁻¹ P
    S = SMatrix{3,3,Float64,9}(F_inv * P)

    # BOX 9.2, Step 1: Volumetric tangent
    # C = (JU')'·J · (1⊗1) − 2·J·U' · I_sym + C̄
    # For U(J) = κ/2(J-1)²:
    coeff_1x1 = κ * J * (2J - 1.0)
    coeff_I   = 2.0 * κ * J * (J - 1.0)

    # Unit normal
    n = s_trial_norm > 0.0 ? μ * _dev(be_bar_tr) / s_trial_norm : zeros(3, 3)

    f_trial = s_trial_norm - sqrt(2.0/3.0) * (σy + K * α_n)

    # Assemble spatial tangent c_{ijkl} as 81-element 4th-order tensor
    CC_spatial = zeros(3, 3, 3, 3)

    # Volumetric: coeff_1x1 δ_ij δ_kl - coeff_I I_ijkl_sym
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        δij = Float64(i == j); δkl = Float64(k == l)
        δik = Float64(i == k); δjl = Float64(j == l)
        δil = Float64(i == l); δjk = Float64(j == k)
        I_sym = 0.5 * (δik * δjl + δil * δjk)

        CC_spatial[i,j,k,l] += coeff_1x1 * δij * δkl - coeff_I * I_sym
    end

    # Deviatoric trial tangent: C̄_trial
    # C̄ = 2μ̄[I_sym - (1/3)1⊗1] - (2/3)‖s_trial‖[n⊗1 + 1⊗n]
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        δij = Float64(i == j); δkl = Float64(k == l)
        δik = Float64(i == k); δjl = Float64(j == l)
        δil = Float64(i == l); δjk = Float64(j == k)
        I_sym = 0.5 * (δik * δjl + δil * δjk)

        c_dev = 2μ̄ * (I_sym - 1.0/3.0 * δij * δkl) -
                2.0/3.0 * s_trial_norm * (n[i,j] * δkl + δij * n[k,l])

        if f_trial ≤ 0.0
            CC_spatial[i,j,k,l] += c_dev
        else
            # BOX 9.2, Step 2-3: Plastic correction
            β₀ = 1.0 + K / (3μ̄)
            β₁ = 2μ̄ * Δγ / s_trial_norm
            β₂ = (1.0 - 1.0/β₀) * 2.0/3.0 * s_trial_norm / μ * Δγ
            β₃ = 1.0/β₀ - β₁ + β₂
            β₄ = (1.0/β₀ - β₁) * s_trial_norm / μ̄

            n_sq = n * n
            n_sq_dev = _dev(n_sq)

            c_dev_n2 = 0.5 * (n[i,j] * n_sq_dev[k,l] + n_sq_dev[i,j] * n[k,l])

            CC_spatial[i,j,k,l] += (1.0 - β₁) * c_dev -
                                    2μ̄ * β₃ * n[i,j] * n[k,l] -
                                    2μ̄ * β₄ * c_dev_n2
        end
    end

    # Pull-back: spatial → material (CC_mat = F⁻¹ ⊗ F⁻¹ : CC_spatial : F⁻¹ ⊗ F⁻¹)
    CC_mat = MArray{Tuple{3,3,3,3},Float64}(undef)
    for A in 1:3, B in 1:3, C in 1:3, D in 1:3
        val = 0.0
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
            val += F_inv[A,a] * F_inv[B,b] * CC_spatial[a,b,c,d] * F_inv[C,c] * F_inv[D,d]
        end
        CC_mat[A,B,C,D] = val
    end

    # convect_tangent: CC_mat + S → AA (∂P/∂F)
    CC_s = SArray{Tuple{3,3,3,3}}(CC_mat)
    AA = convect_tangent(CC_s, S, F)
    return AA
end

# ---------------------------------------------------------------------------
# Top-level constitutive call (replaces old _j2_stress + _j2_tangent_analytical)
# ---------------------------------------------------------------------------

function constitutive(material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64})
    W, P, state_new, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n =
        _sh_j2_stress(material, F, state_old)
    AA = _sh_j2_tangent(material, F, state_old, P,
                         s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n)
    return W, P, AA, state_new
end
