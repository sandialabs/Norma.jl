# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using StaticArrays
using ForwardDiff

function elastic_constants(params::Parameters)
    E = 0.0
    ν = 0.0
    κ = 0.0
    λ = 0.0
    μ = 0.0
    if haskey(params, "elastic modulus") == true
        E = params["elastic modulus"]
        if haskey(params, "Poisson's ratio") == true
            ν = params["Poisson's ratio"]
            κ = E / 3(1 - 2ν)
            λ = E * ν / (1 + ν) / (1 - 2ν)
            μ = E / 2(1 + ν)
        elseif haskey(params, "bulk modulus") == true
            κ = params["bulk modulus"]
            ν = (3κ - E) / 6κ
            λ = (3κ * (3κ - E)) / (9κ - E)
            μ = 3κ * E / (9κ - E)
        elseif haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            R = sqrt(E^2 + 9λ^2 + 2E * λ)
            ν = 2λ / (E + λ + R)
            κ = (E + 3λ + R) / 6
            μ = (E - 3λ + R) / 4
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            ν = E / 2μ - 1
            κ = E * μ / 3(3μ - E)
            λ = μ * (E - 2μ) / (3μ - E)
        else
            norma_abort("Two elastic constants are required but only elastic modulus found")
        end
    elseif haskey(params, "Poisson's ratio") == true
        ν = params["Poisson's ratio"]
        if haskey(params, "bulk modulus") == true
            κ = params["bulk modulus"]
            E = 3κ * (1 - 2ν)
            λ = 3κ * ν / (1 + ν)
            μ = 3κ * (1 - 2ν) / 2(1 + ν)
        elseif haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            E = λ * (1 + ν) * (1 - 2ν) / ν
            κ = λ * (1 + ν) / 3ν
            μ = λ * (1 - 2ν) / 2ν
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = 2μ * (1 + ν)
            κ = 2μ * (1 + ν) / 3(1 - 2ν)
            λ = 2μ * ν / (1 - 2ν)
        else
            norma_abort("Two elastic constants are required but only Poisson's ratio found")
        end
    elseif haskey(params, "bulk modulus") == true
        κ = params["bulk modulus"]
        if haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            E = 9κ * (κ - λ) / (3κ - λ)
            ν = λ / (3κ - λ)
            μ = 3(κ - λ) / 2
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = 9κ * μ / (3κ + μ)
            ν = (3κ - 2μ) / 2(3κ + μ)
            λ = κ - 2μ / 3
        else
            norma_abort("Two elastic constants are required but only bulk modulus found")
        end
    elseif haskey(params, "Lamé's first constant") == true
        λ = params["Lamé's first constant"]
        if haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = μ * (3λ + 2μ) / (λ + μ)
            ν = λ / 2(λ + μ)
            κ = λ + 2μ / 3
        else
            norma_abort("Two elastic constants are required but only Lamé's first constant found")
        end
    elseif haskey(params, "shear modulus") == true
        norma_abort("Two elastic constants are required but only shear modulus found")
    else
        norma_abort("Two elastic constants are required but none found")
    end
    return E, ν, κ, λ, μ
end

@inline number_states(::Elastic) = 0
@inline internal_variable_names(::Solid) = String[]

mutable struct SaintVenant_Kirchhoff <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function SaintVenant_Kirchhoff(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Linear_Elastic <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Linear_Elastic(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Neohookean <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Neohookean(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Reciprocal_Neohookean <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Reciprocal_Neohookean(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct SethHill <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    m::Int
    n::Int
    function SethHill(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ, params["m"], params["n"])
    end
end

# Hencky (logarithmic strain) hyperelasticity (issue #71):
#   ψ = κ/2 (tr E)² + μ devE : devE,  E = ½ log C
# The quadratic energy in logarithmic strain makes uniaxial stress exactly
# linear in log stretch at any deformation, which is what makes this model
# useful for testing strong nonlinearities.  Formulation cross-validated
# against the Hencky model in ConstitutiveModels (Carina).
mutable struct Hencky <: Elastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Hencky(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        return new(E, ν, κ, λ, μ, ρ)
    end
end

# ---------------------------------------------------------------------------
# J2 plasticity — finite deformation, multiplicative split F = Fᵉ * Fᵖ,
# Hencky (logarithmic) elasticity, radial-return in Mandel-stress space,
# linear isotropic hardening.
#
# State vector (10 entries):
#   state[1:9]  = vec(Fᵖ)  column-major 3×3 plastic deformation gradient
#   state[10]   = εᵖ        accumulated equivalent plastic strain
# ---------------------------------------------------------------------------

mutable struct J2Plasticity <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    σy::Float64   # initial yield stress
    H::Float64    # linear isotropic hardening modulus
    function J2Plasticity(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        σy = get(params, "yield stress", 0.0)
        H = get(params, "hardening modulus", 0.0)
        return new(E, ν, κ, λ, μ, ρ, σy, H)
    end
end

@inline number_states(::J2Plasticity) = 10
@inline initial_state(::J2Plasticity) = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

# Names for each entry of the state vector, in storage order.
# state[1:9] = vec(Fᵖ) column-major: Fp_11, Fp_21, Fp_31, Fp_12, ...
# state[10]  = εᵖ (equivalent plastic strain)
@inline internal_variable_names(::J2Plasticity) = [
    "Fp_11", "Fp_21", "Fp_31",
    "Fp_12", "Fp_22", "Fp_32",
    "Fp_13", "Fp_23", "Fp_33",
    "eqps",
]

function second_from_fourth(AA::SArray{Tuple{3,3,3,3},Float64,4})
    # Row-major pair packing: M[3(i-1)+j, 3(k-1)+l] = AA[i,j,k,l].
    # Matches the row-major layout of the gradient operator so that
    # grad_op' * M * grad_op yields the standard element stiffness.
    M = MMatrix{9,9,Float64,81}(undef)
    @inbounds for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        M[3 * (i - 1) + j, 3 * (k - 1) + l] = AA[i, j, k, l]
    end
    return SMatrix(M)
end

const I3 = @SMatrix [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
]

# Placeholder returned in lieu of the material tangent when a matrix-free
# solver (e.g. steepest descent, L-BFGS) is in use and the tangent would only
# be discarded.  Constructing the real tangent (AD Hessian for SethHill, or the
# convected/finite-difference tangents for the other models) dominates the
# per-element cost, so skipping it speeds up every residual-only evaluation.
const ZERO_TANGENT = zero(SArray{Tuple{3,3,3,3},Float64,4,81})

function odot(A::SMatrix{3,3,Float64,9}, B::SMatrix{3,3,Float64,9})
    C = MArray{Tuple{3,3,3,3},Float64}(undef)
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for d in 1:3
                    C[a, b, c, d] = 0.5 * (A[a, c] * B[b, d] + A[a, d] * B[b, c])
                end
            end
        end
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function ox(A::SMatrix{3,3,Float64,9}, B::SMatrix{3,3,Float64,9})
    C = MArray{Tuple{3,3,3,3},Float64}(undef)
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for d in 1:3
                    C[a, b, c, d] = A[a, b] * B[c, d]
                end
            end
        end
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function oxI(A::SMatrix{3,3,Float64,9})
    C = zeros(MArray{Tuple{3,3,3,3},Float64})
    for a in 1:3
        C[:, :, a, a] .= A
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function Iox(B::SMatrix{3,3,Float64,9})
    C = zeros(MArray{Tuple{3,3,3,3},Float64})
    for a in 1:3
        C[a, a, :, :] .= B
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function convect_tangent(CC::SArray{Tuple{3,3,3,3},Float64}, S::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9})
    AA = MArray{Tuple{3,3,3,3},Float64}(undef)
    I_n = @SMatrix [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    for j in 1:3
        for l in 1:3
            M = @SMatrix [
                CC[1, j, l, 1] CC[1, j, l, 2] CC[1, j, l, 3]
                CC[2, j, l, 1] CC[2, j, l, 2] CC[2, j, l, 3]
                CC[3, j, l, 1] CC[3, j, l, 2] CC[3, j, l, 3]
            ]
            G = F * M * F'
            AA[:, j, :, l] .= S[l, j] .* I_n .+ G
        end
    end
    return SArray{Tuple{3,3,3,3}}(AA)
end

# ---------------------------------------------------------------------------
# Strain energy density functions — generic in element type so that
# ForwardDiff can propagate dual numbers through them.
# ---------------------------------------------------------------------------

function strain_energy(material::Linear_Elastic, F::SMatrix{3,3,T,9}) where {T<:Number}
    ∇u = F - I3
    ϵ = 0.5 .* (∇u + ∇u')
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    return 0.5 * λ * trϵ^2 + μ * tr(ϵ * ϵ)
end

function strain_energy(material::SaintVenant_Kirchhoff, F::SMatrix{3,3,T,9}) where {T<:Number}
    C = F' * F
    E = 0.5 .* (C - I3)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    return 0.5 * λ * trE^2 + μ * tr(E * E)
end

function strain_energy(material::Neohookean, F::SMatrix{3,3,T,9}) where {T<:Number}
    C = F' * F
    J2 = det(C)
    Jm23 = inv(cbrt(J2))
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1.0)
    Wdev = 0.5 * μ * (Jm23 * trC - 3.0)
    return Wvol + Wdev
end

function strain_energy(material::Reciprocal_Neohookean, F::SMatrix{3,3,T,9}) where {T<:Number}
    C = F' * F
    J = det(F)
    J² = J * J
    J⁻¹ = 1.0 / J
    Jm23 = inv(cbrt(J²))
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * ((J - 1.0)^2 + (J⁻¹ - 1)^2)
    Wdev = 0.5 * μ * (Jm23 * trC - 3.0)
    return Wvol + Wdev
end

function strain_energy(material::SethHill, F::SMatrix{3,3,T,9}) where {T<:Number}
    C = F' * F
    F⁻¹ = inv(F)
    F⁻ᵀ = F⁻¹'
    J = det(F)
    Jᵐ = J^material.m
    J⁻ᵐ = 1 / Jᵐ
    Cbar = J^(-2 / 3) * C
    Cbar⁻¹ = J^(2 / 3) * F⁻¹ * F⁻ᵀ
    Cbarⁿ = Cbar^material.n
    Cbar⁻ⁿ = Cbar⁻¹^material.n
    trCbarⁿ = tr(Cbarⁿ)
    trCbar⁻ⁿ = tr(Cbar⁻ⁿ)
    trCbar²ⁿ = tr(Cbarⁿ * Cbarⁿ)
    trCbar⁻²ⁿ = tr(Cbar⁻ⁿ * Cbar⁻ⁿ)
    Wbulk = material.κ / 4 / material.m^2 * ((Jᵐ - 1)^2 + (J⁻ᵐ - 1)^2)
    Wshear = material.μ / 4 / material.n^2 * (trCbar²ⁿ + trCbar⁻²ⁿ - 2 * trCbarⁿ - 2 * trCbar⁻ⁿ + 6)
    return Wbulk + Wshear
end

# Principal square root of a symmetric positive definite matrix by the
# Denman-Beavers iteration.  Pure matrix arithmetic (no eigendecomposition),
# quadratically convergent from any SPD start, so a fixed iteration count
# reaches machine precision with a wide margin; extra iterations are
# stationary at the fixed point.  Forward-mode AD passes through: the dual
# parts converge along with the values.
function sqrt_spd(A::SMatrix{3,3,T,9}) where {T<:Number}
    Y = A
    Z = SMatrix{3,3,T,9}(I)
    # The iteration averages exponents at rate one bit per step before its
    # quadratic phase, so the count must cover log2 of the condition number.
    for _ in 1:32
        Y, Z = 0.5 * (Y + inv(Z)), 0.5 * (Z + inv(Y))
    end
    return 0.5 * (Y + Y')
end

# Matrix logarithm of a symmetric positive definite matrix by inverse scaling
# and squaring: square roots bring the spectrum near one, an atanh series
# evaluates the logarithm there, and the roots are undone by doubling.  With
# ‖A - I‖_F ≤ 1/4 every eigenvalue b of B = (A - I)(A + I)⁻¹ satisfies
# |b| ≤ (1/4)/(2 - 1/4) = 1/7, so the series through B¹⁷ truncates below
# machine precision.  No eigendecomposition and no branches on the spectrum,
# so the result is smooth at repeated eigenvalues, including C = I, the
# reference state of every simulation, where it returns exactly zero.
# LinearAlgebra's matrix logarithm is LAPACK-bound and Float64-only; this one
# is generic in the element type so ForwardDiff can differentiate through it
# (stress_ad and tangent_ad).  Comparisons on dual numbers use their value
# parts, which is exactly the branch semantics wanted here.
function log_spd(C::SMatrix{3,3,T,9}) where {T<:Number}
    I3 = SMatrix{3,3,T,9}(I)
    A = 0.5 * (C + C')
    # Factor out the determinant: log C = log Ĉ + ⅓ log(det C) I with Ĉ
    # unimodular.  The scalar part is exact, pure scalings need no matrix
    # work at all, and the remaining spectrum is centered on one, which
    # shortens both the square-root loop and each Denman-Beavers run.
    d = det(A)
    d > 0 || norma_abort("Matrix logarithm requires a positive definite input.")
    A = A / cbrt(d)
    k = 0
    while norm(A - I3) > 0.25
        k == 32 && norma_abort("Matrix logarithm did not converge; the input is not positive definite.")
        A = sqrt_spd(A)
        k += 1
    end
    B = (A - I3) * inv(A + I3)
    B² = B * B
    # log A = 2 (B + B³/3 + B⁵/5 + ... + B¹⁷/17), Horner in B².
    S = B² / 17
    for c in (15, 13, 11, 9, 7, 5, 3)
        S = B² * (S + I3 / c)
    end
    return 2.0^(k + 1) * (B * (S + I3)) + log(d) / 3 * I3
end

function strain_energy(material::Hencky, F::SMatrix{3,3,T,9}) where {T<:Number}
    C = F' * F
    E = 0.5 * log_spd(C)
    trE = tr(E)
    devE = E - trE / 3 * SMatrix{3,3,T,9}(I)
    return 0.5 * material.κ * trE * trE + material.μ * sum(devE .* devE)
end

function constitutive(material::SaintVenant_Kirchhoff, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    C = F' * F
    E = 0.5 .* (C - I3)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * (trE^2) + μ * tr(E * E)
    S = λ * trE .* I3 .+ 2.0 .* μ .* E
    P = F * S
    need_tangent || return W, P, ZERO_TANGENT
    CC_m = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    CC_m[i, j, k, l] = λ * I3[i, j] * I3[k, l] + μ * (I3[i, k] * I3[j, l] + I3[i, l] * I3[j, k])
                end
            end
        end
    end
    CC_s = SArray{Tuple{3,3,3,3}}(CC_m)
    AA = convect_tangent(CC_s, S, F)
    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    ∇u = F - I3
    ϵ = 0.5 .* (∇u + ∇u')
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * (trϵ^2) + μ * tr(ϵ * ϵ)
    σ = λ * trϵ .* I3 .+ 2.0 .* μ .* ϵ
    need_tangent || return W, σ, ZERO_TANGENT
    CC_m = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    CC_m[i, j, k, l] = λ * I3[i, j] * I3[k, l] + μ * (I3[i, k] * I3[j, l] + I3[i, l] * I3[j, k])
                end
            end
        end
    end
    CC_s = SArray{Tuple{3,3,3,3}}(CC_m)
    return W, σ, CC_s
end

function constitutive(material::Neohookean, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    C = F' * F
    J2 = det(C)
    Jm23 = inv(cbrt(J2))
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1.0)
    Wdev = 0.5 * μ * (Jm23 * trC - 3.0)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J2 - 1.0) .* IC
    Sdev = μ .* Jm23 .* (I3 .- (IC .* (trC / 3.0)))
    S = Svol .+ Sdev
    P = F * S
    need_tangent || return W, P, ZERO_TANGENT
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3.0
    CCvol = κ .* (J2 .* ICxIC .- (J2 - 1.0) .* ICoIC)
    CCdev = μJ2n .* (trC .* (ICxIC ./ 3 .+ ICoIC) .- oxI(IC) .- Iox(IC))
    CC = CCvol .+ CCdev
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function constitutive(material::Reciprocal_Neohookean, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    C = F' * F
    J = det(F)
    J² = J * J
    J⁻¹ = 1.0 / J
    J⁻² = J⁻¹ * J⁻¹
    Jm23 = inv(cbrt(J²))
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * ((J - 1.0)^2 + (J⁻¹ - 1)^2)
    Wdev = 0.5 * μ * (Jm23 * trC - 3.0)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J² - J - J⁻² + J⁻¹) .* IC
    Sdev = μ .* Jm23 .* (I3 .- (IC .* (trC / 3.0)))
    S = Svol .+ Sdev
    P = F * S
    need_tangent || return W, P, ZERO_TANGENT
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3.0
    CCvol = 0.5 * κ .* (2.0 * J² - J + 2.0 * J⁻² - J⁻¹) .* ICxIC .- κ .* (J² - J - J⁻² + J⁻¹) .* ICoIC
    CCdev = μJ2n .* (trC .* (ICxIC ./ 3 .+ ICoIC) .- oxI(IC) .- Iox(IC))
    CC = CCvol .+ CCdev
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

# ---------------------------------------------------------------------------
# AD primitives — composable building blocks for mixed manual/AD constitutive
# implementations.  Any Elastic model with a strain_energy method can use
# these independently, e.g. manual stress + AD tangent, or vice-versa.
# ---------------------------------------------------------------------------

function stress_ad(material::Elastic, F::SMatrix{3,3,Float64,9})
    Fv = Vector{Float64}(undef, 9)
    Fv .= vec(F)
    W_func = x -> strain_energy(material, SMatrix{3,3,eltype(x),9}(x))
    W = W_func(Fv)
    Pv = ForwardDiff.gradient(W_func, Fv)
    P = SMatrix{3,3,Float64,9}(reshape(Pv, 3, 3))
    return W, P
end

function tangent_ad(material::Elastic, F::SMatrix{3,3,Float64,9})
    Fv = Vector{Float64}(undef, 9)
    Fv .= vec(F)
    W_func = x -> strain_energy(material, SMatrix{3,3,eltype(x),9}(x))
    Hmat = ForwardDiff.hessian(W_func, Fv)
    return SArray{Tuple{3,3,3,3},Float64,4,81}(reshape(Hmat, 3, 3, 3, 3))
end

# SethHill: manual (W, P) — shares intermediates — plus AD tangent.
function constitutive(material::SethHill, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    C = F' * F
    F⁻¹ = inv(F)
    F⁻ᵀ = F⁻¹'
    J = det(F)
    Jᵐ = J^material.m
    J⁻ᵐ = 1.0 / Jᵐ
    J²ᵐ = Jᵐ * Jᵐ
    J⁻²ᵐ = 1.0 / J²ᵐ
    Cbar = J^(-2 / 3) * C
    Cbar⁻¹ = J^(2 / 3) * F⁻¹ * F⁻ᵀ
    Cbarⁿ = Cbar^material.n
    Cbar⁻ⁿ = Cbar⁻¹^material.n
    Cbar²ⁿ = Cbarⁿ * Cbarⁿ
    Cbar⁻²ⁿ = Cbar⁻ⁿ * Cbar⁻ⁿ
    trCbarⁿ = tr(Cbarⁿ)
    trCbar⁻ⁿ = tr(Cbar⁻ⁿ)
    trCbar²ⁿ = tr(Cbar²ⁿ)
    trCbar⁻²ⁿ = tr(Cbar⁻²ⁿ)
    Wbulk = material.κ / 4 / material.m^2 * ((Jᵐ - 1)^2 + (J⁻ᵐ - 1)^2)
    Wshear = material.μ / 4 / material.n^2 * (trCbar²ⁿ + trCbar⁻²ⁿ - 2 * trCbarⁿ - 2 * trCbar⁻ⁿ + 6)
    W = Wbulk + Wshear
    Pbulk = material.κ / 2 / material.m * (J²ᵐ - Jᵐ - J⁻²ᵐ + J⁻ᵐ) * F⁻ᵀ
    Pshear =
        material.μ / material.n *
        (1 / 3 * (-trCbar²ⁿ + trCbarⁿ + trCbar⁻²ⁿ - trCbar⁻ⁿ) * F⁻ᵀ + F⁻ᵀ * (Cbar²ⁿ - Cbarⁿ - Cbar⁻²ⁿ + Cbar⁻ⁿ))
    P = Pbulk + Pshear
    need_tangent || return W, P, ZERO_TANGENT
    AA = tangent_ad(material, F)
    return W, P, AA
end

function constitutive(material::Hencky, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    κ = material.κ
    μ = material.μ
    I3 = SMatrix{3,3,Float64,9}(I)
    C = F' * F
    E = 0.5 * log_spd(C)
    trE = tr(E)
    devE = E - trE / 3 * I3
    W = 0.5 * κ * trE * trE + μ * sum(devE .* devE)
    # For isotropic ψ(E) with E = ½ log C, the work-conjugate stress
    # M = ∂ψ/∂E commutes with C, so S = 2 ∂ψ/∂C = C⁻¹ M exactly.
    M = κ * trE * I3 + 2.0 * μ * devE
    S = inv(C) * M
    P = F * S
    need_tangent || return W, P, ZERO_TANGENT
    AA = tangent_ad(material, F)
    return W, P, AA
end

# ---------------------------------------------------------------------------
# Full AD fallback for Elastic models without any manual implementation.
# Specific methods above take dispatch priority; this catches anything else
# (e.g. new models added by a developer who only defines strain_energy).
# ---------------------------------------------------------------------------

function constitutive(material::Elastic, F::SMatrix{3,3,Float64,9}; need_tangent::Bool=true)
    Fv = Vector{Float64}(undef, 9)
    Fv .= vec(F)
    W_func = x -> strain_energy(material, SMatrix{3,3,eltype(x),9}(x))
    W = W_func(Fv)
    Pv = ForwardDiff.gradient(W_func, Fv)
    P = SMatrix{3,3,Float64,9}(reshape(Pv, 3, 3))
    need_tangent || return W, P, ZERO_TANGENT
    Hmat = ForwardDiff.hessian(W_func, Fv)
    AA = SArray{Tuple{3,3,3,3},Float64,4,81}(reshape(Hmat, 3, 3, 3, 3))
    return W, P, AA
end

# Deviatoric part of a matrix.  Used by the Simo-Hughes J2 model below.
# Subtracts through UniformScaling rather than I(3): the latter builds a
# Diagonal{Bool,Vector{Bool}}, which allocates, and this sits on the hot path
# of every plastic quadrature point.
function _dev(A::AbstractMatrix{T}) where {T}
    return A - (tr(A) / 3) * I
end

# ---------------------------------------------------------------------------
# Material factory and kinematics
# ---------------------------------------------------------------------------

function create_material(params::Parameters)
    model_name = params["model"]
    if model_name == "linear elastic"
        return Linear_Elastic(params)
    elseif model_name == "Saint-Venant Kirchhoff"
        return SaintVenant_Kirchhoff(params)
    elseif model_name == "neohookean"
        return Neohookean(params)
    elseif model_name == "r-neohookean"
        return Reciprocal_Neohookean(params)
    elseif model_name == "seth-hill"
        return SethHill(params)
    elseif model_name == "hencky"
        return Hencky(params)
    elseif model_name == "j2 plasticity"
        return J2Plasticity(params)
    else
        norma_abort("Unknown material model : $model_name")
    end
    return nothing
end

function get_kinematics(material::Solid)
    if material isa Linear_Elastic
        return Infinitesimal
    elseif material isa SaintVenant_Kirchhoff
        return Finite
    elseif material isa Neohookean
        return Finite
    elseif material isa Reciprocal_Neohookean
        return Finite
    elseif material isa SethHill
        return Finite
    elseif material isa Hencky
        return Finite
    elseif material isa J2Plasticity
        return Finite
    end
    norma_abort("Unknown material model : $(typeof(material))")
    return nothing
end

function get_p_wave_modulus(material::Solid)
    return material.λ + 2.0 * material.μ
end

# ---------------------------------------------------------------------------
# Simo-Hughes J2 finite-deformation plasticity (BOX 9.1 + 9.2)
# ---------------------------------------------------------------------------
#
# Reference: Simo & Hughes, Computational Inelasticity, pp 317-321.
#
#   Kinematics:    F = Fᵉ Fᵖ  (multiplicative decomposition)
#   Elasticity:    neo-Hookean vol/dev split, U(J) + μ/2 (tr b̄ᵉ - 3),
#                  with U(J) = κ/2 (J - 1)²
#   Yield surface: von Mises, f = ‖s‖ - √(2/3)(σy + K α)
#   Flow rule:     associated, isochoric
#   Hardening:     linear isotropic
#
# This is the same model, and the same BOX 9.1/9.2 algorithm, that
# ConstitutiveModels.jl provides to Carina, so the two codes agree on J2.
#
# On the tangent: BOX 9.2 is the MAJOR-SYMMETRIC PART of the exact Jacobian of
# the stress update, not the Jacobian itself.  The return map evaluates the
# effective shear modulus μ̄ = μ tr(b̄ᵉ_trial)/3 at the trial state, so the
# discrete update is not exactly variational and its true Jacobian carries a
# small antisymmetric part, which BOX 9.2 drops -- deliberately, since the
# solver assembles a symmetric operator.  Measured in test/constitutive.jl
# ("J2Plasticity Consistent Tangent"): on steps that yield, the symmetric part
# is recovered to ~1e-8 while the discarded antisymmetric part is O(1e-4) of
# the tangent norm, so Newton stays fast but is not fully quadratic on those
# steps.  The elastic branch is the exact Jacobian.

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

    # Unit normal.  Both branches must yield the same concrete type; a bare
    # zeros(3,3) here widens n to a Union and sends every use below through a
    # dynamic dispatch.
    n = s_trial_norm > 0.0 ?
        SMatrix{3,3,Float64,9}(μ * _dev(be_bar_tr) / s_trial_norm) :
        zero(SMatrix{3,3,Float64,9})

    f_trial = s_trial_norm - sqrt(2.0/3.0) * (σy + K * α_n)
    plastic = f_trial > 0.0

    # BOX 9.2, Steps 2-3: plastic correction coefficients, and dev(n²).  None
    # of these depend on (i,j,k,l), so they are computed once here rather than
    # 81 times inside the assembly loop.
    β₁ = 0.0
    β₃ = 0.0
    β₄ = 0.0
    n_sq_dev = zero(SMatrix{3,3,Float64,9})
    if plastic
        β₀ = 1.0 + K / (3μ̄)
        β₁ = 2μ̄ * Δγ / s_trial_norm
        β₂ = (1.0 - 1.0/β₀) * 2.0/3.0 * s_trial_norm / μ * Δγ
        β₃ = 1.0/β₀ - β₁ + β₂
        β₄ = (1.0/β₀ - β₁) * s_trial_norm / μ̄
        n_sq_dev = SMatrix{3,3,Float64,9}(_dev(n * n))
    end

    # Spatial tangent c_{ijkl}: volumetric (BOX 9.2 Step 1) + deviatoric trial,
    # plus the plastic correction when the step yielded.  Assembled in a single
    # pass into a stack-allocated MArray; the previous form made two passes over
    # a heap zeros(3,3,3,3).
    CC_spatial = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        δij = Float64(i == j); δkl = Float64(k == l)
        δik = Float64(i == k); δjl = Float64(j == l)
        δil = Float64(i == l); δjk = Float64(j == k)
        I_sym = 0.5 * (δik * δjl + δil * δjk)

        # C = (JU')'·J · (1⊗1) − 2·J·U' · I_sym,  for U(J) = κ/2 (J-1)²
        c = coeff_1x1 * δij * δkl - coeff_I * I_sym

        # C̄_trial = 2μ̄[I_sym - (1/3)1⊗1] - (2/3)‖s_trial‖[n⊗1 + 1⊗n]
        c_dev = 2μ̄ * (I_sym - 1.0/3.0 * δij * δkl) -
                2.0/3.0 * s_trial_norm * (n[i,j] * δkl + δij * n[k,l])

        if plastic
            c_dev_n2 = 0.5 * (n[i,j] * n_sq_dev[k,l] + n_sq_dev[i,j] * n[k,l])
            c += (1.0 - β₁) * c_dev -
                 2μ̄ * β₃ * n[i,j] * n[k,l] -
                 2μ̄ * β₄ * c_dev_n2
        else
            c += c_dev
        end
        CC_spatial[i,j,k,l] = c
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
# Top-level constitutive call
# ---------------------------------------------------------------------------

function constitutive(material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64}; need_tangent::Bool=true)
    W, P, state_new, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n =
        _sh_j2_stress(material, F, state_old)
    need_tangent || return W, P, ZERO_TANGENT, state_new
    AA = _sh_j2_tangent(material, F, state_old, P,
                         s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n)
    return W, P, AA, state_new
end
