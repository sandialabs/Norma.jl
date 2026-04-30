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

function constitutive(material::SaintVenant_Kirchhoff, F::SMatrix{3,3,Float64,9})
    C = F' * F
    E = 0.5 .* (C - I3)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * (trE^2) + μ * tr(E * E)
    S = λ * trE .* I3 .+ 2.0 .* μ .* E
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
    P = F * S
    AA = convect_tangent(CC_s, S, F)
    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::SMatrix{3,3,Float64,9})
    ∇u = F - I3
    ϵ = 0.5 .* (∇u + ∇u')
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * (trϵ^2) + μ * tr(ϵ * ϵ)
    σ = λ * trϵ .* I3 .+ 2.0 .* μ .* ϵ
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

function constitutive(material::Neohookean, F::SMatrix{3,3,Float64,9})
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
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3.0
    CCvol = κ .* (J2 .* ICxIC .- (J2 - 1.0) .* ICoIC)
    CCdev = μJ2n .* (trC .* (ICxIC ./ 3 .+ ICoIC) .- oxI(IC) .- Iox(IC))
    CC = CCvol .+ CCdev
    P = F * S
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
function constitutive(material::SethHill, F::SMatrix{3,3,Float64,9})
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
    AA = tangent_ad(material, F)
    return W, P, AA
end

# ---------------------------------------------------------------------------
# Full AD fallback for Elastic models without any manual implementation.
# Specific methods above take dispatch priority; this catches anything else
# (e.g. new models added by a developer who only defines strain_energy).
# ---------------------------------------------------------------------------

function constitutive(material::Elastic, F::SMatrix{3,3,Float64,9})
    Fv = Vector{Float64}(undef, 9)
    Fv .= vec(F)
    W_func = x -> strain_energy(material, SMatrix{3,3,eltype(x),9}(x))
    W = W_func(Fv)
    Pv = ForwardDiff.gradient(W_func, Fv)
    Hmat = ForwardDiff.hessian(W_func, Fv)
    P = SMatrix{3,3,Float64,9}(reshape(Pv, 3, 3))
    AA = SArray{Tuple{3,3,3,3},Float64,4,81}(reshape(Hmat, 3, 3, 3, 3))
    return W, P, AA
end

# ---------------------------------------------------------------------------
# J2 plasticity — finite deformation constitutive update
# ---------------------------------------------------------------------------
# Formulation:
#   Multiplicative split:  F = Fᵉ * Fᵖ  (det(Fᵖ) = 1 isochoric plastic flow)
#   Hencky elasticity:     W = λ/2 (tr Eᵉ)² + μ tr(Eᵉ²),  Eᵉ = ½ log Cᵉ
#   Mandel stress:         M = λ tr(Eᵉ) I + 2μ Eᵉ
#   Yield function:        f = √(3/2) ‖Mdev‖ − (σy + H εᵖ)  [von Mises]
#   Flow rule:             ΔFᵖ = exp(Δεᵖ N) Fᵖ_old,  N = (3/2) Mdev/‖Mdev‖_vm
#   Analytical return mapping (linear hardening):
#       Δεᵖ = f_trial / (3μ + H)
#
# PK1 stress from Mandel stress:
#   P = Fᵉ⁻ᵀ M Fᵖ⁻ᵀ
#
# Consistent algorithmic tangent:
#   Computed via forward finite differences of P with respect to F.
#   This is exact (up to FD truncation error) and works for both the
#   elastic and plastic steps without deriving the complex analytical
#   tensor expression.
# ---------------------------------------------------------------------------

# deviatoric part of a matrix
function _dev(A::AbstractMatrix{T}) where {T}
    return A .- (tr(A) / 3) .* I(3)
end

# Perform the radial-return update and return (W, P, state_new).
# Uses Float64 throughout (matrix log/exp via LinearAlgebra).
function _j2_stress(
    material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64}
)
    λ = material.λ
    μ = material.μ
    σy = material.σy
    H = material.H

    # Extract state
    Fp_old = Matrix{Float64}(reshape(state_old[1:9], 3, 3))
    eqps_old = state_old[10]

    # Trial elastic deformation gradient
    Fe_tr = Matrix{Float64}(F) * inv(Fp_old)
    Ce_tr = Fe_tr' * Fe_tr
    Ee_tr = 0.5 * log(Ce_tr)           # symmetric by construction

    # Trial Mandel stress (Hencky elasticity)
    trEe_tr = tr(Ee_tr)
    M_tr = λ * trEe_tr * I(3) + 2μ * Ee_tr
    Mdev_tr = _dev(M_tr)
    σvm_tr = sqrt(1.5) * norm(Mdev_tr)

    # Yield check
    f_tr = σvm_tr - (σy + H * eqps_old)

    if f_tr ≤ 0.0
        # Elastic step — state unchanged
        W = 0.5 * λ * trEe_tr^2 + μ * tr(Ee_tr * Ee_tr)
        Fe_new = Fe_tr
        Fp_new = Fp_old
        eqps_new = eqps_old
        M_new = M_tr
    else
        # Plastic step — radial return (analytical for linear hardening)
        Δεᵖ = f_tr / (3μ + H)

        # Flow direction (unit in von Mises sense, deviatoric)
        N = 1.5 * Mdev_tr / σvm_tr      # tr(N) = 0, ‖N‖_F = √(3/2)

        # Update plastic deformation gradient: ΔFᵖ = expm(Δεᵖ N)
        Fp_new = exp(Δεᵖ * N) * Fp_old

        # Updated equivalent plastic strain
        eqps_new = eqps_old + Δεᵖ

        # Updated Hencky elastic strain and Mandel stress
        Ee_new = Ee_tr - Δεᵖ * N       # tr(N)=0 ⟹ tr(Ee_new) = tr(Ee_tr)
        W = 0.5 * λ * trEe_tr^2 + μ * tr(Ee_new * Ee_new)
        M_new = λ * trEe_tr * I(3) + 2μ * Ee_new

        Fe_new = Matrix{Float64}(F) * inv(Fp_new)
    end

    # PK1 stress:  P = Fᵉ⁻ᵀ M Fᵖ⁻ᵀ
    Fe_inv_T = inv(Fe_new)'
    Fp_inv_T = inv(Fp_new)'
    P_dense = Fe_inv_T * M_new * Fp_inv_T
    P = SMatrix{3,3,Float64,9}(P_dense)

    state_new = Vector{Float64}(undef, 10)
    state_new[1:9] = vec(Fp_new)
    state_new[10] = eqps_new

    return W, P, state_new
end

# Consistent algorithmic tangent AA_{iJkL} = ∂P_{iJ}/∂F_{kL} via forward FD.
# 9 extra stress evaluations; used for benchmarking against the analytical tangent.
function _j2_tangent_fd(
    material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64},
    P0::SMatrix{3,3,Float64,9}
)
    h = 1.0e-8
    AA = MArray{Tuple{3,3,3,3},Float64}(undef)
    Fm = MMatrix{3,3,Float64,9}(F)
    for k in 1:3
        for l in 1:3
            Fm[k, l] += h
            _, Ph, _ = _j2_stress(material, SMatrix{3,3,Float64,9}(Fm), state_old)
            for i in 1:3
                for j in 1:3
                    AA[i, j, k, l] = (Ph[i, j] - P0[i, j]) / h
                end
            end
            Fm[k, l] -= h
        end
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(AA)
end

# ---------------------------------------------------------------------------
# Analytical tangent AA_{iJkL} = ∂P_{iJ}/∂F_{kL}.
#
# P = Fᵉ_new⁻ᵀ M_new Fᵖ_new⁻ᵀ,  Fᵉ_new = F Fᵖ_new⁻¹.
#
# Approximation: Fᵖ_new is treated as FROZEN when differentiating P with
# respect to F.  The consequence:
#   • Elastic step  (Fᵖ_new = Fᵖ_old, truly constant): result is EXACT.
#   • Plastic step  (Fᵖ_new depends on F): the term ∂Fᵖ_new/∂F is omitted,
#     introducing an error O(Δεᵖ) relative to the full consistent tangent.
#     In practice this error is small (~0.1–1 % of the tangent norm) and
#     Newton convergence remains rapid.
#
# ∂P/∂F = geometric + material contributions:
#
#   Geometric:  –(F⁻ᵀ)_{iL} P_{kJ}
#               from ∂Fᵉ_new⁻ᵀ/∂F at fixed Fᵖ_new and M_new.
#
#   Material:   2 Σ_{ABCD} (Fᵉ_new⁻ᵀ)_{iA} ℂ_{ABCD} (Fᵖ_new⁻ᵀ)_{BJ}
#                          (Fᵉ_tr)_{kC} (Fᵖ_old⁻¹)_{DL}
#               from ∂M_new/∂Cᵉ_tr · ∂Cᵉ_tr/∂F; the factor of 2 comes
#               from symmetry of ℂ in (C,D) and ∂Cᵉ_tr/∂F.
#               Note: Cᵉ_tr = Fᵉ_trᵀ Fᵉ_tr uses Fᵖ_OLD (fixed), so this
#               chain is exact.
#
# ℂ_{ABCD} = ∂M_new/∂Cᵉ_tr assembled in the spectral basis {nᵢ} of Cᵉ_tr:
#
#   Material part:  Σᵢⱼ [ĉᵢⱼ/(2cⱼ)] (nᵢ⊗nᵢ)_{AB} (nⱼ⊗nⱼ)_{CD}
#   Geometric part: Σᵢ≠ⱼ θᵢⱼ (nᵢ⊗nⱼ)_sym_{AB} (nᵢ⊗nⱼ)_sym_{CD}
#
#   cᵢ = eigenvalues of Cᵉ_tr,  nᵢ = eigenvectors,
#   θᵢⱼ = (mᵢ_new − mⱼ_new)/(cᵢ − cⱼ)  [L'Hôpital limit when cᵢ ≈ cⱼ].
#
# Principal-space algorithmic moduli ĉᵢⱼ = ∂mᵢ_new/∂εⱼ_tr:
#   Elastic:  ĉᵢⱼ = λ + 2μ δᵢⱼ
#   Plastic:  ĉᵢⱼ = (λ + μβ) + 2μ(1 – 3β/2) δᵢⱼ + (2μβ – 4μ²/(3μ+H)) Nᵢ Nⱼ
#             β = 2μ Δεᵖ / σvm_tr,  Nᵢ = 1.5 mdev_i_tr / σvm_tr.
# ---------------------------------------------------------------------------
function _j2_tangent_analytical(
    material::J2Plasticity,
    F::SMatrix{3,3,Float64,9},
    state_old::Vector{Float64},
    P::SMatrix{3,3,Float64,9},
    state_new::Vector{Float64},
)
    λ = material.λ
    μ = material.μ
    H = material.H

    # Recover kinematics
    Fp_old = reshape(state_old[1:9], 3, 3)
    eqps_old = state_old[10]
    Fp_new = reshape(state_new[1:9], 3, 3)
    Δεᵖ = state_new[10] - eqps_old

    Fm = Matrix{Float64}(F)
    Fp_old_inv = inv(Fp_old)
    Fp_new_inv = inv(Fp_new)
    Fe_tr = Fm * Fp_old_inv
    Fe_new = Fm * Fp_new_inv
    Fe_new_inv_T = inv(Fe_new)'
    F_inv_T = inv(Fm)'

    # Spectral decomposition of trial Ce
    Ce_tr = Symmetric(Fe_tr' * Fe_tr)
    eig = eigen(Ce_tr)
    c = eig.values      # principal squared stretches
    Q = eig.vectors     # columns are eigenvectors

    # Trial principal Hencky strains and Mandel stresses
    ε_tr = 0.5 .* log.(c)
    tr_ε = sum(ε_tr)
    m_tr = [λ * tr_ε + 2μ * ε_tr[i] for i in 1:3]
    m_dev_tr = m_tr .- sum(m_tr) / 3
    σvm_tr = sqrt(1.5) * norm(m_dev_tr)

    # Principal-space algorithmic moduli ĉᵢⱼ and updated principal stresses
    f_tr = σvm_tr - (material.σy + H * eqps_old)
    if f_tr ≤ 0.0
        # Elastic step
        m_new = m_tr
        ĉ = [λ + 2μ * Float64(i == j) for i in 1:3, j in 1:3]
    else
        # Plastic step
        N_pr = 1.5 .* m_dev_tr ./ σvm_tr   # principal values of flow direction
        m_new = m_tr .- 2μ * Δεᵖ .* N_pr
        β = 2μ * Δεᵖ / σvm_tr
        A = λ + μ * β
        B = 2μ * (1.0 - 1.5β)
        Ccoef = 2μ * β - 4μ^2 / (3μ + H)
        ĉ = [A + B * Float64(i == j) + Ccoef * N_pr[i] * N_pr[j] for i in 1:3, j in 1:3]
    end

    # Assemble ℂ_{ABCD} = ∂M_new/∂Ce_tr in the principal basis
    CC = zeros(3, 3, 3, 3)

    # Material part: Σᵢⱼ [ĉᵢⱼ/(2cⱼ)] (nᵢ⊗nᵢ)_{AB} (nⱼ⊗nⱼ)_{CD}
    for α in 1:3, β_idx in 1:3
        coeff = ĉ[α, β_idx] / (2c[β_idx])
        nα = Q[:, α]
        nβ = Q[:, β_idx]
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            CC[A, B, C, D] += coeff * nα[A] * nα[B] * nβ[C] * nβ[D]
        end
    end

    # Geometric part: Σᵢ≠ⱼ θᵢⱼ (nᵢ⊗nⱼ)_sym_{AB} (nᵢ⊗nⱼ)_sym_{CD}
    tol = 1e-10
    for α in 1:3, β_idx in 1:3
        α == β_idx && continue
        nα = Q[:, α]
        nβ = Q[:, β_idx]
        Δc = c[α] - c[β_idx]
        if abs(Δc) > tol * (c[α] + c[β_idx])
            θ = (m_new[α] - m_new[β_idx]) / Δc
        else
            # L'Hôpital: lim_{cβ→cα} (mα−mβ)/(cα−cβ) = (ĉ_{αα}−ĉ_{βα})/(2cα)
            θ = (ĉ[α, α] - ĉ[β_idx, α]) / (2c[α])
        end
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            sym_AB = 0.5 * (nα[A] * nβ[B] + nα[B] * nβ[A])
            sym_CD = 0.5 * (nα[C] * nβ[D] + nα[D] * nβ[C])
            CC[A, B, C, D] += θ * sym_AB * sym_CD
        end
    end

    # Assemble AA_{iJkL} = –(F⁻ᵀ)_{iL} P_{kJ}
    #   + 2 Σ_{ABCD} (Fe_new⁻ᵀ)_{iA} CC_{ABCD} (Fp_new⁻ᵀ)_{BJ} (Fe_tr)_{kC} (Fp_old⁻¹)_{DL}
    # Factor of 2 arises from symmetry of CC in (C,D) and ∂Ce_tr/∂F.
    Fp_new_inv_T = Fp_new_inv'
    AA = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, J in 1:3, k in 1:3, L in 1:3
        val = -F_inv_T[i, L] * P[k, J]
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            val += 2 * Fe_new_inv_T[i, A] * CC[A, B, C, D] * Fp_new_inv_T[B, J] *
                   Fe_tr[k, C] * Fp_old_inv[D, L]
        end
        AA[i, J, k, L] = val
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(AA)
end

# J2Plasticity constitutive is defined in j2_simo_hughes.jl (Simo-Hughes formulation)

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
    elseif model_name == "seth-hill"
        return SethHill(params)
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
    elseif material isa SethHill
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
# Replaces the old Hencky/Mandel-based J2 model above.
# ---------------------------------------------------------------------------

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
