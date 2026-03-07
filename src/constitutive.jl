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

function second_from_fourth(AA::SArray{Tuple{3,3,3,3},Float64,4})
    # Reshape the 3x3x3x3 tensor to 9x9 directly
    return SMatrix{9,9,Float64,81}(reshape(AA, 9, 9)')
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
function _j2_tangent(
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

function constitutive(material::J2Plasticity, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64})
    W, P, state_new = _j2_stress(material, F, state_old)
    AA = _j2_tangent(material, F, state_old, P)
    return W, P, AA, state_new
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
