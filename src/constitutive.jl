# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using StaticArrays

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

mutable struct PlasticLinearHardeningInfinitessimal <: Inelastic
    # Slow, but simple J2 plasticity with linear hardening and isotropic elasticity
    # Requires Lame Constants, Initial Yield, Hardening Modulus
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    σy::Float64
    H::Float64
    function PlasticLinearHardeningInfinitessimal(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        σy = get(params, "initial yield", 0.0)
        H = get(params, "hardening modulus", 0.0)
        return new(E, ν, κ, λ, μ, ρ, σy, H)
    end
end

@inline number_states(::PlasticLinearHardeningInfinitessimal) = 7
@inline initial_state(::PlasticLinearHardeningInfinitessimal) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

mutable struct PlasticLinearHardening <: Inelastic
    # Slow, but simple J2 plasticity with linear hardening and isotropic elasticity
    # Requires Lame Constants, Initial Yield, Hardening Modulus
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    σy::Float64
    H::Float64
    function PlasticLinearHardening(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        σy = get(params, "initial yield", 0.0)
        H = get(params, "hardening modulus", 0.0)
        return new(E, ν, κ, λ, μ, ρ, σy, H)
    end
end

@inline number_states(::PlasticLinearHardening) = 7
@inline initial_state(::PlasticLinearHardening) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,               # 1-7: Voigt Ep and EQPS
                                                   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,     # 8-16: F_(n-1)
                                                   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,     # 17-25: Elastic left Cauchy-Green
                                                    ] 

mutable struct J2 <: Inelastic
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    Y₀::Float64
    n::Float64
    ε₀::Float64
    Sᵥᵢₛ₀::Float64
    m::Float64
    ∂ε∂t₀::Float64
    Cₚ::Float64
    β::Float64
    T₀::Float64
    Tₘ::Float64
    M::Float64
    function J2(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = get(params, "density", 0.0)
        Y₀ = get(params, "yield stress", 0.0)
        n = get(params, "hardening exponent", 0.0)
        ε₀ = get(params, "reference plastic strain", 0.0)
        Sᵥᵢₛ₀ = get(params, "reference viscoplastic stress", 0.0)
        m = get(params, "rate dependence exponent", 0.0)
        ∂ε∂t₀ = get(params, "reference plastic strain rate", 0.0)
        Cₚ = get(params, "specific heat capacity", 0.0)
        β = get(params, "Taylor-Quinney coefficient", 0.0)
        T₀ = get(params, "reference temperature", 0.0)
        Tₘ = get(params, "melting temperature", 0.0)
        M = get(params, "thermal softening exponent", 0.0)
        κ = E / (1.0 - 2.0 * ν) / 3.0
        μ = E / (1.0 + ν) / 2.0
        return new(E, ν, κ, λ, μ, ρ, Y₀, n, ε₀, Sᵥᵢₛ₀, m, ∂ε∂t₀, Cₚ, β, T₀, Tₘ, M)
    end
end

@inline number_states(::J2) = 10
@inline initial_state(::J2) = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

function temperature_multiplier(material::J2, T::Float64)
    T₀ = material.T₀
    Tₘ = material.Tₘ
    M = material.M
    return M > 0.0 ? 1.0 - ((T - T₀) / (Tₘ - T₀))^M : 1.0
end

function hardening_potential(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    exponent = (1.0 + n) / n
    return n > 0.0 ? Y₀ * ε₀ / exponent * ((1.0 + ε / ε₀)^exponent - 1.0) : Y₀ * ε
end

function hardening_rate(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    exponent = (1.0 - n) / n
    return n > 0.0 ? Y₀ / ε₀ / n * (1.0 + ε / ε₀)^exponent : 0.0
end

function flow_strength(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    return n > 0.0 ? Y₀ * (1.0 + ε / ε₀)^(1.0 / n) : Y₀
end

function viscoplastic_dual_kinetic_potential(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    exponent = (1.0 + m) / m
    return if Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0
        Δt * Sᵥᵢₛ₀ * ∂ε∂t₀ / exponent * (Δε / Δt / ∂ε∂t₀)^exponent
    else
        0.0
    end
end

function viscoplastic_stress(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    return if Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0
        Sᵥᵢₛ₀ / ∂ε∂t₀ / Δt / m * (Δε / Δt / ∂ε∂t₀)^((1.0 - m) / m)
    else
        0.0
    end
end

function viscoplastic_hardening_rate(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    return Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0 ? Sᵥᵢₛ₀ * (Δε / Δt / ∂ε∂t₀)^(1.0 / m) : 0.0
end

function vol(A::Matrix{Float64})
    return tr(A) * I(3) / 3.0
end

function dev(A::Matrix{Float64})
    return A - vol(A)
end

function stress_update(material::J2, F::Matrix{Float64}, Fᵖ::Matrix{Float64}, εᵖ::Float64, Δt::Float64)
    max_rma_iter = 64
    max_ls_iter = 64

    κ = material.κ
    μ = material.μ
    λ = material.λ
    J = det(F)

    Fᵉ = F * inv(Fᵖ)
    Cᵉ = Fᵉ' * Fᵉ
    Eᵉ = 0.5 * log(Cᵉ)
    M = λ * tr(Eᵉ) * I(3) + 2.0 * μ * Eᵉ
    Mᵈᵉᵛ = dev(M)
    σᵛᵐ = sqrt(1.5) * norm(Mᵈᵉᵛ)
    σᵛᵒˡ = κ * vol(Eᵉ)

    Y = flow_strength(material, εᵖ)
    r = σᵛᵐ - Y
    r0 = r

    Δεᵖ = 0.0
    r_tol = 1e-10
    Δεᵖ_tol = 1e-10

    rma_iter = 0
    rma_converged = r ≤ r_tol
    while rma_converged == false
        if rma_iter == max_rma_iter
            break
        end
        Δεᵖ₀ = Δεᵖ
        merit_old = r * r
        H = hardening_rate(material, εᵖ + Δεᵖ) + viscoplastic_hardening_rate(material, Δεᵖ, Δt)
        ∂r = -3.0 * μ - H
        δεᵖ = -r / ∂r

        # line search
        ls_iter = 0
        α = 1.0
        backtrack_factor = 0.1
        decrease_factor = 1.0e-05
        ls_converged = false
        while ls_converged == false
            if ls_iter == max_ls_iter
                # line search has failed to satisfactorily improve newton step
                # just take the full newton step and hope for the best
                α = 1
                break
            end
            ls_iter += 1
            Δεᵖ = max(Δεᵖ₀ + α * δεᵖ, 0.0)
            Y = flow_strength(material, εᵖ + Δεᵖ) + viscoplastic_stress(material, Δεᵖ, Δt)
            r = σᵛᵐ - 3.0 * μ * Δεᵖ - Y

            merit_new = r * r
            decrease_tol = 1.0 - 2.0 * α * decrease_factor
            if merit_new <= decrease_tol * merit_old
                merit_old = merit_new
                ls_converged = true
            else
                α₀ = α
                α = α₀ * α₀ * merit_old / (merit_new - merit_old + 2.0 * α₀ * merit_old)
                if backtrack_factor * α₀ > α
                    α = backtrack_factor * α₀
                end
            end
        end
        rma_converged = abs(r / r0) < r_tol || Δεᵖ < Δεᵖ_tol
        rma_iter += 1
    end
    if rma_converged == false
        norma_log(2, :warning, "J2 stress update did not converge to specified tolerance")
    end

    Nᵖ = σᵛᵐ > 0.0 ? 1.5 * Mᵈᵉᵛ / σᵛᵐ : zeros(3, 3)
    ΔFᵖ = exp(Δεᵖ * Nᵖ)
    Fᵖ = ΔFᵖ * Fᵖ
    εᵖ += Δεᵖ

    ΔEᵉ = Δεᵖ * Nᵖ
    Mᵈᵉᵛ -= 2.0 * μ * ΔEᵉ
    σᵛᵐ = sqrt(1.5) * norm(Mᵈᵉᵛ)
    σᵈᵉᵛ = inv(Fᵉ)' * Mᵈᵉᵛ * Fᵉ' / J
    σ = σᵈᵉᵛ + σᵛᵒˡ

    eʸ = (σᵛᵐ - Y) / Y
    Fᵉ = F * inv(Fᵖ)
    Cᵉ = Fᵉ' * Fᵉ
    Eᵉ = 0.5 * log(Cᵉ)
    M = λ * tr(Eᵉ) * I(3) + 2.0 * μ * Eᵉ
    eᴹ = norm(Mᵈᵉᵛ - dev(M)) / norm(Mᵈᵉᵛ)
    return Fᵉ, Fᵖ, εᵖ, σ
end

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
        C[:, :, a, a] .= A  # Fill diagonal blocks directly
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function Iox(B::SMatrix{3,3,Float64,9})
    C = zeros(MArray{Tuple{3,3,3,3},Float64})
    for a in 1:3
        C[a, a, :, :] .= B  # Fill diagonal blocks directly
    end
    return SArray{Tuple{3,3,3,3}}(C)
end

function convect_tangent(CC::SArray{Tuple{3,3,3,3},Float64}, S::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9})
    # Pre-allocate the 4D output as mutable static array
    AA = MArray{Tuple{3,3,3,3},Float64}(undef)

    # Identity matrix for 3D
    I_n = @SMatrix [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]

    for j in 1:3
        for l in 1:3
            # Extract slice M[p,q] = CC[p, j, l, q]
            M = @SMatrix [
                CC[1, j, l, 1] CC[1, j, l, 2] CC[1, j, l, 3]
                CC[2, j, l, 1] CC[2, j, l, 2] CC[2, j, l, 3]
                CC[3, j, l, 1] CC[3, j, l, 2] CC[3, j, l, 3]
            ]

            # Compute G = F * M * Fᵀ
            G = F * M * F'

            # Fill the result tensor
            AA[:, j, :, l] .= S[l, j] .* I_n .+ G
        end
    end
    return SArray{Tuple{3,3,3,3}}(AA)  # Convert to immutable for better efficiency
end

function second_from_fourth(AA::SArray{Tuple{3,3,3,3},Float64,4})
    # Reshape the 3x3x3x3 tensor to 9x9 directly
    return SMatrix{9,9,Float64,81}(reshape(AA, 9, 9)')
end

const I3 = @SMatrix [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
]

function constitutive(material::SaintVenant_Kirchhoff, F::SMatrix{3,3,Float64,9})
    C = F' * F
    E = 0.5 .* (C - I3)

    λ = material.λ
    μ = material.μ

    trE = tr(E)
    # Strain energy
    W = 0.5 * λ * (trE^2) + μ * tr(E * E)

    # 2nd Piola-Kirchhoff stress
    S = λ * trE .* I3 .+ 2.0 .* μ .* E

    # 4th-order elasticity tensor CC
    # Build it in an MArray, then convert to SArray
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

    # 1st Piola-Kirchhoff stress
    P = F * S

    # Convert the 4th-order tensor for large-deformation convect_tangent
    AA = convect_tangent(CC_s, S, F)

    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::SMatrix{3,3,Float64,9})
    ∇u = F - I3
    ϵ = 0.5 .* (∇u + ∇u')

    λ = material.λ
    μ = material.μ

    trϵ = tr(ϵ)
    # Strain energy
    W = 0.5 * λ * (trϵ^2) + μ * tr(ϵ * ϵ)

    σ = λ * trϵ .* I3 .+ 2.0 .* μ .* ϵ

    # 4th-order elasticity tensor
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

    # Decompose energy: volumetric + deviatoric
    Wvol = 0.25 * κ * (J2 - log(J2) - 1.0)
    Wdev = 0.5 * μ * (Jm23 * trC - 3.0)
    W = Wvol + Wdev

    # Inverse of C
    IC = inv(C)

    # S = Svol + Sdev
    Svol = 0.5 * κ * (J2 - 1.0) .* IC
    Sdev = μ .* Jm23 .* (I3 .- (IC .* (trC / 3.0)))
    S = Svol .+ Sdev

    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3.0

    CCvol = κ .* (J2 .* ICxIC .- (J2 - 1.0) .* ICoIC)
    CCdev = μJ2n .* (trC .* (ICxIC ./ 3 .+ ICoIC) .- oxI(IC) .- Iox(IC))

    CC = CCvol .+ CCdev

    # 1st Piola stress
    P = F * S
    # Large-deformation tangent
    AA = convect_tangent(CC, S, F)

    return W, P, AA
end

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
    AA_m = zeros(MArray{Tuple{3,3,3,3},Float64})
    AA = SArray{Tuple{3,3,3,3}}(AA_m)
    return W, P, AA
end

function major_symmetric_part_fourth_order(tensor::MArray{Tuple{3,3,3,3},Float64})

    n1, n2, n3, n4 = size(tensor)
    symmetric_tensor = zeros(n1, n2, n3, n4)

    # Loop over all indices
    for i in 1:n1, j in 1:n2, k in 1:n3, l in 1:n4
        # Average over major symmetry (swap ij and kl)
        symmetric_tensor[i, j, k, l] = (tensor[i, j, k, l] + tensor[k, l, i, j]) / 2
    end

    return symmetric_tensor
end


function constitutive(material::PlasticLinearHardening, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64})
    state_new = copy(state_old)
    state_new[8:16] = vec(F)
    F⁻¹ = inv(F)
    F⁻ᵀ = transpose(F⁻¹)

    # Simo/Hughes Box 9.1-9.2
    ROOT23 = sqrt(2/3)

    λ = material.λ
    μ = material.μ
    κ = material.κ
    H = material.H

    Ep_old_voigt = state_old[1:6]
    Ep_old = zeros(Float64, 3, 3)
    
    # Map Voigt components to the strain tensor
    Ep_old[1, 1] = Ep_old_voigt[1]  # ε_xx
    Ep_old[2, 2] = Ep_old_voigt[2]  # ε_yy
    Ep_old[3, 3] = Ep_old_voigt[3]  # ε_zz
    Ep_old[2, 3] = Ep_old_voigt[4] / 2  # ε_yz
    Ep_old[3, 2] = Ep_old_voigt[4] / 2  # ε_yz
    Ep_old[1, 3] = Ep_old_voigt[5] / 2  # ε_zx
    Ep_old[3, 1] = Ep_old_voigt[5] / 2  # ε_zx
    Ep_old[1, 2] = Ep_old_voigt[6] / 2  # ε_xy
    Ep_old[2, 1] = Ep_old_voigt[6] / 2  # ε_xy

    eqps_old = state_old[7]

    # Map Deformation Gradient from Last step
    Fn_old = reshape(state_old[8:16], 3, 3)
    #println("Last Def Grad", Fn_old)
    #println("Attempt F", F)
    # Relative def grad
    f_n = F * inv(Fn_old)
    #println("Relative Def Grad", f_n)
    
    # Map elastic left Cauchy-Green from last step
    bᵉ = reshape(state_old[17:25], 3, 3)
    println("Last left CG tensor", bᵉ)

    # Compute elastic predictors
    f̄ₙ = det(f_n)^(-1/3) * f_n
    b̄ᵉₙ =  f̄ₙ * bᵉ * transpose(f̄ₙ)
    b̄ᵉₙ_dev = b̄ᵉₙ - 1/3*tr(b̄ᵉₙ)*I3
    s_trial = μ * b̄ᵉₙ_dev

    s_norm = norm(s_trial)
    f_trial = s_norm - ROOT23 * (material.σy + H * eqps_old)
    N = s_trial / s_norm
    NN = N * N'
    dev_NN = NN - 1/3*tr(NN)*I3

    Īᵉ = 1/3*tr(b̄ᵉₙ)
    μ̄ = Īᵉ * μ
    Jₙ = det(F)
    # Addition of the elastic mean stress
    p = κ/2 * (Jₙ^2 - 1)/Jₙ
    println("pressure", p)


    #euler_strain = 1/2*(I3 - inv(b̄ᵉₙ))
    #kirchhoff_stress = zeros((3,3))
    c_trial = MArray{Tuple{3,3,3,3},Float64}(undef)
    c̄_trial = MArray{Tuple{3,3,3,3},Float64}(undef)

    NxdevN2 = MArray{Tuple{3,3,3,3},Float64}(undef)

    delta_0 = 1 
    f_1 = (1/delta_0)
    delta_1 = 2 * μ̄  
    delta_2 = 2 * s_norm 

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    NxdevN2[i,j,k,l] = N[i,j] * dev_NN[k,l]
                end
            end
        end
    end
    symm_part_C = major_symmetric_part_fourth_order(NxdevN2)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    #c̄_trial[i, j, k, l] =  2*μ̄*(I3[i,k] * I3[j,l] - I3[i,j]*I3[k,l]/3.0) - 2/3*s_norm * ( N[i,j] * I3[k,l] + I3[i,j]*N[k,l] )
                    c_trial[i, j, k, l] =  2*μ̄*(I3[i,j]*I3[k,l] - I3[i,j]*I3[k,l]/3.0 ) - 2/3*( s_trial[i,j] * I3[k,l] + I3[i,j]*s_trial[k,l] ) #+ c̄_trial[i, j, k, l]
                    c_trial[i,j,k,l] -= delta_1 * N[i,j] * N[k,l]

                   
                    c_trial[i,j,k,l] -= delta_2 * symm_part_C[i,j,k,l]
                    # c_trial[i, j, k, l] -= (2*μ̄*(I3[i,j]*I3[k,l] - I3[i,j]*I3[k,l]/3.0 ) - 2/3*( s_trial[i,j] * I3[k,l] + I3[i,j]*s_trial[k,l] ))*2*μ̄*0
                    # From the Book: κ*Jₙ *Jₙ * I3[i,j]*I3[k,l] - 2*κ/2*(Jₙ^2 - 1) + c̄_trial[i, j, k, l]
                    #kirchhoff_stress[i,j] += c̄_trial[i, j, k, l] * euler_strain[k, l]

                    end
            end
        end
    end


    if f_trial <= 0
        τ = Jₙ*p*I3 + s_trial

        state_new[1:6] .= Ep_old_voigt
        state_new[7] = eqps_old

        state_new[17:25] = vec(b̄ᵉₙ_dev + Īᵉ*I3)

 
        

        C = F' * F
        E = 0.5 .* (C - I3)

        trE = tr(E)
        # Strain energy
        W = 0.5 * λ * (trE^2) + μ * tr(E * E)

        S = F⁻¹ * τ * F⁻ᵀ
        P = τ * F⁻ᵀ
        P_from_PK2 = F *S
        println("Alg Stress", P)
        println("Real Stress", P_from_PK2 )
        println("Deviatoric stress", s_trial * F⁻ᵀ)
        println("Dev from PK2", P_from_PK2 - 1/3*tr(P_from_PK2)*I3)
        println("Pressure from PK2 ", 1/3*tr(P_from_PK2))
        C_trial = MArray{Tuple{3,3,3,3},Float64}(undef)

        for PP in 1:3
            for QQ in 1:3
                for RR in 1:3
                    for SS in 1:3
                        C_trial[PP, QQ, RR, SS] =  sum(
                             c_trial[i,j,k,l] * F⁻¹[PP, i] * F⁻¹[QQ, j] * F⁻¹[RR, k] * F⁻¹[SS, l]
                             for i in 1:3, j in 1:3, k in 1:3, l in 1:3
                         )
                    end
                end
            end
        end
        
        CC = SArray{Tuple{3,3,3,3}}(C_trial)
        AA = convect_tangent(CC, S, F)
        return W, P, AA, state_new
    end
    norma_abort("No plasticity for now.")
    # Linear hardening allows for analytical solution of delta gamma
    Δγ = (f_trial/2/μ̄)/(1 + κ/3/μ̄)
    

    s_n = s_trial - 2*μ̄*Δγ*N
    eqps = eqps_old + ROOT23*Δγ
    Ep = Ep_old + Δγ * N

    
    # Kirchhoff stress
    τ = Jₙ*p*I3 + s_n
    P = τ * transpose(F⁻¹)
    b̄ᵉₙ = s_n/μ + Īᵉ*I3 
    # Update intermediate configuration 
    state_new[17:25] = vec(b̄ᵉₙ)

    # Break Ep into voigt notation
    EpV = zeros(Float64, 6)
    EpV[1] = Ep[1, 1]
    EpV[2] = Ep[2, 2]
    EpV[3] = Ep[3, 3]
    EpV[4] = Ep[2, 3] * 2
    EpV[5] = Ep[1, 3] * 2
    EpV[6] = Ep[1, 2] * 2
    state_new[1:6] .= EpV
    state_new[7] = eqps

    W = 0.5*κ*(0.5*(Jₙ^2 -1) - log(Jₙ)) +  0.5 * μ * (tr(b̄ᵉₙ) - 3) 

    # tangent modulus
    β₀ = 1 + H/3/μ̄
    β₁ = 2*μ̄*Δγ / s_norm
    β₂ = (1 - 1/β₀)*2/3*s_norm/μ̄*Δγ
    β₃ = 1/β₀ - β₁ + β₂
    β₄ = (1/β₀ - β₁) * s_norm/μ̄

    

    CC_m = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    comp_1 = N[i,j] * dev_NN[k,l]
                    comp_2 = N[j,i] * dev_NN[k,l]
                    comp_3 = N[i,j] * dev_NN[l,k]
                    comp_4 = N[k,l] * dev_NN[i,j]
                    symm_comp = 1/4*(comp_1 + comp_2 + comp_3 + comp_4)
                    CC_m[i, j, k, l] = C_trial[i,j,k,l] - β₁*C̄_trial[i, j, k, l] - 2*μ̄*β₃*N[i,j]*N[k,l] - 2*μ̄*β₄*symm_comp
                end
            end
        end
    end
    CC_s = SArray{Tuple{3,3,3,3}}(CC_m)

    return W, τ, CC_s, state_new
end


function constitutive(material::PlasticLinearHardeningInfinitessimal, F::SMatrix{3,3,Float64,9}, state_old::Vector{Float64})
    # Basic Implementation of Simo Hughes Radial Return 3.3.1
    # Simplified to only support linear isotropic hardening
    # No Kinematic Hardening
    ROOT23 = sqrt(2/3)

    λ = material.λ
    μ = material.μ
    κ = material.κ

    Ep_old_voigt = state_old[1:6]
    Ep_old = zeros(Float64, 3, 3)

    # Map Voigt components to the strain tensor
    Ep_old[1, 1] = Ep_old_voigt[1]  # ε_xx
    Ep_old[2, 2] = Ep_old_voigt[2]  # ε_yy
    Ep_old[3, 3] = Ep_old_voigt[3]  # ε_zz
    Ep_old[2, 3] = Ep_old_voigt[4] / 2  # ε_yz
    Ep_old[3, 2] = Ep_old_voigt[4] / 2  # ε_yz
    Ep_old[1, 3] = Ep_old_voigt[5] / 2  # ε_zx
    Ep_old[3, 1] = Ep_old_voigt[5] / 2  # ε_zx
    Ep_old[1, 2] = Ep_old_voigt[6] / 2  # ε_xy
    Ep_old[2, 1] = Ep_old_voigt[6] / 2  # ε_xy

    eqps_old = state_old[7]

    # Calculate the current step's strain
    ∇u = F - I3
    ϵ = 0.5 .* (∇u + ∇u')
    trϵ = tr(ϵ)

    iso_ϵ = tr(ϵ) / 3 * I3
    # Deviatoric strain increment
    dev_ϵ = ϵ - iso_ϵ
    dev_ϵᵖ = Ep_old - tr(Ep_old)/3 * I3

    # Compute deviatoric trial stress
    trial_stress = 2 * material.μ * (dev_ϵ - dev_ϵᵖ)

    f_trial = norm(trial_stress) - ROOT23 * (material.σy + material.H * eqps_old)

    if f_trial <= 0
        σ = λ * trϵ .* I3 .+ 2.0 .* μ .* ϵ
        W = 0.5 * λ * (trϵ^2) + μ * tr(ϵ * ϵ)
        # 4th-order elasticity tensor
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
        # println("Elastic step", σ)
        state_old[1:6] .= Ep_old_voigt
        state_old[7] = eqps_old
        return W, σ, CC_s, state_old
    end

    # Linear hardening allows for analytical solution of delta gamma
    deltaGamma = (f_trial/(1 + material.H/3/material.μ))/2/material.μ
    eqps = eqps_old + ROOT23*deltaGamma
    N = trial_stress / norm(trial_stress)
    Ep = Ep_old + deltaGamma * N
    σ = κ * tr(ϵ) .* I3 + trial_stress - 2 * material.μ * deltaGamma * N

    state_new = copy(state_old)
    # Break Ep into voigt notation
    EpV = zeros(Float64, 6)
    EpV[1] = Ep[1, 1]
    EpV[2] = Ep[2, 2]
    EpV[3] = Ep[3, 3]
    EpV[4] = Ep[2, 3] * 2
    EpV[5] = Ep[1, 3] * 2
    EpV[6] = Ep[1, 2] * 2
    state_new[1:6] .= EpV
    state_new[7] = eqps

    ϵe = ϵ - Ep
    W = 0.5 * material.λ * (tr(ϵe)^2) + material.μ * tr(ϵe * ϵe)

    # tangent modulus
    θ = 1 - 2 * deltaGamma * μ / norm(trial_stress)
    θ_bar = 1/(1 + material.H/3/μ) - (1 - θ)
    CC_m = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    CC_m[i, j, k, l] =
                        κ * I3[i, j] * I3[k, l] + 2 * μ * θ * (I3[i, k] * I3[j, l] - I3[i, j] * I3[k, l]/3.0) -
                        2 * μ * θ_bar * N[i, j] * N[k, l]
                end
            end
        end
    end
    CC_s = SArray{Tuple{3,3,3,3}}(CC_m)

    return W, σ, CC_s, state_new

end

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
    elseif model_name == "PlasticLinearHardeningInfinitessimal"
        return PlasticLinearHardeningInfinitessimal(params)
    elseif model_name == "PlasticLinearHardening"
        return PlasticLinearHardening(params)
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
    elseif material isa PlasticLinearHardeningInfinitessimal
        return Infinitesimal
    elseif material isa PlasticLinearHardening
        return Finite
    end
    norma_abort("Unknown material model : $(typeof(material))")
    return nothing
end

function get_p_wave_modulus(material::Solid)
    return material.λ + 2.0 * material.μ
end
