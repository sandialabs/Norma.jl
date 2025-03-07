# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

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
            error("Two elastic constants are required but only elastic modulus found")
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
            error("Two elastic constants are required but only Poisson's ratio found")
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
            error("Two elastic constants are required but only bulk modulus found")
        end
    elseif haskey(params, "Lamé's first constant") == true
        λ = params["Lamé's first constant"]
        if haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = μ * (3λ + 2μ) / (λ + μ)
            ν = λ / 2(λ + μ)
            κ = λ + 2μ / 3
        else
            error("Two elastic constants are required but only Lamé's first constant found")
        end
    elseif haskey(params, "shear modulus") == true
        error("Two elastic constants are required but only shear modulus found")
    else
        error("Two elastic constants are required but none found")
    end
    return E, ν, κ, λ, μ
end

mutable struct SaintVenant_Kirchhoff <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function SaintVenant_Kirchhoff(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Linear_Elastic <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Linear_Elastic(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Neohookean <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Neohookean(params::Parameters)
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        return new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct SethHill <: Solid
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
        ρ = params["density"]
        return new(E, ν, κ, λ, μ, ρ, params["m"], params["n"])
    end
end

mutable struct J2 <: Solid
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
        ρ = params["density"]
        Y₀ = params["yield stress"]
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

function stress_update(
    material::J2, F::Matrix{Float64}, Fᵖ::Matrix{Float64}, εᵖ::Float64, Δt::Float64
)
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
        H =
            hardening_rate(material, εᵖ + Δεᵖ) +
            viscoplastic_hardening_rate(material, Δεᵖ, Δt)
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
        println("J2 stress update did not converge to specified tolerance")
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

struct Linear_Isotropic <: Thermal
    κ::Float64
    function Linear_Isotropic(params::Parameters)
        κ = params["thermal conductivity"]
        return new(κ)
    end
end

function odot(A::Matrix{Float64}, B::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a in 1:n
        for b in 1:n
            for c in 1:n
                for d in 1:n
                    C[a, b, c, d] = A[a, c] * B[b, d] + A[a, d] * B[b, c]
                end
            end
        end
    end
    C = 0.5 * C
    return C
end

function ox(A::Matrix{Float64}, B::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a in 1:n
        for b in 1:n
            for c in 1:n
                for d in 1:n
                    C[a, b, c, d] = A[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function oxI(A::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a in 1:n
        for b in 1:n
            for c in 1:n
                for d in 1:n
                    C[a, b, c, d] = A[a, b] * I[c, d]
                end
            end
        end
    end
    return C
end

function Iox(B::Matrix{Float64})
    n, _ = size(B)
    C = zeros(n, n, n, n)
    for a in 1:n
        for b in 1:n
            for c in 1:n
                for d in 1:n
                    C[a, b, c, d] = I[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function convect_tangent(CC::Array{Float64}, S::Matrix{Float64}, F::Matrix{Float64})
    n, _ = size(F)
    AA = zeros(n, n, n, n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                for l in 1:n
                    s = 0.0
                    for p in 1:n
                        for q in 1:n
                            s = s + F[i, p] * CC[p, j, l, q] * F[k, q]
                        end
                    end
                    AA[i, j, k, l] = S[l, j] * I[i, k] + s
                end
            end
        end
    end
    return AA
end

function second_from_fourth(AA::Array{Float64})
    n, _, _, _ = size(AA)
    A = zeros(n * n, n * n)
    for i in 1:n
        for j in 1:n
            p = n * (i - 1) + j
            for k in 1:n
                for l in 1:n
                    q = n * (k - 1) + l
                    A[p, q] = AA[i, j, k, l]
                end
            end
        end
    end
    return A
end

function constitutive(material::SaintVenant_Kirchhoff, F::Matrix{Float64})
    C = F' * F
    E = 0.5 * (C - I)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * trE * trE + μ * tr(E * E)
    S = λ * trE * I + 2.0 * μ * E
    CC = zeros(3, 3, 3, 3)
    for i in 1:3
        for j in 1:3
            δᵢⱼ = I[i, j]
            for k in 1:3
                δᵢₖ = I[i, k]
                δⱼₖ = I[j, k]
                for l in 1:3
                    δᵢₗ = I[i, l]
                    δⱼₗ = I[j, l]
                    δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::Matrix{Float64})
    ∇u = F - I
    ϵ = MiniTensor.symm(∇u)
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * trϵ * trϵ + μ * tr(ϵ * ϵ)
    σ = λ * trϵ * I + 2.0 * μ * ϵ
    CC = zeros(3, 3, 3, 3)
    for i in 1:3
        for j in 1:3
            δᵢⱼ = I[i, j]
            for k in 1:3
                δᵢₖ = I[i, k]
                δⱼₖ = I[j, k]
                for l in 1:3
                    δᵢₗ = I[i, l]
                    δⱼₗ = I[j, l]
                    δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    return W, σ, CC
end

function constitutive(material::Neohookean, F::Matrix{Float64})
    C = F' * F
    J2 = det(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J2 - 1) * IC
    Sdev = μ * Jm23 * (I - IC * trC / 3)
    S = Svol + Sdev
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3
    CCvol = κ * (J2 * ICxIC - (J2 - 1) * ICoIC)
    CCdev = μJ2n * (trC * (ICxIC / 3 + ICoIC) - oxI(IC) - Iox(IC))
    CC = CCvol + CCdev
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function constitutive(material::SethHill, F::Matrix{Float64})
    C = F' * F
    F⁻ᵀ = inv(F)'
    J = det(F)
    Jᵐ = J^material.m
    J⁻ᵐ = 1.0 / Jᵐ
    J²ᵐ = Jᵐ * Jᵐ
    J⁻²ᵐ = 1.0 / J²ᵐ
    Cbar = J^(-2 / 3) * C
    Cbarⁿ = Cbar^material.n
    Cbar⁻ⁿ = inv(Cbarⁿ)
    Cbar²ⁿ = Cbarⁿ * Cbarⁿ
    Cbar⁻²ⁿ = Cbar⁻ⁿ * Cbar⁻ⁿ
    trCbarⁿ = tr(Cbarⁿ)
    trCbar⁻ⁿ = tr(Cbar⁻ⁿ)
    trCbar²ⁿ = tr(Cbar²ⁿ)
    trCbar⁻²ⁿ = tr(Cbar⁻²ⁿ)
    Wbulk = material.κ / 4 / material.m^2 * ((Jᵐ - 1)^2 + (J⁻ᵐ - 1)^2)
    Wshear =
        material.μ / 4 / material.n^2 *
        (trCbar²ⁿ + trCbar⁻²ⁿ - 2 * trCbarⁿ - 2 * trCbar⁻ⁿ + 6)
    W = Wbulk + Wshear
    Pbulk = material.κ / 2 / material.m * (J²ᵐ - Jᵐ - J⁻²ᵐ + J⁻ᵐ) * F⁻ᵀ
    Pshear =
        material.μ / material.n * (
            1 / 3 * (-trCbar²ⁿ + trCbarⁿ + trCbar⁻²ⁿ - trCbar⁻ⁿ) * F⁻ᵀ +
            F⁻ᵀ * (Cbar²ⁿ - Cbarⁿ - Cbar⁻²ⁿ + Cbar⁻ⁿ)
        )
    P = Pbulk + Pshear
    AA = zeros(3, 3, 3, 3)
    return W, P, AA
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
    else
        error("Unknown material model : ", model_name)
    end
end

function get_kinematics(material::Solid)
    if typeof(material) == Linear_Elastic
        return Infinitesimal
    elseif typeof(material) == SaintVenant_Kirchhoff
        return Finite
    elseif typeof(material) == Neohookean
        return Finite
    elseif typeof(material) == SethHill
        return Finite
    end
    return error("Unknown material model : ", typeof(material))
end

function get_p_wave_modulus(material::Solid)
    return material.λ + 2.0 * material.μ
end
