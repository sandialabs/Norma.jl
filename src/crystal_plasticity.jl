# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Crystal plasticity, ported from LCM's CrystalPlasticity model
# (src/LCM/models/CrystalPlasticity in sandialabs/LCM), including the
# corrections of sandialabs/LCM#127.
#
# Kinematics. F = Fe Fp. Each slip system has a unit direction s and normal
# n in the lattice frame; with the lattice orientation Q its Schmid tensor in
# the sample frame is P_α = (Q s_α) ⊗ (Q n_α). Over a step the plastic part
# follows the exponential map
#
#     Fp_{n+1} = exp(Σ_α Δγ_α P_α) Fp_n,
#
# which keeps det Fp = 1. The elastic response is St. Venant-Kirchhoff on
# Fe, with the lattice elasticity tensor rotated by Q,
#
#     Ee = ½ (Feᵀ Fe - I),   S = ℂ : Ee,   σ = Fe S Feᵀ / det Fe,
#
# and the resolved shear stress of a system is the Mandel stress projected
# on it, τ_α = (Feᵀ Fe S) : P_α.
#
# Slip. Every system belongs to a slip family with a flow rule, which gives
# the slip rate γ̇_α from τ_α and the slip resistance g_α, and a hardening
# law, which gives the hardening state h_α (the resistance itself, or a
# dislocation density) from the slip rates of the family.
#
# Integration over a step of length Δt. The predictor takes the slip rates
# of the previous step, γ_{n+1} = γ_n + Δt γ̇_n, and hardens with them.
#   explicit:  τ from Fp_n, g from the predicted rates,
#              γ_{n+1} = γ_n + Δt γ̇(τ, g);
#   implicit:  Newton on the slip residual
#              r(γ) = γ - γ_n - Δt γ̇(τ(γ), g),
#              with g hardened from h_n with the rates (γ - γ_n)/Δt, or on
#              the slip and hardness residual, with unknowns (γ, h) and the
#              additional r_h = h - h(γ̇, h). After convergence the hardening
#              is evaluated once more with the converged rates.
# A local failure (no convergence, a slip rate beyond the family's maximum,
# a singular deformation) throws a DomainError, which the model turns into
# a cut of the global step. The tangent is the exact algorithmic one, from
# forward differentiation through the explicit update or through the
# implicit function γ(F) of the implicit schemes.
#
# State of an integration point, for N slip systems:
#     Fp (9, column-major), γ (N), γ̇ (N), h (N), τ (N), eqps (1).

using ForwardDiff
using LinearAlgebra
using StaticArrays

# ---------------------------------------------------------------------------
# Flow rules
# ---------------------------------------------------------------------------

abstract type FlowRule end

struct NoFlow <: FlowRule end

struct PowerLawFlow <: FlowRule
    γ̇0::Float64        # reference slip rate
    m::Float64         # rate exponent
    min_ratio::Float64 # below this τ/g the rate underflows to zero
    max_ratio::Float64 # beyond this τ/g the rate overflows
end

struct PowerLawDragFlow <: FlowRule
    γ̇0::Float64
    m::Float64
    drag::Float64
    min_ratio::Float64
    max_ratio::Float64
end

struct ThermalActivationFlow <: FlowRule
    γ̇0::Float64 # reference slip rate
    F0::Float64 # activation energy, in units of kT
    s_t::Float64 # thermal resistance
    p::Float64
    q::Float64
end

power_law_ratio_range(m::Real) = ((2 * floatmin(Float64))^(0.5 / m), (0.5 * floatmax(Float64))^(0.5 / m))

slip_rate(::NoFlow, τ, g) = zero(τ / g)

# γ̇ = γ̇0 |τ/g|^(m-1) τ/g
function slip_rate(rule::PowerLawFlow, τ, g)
    x = τ / g
    abs(x) > rule.max_ratio && throw(DomainError(x, "Crystal plasticity: τ/g is beyond the range of the power law."))
    abs(x) > rule.min_ratio || return zero(x)
    return rule.γ̇0 * abs(x)^(rule.m - 1) * x
end

# The power law and a linear viscous drag γ̇0 (τ/g)/d combined in series.
function slip_rate(rule::PowerLawDragFlow, τ, g)
    x = τ / g
    abs(x) < rule.min_ratio && return zero(x)
    viscous = x / rule.drag
    abs(x) > rule.max_ratio && return rule.γ̇0 * viscous
    power = abs(x)^(rule.m - 1)
    rate = power * x
    if rule.drag * power > eps(Float64)
        rate = 1 / (1 / rate + 1 / viscous)
    end
    return rule.γ̇0 * rate
end

# γ̇ = γ̇0 exp(-F0 (1 - ξ^p)^q) sign(τ), ξ = (|τ| - g)/s_t. Beyond ξ = 1 the
# barrier vanishes and the rate stays at γ̇0 (LCM returns NaN there).
function slip_rate(rule::ThermalActivationFlow, τ, g)
    ξ = (abs(τ) - g) / rule.s_t
    ξ > 2 * floatmin(Float64) || return zero(ξ)
    ξ = min(ξ, one(ξ))
    return rule.γ̇0 * exp(-rule.F0 * (1 - ξ^rule.p)^rule.q) * sign(τ)
end

# ---------------------------------------------------------------------------
# Hardening laws
# ---------------------------------------------------------------------------

abstract type HardeningLaw end

struct NoHardening <: HardeningLaw
    h0::Float64 # constant slip resistance
end

# ġ = (H - R (g - g0)) Σ|γ̇|, integrated exactly over the step.
struct LinearMinusRecoveryHardening <: HardeningLaw
    H::Float64
    R::Float64
    h0::Float64
end

# ġ_α = rate · 2 Σ_β q_αβ |γ̇_β| (gs - g_α)/(gs - g0), gs = gs0 (Σ|γ̇|/γ̇s)^ω.
struct SaturationHardening <: HardeningLaw
    rate::Float64
    gs0::Float64
    ω::Float64
    γ̇s::Float64
    h0::Float64
end

# ρ̇_α = max(k1 √ρF_α - k2 ρ_α, 0) |γ̇_α|, g_α = c μ b √ρP_α, with the forest
# and parallel densities ρF = A ρ and ρP = B ρ of the family.
struct DislocationDensityHardening <: HardeningLaw
    c::Float64  # geometric factor
    k1::Float64 # generation factor
    k2::Float64 # annihilation factor
    μ::Float64  # shear modulus
    b::Float64  # Burgers vector magnitude
    h0::Float64 # initial dislocation density
end

initial_hardening_state(law::HardeningLaw) = law.h0

struct SlipFamily{FR<:FlowRule,HL<:HardeningLaw}
    flow::FR
    hardening::HL
    systems::Vector{Int}     # indices of the family's systems
    latent::Matrix{Float64}  # interaction matrix, over the family's systems
    aux::Matrix{Float64}     # parallel matrix of the dislocation density law
    max_rate::Float64        # largest admissible |γ̇| (LCM's Max Slip Increment)
end

# Each law writes the new hardening state h_np1 and the slip resistance g of
# its family's systems. h_in is the hardening state the semi-implicit laws
# evaluate their rates with: h_n for the slip residual, the unknown for the
# slip and hardness residual.
function harden!(h_np1, g, ::NoHardening, family::SlipFamily, Δt, γ̇, h_n, h_in; explicit=false)
    for i in family.systems
        h_np1[i] = h_n[i]
        g[i] = h_np1[i]
    end
    return nothing
end

function harden!(h_np1, g, law::LinearMinusRecoveryHardening, family::SlipFamily, Δt, γ̇, h_n, h_in; explicit=false)
    ΔΓ = Δt * sum(abs(γ̇[i]) for i in family.systems)
    for i in family.systems
        if law.R > 0.0
            h_∞ = law.h0 + law.H / law.R
            h_np1[i] = h_∞ + (h_n[i] - h_∞) * exp(-law.R * ΔΓ)
        else
            h_np1[i] = h_n[i] + law.H * ΔΓ
        end
        g[i] = h_np1[i]
    end
    return nothing
end

# With `explicit`, the update is the final forward Euler step of the slip
# residual or of the explicit scheme. It overshoots the saturation stress
# when the step slips more than about (gs - g0)/(2 rate), and diverges beyond
# twice that, so a step that crosses gs is cut.
function harden!(h_np1, g, law::SaturationHardening, family::SlipFamily, Δt, γ̇, h_n, h_in; explicit=false)
    systems = family.systems
    rates = [abs(γ̇[i]) for i in systems]
    driver = 2 * (family.latent * rates)
    Γ̇ = sum(rates)
    if Γ̇ < floatmin(Float64)
        for (k, i) in enumerate(systems)
            h_np1[i] = h_n[i] + Δt * driver[k]
            g[i] = h_np1[i]
        end
        return nothing
    end
    gs = law.ω > 0.0 ? law.gs0 * (Γ̇ / law.γ̇s)^law.ω : law.gs0 + zero(Γ̇)
    for (k, i) in enumerate(systems)
        g[i] = h_in[i]
        # At a saturation stress equal to the initial hardening state the
        # crystal starts saturated and does not harden (the law is 0/0 there).
        h_np1[i] =
            gs == law.h0 ? h_n[i] + zero(gs) : h_n[i] + Δt * law.rate * driver[k] * (gs - h_in[i]) / (gs - law.h0)
        explicit &&
            (gs - h_np1[i]) * (gs - h_n[i]) < 0 &&
            throw(DomainError(h_np1[i], "Crystal plasticity: the hardening overshoots the saturation stress."))
    end
    return nothing
end

function harden!(h_np1, g, law::DislocationDensityHardening, family::SlipFamily, Δt, γ̇, h_n, h_in; explicit=false)
    systems = family.systems
    forest = family.latent * [h_in[i] for i in systems]
    for (k, i) in enumerate(systems)
        forest[k] < 0.0 && throw(DomainError(forest[k], "Crystal plasticity: negative forest dislocation density."))
        generation = law.k1 * sqrt(forest[k])
        annihilation = law.k2 * h_in[i]
        driver = generation > annihilation ? generation - annihilation : zero(generation - annihilation)
        h_np1[i] = h_n[i] + Δt * driver * abs(γ̇[i])
    end
    parallel = family.aux * [h_np1[i] for i in systems]
    for (k, i) in enumerate(systems)
        parallel[k] < 0.0 &&
            throw(DomainError(parallel[k], "Crystal plasticity: negative parallel dislocation density."))
        g[i] = law.c * law.μ * law.b * sqrt(parallel[k])
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Slip systems
# ---------------------------------------------------------------------------

const SQRT3 = sqrt(3.0)

# LCM's tables (src/LCM/utils/python/lcm_preprocess/create_slip_systems.py),
# as (direction, normal) in the lattice frame, unnormalized.
const CUBIC_SLIP_SYSTEMS = Dict{String,Vector{NTuple{2,NTuple{3,Float64}}}}(
    "fcc" => [
        ((-1, 1, 0), (1, 1, 1)),
        ((0, -1, 1), (1, 1, 1)),
        ((1, 0, -1), (1, 1, 1)),
        ((-1, -1, 0), (-1, 1, 1)),
        ((1, 0, 1), (-1, 1, 1)),
        ((0, 1, -1), (-1, 1, 1)),
        ((1, -1, 0), (-1, -1, 1)),
        ((0, 1, 1), (-1, -1, 1)),
        ((-1, 0, -1), (-1, -1, 1)),
        ((1, 1, 0), (1, -1, 1)),
        ((-1, 0, 1), (1, -1, 1)),
        ((0, -1, -1), (1, -1, 1)),
    ],
    "bcc 110" => [
        ((1, 1, -1), (0, 1, 1)),
        ((1, -1, 1), (0, 1, 1)),
        ((-1, 1, 1), (1, 1, 0)),
        ((1, -1, 1), (1, 1, 0)),
        ((-1, 1, 1), (1, 0, 1)),
        ((1, 1, -1), (1, 0, 1)),
        ((1, 1, 1), (0, 1, -1)),
        ((-1, 1, 1), (0, 1, -1)),
        ((1, 1, 1), (-1, 1, 0)),
        ((1, 1, -1), (-1, 1, 0)),
        ((1, 1, 1), (1, 0, -1)),
        ((1, -1, 1), (1, 0, -1)),
    ],
    "bcc 112" => [
        ((1, -1, -1), (2, 1, 1)),
        ((-1, -1, 1), (1, 1, 2)),
        ((-1, 1, -1), (1, 2, 1)),
        ((-1, -1, -1), (-2, 1, 1)),
        ((1, -1, 1), (-1, 1, 2)),
        ((1, 1, -1), (-1, 2, 1)),
        ((1, -1, 1), (2, 1, -1)),
        ((-1, -1, -1), (1, 1, -2)),
        ((-1, 1, 1), (1, 2, -1)),
        ((1, 1, -1), (2, -1, 1)),
        ((-1, 1, 1), (1, -1, 2)),
        ((-1, -1, -1), (1, -2, 1)),
    ],
    "bcc 123" => [
        ((1, -1, -1), (3, 1, 2)),
        ((-1, -1, 1), (1, 2, 3)),
        ((-1, 1, -1), (2, 3, 1)),
        ((-1, -1, 1), (2, 1, 3)),
        ((-1, 1, -1), (1, 3, 2)),
        ((1, -1, -1), (3, 2, 1)),
        ((-1, -1, -1), (-3, 1, 2)),
        ((1, -1, 1), (-1, 2, 3)),
        ((1, 1, -1), (-2, 3, 1)),
        ((1, -1, 1), (-2, 1, 3)),
        ((1, 1, -1), (-1, 3, 2)),
        ((-1, -1, -1), (-3, 2, 1)),
        ((1, -1, 1), (3, 1, -2)),
        ((-1, -1, -1), (1, 2, -3)),
        ((-1, 1, 1), (2, 3, -1)),
        ((-1, -1, -1), (2, 1, -3)),
        ((-1, 1, 1), (1, 3, -2)),
        ((1, -1, 1), (3, 2, -1)),
        ((1, 1, -1), (3, -1, 2)),
        ((-1, 1, 1), (1, -2, 3)),
        ((-1, -1, -1), (2, -3, 1)),
        ((-1, 1, 1), (2, -1, 3)),
        ((-1, -1, -1), (1, -3, 2)),
        ((1, 1, -1), (3, -2, 1)),
    ],
)

# The hexagonal sets, with lattice axis 1 along a1 and axis 3 along c. LCM
# lists the pyramidal <a> systems with the prismatic ones; here they are a
# set of their own. LCM's sixth first-order <c+a> system (its HCP system 18)
# has the direction (1, 0, 2 c/a), which is not orthogonal to its normal;
# (1, 0, c/a) is.
function hexagonal_slip_systems(name::AbstractString, ca::Real)
    s3 = SQRT3
    name == "hcp basal" && return [
        ((1.0, -s3, 0.0), (0.0, 0.0, 1.0)), ((1.0, s3, 0.0), (0.0, 0.0, 1.0)), ((-1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    ]
    name == "hcp prismatic" && return [
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)), ((1.0, s3, 0.0), (-s3, 1.0, 0.0)), ((-1.0, s3, 0.0), (-s3, -1.0, 0.0))
    ]
    name == "hcp pyramidal a" && return [
        ((1.0, 0.0, 0.0), (0.0, -2ca, s3)),
        ((1.0, s3, 0.0), (s3 * ca, -ca, s3)),
        ((-1.0, s3, 0.0), (s3 * ca, ca, s3)),
        ((-1.0, 0.0, 0.0), (0.0, 2ca, s3)),
        ((-1.0, -s3, 0.0), (-s3 * ca, ca, s3)),
        ((1.0, -s3, 0.0), (-s3 * ca, -ca, s3)),
    ]
    name == "hcp pyramidal c+a first order" && return [
        ((1.0, s3, 2ca), (0.0, -2ca, s3)),
        ((-1.0, s3, 2ca), (s3 * ca, -ca, s3)),
        ((-1.0, 0.0, ca), (s3 * ca, ca, s3)),
        ((-1.0, -s3, 2ca), (0.0, 2ca, s3)),
        ((1.0, -s3, 2ca), (-s3 * ca, ca, s3)),
        ((1.0, 0.0, ca), (-s3 * ca, -ca, s3)),
        ((-1.0, s3, 2ca), (0.0, -2ca, s3)),
        ((-1.0, 0.0, ca), (s3 * ca, -ca, s3)),
        ((-1.0, -s3, 2ca), (s3 * ca, ca, s3)),
        ((1.0, -s3, 2ca), (0.0, 2ca, s3)),
        ((1.0, 0.0, ca), (-s3 * ca, ca, s3)),
        ((1.0, s3, 2ca), (-s3 * ca, -ca, s3)),
    ]
    name == "hcp pyramidal c+a second order" && return [
        ((-1.0, s3, 2ca), (ca, -s3 * ca, 2.0)),
        ((-1.0, 0.0, ca), (ca, 0.0, 1.0)),
        ((-1.0, -s3, 2ca), (ca, s3 * ca, 2.0)),
        ((1.0, -s3, 2ca), (-ca, s3 * ca, 2.0)),
        ((1.0, 0.0, ca), (-ca, 0.0, 1.0)),
        ((1.0, s3, 2ca), (-ca, -s3 * ca, 2.0)),
    ]
    return nothing
end

const SLIP_SYSTEM_SET_NAMES = vcat(
    sort(collect(keys(CUBIC_SLIP_SYSTEMS))),
    [
        "hcp basal",
        "hcp prismatic",
        "hcp pyramidal a",
        "hcp pyramidal c+a first order",
        "hcp pyramidal c+a second order",
    ],
)

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

mutable struct CrystalPlasticity <: Solid
    CC::SArray{Tuple{3,3,3,3},Float64,4,81}  # lattice-frame elasticity tensor
    M::Float64                                # largest wave modulus ρ c²
    ρ::Float64
    s::Vector{SVector{3,Float64}}            # unit slip directions, lattice frame
    n::Vector{SVector{3,Float64}}            # unit slip normals, lattice frame
    families::Vector{SlipFamily}
    h_initial::Vector{Float64}               # initial hardening state of each system
    explicit::Bool
    slip_hardness::Bool                      # implicit with the slip and hardness residual
    rate_predictor::Bool
    absolute_tolerance::Float64
    relative_tolerance::Float64
    minimum_iterations::Int
    maximum_iterations::Int
end

function required_parameter(params::AbstractDict, key::String, context::String)
    haskey(params, key) || norma_abort("Crystal plasticity: $context requires `$key`.")
    return Float64(params[key])
end

# Parameters outside their physical range gave misleading failures far from
# the input (a rate exponent m <= 0 inverts the range of the power law, a
# zero slip resistance makes τ/g infinite) or were silently ignored (a
# negative recovery modulus or saturation exponent).
function checked_parameter(value::Float64, key::String, context::String; positive::Bool=true)
    (positive ? value > 0.0 : value >= 0.0) || norma_abortf(
        "Crystal plasticity: `%s` of %s must be %s (got %g).",
        key,
        context,
        positive ? "positive" : "nonnegative",
        value,
    )
    return value
end
positive_parameter(params, key, context) = checked_parameter(required_parameter(params, key, context), key, context)

# A family without a flow rule, or a flow rule or hardening law without a
# type, used to fall back to `none` and silently made the crystal elastic or
# not harden; both are required, and `none` must be said.
function create_flow_rule(params)
    params === nothing && norma_abort(
        "Crystal plasticity: every slip family requires a `flow rule` (with `type: none` for a family that does not slip).",
    )
    haskey(params, "type") || norma_abort(
        "Crystal plasticity: a `flow rule` requires a `type`: power law, power law with drag, thermal activation, or none.",
    )
    type = lowercase(string(params["type"]))
    context = "the $type flow rule"
    if type == "power law"
        m = positive_parameter(params, "rate exponent", context)
        return PowerLawFlow(positive_parameter(params, "reference slip rate", context), m, power_law_ratio_range(m)...)
    elseif type == "power law with drag"
        m = positive_parameter(params, "rate exponent", context)
        return PowerLawDragFlow(
            positive_parameter(params, "reference slip rate", context),
            m,
            positive_parameter(params, "drag coefficient", context),
            power_law_ratio_range(m)...,
        )
    elseif type == "thermal activation"
        return ThermalActivationFlow(
            positive_parameter(params, "reference slip rate", context),
            checked_parameter(
                required_parameter(params, "activation energy", context), "activation energy", context; positive=false
            ),
            positive_parameter(params, "thermal resistance", context),
            positive_parameter(params, "p exponent", context),
            positive_parameter(params, "q exponent", context),
        )
    elseif type == "none"
        return NoFlow()
    end
    return norma_abort(
        "Crystal plasticity: unknown flow rule '$type'; expected power law, power law with drag, thermal activation, or none.",
    )
end

function create_hardening_law(params)
    params === nothing && norma_abort("Crystal plasticity: every slip family requires a `hardening law`.")
    haskey(params, "type") || norma_abort(
        "Crystal plasticity: a `hardening law` requires a `type`: linear minus recovery, saturation, " *
        "dislocation density, or none.",
    )
    type = lowercase(string(params["type"]))
    context = "the $type hardening law"
    h0 = positive_parameter(params, "initial hardening state", context)
    optional(key, default) = checked_parameter(Float64(get(params, key, default)), key, context; positive=false)
    if type == "linear minus recovery"
        return LinearMinusRecoveryHardening(
            required_parameter(params, "hardening modulus", context), optional("recovery modulus", 0.0), h0
        )
    elseif type == "saturation"
        gs0 = required_parameter(params, "initial saturation stress", context)
        gs0 >= h0 || norma_abortf(
            "Crystal plasticity: the initial saturation stress (%g) of a saturation hardening law must be at " *
            "least its initial hardening state (%g); below it the law hardens without bound.",
            gs0,
            h0,
        )
        return SaturationHardening(
            checked_parameter(
                required_parameter(params, "hardening rate", context), "hardening rate", context; positive=false
            ),
            gs0,
            optional("saturation exponent", 0.0),
            checked_parameter(Float64(get(params, "reference slip rate", 1.0)), "reference slip rate", context),
            h0,
        )
    elseif type == "dislocation density"
        nonnegative(key) = checked_parameter(required_parameter(params, key, context), key, context; positive=false)
        return DislocationDensityHardening(
            positive_parameter(params, "geometric factor", context),
            nonnegative("generation factor"),
            nonnegative("annihilation factor"),
            positive_parameter(params, "shear modulus", context),
            positive_parameter(params, "burgers vector magnitude", context),
            h0,
        )
    elseif type == "none"
        return NoHardening(h0)
    end
    return norma_abort(
        "Crystal plasticity: unknown hardening law '$type'; expected linear minus recovery, saturation, " *
        "dislocation density, or none.",
    )
end

# The slip systems of a family as (direction, normal, initial hardening
# state or nothing), from the name of a built-in set or an explicit list.
function family_slip_systems(params::AbstractDict)
    haskey(params, "slip systems") || norma_abort("Crystal plasticity: every slip family requires `slip systems`.")
    systems = params["slip systems"]
    systems isa AbstractString ||
        systems isa AbstractVector ||
        norma_abort(
            "Crystal plasticity: `slip systems` must be the name of a set or a list of systems, each with a `direction` and a `normal`.",
        )
    if systems isa AbstractString
        name = lowercase(systems)
        table = if haskey(CUBIC_SLIP_SYSTEMS, name)
            CUBIC_SLIP_SYSTEMS[name]
        else
            ca = checked_parameter(Float64(get(params, "c/a ratio", sqrt(8 / 3))), "c/a ratio", "a slip family")
            hexagonal_slip_systems(name, ca)
        end
        table === nothing && norma_abort(
            "Crystal plasticity: unknown slip system set '$systems'; expected one of " *
            join(SLIP_SYSTEM_SET_NAMES, ", ") *
            ".",
        )
        return [(SVector{3,Float64}(d), SVector{3,Float64}(n), nothing) for (d, n) in table]
    end
    return map(systems) do system
        system isa AbstractDict && haskey(system, "direction") && haskey(system, "normal") ||
            norma_abort("Crystal plasticity: every slip system requires a `direction` and a `normal`.")
        h = if haskey(system, "initial hardening state")
            checked_parameter(Float64(system["initial hardening state"]), "initial hardening state", "a slip system")
        else
            nothing
        end
        return (SVector{3,Float64}(Float64.(system["direction"])), SVector{3,Float64}(Float64.(system["normal"])), h)
    end
end

# Interaction matrices of a family, from the unit lattice-frame directions
# and normals of its systems; all are invariant under the orientation.
function interaction_matrices(law::HardeningLaw, s::Vector{SVector{3,Float64}}, n::Vector{SVector{3,Float64}})
    num = length(s)
    if law isa LinearMinusRecoveryHardening
        return ones(num, num), zeros(num, num)
    elseif law isa SaturationHardening
        sym_P = [0.5 * (s[i] * n[i]' + n[i] * s[i]') for i in 1:num]
        return [abs(sum(sym_P[i] .* sym_P[j])) for i in 1:num, j in 1:num], zeros(num, num)
    elseif law isa DislocationDensityHardening
        A = [abs(dot(n[i], normalize(cross(n[j], s[j])))) for i in 1:num, j in 1:num]
        return A, sqrt.(max.(1.0 .- A .^ 2, 0.0))
    end
    return zeros(num, num), zeros(num, num)
end

function CrystalPlasticity(params::Parameters)
    CC = voigt_to_tensor(anisotropic_elastic_constants(params))
    ρ = get(params, "density", 0.0)
    scheme = lowercase(string(get(params, "integration scheme", "implicit")))
    scheme in ("implicit", "explicit") ||
        norma_abort("Crystal plasticity: `integration scheme` must be implicit or explicit (got '$scheme').")
    residual = lowercase(string(get(params, "residual", "slip")))
    residual in ("slip", "slip and hardness") ||
        norma_abort("Crystal plasticity: `residual` must be slip or slip and hardness (got '$residual').")
    predictor = lowercase(string(get(params, "slip predictor", "rate")))
    predictor in ("rate", "none") ||
        norma_abort("Crystal plasticity: `slip predictor` must be rate or none (got '$predictor').")
    # The explicit update hardens with the predicted rates, which `none` sets
    # to zero, so the hardening state would never change (as in LCM).
    scheme == "explicit" &&
        predictor == "none" &&
        norma_abort(
            "Crystal plasticity: the explicit scheme hardens with the slip rates of the previous step, " *
            "which `slip predictor: none` sets to zero, so its hardening would never change; " *
            "use `slip predictor: rate` (the default) or the implicit scheme.",
        )
    solver = get(params, "local solver", Dict{String,Any}())
    s = SVector{3,Float64}[]
    n = SVector{3,Float64}[]
    h_initial = Float64[]
    families = SlipFamily[]
    families_params = get(params, "slip families", [])
    families_params isa AbstractVector && all(f -> f isa AbstractDict, families_params) ||
        norma_abort("Crystal plasticity: `slip families` must be a list of families.")
    for family_params in families_params
        flow = create_flow_rule(get(family_params, "flow rule", nothing))
        hardening = create_hardening_law(get(family_params, "hardening law", nothing))
        first_system = length(s) + 1
        family_systems = family_slip_systems(family_params)
        isempty(family_systems) && norma_abort("Crystal plasticity: a slip family has no slip systems.")
        for (direction, normal, h) in family_systems
            norm(direction) > 0.0 && norm(normal) > 0.0 ||
                norma_abort("Crystal plasticity: slip directions and normals must be nonzero.")
            d, m = normalize(direction), normalize(normal)
            abs(dot(d, m)) <= 1.0e-8 || norma_abortf(
                "Crystal plasticity: the slip direction %s is not orthogonal to the normal %s.",
                string(Tuple(direction)),
                string(Tuple(normal)),
            )
            push!(s, d)
            push!(n, m)
            push!(h_initial, h === nothing ? initial_hardening_state(hardening) : h)
        end
        systems = collect(first_system:length(s))
        latent, aux = interaction_matrices(hardening, s[systems], n[systems])
        max_rate = checked_parameter(
            Float64(get(family_params, "maximum slip rate", log(floatmax(Float64)))),
            "maximum slip rate",
            "a slip family",
        )
        push!(families, SlipFamily(flow, hardening, systems, latent, aux, max_rate))
    end
    return CrystalPlasticity(
        CC,
        maximum_wave_modulus(CC),
        ρ,
        s,
        n,
        families,
        h_initial,
        scheme == "explicit",
        residual == "slip and hardness",
        predictor == "rate",
        Float64(get(solver, "absolute tolerance", 1.0e-12)),
        Float64(get(solver, "relative tolerance", 1.0e-14)),
        Int(get(solver, "minimum iterations", 0)),
        Int(get(solver, "maximum iterations", 50)),
    )
end

has_orientation(::CrystalPlasticity) = true

number_slip_systems(material::CrystalPlasticity) = length(material.s)

number_states(material::CrystalPlasticity) = 10 + 4 * number_slip_systems(material)

function initial_state(material::CrystalPlasticity)
    N = number_slip_systems(material)
    return vcat(vec(Matrix{Float64}(I, 3, 3)), zeros(N), zeros(N), material.h_initial, zeros(N), [0.0])
end

function internal_variable_names(material::CrystalPlasticity)
    systems = 1:number_slip_systems(material)
    return vcat(
        ["Fp_$i$j" for j in 1:3 for i in 1:3],
        ["gamma_$α" for α in systems],
        ["gamma_dot_$α" for α in systems],
        ["tau_hard_$α" for α in systems],
        ["tau_$α" for α in systems],
        ["eqps"],
    )
end

get_p_wave_modulus(material::CrystalPlasticity) = material.M

# ---------------------------------------------------------------------------
# Constitutive update
# ---------------------------------------------------------------------------

function slip_rates(material::CrystalPlasticity, τ::AbstractVector, g::AbstractVector)
    γ̇ = Vector{promote_type(eltype(τ), eltype(g))}(undef, length(τ))
    for family in material.families
        family_slip_rates!(γ̇, family, τ, g)
    end
    return γ̇
end

# The rates of one family's systems. A function of its own, so that the flow
# rule, whose type the family carries, is known inside the loop.
function family_slip_rates!(γ̇::AbstractVector, family::SlipFamily, τ::AbstractVector, g::AbstractVector)
    for i in family.systems
        # A NaN slip resistance (or shear) would give a zero rate, since
        # every comparison of the flow rules is false for it.
        isfinite(τ[i]) && isfinite(g[i]) ||
            throw(DomainError(g[i], "Crystal plasticity: a resolved shear or slip resistance is not finite."))
        rate = slip_rate(family.flow, τ[i], g[i])
        abs(rate) > family.max_rate &&
            throw(DomainError(rate, "Crystal plasticity: a slip rate exceeds the maximum of its family."))
        γ̇[i] = rate
    end
    return nothing
end

function harden(
    material::CrystalPlasticity,
    Δt::Float64,
    γ̇::AbstractVector,
    h_n::AbstractVector,
    h_in::AbstractVector;
    explicit::Bool=false,
)
    T = promote_type(eltype(γ̇), eltype(h_in), Float64)
    N = length(h_n)
    h_np1 = Vector{T}(undef, N)
    g = Vector{T}(undef, N)
    for family in material.families
        harden!(h_np1, g, family.hardening, family, Δt, γ̇, h_n, h_in; explicit=explicit)
    end
    return h_np1, g
end

function plastic_deformation(Fp_n::SMatrix{3,3,Float64,9}, Δγ::AbstractVector, P::Vector{SMatrix{3,3,Float64,9}})
    isempty(P) && return Fp_n
    return exp(sum(Δγ[α] * P[α] for α in eachindex(P))) * Fp_n
end

# Elastic response on Fe = F Fp⁻¹: Cauchy stress, second Piola-Kirchhoff
# stress in the intermediate configuration, resolved shear stresses and Fe.
function crystal_stress(
    CC::SArray{Tuple{3,3,3,3},Float64,4,81}, F::SMatrix{3,3}, Fp::SMatrix{3,3}, P::Vector{SMatrix{3,3,Float64,9}}
)
    # Written so that NaN fails each check, as every comparison with it is false.
    all(isfinite, Fp) && maximum(abs, Fp) <= 1.0e100 ||
        throw(DomainError(maximum(abs, Fp), "Crystal plasticity: unbounded Fp."))
    det(Fp) != 0.0 || throw(DomainError(0.0, "Crystal plasticity: singular Fp."))
    Fe = F * inv(Fp)
    Je = det(Fe)
    abs(Je) >= sqrt(floatmin(Float64)) || throw(DomainError(Je, "Crystal plasticity: singular Fe."))
    Ce = Fe' * Fe
    S = double_contraction(CC, 0.5 * (Ce - I))
    σ = Fe * S * Fe' / Je
    mandel = Ce * S
    τ = [sum(Pα .* mandel) for Pα in P]
    return σ, S, τ, Fe
end

# Everything a step needs besides the deformation gradient.
struct CrystalStep
    CC::SArray{Tuple{3,3,3,3},Float64,4,81}  # sample-frame elasticity tensor
    P::Vector{SMatrix{3,3,Float64,9}}         # sample-frame Schmid tensors
    Fp_n::SMatrix{3,3,Float64,9}
    γ_n::Vector{Float64}
    h_n::Vector{Float64}
    Δt::Float64
end

slip_rates_of(step::CrystalStep, γ::AbstractVector) = step.Δt > 0.0 ? (γ .- step.γ_n) ./ step.Δt : zero(γ)

# Residual of the implicit schemes, generic in the unknowns x and in F.
function crystal_residual(material::CrystalPlasticity, step::CrystalStep, x::AbstractVector, F::SMatrix{3,3})
    N = length(step.γ_n)
    γ = x[1:N]
    Δγ = γ .- step.γ_n
    norm(Δγ) > log(floatmax(Float64)) && throw(DomainError(norm(Δγ), "Crystal plasticity: slip increment overflow."))
    Fp = plastic_deformation(step.Fp_n, Δγ, step.P)
    _, _, τ, _ = crystal_stress(step.CC, F, Fp, step.P)
    h_in = material.slip_hardness ? x[(N + 1):(2N)] : step.h_n
    h, g = harden(material, step.Δt, slip_rates_of(step, γ), step.h_n, h_in)
    # With Δt = 0 the slips cannot change, and the slip rates, whose cap
    # could reject a purely elastic state, are not needed.
    r_γ = step.Δt > 0.0 ? Δγ .- step.Δt .* slip_rates(material, τ, g) : Δγ .+ zero.(τ)
    material.slip_hardness || return r_γ
    return vcat(r_γ, x[(N + 1):(2N)] .- h)
end

# Response at the converged unknowns of an implicit scheme: the hardening is
# evaluated once more with the converged slip rates.
function crystal_response(material::CrystalPlasticity, step::CrystalStep, x::AbstractVector, F::SMatrix{3,3})
    N = length(step.γ_n)
    γ = x[1:N]
    γ̇ = slip_rates_of(step, γ)
    h_in = material.slip_hardness ? x[(N + 1):(2N)] : step.h_n
    h, _ = harden(material, step.Δt, γ̇, step.h_n, h_in; explicit=(!material.slip_hardness))
    Fp = plastic_deformation(step.Fp_n, γ .- step.γ_n, step.P)
    σ, S, τ, Fe = crystal_stress(step.CC, F, Fp, step.P)
    return (det(F) * σ * inv(F)', S, Fe, Fp, γ, γ̇, h, τ)
end

# The explicit update: τ from Fp_n, the hardening from the predicted rates.
function crystal_explicit_response(
    material::CrystalPlasticity, step::CrystalStep, γ̇_pred::Vector{Float64}, h_pred::Vector{Float64}, F::SMatrix{3,3}
)
    _, _, τ_trial, _ = crystal_stress(step.CC, F, step.Fp_n, step.P)
    h, g = harden(material, step.Δt, γ̇_pred, step.h_n, h_pred; explicit=true)
    γ = step.Δt > 0.0 ? step.γ_n .+ step.Δt .* slip_rates(material, τ_trial, g) : step.γ_n .+ zero.(τ_trial)
    # The bound of the implicit residual: beyond it the exponential map of Fp
    # overflows.
    norm(γ .- step.γ_n) > log(floatmax(Float64)) &&
        throw(DomainError(norm(γ .- step.γ_n), "Crystal plasticity: slip increment overflow."))
    γ̇ = step.Δt > 0.0 ? (γ .- step.γ_n) ./ step.Δt : γ̇_pred .+ zero(γ)
    Fp = plastic_deformation(step.Fp_n, γ .- step.γ_n, step.P)
    σ, S, τ, Fe = crystal_stress(step.CC, F, Fp, step.P)
    return (det(F) * σ * inv(F)', S, Fe, Fp, γ, γ̇, h, τ)
end

# Norm of the residual for the convergence test: slips are dimensionless,
# and the hardness residual is scaled by the hardness at the start of the
# step. The scaling does not change the Newton iterates.
function residual_norm(material::CrystalPlasticity, step::CrystalStep, r::AbstractVector)
    material.slip_hardness || return norm(r)
    N = length(step.γ_n)
    scale = [abs(h) > 0.0 ? abs(h) : 1.0 for h in step.h_n]
    return sqrt(sum(abs2, r[1:N]) + sum(abs2, r[(N + 1):(2N)] ./ scale))
end

# Newton's method with backtracking on the residual norm.
function solve_crystal_residual(
    material::CrystalPlasticity, step::CrystalStep, x0::Vector{Float64}, F::SMatrix{3,3,Float64,9}
)
    residual = x -> crystal_residual(material, step, x, F)
    x = copy(x0)
    r = residual(x)
    e0 = residual_norm(material, step, r)
    e = e0
    iteration = 0
    while true
        all(isfinite, r) || throw(DomainError(e, "Crystal plasticity: non-finite local residual."))
        converged = e <= material.absolute_tolerance || e <= material.relative_tolerance * e0
        converged && iteration >= material.minimum_iterations && return x
        iteration == material.maximum_iterations &&
            throw(DomainError(e, "Crystal plasticity: the local Newton iteration did not converge."))
        # The Jacobian of a Vector has a chunk size known only at run time, so
        # its type is asserted to keep the iteration type-stable.
        δ = ForwardDiff.jacobian(residual, x)::Matrix{Float64} \ r
        all(isfinite, δ) || throw(DomainError(e, "Crystal plasticity: singular local Jacobian."))
        α = 1.0
        while true
            x_trial = x .- α .* δ
            r_trial = try
                residual(x_trial)
            catch error
                error isa _MATH_ERRORS || rethrow()
                nothing
            end
            if r_trial !== nothing && all(isfinite, r_trial)
                e_trial = residual_norm(material, step, r_trial)
                if e_trial <= (1.0 - 1.0e-4 * α) * e || α < 2.0^-10
                    x, r, e = x_trial, r_trial, e_trial
                    break
                end
            end
            α /= 2
            α < 2.0^-20 && throw(DomainError(e, "Crystal plasticity: the local line search failed."))
        end
        iteration += 1
    end
end

# Predictor: the rates of the previous step, and the hardening they give
# (LCM's rate predictor), or no slip. A function of its own, so that the
# closures of the update capture variables assigned once, which keeps it
# type-stable.
function crystal_prediction(
    material::CrystalPlasticity, Δt::Float64, γ_n::Vector{Float64}, γ̇_n::Vector{Float64}, h_n::Vector{Float64}
)
    N = length(γ_n)
    if Δt > 0.0 && material.rate_predictor && N > 0
        h_pred, _ = harden(material, Δt, γ̇_n, h_n, h_n)
        return γ̇_n, γ_n .+ Δt .* γ̇_n, h_pred
    end
    return zeros(N), copy(γ_n), copy(h_n)
end

function plastic_velocity_gradient(Δt::Float64, γ̇::Vector{Float64}, P::Vector{SMatrix{3,3,Float64,9}})
    Δt > 0.0 && !isempty(P) || return zero(SMatrix{3,3,Float64,9})
    return sum(γ̇[α] * P[α] for α in eachindex(P))
end

function constitutive(
    material::CrystalPlasticity,
    F::SMatrix{3,3,Float64,9},
    Q::SMatrix{3,3,Float64,9},
    state_old::Vector{Float64},
    Δt::Float64;
    need_tangent::Bool=true,
)
    N = number_slip_systems(material)
    Fp_n = SMatrix{3,3,Float64,9}(state_old[1:9])
    γ_n = state_old[10:(9 + N)]
    γ̇_n = state_old[(10 + N):(9 + 2N)]
    h_n = state_old[(10 + 2N):(9 + 3N)]
    eqps_n = state_old[10 + 4N]
    P = [SMatrix{3,3,Float64,9}((Q * material.s[α]) * (Q * material.n[α])') for α in 1:N]
    step = CrystalStep(rotate_tensor4(Q, material.CC), P, Fp_n, γ_n, h_n, Δt)
    γ̇_pred, γ_pred, h_pred = crystal_prediction(material, Δt, γ_n, γ̇_n, h_n)
    if material.explicit || N == 0
        response = f -> crystal_explicit_response(material, step, γ̇_pred, h_pred, SMatrix{3,3,eltype(f),9}(f))
        P1, S, Fe, Fp, γ, γ̇, h, τ = response(F)
        AA = if need_tangent
            SArray{Tuple{3,3,3,3},Float64,4,81}(ForwardDiff.jacobian(f -> vec(response(f)[1]), SVector{9}(F)))
        else
            ZERO_TANGENT
        end
    else
        x0 = material.slip_hardness ? vcat(γ_pred, h_pred) : γ_pred
        x = solve_crystal_residual(material, step, x0, F)
        P1, S, Fe, Fp, γ, γ̇, h, τ = crystal_response(material, step, x, F)
        if need_tangent
            # Implicit function theorem: dx/dF = -(∂r/∂x)⁻¹ ∂r/∂F, then the
            # response is differentiated along the linearized x(F).
            f0 = SVector{9,Float64}(F)
            J_x = ForwardDiff.jacobian(y -> crystal_residual(material, step, y, F), x)::Matrix{Float64}
            J_F = ForwardDiff.jacobian(
                f -> crystal_residual(material, step, x, SMatrix{3,3,eltype(f),9}(f)), f0
            )::Matrix{Float64}
            G = -(J_x \ J_F)
            linearized = f -> vec(crystal_response(material, step, x .+ G * (f .- f0), SMatrix{3,3,eltype(f),9}(f))[1])
            AA = SArray{Tuple{3,3,3,3},Float64,4,81}(ForwardDiff.jacobian(linearized, f0))
        else
            AA = ZERO_TANGENT
        end
    end
    Lp = plastic_velocity_gradient(Δt, γ̇, P)
    Dp = 0.5 * (Lp + Lp')
    eqps = eqps_n + Δt * sqrt(2.0 / 3.0 * sum(Dp .* Dp))
    Ee = 0.5 * (Fe' * Fe - I)
    W = 0.5 * sum(Ee .* S)
    state_new = Vector{Float64}(vcat(vec(Fp), γ, γ̇, h, τ, [eqps]))
    return W, P1, AA, state_new
end
