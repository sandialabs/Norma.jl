# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Three-field Hu-Washizu element with a projected volumetric strain.
#
# The fields are the motion, an independent volumetric logarithmic strain θ̄,
# and a pressure p̄.  Both auxiliary fields live in the same element-local
# discontinuous space, which is what lets them be eliminated element by element
# rather than globally.
#
# With θ̄ and p̄ in one space Q_e, stationarity gives p̄ = κθ̄ exactly (not merely
# after a projection) and θ̄ = Π_e θ(u), the L² projection of θ = log J onto
# Q_e.  The constraint term of the functional then vanishes identically:
#
#     ∫ p̄ (θ - θ̄) = κ ∫ θ̄ (θ - θ̄) = 0
#
# by orthogonality, since κθ̄ ∈ Q_e and (θ - θ̄) ⊥ Q_e.  So the three-field
# functional collapses to a displacement functional carrying a projected
# volumetric strain,
#
#     Π[u] = ∫ W_dev(dev ε(u)) + ∫ (κ/2) (Π_e log J)² ,
#
# whose stationarity conditions are exactly the three-field equations.  This is
# the mean-dilatation mechanism of Simo, Taylor and Pister (1985); what is
# open, and what this element exists to explore, is the choice of Q_e.  A
# constant Q_e is the classical element; a linear one constrains all four
# volumetric directions a quadratic tetrahedron can produce.
#
# Residual and tangent are taken from the element energy by automatic
# differentiation rather than derived by hand.  The projection couples the
# quadrature points of an element, so the tangent is not a sum of per-point
# contributions and cannot be assembled the way the displacement element's is.
# Differentiating the energy is exact, and it keeps this file honest: there is
# one definition of the formulation, not a residual and a tangent that can
# drift apart.

using ForwardDiff
using LinearAlgebra
using StaticArrays

# Element-local pressure space, evaluated at a parametric point.  `order` 0 is
# the classical element-constant dilatation; order 1 is linear, which on a
# quadratic tetrahedron matches the dimension of the volumetric strain the
# displacement space can actually produce.
@inline pressure_shape(::Val{0}, ξ::AbstractVector) = SVector{1,Float64}(1.0)
@inline pressure_shape(::Val{1}, ξ::AbstractVector) = SVector{4,Float64}(1.0, ξ[1], ξ[2], ξ[3])

# Integer entry point.  Its return type is a Union of the two SVector sizes, so
# anything that calls it in a loop must go through `Val` instead -- see
# `projected_theta`.
@inline function pressure_shape(order::Int, ξ::AbstractVector)
    if order == 0
        return pressure_shape(Val(0), ξ)
    elseif order == 1
        return pressure_shape(Val(1), ξ)
    end
    return error("three-field pressure order must be 0 or 1, got $order")
end

@inline pressure_dim(order::Int) = order == 0 ? 1 : 4

"""
    check_three_field_quadrature(order, num_points)

Refuse a quadrature rule with no more points than the pressure space has
dimensions.

With `num_points == pressure_dim(order)` the L² projection is not a projection
at all: it interpolates θ exactly at the quadrature points, so θ̄ == θ wherever
the energy is sampled and the element is silently identical to a displacement
element.  Measured on one TETRA10 with a non-affine field: the four-point rule
gives a relative energy difference of exactly 0.0 at order 1, against 0.372 at
order 0.  That is a formulation that appears to run and does nothing, which is
worse than one that fails.
"""
function check_three_field_quadrature(order::Int, num_points::Int)
    m = pressure_dim(order)
    num_points > m || error(
        "three-field pressure order $order needs a quadrature rule with more " *
        "than $m points to be a projection rather than interpolation; got " *
        "$num_points. The element would be identical to a displacement " *
        "element. Use a richer rule (TETRA10 has 5- and 14-point rules).")
    return nothing
end

"""
    projected_theta(order, ξ, θs, dvols)

`L²` projection of the pointwise volumetric strain θ = log J onto the
element-local pressure space, evaluated back at each quadrature point.

Returns the projected values in the same order as the inputs.  The projection
matrix is the element pressure mass matrix, which for order 0 is a scalar and
for order 1 a 4×4 solve -- small enough that forming it per element per
evaluation costs nothing measurable next to the constitutive work.
"""
function projected_theta(order::Int, ξ, θs::AbstractVector{T}, dvols) where {T}
    # Branch on the order ONCE, here, and hand a Val to the kernel.  Calling
    # pressure_shape(::Int, ...) inside the loop would infer φ as a Union of
    # SVector{1} and SVector{4} and send every use of it through a dynamic
    # dispatch.
    if order == 0
        return _projected_theta(Val(0), ξ, θs, dvols)
    elseif order == 1
        return _projected_theta(Val(1), ξ, θs, dvols)
    end
    return error("three-field pressure order must be 0 or 1, got $order")
end

function _projected_theta(v::Val{K}, ξ, θs::AbstractVector{T}, dvols) where {K,T}
    m = pressure_dim(K)
    M = zeros(T, m, m)
    b = zeros(T, m)
    for q in eachindex(θs)
        φ = pressure_shape(v, view(ξ, :, q))
        w = dvols[q]
        for i in 1:m
            b[i] += w * φ[i] * θs[q]
            for j in 1:m
                M[i, j] += w * φ[i] * φ[j]
            end
        end
    end
    c = M \ b
    return [dot(pressure_shape(v, view(ξ, :, q)), c) for q in eachindex(θs)]
end

"""
    three_field_energy(material, order, Np, dN, ip_weights, X, u_flat, states)

Strain energy of one element under the three-field formulation, as a function
of the element displacement degrees of freedom.

The deviatoric response is evaluated pointwise, at the quadrature point where
the strain is computed, so the constitutive law never sees a quantity from
another discretization; only the volumetric strain is projected.  Written
generically in the element type of `u_flat` so `ForwardDiff` can carry dual
numbers through it.
"""
function three_field_energy(material, order::Int, dN, ip_weights, ξ, X,
                            u_flat::AbstractVector{T}) where {T}
    num_nodes = size(X, 2)
    num_points = length(ip_weights)
    u = reshape(u_flat, 3, num_nodes)
    x = X + u

    Fs = Vector{SMatrix{3,3,T,9}}(undef, num_points)
    θs = Vector{T}(undef, num_points)
    dvols = Vector{Float64}(undef, num_points)
    for q in 1:num_points
        dNdξ = dN[:, :, q]
        dXdξ = SMatrix{3,3,Float64,9}(dNdξ * X')
        dNdX = dXdξ \ dNdξ
        F = SMatrix{3,3,T,9}(x * dNdX')
        Fs[q] = F
        θs[q] = log(det(F))
        dvols[q] = det(dXdξ) * ip_weights[q]
    end

    θbars = projected_theta(order, ξ, θs, dvols)

    # W_dev at the pointwise strain, plus the volumetric energy at the projected
    # strain.  The subtraction is exact rather than approximate BECAUSE the
    # energy splits additively in the logarithmic measure:
    #
    #     W(F) = (κ/2)(tr E)² + μ‖dev E‖²,   E = ½ log C,   tr E = log J = θ,
    #
    # so W(F) - (κ/2)θ² is exactly the deviatoric part, with nothing left over.
    # That identity is why the formulation is built on the Hencky measure and
    # is checked by `three_field_supported`; a material without the additive
    # split would leave a volumetric remainder in the "deviatoric" term and the
    # element would be silently wrong.
    κ = material.κ
    energy = zero(T)
    for q in 1:num_points
        W_total = strain_energy(material, Fs[q])
        W_vol_pointwise = 0.5 * κ * θs[q]^2
        W_vol_projected = 0.5 * κ * θbars[q]^2
        energy += (W_total - W_vol_pointwise + W_vol_projected) * dvols[q]
    end
    return energy
end

# Materials whose stored energy splits additively into volumetric and
# deviatoric parts in the logarithmic strain measure, which is what makes the
# subtraction in `three_field_energy` exact.  Anything else must be refused
# rather than approximated: the error would be a quiet inaccuracy in the
# deviatoric response, not a visible failure.
three_field_supported(::Hencky) = true
three_field_supported(::Any) = false

"""
    three_field_element(material, order, dN, ip_weights, ξ, X, u_flat; need_tangent)

Energy, internal force and (optionally) stiffness of one element, all from the
same energy function by automatic differentiation.
"""
function three_field_element(material, order::Int, dN, ip_weights, ξ, X,
                             u_flat::Vector{Float64}; need_tangent::Bool=true)
    check_three_field_quadrature(order, length(ip_weights))
    three_field_supported(material) || error(
        "three-field formulation requires a material whose stored energy splits " *
        "additively in the logarithmic strain measure; got $(typeof(material)). " *
        "Hencky is supported.")
    f = uu -> three_field_energy(material, order, dN, ip_weights, ξ, X, uu)
    energy = f(u_flat)
    force = ForwardDiff.gradient(f, u_flat)
    stiff = need_tangent ? ForwardDiff.hessian(f, u_flat) : zeros(0, 0)
    return energy, force, stiff
end
