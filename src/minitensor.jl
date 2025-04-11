# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using LinearAlgebra
using LinearAlgebra: norm, tr, dot, normalize
using StaticArrays

#
# Lie groups and Lie algebra utilities, mostly algebra of rotations
# SO(3) and so(3).
#

"""
    skew(A::AbstractMatrix)

Returns the skew-symmetric part of a square matrix: (A - A') / 2
"""
function skew(A::AbstractMatrix{T}) where {T<:Number}
    @assert size(A, 1) == size(A, 2) "Matrix must be square"
    return 0.5 * (A - A')
end

"""
    symm(A::AbstractMatrix)

Returns the symmetric part of a square matrix: (A + A') / 2
"""
function symm(A::AbstractMatrix{T}) where {T<:Number}
    @assert size(A, 1) == size(A, 2) "Matrix must be square"
    return 0.5 * (A + A')
end

"""
    hat(w::SVector{3, Float64}) -> SMatrix{3,3,Float64}

Convert a 3-vector into a 3x3 skew-symmetric matrix.
"""
function hat(w::SVector{3,Float64})::SMatrix{3,3,Float64}
    return @SMatrix [
        0.0 -w[3] w[2]
        w[3] 0.0 -w[1]
        -w[2] w[1] 0.0
    ]
end

"""
    vee(W::SMatrix{3,3,Float64}) -> SVector{3, Float64}

Convert a 3x3 skew-symmetric matrix into a 3-vector.
"""
function vee(W::SMatrix{3,3,Float64})::SVector{3,Float64}
    return @SVector [W[3, 2], W[1, 3], W[2, 1]]
end

"""
    bch(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where T <: Number -> Matrix{T}

Computes the Baker–Campbell–Hausdorff (BCH) series for two square matrices `x` and `y` up to 8th-order terms.

The BCH formula gives an expression for the matrix `z` such that:

    exp(z) = exp(x) * exp(y)

This expression allows the composition of Lie group elements to be approximated using only Lie algebraic operations.
It is particularly relevant for matrix Lie groups where `x` and `y` are elements of the Lie algebra
(e.g., skew-symmetric matrices for SO(3)).

# Arguments
- `x::AbstractMatrix{T}`: A square matrix (element of a Lie algebra).
- `y::AbstractMatrix{T}`: Another square matrix of the same size and type as `x`.

# Returns
- `z::Matrix{T}`: The BCH approximation of `log(exp(x) * exp(y))` up to 8th order.

# Notes
- The matrices `x` and `y` must be square and of the same size.
- Coefficients are based on Goldberg’s algorithm, as implemented in Mathematica following:
  Weyrauch & Scholz, *Computing the BCH series and the Zassenhaus product*, CPC, 2009.
- This implementation is symbolic and is not optimized for runtime performance.

# References
- Weyrauch & Scholz, *Computer Physics Communications*, 180(9), 1558–1565, 2009.
- https://en.wikipedia.org/wiki/Baker–Campbell–Hausdorff_formula

# Example
```julia
x = randn(3, 3)
x = 0.5 * (x - x')  # make x skew-symmetric

y = randn(3, 3)
y = 0.5 * (y - y')

z = bch(x, y)
```
"""
function bch(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    @assert size(x) == size(y) "Input matrices must be the same size"
    @assert size(x, 1) == size(x, 2) "Matrices must be square"
    z1 = bch_term1(x, y)
    z2 = bch_term2(x, y)
    z3 = bch_term3(x, y)
    z4 = bch_term4(x, y)
    z5 = bch_term5(x, y)
    z6 = bch_term6(x, y)
    z7 = bch_term7(x, y)
    z8 = bch_term8(x, y)
    z = z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8
    return z
end

function bch_term1(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return x + y
end

function bch_term2(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return 0.5 * (x * y - y * x)
end

function bch_term3(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return x * x * y / 12 - x * y * x / 6 + x * y * y / 12 + y * x * x / 12 - y * x * y / 6 + y * y * x / 12
end

function bch_term4(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return x * x * y * y / 24 - x * y * x * y / 12 + y * x * y * x / 12 - y * y * x * x / 24
end

function bch_term5(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return -x * x * x * x * y / 720 + x * x * x * y * x / 180 + x * x * x * y * y / 180 - x * x * y * x * x / 120 -
           x * x * y * x * y / 120 - x * x * y * y * x / 120 +
           x * x * y * y * y / 180 +
           x * y * x * x * x / 180 - x * y * x * x * y / 120 + x * y * x * y * x / 30 - x * y * x * y * y / 120 -
           x * y * y * x * x / 120 - x * y * y * x * y / 120 + x * y * y * y * x / 180 - x * y * y * y * y / 720 -
           y * x * x * x * x / 720 + y * x * x * x * y / 180 - y * x * x * y * x / 120 - y * x * x * y * y / 120 -
           y * x * y * x * x / 120 + y * x * y * x * y / 30 - y * x * y * y * x / 120 +
           y * x * y * y * y / 180 +
           y * y * x * x * x / 180 - y * y * x * x * y / 120 - y * y * x * y * x / 120 - y * y * x * y * y / 120 +
           y * y * y * x * x / 180 +
           y * y * y * x * y / 180 - y * y * y * y * x / 720
end

function bch_term6(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return -x * x * x * x * y * y / 1440 + x * x * x * y * x * y / 360 + x * x * x * y * y * y / 360 -
           x * x * y * x * x * y / 240 - x * x * y * x * y * y / 240 - x * x * y * y * x * y / 240 -
           x * x * y * y * y * y / 1440 + x * y * x * x * x * y / 360 - x * y * x * x * y * y / 240 +
           x * y * x * y * x * y / 60 +
           x * y * x * y * y * y / 360 - x * y * y * x * x * y / 240 - x * y * y * x * y * y / 240 +
           x * y * y * y * x * y / 360 - y * x * x * x * y * x / 360 +
           y * x * x * y * x * x / 240 +
           y * x * x * y * y * x / 240 - y * x * y * x * x * x / 360 - y * x * y * x * y * x / 60 +
           y * x * y * y * x * x / 240 - y * x * y * y * y * x / 360 +
           y * y * x * x * x * x / 1440 +
           y * y * x * x * y * x / 240 +
           y * y * x * y * x * x / 240 +
           y * y * x * y * y * x / 240 - y * y * y * x * x * x / 360 - y * y * y * x * y * x / 360 +
           y * y * y * y * x * x / 1440
end

function bch_term7(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return x * x * x * x * x * x * y / 30240 - x * x * x * x * x * y * x / 5040 - x * x * x * x * x * y * y / 5040 +
           x * x * x * x * y * x * x / 2016 +
           x * x * x * x * y * x * y / 2016 +
           x * x * x * x * y * y * x / 2016 +
           x * x * x * x * y * y * y / 3780 - x * x * x * y * x * x * x / 1512 - x * x * x * y * x * x * y / 5040 -
           x * x * x * y * x * y * x / 630 - x * x * x * y * x * y * y / 5040 - x * x * x * y * y * x * x / 5040 -
           x * x * x * y * y * x * y / 5040 - x * x * x * y * y * y * x / 1512 +
           x * x * x * y * y * y * y / 3780 +
           x * x * y * x * x * x * x / 2016 - x * x * y * x * x * x * y / 5040 + x * x * y * x * x * y * x / 840 -
           x * x * y * x * x * y * y / 1120 +
           x * x * y * x * y * x * x / 840 +
           x * x * y * x * y * x * y / 840 +
           x * x * y * x * y * y * x / 840 - x * x * y * x * y * y * y / 5040 - x * x * y * y * x * x * x / 5040 -
           x * x * y * y * x * x * y / 1120 + x * x * y * y * x * y * x / 840 - x * x * y * y * x * y * y / 1120 -
           x * x * y * y * y * x * x / 5040 - x * x * y * y * y * x * y / 5040 + x * x * y * y * y * y * x / 2016 -
           x * x * y * y * y * y * y / 5040 - x * y * x * x * x * x * x / 5040 + x * y * x * x * x * x * y / 2016 -
           x * y * x * x * x * y * x / 630 - x * y * x * x * x * y * y / 5040 +
           x * y * x * x * y * x * x / 840 +
           x * y * x * x * y * x * y / 840 +
           x * y * x * x * y * y * x / 840 - x * y * x * x * y * y * y / 5040 - x * y * x * y * x * x * x / 630 +
           x * y * x * y * x * x * y / 840 - x * y * x * y * x * y * x / 140 +
           x * y * x * y * x * y * y / 840 +
           x * y * x * y * y * x * x / 840 +
           x * y * x * y * y * x * y / 840 - x * y * x * y * y * y * x / 630 +
           x * y * x * y * y * y * y / 2016 +
           x * y * y * x * x * x * x / 2016 - x * y * y * x * x * x * y / 5040 + x * y * y * x * x * y * x / 840 -
           x * y * y * x * x * y * y / 1120 +
           x * y * y * x * y * x * x / 840 +
           x * y * y * x * y * x * y / 840 +
           x * y * y * x * y * y * x / 840 - x * y * y * x * y * y * y / 5040 - x * y * y * y * x * x * x / 1512 -
           x * y * y * y * x * x * y / 5040 - x * y * y * y * x * y * x / 630 - x * y * y * y * x * y * y / 5040 +
           x * y * y * y * y * x * x / 2016 +
           x * y * y * y * y * x * y / 2016 - x * y * y * y * y * y * x / 5040 +
           x * y * y * y * y * y * y / 30240 +
           y * x * x * x * x * x * x / 30240 - y * x * x * x * x * x * y / 5040 +
           y * x * x * x * x * y * x / 2016 +
           y * x * x * x * x * y * y / 2016 - y * x * x * x * y * x * x / 5040 - y * x * x * x * y * x * y / 630 -
           y * x * x * x * y * y * x / 5040 - y * x * x * x * y * y * y / 1512 - y * x * x * y * x * x * x / 5040 +
           y * x * x * y * x * x * y / 840 +
           y * x * x * y * x * y * x / 840 +
           y * x * x * y * x * y * y / 840 - y * x * x * y * y * x * x / 1120 + y * x * x * y * y * x * y / 840 -
           y * x * x * y * y * y * x / 5040 +
           y * x * x * y * y * y * y / 2016 +
           y * x * y * x * x * x * x / 2016 - y * x * y * x * x * x * y / 630 +
           y * x * y * x * x * y * x / 840 +
           y * x * y * x * x * y * y / 840 +
           y * x * y * x * y * x * x / 840 - y * x * y * x * y * x * y / 140 + y * x * y * x * y * y * x / 840 -
           y * x * y * x * y * y * y / 630 - y * x * y * y * x * x * x / 5040 +
           y * x * y * y * x * x * y / 840 +
           y * x * y * y * x * y * x / 840 +
           y * x * y * y * x * y * y / 840 - y * x * y * y * y * x * x / 5040 - y * x * y * y * y * x * y / 630 +
           y * x * y * y * y * y * x / 2016 - y * x * y * y * y * y * y / 5040 - y * y * x * x * x * x * x / 5040 +
           y * y * x * x * x * x * y / 2016 - y * y * x * x * x * y * x / 5040 - y * y * x * x * x * y * y / 5040 -
           y * y * x * x * y * x * x / 1120 + y * y * x * x * y * x * y / 840 - y * y * x * x * y * y * x / 1120 -
           y * y * x * x * y * y * y / 5040 - y * y * x * y * x * x * x / 5040 +
           y * y * x * y * x * x * y / 840 +
           y * y * x * y * x * y * x / 840 +
           y * y * x * y * x * y * y / 840 - y * y * x * y * y * x * x / 1120 + y * y * x * y * y * x * y / 840 -
           y * y * x * y * y * y * x / 5040 +
           y * y * x * y * y * y * y / 2016 +
           y * y * y * x * x * x * x / 3780 - y * y * y * x * x * x * y / 1512 - y * y * y * x * x * y * x / 5040 -
           y * y * y * x * x * y * y / 5040 - y * y * y * x * y * x * x / 5040 - y * y * y * x * y * x * y / 630 -
           y * y * y * x * y * y * x / 5040 - y * y * y * x * y * y * y / 1512 +
           y * y * y * y * x * x * x / 3780 +
           y * y * y * y * x * x * y / 2016 +
           y * y * y * y * x * y * x / 2016 +
           y * y * y * y * x * y * y / 2016 - y * y * y * y * y * x * x / 5040 - y * y * y * y * y * x * y / 5040 +
           y * y * y * y * y * y * x / 30240
end

function bch_term8(x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T<:Number}
    return x * x * x * x * x * x * y * y / 60480 - x * x * x * x * x * y * x * y / 10080 -
           x * x * x * x * x * y * y * y / 10080 +
           x * x * x * x * y * x * x * y / 4032 +
           x * x * x * x * y * x * y * y / 4032 +
           x * x * x * x * y * y * x * y / 4032 +
           23 * x * x * x * x * y * y * y * y / 120960 - x * x * x * y * x * x * x * y / 3024 -
           x * x * x * y * x * x * y * y / 10080 - x * x * x * y * x * y * x * y / 1260 -
           x * x * x * y * x * y * y * y / 3024 - x * x * x * y * y * x * x * y / 10080 -
           x * x * x * y * y * x * y * y / 10080 - x * x * x * y * y * y * x * y / 3024 -
           x * x * x * y * y * y * y * y / 10080 + x * x * y * x * x * x * x * y / 4032 -
           x * x * y * x * x * x * y * y / 10080 + x * x * y * x * x * y * x * y / 1680 -
           x * x * y * x * x * y * y * y / 10080 +
           x * x * y * x * y * x * x * y / 1680 +
           x * x * y * x * y * x * y * y / 1680 +
           x * x * y * x * y * y * x * y / 1680 +
           x * x * y * x * y * y * y * y / 4032 - x * x * y * y * x * x * x * y / 10080 -
           x * x * y * y * x * x * y * y / 2240 + x * x * y * y * x * y * x * y / 1680 -
           x * x * y * y * x * y * y * y / 10080 - x * x * y * y * y * x * x * y / 10080 -
           x * x * y * y * y * x * y * y / 10080 +
           x * x * y * y * y * y * x * y / 4032 +
           x * x * y * y * y * y * y * y / 60480 - x * y * x * x * x * x * x * y / 10080 +
           x * y * x * x * x * x * y * y / 4032 - x * y * x * x * x * y * x * y / 1260 -
           x * y * x * x * x * y * y * y / 3024 +
           x * y * x * x * y * x * x * y / 1680 +
           x * y * x * x * y * x * y * y / 1680 +
           x * y * x * x * y * y * x * y / 1680 +
           x * y * x * x * y * y * y * y / 4032 - x * y * x * y * x * x * x * y / 1260 +
           x * y * x * y * x * x * y * y / 1680 - x * y * x * y * x * y * x * y / 280 -
           x * y * x * y * x * y * y * y / 1260 +
           x * y * x * y * y * x * x * y / 1680 +
           x * y * x * y * y * x * y * y / 1680 - x * y * x * y * y * y * x * y / 1260 -
           x * y * x * y * y * y * y * y / 10080 + x * y * y * x * x * x * x * y / 4032 -
           x * y * y * x * x * x * y * y / 10080 + x * y * y * x * x * y * x * y / 1680 -
           x * y * y * x * x * y * y * y / 10080 +
           x * y * y * x * y * x * x * y / 1680 +
           x * y * y * x * y * x * y * y / 1680 +
           x * y * y * x * y * y * x * y / 1680 +
           x * y * y * x * y * y * y * y / 4032 - x * y * y * y * x * x * x * y / 3024 -
           x * y * y * y * x * x * y * y / 10080 - x * y * y * y * x * y * x * y / 1260 -
           x * y * y * y * x * y * y * y / 3024 +
           x * y * y * y * y * x * x * y / 4032 +
           x * y * y * y * y * x * y * y / 4032 - x * y * y * y * y * y * x * y / 10080 +
           y * x * x * x * x * x * y * x / 10080 - y * x * x * x * x * y * x * x / 4032 -
           y * x * x * x * x * y * y * x / 4032 +
           y * x * x * x * y * x * x * x / 3024 +
           y * x * x * x * y * x * y * x / 1260 +
           y * x * x * x * y * y * x * x / 10080 +
           y * x * x * x * y * y * y * x / 3024 - y * x * x * y * x * x * x * x / 4032 -
           y * x * x * y * x * x * y * x / 1680 - y * x * x * y * x * y * x * x / 1680 -
           y * x * x * y * x * y * y * x / 1680 + y * x * x * y * y * x * x * x / 10080 -
           y * x * x * y * y * x * y * x / 1680 + y * x * x * y * y * y * x * x / 10080 -
           y * x * x * y * y * y * y * x / 4032 +
           y * x * y * x * x * x * x * x / 10080 +
           y * x * y * x * x * x * y * x / 1260 - y * x * y * x * x * y * x * x / 1680 -
           y * x * y * x * x * y * y * x / 1680 +
           y * x * y * x * y * x * x * x / 1260 +
           y * x * y * x * y * x * y * x / 280 - y * x * y * x * y * y * x * x / 1680 +
           y * x * y * x * y * y * y * x / 1260 - y * x * y * y * x * x * x * x / 4032 -
           y * x * y * y * x * x * y * x / 1680 - y * x * y * y * x * y * x * x / 1680 -
           y * x * y * y * x * y * y * x / 1680 +
           y * x * y * y * y * x * x * x / 3024 +
           y * x * y * y * y * x * y * x / 1260 - y * x * y * y * y * y * x * x / 4032 +
           y * x * y * y * y * y * y * x / 10080 - y * y * x * x * x * x * x * x / 60480 -
           y * y * x * x * x * x * y * x / 4032 +
           y * y * x * x * x * y * x * x / 10080 +
           y * y * x * x * x * y * y * x / 10080 +
           y * y * x * x * y * x * x * x / 10080 - y * y * x * x * y * x * y * x / 1680 +
           y * y * x * x * y * y * x * x / 2240 +
           y * y * x * x * y * y * y * x / 10080 - y * y * x * y * x * x * x * x / 4032 -
           y * y * x * y * x * x * y * x / 1680 - y * y * x * y * x * y * x * x / 1680 -
           y * y * x * y * x * y * y * x / 1680 + y * y * x * y * y * x * x * x / 10080 -
           y * y * x * y * y * x * y * x / 1680 + y * y * x * y * y * y * x * x / 10080 -
           y * y * x * y * y * y * y * x / 4032 +
           y * y * y * x * x * x * x * x / 10080 +
           y * y * y * x * x * x * y * x / 3024 +
           y * y * y * x * x * y * x * x / 10080 +
           y * y * y * x * x * y * y * x / 10080 +
           y * y * y * x * y * x * x * x / 3024 +
           y * y * y * x * y * x * y * x / 1260 +
           y * y * y * x * y * y * x * x / 10080 +
           y * y * y * x * y * y * y * x / 3024 - 23 * y * y * y * y * x * x * x * x / 120960 -
           y * y * y * y * x * x * y * x / 4032 - y * y * y * y * x * y * x * x / 4032 -
           y * y * y * y * x * y * y * x / 4032 +
           y * y * y * y * y * x * x * x / 10080 +
           y * y * y * y * y * x * y * x / 10080 - y * y * y * y * y * y * x * x / 60480
end

#
# Most of the functions that follow were adopted and adapted from the
# ONDAP FEM code.  See: Object-oriented finite-element dynamic
# simulation of geometrically nonlinear space structures, Victor
# Balopoulos, Ph.D. dissertation, Cornell University, 1997
#

"""
    q_of_rt(R::SMatrix{3,3,Float64}) -> SVector{4,Float64}

Compute the quaternion corresponding to a rotation matrix `R`, using Spurrier's singularity-free algorithm.

This method avoids singularities by selecting the largest of the trace and diagonal elements of the matrix,
and computing the quaternion components accordingly. This ensures a stable evaluation of the quaternion,
which is otherwise a quadratic function of the rotation matrix.
"""
function q_of_rt(R::SMatrix{3,3,Float64})::SVector{4,Float64}
    trR = LinearAlgebra.tr(R)
    maxm = trR
    maxi = 4
    for i in 1:3
        if R[i, i] > maxm
            maxm = R[i, i]
            maxi = i
        end
    end
    if maxi == 4
        root = sqrt(maxm + 1.0)
        factor = 0.5 / root
        return @SVector [
            0.5 * root, factor * (R[3, 2] - R[2, 3]), factor * (R[1, 3] - R[3, 1]), factor * (R[2, 1] - R[1, 2])
        ]
    elseif maxi == 3
        root = sqrt(2.0 * maxm + 1.0 - trR)
        factor = 0.5 / root
        return @SVector [
            factor * (R[2, 1] - R[1, 2]), factor * (R[1, 3] + R[3, 1]), factor * (R[2, 3] + R[3, 2]), 0.5 * root
        ]
    elseif maxi == 2
        root = sqrt(2.0 * maxm + 1.0 - trR)
        factor = 0.5 / root
        return @SVector [
            factor * (R[1, 3] - R[3, 1]), factor * (R[1, 2] + R[2, 1]), 0.5 * root, factor * (R[2, 3] + R[3, 2])
        ]
    else
        root = sqrt(2.0 * maxm + 1.0 - trR)
        factor = 0.5 / root
        return @SVector [
            factor * (R[3, 2] - R[2, 3]), 0.5 * root, factor * (R[1, 2] + R[2, 1]), factor * (R[1, 3] + R[3, 1])
        ]
    end
end

"""
    rv_of_q(q::SVector{4,Float64}) -> SVector{3,Float64}

Return the principal rotation vector associated with the quaternion `q`.

Maps a quaternion `qq = (qs, qv)` to its corresponding principal rotation pseudo-vector `aa` with |aa| ≤ π.
If qs ≈ 1 (i.e., small rotations), uses an alternate evaluation to avoid numerical instability.
"""
function rv_of_q(qq::SVector{4,Float64})::SVector{3,Float64}
    q = qq[1] >= 0.0 ? qq : -qq
    qs = q[1]
    qv = SVector{3,Float64}(q[2:4])
    qvnorm = LinearAlgebra.norm(qv)
    aanorm = 2.0 * (qvnorm < sqrt(0.5) ? asin(qvnorm) : acos(qs))
    coef = qvnorm < sqrt(eps()) ? 2.0 : aanorm / qvnorm
    return coef * qv
end

"""
    psi(x::Float64) -> Float64

Computes sin(x)/x with asymptotic expansions near zero to avoid division by small values.
This function is frequently encountered in rotational algebra when evaluating exp/log maps of skew-symmetric matrices.
"""
function psi(x::Float64)::Float64
    y = abs(x)
    e2 = sqrt(eps())
    e4 = sqrt(e2)
    if y > e4
        return sin(y) / y
    elseif y > e2
        return 1.0 - y * y / 6.0
    else
        return 1.0
    end
end

"""
    q_of_rv(aa::SVector{3,Float64}) -> SVector{4,Float64}

Convert a rotation vector `aa` into a quaternion representation.

Implements the formula:
    qv = sin(|aa| / 2) * aa / |aa|
    qs = cos(|aa| / 2)
with asymptotic handling for small angles using `psi`.
"""
function q_of_rv(aa::SVector{3,Float64})::SVector{4,Float64}
    halfnorm = 0.5 * LinearAlgebra.norm(aa)
    temp = 0.5 * psi(halfnorm)
    return @SVector [cos(halfnorm), temp * aa[1], temp * aa[2], temp * aa[3]]
end

"""
    rt_of_q(q::SVector{4,Float64}) -> SMatrix{3,3,Float64}

Return the rotation matrix associated with the quaternion `q`.

Based on the identity:
    R = 2 * qv⊗qv + 2 * qs * hat(qv) + (2 * qs^2 - 1) * I
where `⊗` denotes the outer product.
"""
function rt_of_q(qq::SVector{4,Float64})::SMatrix{3,3,Float64}
    qs = qq[1]
    qv = SVector{3,Float64}(qq[2:4])
    I = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    return 2.0 * qv * qv' + 2.0 * qs * hat(qv) + (2.0 * qs * qs - 1.0) * I
end

"""
    rv_of_rt(R::SMatrix{3,3,Float64}) -> SVector{3,Float64}

Return the rotation vector corresponding to the rotation matrix `R`.
"""
function rv_of_rt(R::SMatrix{3,3,Float64})::SVector{3,Float64}
    return rv_of_q(q_of_rt(R))
end

"""
    rt_of_rv(w::SVector{3,Float64}) -> SMatrix{3,3,Float64}

Return the rotation matrix associated with a rotation vector `w`.
"""
function rt_of_rv(w::SVector{3,Float64})::SMatrix{3,3,Float64}
    return rt_of_q(q_of_rv(w))
end

"""
    rv_continue(old::SVector{3,Float64}, prev::SVector{3,Float64}) -> SVector{3,Float64}

Returns a rotation pseudo-vector that is equivalent to `old` (represents the same rotation), but is as close as possible to `prev`.

This helps enforce continuity in incremental rotation updates by accounting for 2π ambiguities in the exponential map.
"""
function rv_continue(old::SVector{3,Float64}, prev::SVector{3,Float64})::SVector{3,Float64}
    norm_old = LinearAlgebra.norm(old)
    if norm_old > 0.0
        unit = LinearAlgebra.normalize(old)
        proj = LinearAlgebra.dot(unit, prev)
    else
        unit = LinearAlgebra.normalize(prev)
        proj = LinearAlgebra.norm(prev)
    end
    if proj == 0.0
        return old
    end
    kk = round(0.5 * proj / π)
    if kk == 0.0
        return old
    end
    proj = 2.0 * kk * π + norm_old
    return proj * unit
end
