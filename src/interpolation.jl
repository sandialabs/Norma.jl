# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using StaticArrays

function barycentric_shape_functions(::Val{2}, ::Val{3}, ξ::SVector{2,T}) where {T}
    N = @SVector [one(T) - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = @SMatrix [
        -one(T) one(T) zero(T)
        -one(T) zero(T) one(T)
    ]
    ddN = zeros(SArray{Tuple{2,2,3},T,3})
    return N, dN, ddN
end

function barycentric_shape_functions(::Val{2}, ::Val{6}, ξ::SVector{2,T}) where {T}
    t0 = one(T) - ξ[1] - ξ[2]
    t1 = ξ[1]
    t2 = ξ[2]
    N = @SVector [
        t0 * (2t0 - one(T)), #node 0 - corner
        t1 * (2t1 - one(T)), #node 1 - corner
        t2 * (2t2 - one(T)), #node 2 - corner
        4t0 * t1,            #node 3 - middle
        4t1 * t2,            #node 4 - middle
        4t2 * t0,            #node 5 - middle
    ]
    dN = @SMatrix [
        one(T)-4t0 4t1-one(T) zero(T) 4t0-4t1 4t2 -4t2
        one(T)-4t0 zero(T) 4t2-one(T) -4t1 4t1 4t0-4t2
    ]
    ddN = MArray{Tuple{2,2,6},Float64,3}(undef)
    ddN[1, 1, :] = @SVector [4, 4, 0, -8, 0, 0]
    ddN[1, 2, :] = @SVector [4, 0, 0, -4, 4, -4]
    ddN[2, 1, :] = @SVector [4, 0, 0, -4, 4, -4]
    ddN[2, 2, :] = @SVector [4, 0, 4, 0, 0, -8]
    return N, dN, SArray(ddN)
end

function barycentric_shape_functions(::Val{3}, ::Val{4}, ξ::SVector{3,T}) where {T}
    N = @SVector [one(T) - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3]]
    dN = @SMatrix [
        -one(T) one(T) zero(T) zero(T)
        -one(T) zero(T) one(T) zero(T)
        -one(T) zero(T) zero(T) one(T)
    ]
    ddN = zeros(SArray{Tuple{3,3,4},T,3})
    return N, dN, ddN
end

function barycentric_shape_functions(::Val{3}, ::Val{10}, ξ::SVector{3,T}) where {T}
    t0 = one(T) - ξ[1] - ξ[2] - ξ[3]
    t1 = ξ[1]
    t2 = ξ[2]
    t3 = ξ[3]
    N = @SVector [
        t0 * (2t0 - one(T)),  # node 1
        t1 * (2t1 - one(T)),  # node 2
        t2 * (2t2 - one(T)),  # node 3
        t3 * (2t3 - one(T)),  # node 4
        4t0 * t1,             # node 5
        4t1 * t2,             # node 6
        4t2 * t0,             # node 7
        4t0 * t3,             # node 8
        4t1 * t3,             # node 9
        4t2 * t3,             # node 10
    ]
    dN = @SMatrix [
        (one(T)-4t0) (4t1-one(T)) zero(T) zero(T) 4(t0 - t1) 4t2 -4t2 -4t3 4t3 zero(T)
        (one(T)-4t0) zero(T) (4t2-one(T)) zero(T) -4t1 4t1 4(t0 - t2) -4t3 zero(T) 4t3
        (one(T)-4t0) zero(T) zero(T) (4t3-one(T)) -4t1 zero(T) -4t2 4(t0 - t3) 4t1 4t2
    ]
    ddN = MArray{Tuple{3,3,10},Float64,3}(undef)
    ddN[1, 1, :] = @SVector [4, 4, 0, 0, -8, 0, 0, 0, 0, 0]
    ddN[1, 2, :] = @SVector [4, 0, 0, 0, -4, 4, -4, 0, 0, 0]
    ddN[1, 3, :] = @SVector [4, 0, 0, 0, -4, 0, 0, -4, 4, 0]
    ddN[2, 1, :] = @SVector [4, 0, 0, 0, -4, 4, -4, 0, 0, 0]
    ddN[2, 2, :] = @SVector [4, 0, 4, 0, 0, 0, -8, 0, 0, 0]
    ddN[2, 3, :] = @SVector [4, 0, 0, 0, 0, 0, -4, -4, 0, 4]
    ddN[3, 1, :] = @SVector [4, 0, 0, 0, -4, 0, 0, -4, 4, 0]
    ddN[3, 2, :] = @SVector [4, 0, 0, 0, 0, 0, -4, -4, 0, 4]
    ddN[3, 3, :] = @SVector [4, 0, 0, 4, 0, 0, 0, -8, 0, 0]
    return N, dN, SArray(ddN)
end

function barycentric_quadrature(::Val{2}, ::Val{3}, ::Val{1})
    ξ = @SMatrix [1 / 3; 1 / 3]
    w = @SVector [1 / 2]
    return ξ, w
end

function barycentric_quadrature(::Val{2}, ::Val{3}, ::Val{3})
    ξ = @SMatrix [
        1/6 4/6 1/6
        1/6 1/6 4/6
    ]
    w = @SVector [1 / 6, 1 / 6, 1 / 6]
    return ξ, w
end

function barycentric_quadrature(::Val{2}, ::Val{6}, ::Val{3})
    ξ = @SMatrix [
        1/6 4/6 1/6
        1/6 1/6 4/6
    ]
    w = @SVector [1 / 6, 1 / 6, 1 / 6]
    return ξ, w
end

function barycentric_quadrature(::Val{2}, ::Val{6}, ::Val{4})
    ξ = @SMatrix [
        1/3 1/5 3/5 1/5
        1/3 3/5 1/5 1/5
    ]
    w = @SVector [-27 / 96, 25 / 96, 25 / 96, 25 / 96]
    return ξ, w
end

# Dunavant degree-5 rule for triangles: 7 points, all positive weights.
# Exact for polynomials up to degree 5 (needed for TRI6 mass matrix, which
# requires degree 4 accuracy for N_i * N_j).  With 7 > 6 nodes this rule
# also ensures the per-element mass matrix Σ_p w_p N_p N_p^T |J_p| is
# positive definite and hence non-singular — unlike the 4-point rule above
# which yields a rank-deficient (rank ≤ 4) 6×6 element matrix.
# Source: D.A. Dunavant, IJNME 21:1129-1148, 1985.
function barycentric_quadrature(::Val{2}, ::Val{6}, ::Val{7})
    a1 = (6 - sqrt(15)) / 21      # ≈ 0.10128650732345633
    b1 = 1 - 2 * a1               # ≈ 0.79742698535308734
    a2 = (6 + sqrt(15)) / 21      # ≈ 0.47014206410511508
    b2 = 1 - 2 * a2               # ≈ 0.05971587178976984
    w0 = 9 / 80                   # centroid weight
    w1 = (155 - sqrt(15)) / 2400  # ≈ 0.06296959027241357
    w2 = (155 + sqrt(15)) / 2400  # ≈ 0.06619707767059810
    ξ = @SMatrix [
        1/3 a1 b1 a1 a2 b2 a2
        1/3 a1 a1 b1 a2 a2 b2
    ]
    w = @SVector [w0, w1, w1, w1, w2, w2, w2]
    return ξ, w
end

function barycentric_quadrature(::Val{3}, ::Val{4}, ::Val{1})
    ξ = @SMatrix [1 / 4; 1 / 4; 1 / 4]
    w = @SVector [1 / 6]
    return ξ, w
end

function barycentric_quadrature(::Val{3}, ::Val{4}, ::Val{4})
    s = sqrt(5)
    a = (5 + 3s) / 20
    b = (5 - s) / 20
    ξ = @SMatrix [
        b a b b
        b b a b
        b b b a
    ]
    w = @SVector [1 / 24, 1 / 24, 1 / 24, 1 / 24]
    return ξ, w
end

function barycentric_quadrature(::Val{3}, ::Val{10}, ::Val{4})
    s = sqrt(5)
    a = (5 + 3s) / 20
    b = (5 - s) / 20
    ξ = @SMatrix [
        b a b b
        b b a b
        b b b a
    ]
    w = @SVector [1 / 24, 1 / 24, 1 / 24, 1 / 24]
    return ξ, w
end

function barycentric_quadrature(::Val{3}, ::Val{10}, ::Val{5})
    ξ = @SMatrix [
        1/4 1/6 1/6 1/6 1/2
        1/4 1/6 1/6 1/2 1/6
        1/4 1/6 1/2 1/6 1/6
    ]
    w = @SVector [-2 / 15, 3 / 40, 3 / 40, 3 / 40, 3 / 40]
    return ξ, w
end

# Keast rule 6: degree 4, 14 points, all positive weights — required for a PD
# consistent mass on TET10 (the 4-pt and 5-pt rules give rank-deficient or
# sign-indefinite element mass).  Source: Patrick Keast, "Moderate Degree
# Tetrahedral Quadrature Formulas," CMAME 55(3):339-348, 1986.  Coefficients
# transcribed from John Burkardt's `tetrahedron_keast_rule` (LGPL).
function barycentric_quadrature(::Val{3}, ::Val{10}, ::Val{14})
    A2 = 0.698419704324386603
    B2 = 0.100526765225204467
    A3 = 0.0568813795204234229
    B3 = 0.314372873493192195
    w1 = 0.00317460317460317450
    w2 = 0.0147649707904967828
    w3 = 0.0221397911142651221
    ξ = @SMatrix [
        0.0 0.5 0.5 0.5 0.0 0.0 A2 B2 B2 B2 A3 B3 B3 B3
        0.5 0.0 0.5 0.0 0.5 0.0 B2 B2 B2 A2 B3 B3 B3 A3
        0.5 0.5 0.0 0.0 0.0 0.5 B2 B2 A2 B2 B3 B3 A3 B3
    ]
    w = @SVector [w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w3, w3, w3, w3]
    return ξ, w
end

function barycentric(::Val{D}, ::Val{N}, ::Val{G}) where {D,N,G}
    ξ, w = barycentric_quadrature(Val(D), Val(N), Val(G))
    T = eltype(ξ)
    Nmat = MMatrix{N,G,T}(undef)
    dNmat = MArray{Tuple{D,N,G},T}(undef)
    for p in 1:G
        ξp = SVector(ξ[:, p])
        Np, dNp, _ = barycentric_shape_functions(Val(D), Val(N), ξp)
        @inbounds for i in 1:N
            Nmat[i, p] = Np[i]
        end
        @inbounds for i in 1:D
            for j in 1:N
                dNmat[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(Nmat), SArray(dNmat), w, ξ
end

function lagrangian_shape_functions(::Val{1}, ::Val{2}, ξ::SVector{1,T}) where {T}
    N = @SVector [0.5 * (1.0 - ξ[1]), 0.5 * (1.0 + ξ[1])]
    dN = @SMatrix [-0.5 0.5]
    ddN = zeros(SArray{Tuple{1,1,2},T,3})
    return N, dN, ddN
end

function lagrangian_shape_functions(::Val{2}, ::Val{4}, ξ::SVector{2,T}) where {T}
    r, s = ξ
    ra = @SVector [-1, 1, 1, -1]
    sa = @SVector [-1, -1, 1, 1]
    N = MVector{4,T}(undef)
    dN = MMatrix{2,4,T}(undef)
    ddN = MArray{Tuple{2,2,4},T,3}(undef)
    for p in 1:4
        N[p] = 0.25 * (1 + ra[p] * r) * (1 + sa[p] * s)
        dN[1, p] = 0.25 * ra[p] * (1 + sa[p] * s)
        dN[2, p] = 0.25 * sa[p] * (1 + ra[p] * r)
        ddN[1, 1, p] = 0
        ddN[1, 2, p] = 0.25 * ra[p] * sa[p]
        ddN[2, 1, p] = ddN[1, 2, p]
        ddN[2, 2, p] = 0
    end
    return SVector(N), SMatrix(dN), SArray(ddN)
end

function lagrangian_shape_functions(::Val{3}, ::Val{8}, ξ::SVector{3,T}) where {T}
    r, s, t = ξ
    ra = @SVector [-1, 1, 1, -1, -1, 1, 1, -1]
    sa = @SVector [-1, -1, 1, 1, -1, -1, 1, 1]
    ta = @SVector [-1, -1, -1, -1, 1, 1, 1, 1]
    N = MVector{8,T}(undef)
    dN = MMatrix{3,8,T}(undef)
    ddN = MArray{Tuple{3,3,8},T,3}(undef)
    for p in 1:8
        r_p, s_p, t_p = ra[p], sa[p], ta[p]
        N[p] = 0.125 * (1 + r_p * r) * (1 + s_p * s) * (1 + t_p * t)
        dN[1, p] = 0.125 * r_p * (1 + s_p * s) * (1 + t_p * t)
        dN[2, p] = 0.125 * (1 + r_p * r) * s_p * (1 + t_p * t)
        dN[3, p] = 0.125 * (1 + r_p * r) * (1 + s_p * s) * t_p
        ddN[1, 1, p] = zero(T)
        ddN[1, 2, p] = 0.125 * r_p * s_p * (1 + t_p * t)
        ddN[1, 3, p] = 0.125 * r_p * t_p * (1 + s_p * s)
        ddN[2, 1, p] = ddN[1, 2, p]
        ddN[2, 2, p] = zero(T)
        ddN[2, 3, p] = 0.125 * s_p * t_p * (1 + r_p * r)
        ddN[3, 1, p] = ddN[1, 3, p]
        ddN[3, 2, p] = ddN[2, 3, p]
        ddN[3, 3, p] = zero(T)
    end
    return SVector(N), SMatrix(dN), SArray(ddN)
end

function lagrangian_quadrature(::Val{1}, ::Val{2}, ::Val{1})
    ξ = @SMatrix [0.0]
    w = @SVector [2.0]
    return ξ, w
end

function lagrangian_quadrature(::Val{1}, ::Val{2}, ::Val{2})
    g = 1 / sqrt(3)
    ξ = @SMatrix [-g g]
    w = @SVector [1.0, 1.0]
    return ξ, w
end

function lagrangian_quadrature(::Val{2}, ::Val{4}, ::Val{4})
    g = 1 / sqrt(3)
    ξ = @SMatrix [
        -g g g -g
        -g -g g g
    ]
    w = @SVector [1.0, 1.0, 1.0, 1.0]
    return ξ, w
end

function lagrangian_quadrature(::Val{2}, ::Val{4}, ::Val{9})
    g = sqrt(3 / 5)
    ξ = @SMatrix [
        -g g g -g 0 g 0 -g 0
        -g -g g g -g 0 g 0 0
    ]
    w = @SVector [25 / 81, 25 / 81, 25 / 81, 25 / 81, 40 / 81, 40 / 81, 40 / 81, 40 / 81, 64 / 81]
    return ξ, w
end

function lagrangian_quadrature(::Val{3}, ::Val{8}, ::Val{8})
    g = 1 / sqrt(3)
    ξ = @SMatrix [
        -g g g -g -g g g -g
        -g -g g g -g -g g g
        -g -g -g -g g g g g
    ]
    w = @SVector [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    return ξ, w
end

function lagrangian(::Val{D}, ::Val{N}, ::Val{G}) where {D,N,G}
    ξ, w = lagrangian_quadrature(Val(D), Val(N), Val(G))
    T = eltype(ξ)
    Nmat = MMatrix{N,G,T}(undef)
    dNmat = MArray{Tuple{D,N,G},T}(undef)
    for p in 1:G
        ξp = SVector(ξ[:, p])
        Np, dNp, _ = lagrangian_shape_functions(Val(D), Val(N), ξp)
        @inbounds for i in 1:N
            Nmat[i, p] = Np[i]
        end
        @inbounds for i in 1:D
            for j in 1:N
                dNmat[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(Nmat), SArray(dNmat), w, ξ
end

function element_type_from_string(s::AbstractString)::ElementType
    if s == "BAR2" || s == "BAR"
        return BAR2
    elseif s == "TRI3" || s == "TRI"
        return TRI3
    elseif s == "TRI6"
        return TRI6
    elseif s == "QUAD4" || s == "QUAD"
        return QUAD4
    elseif s == "TETRA4" || s == "TETRA"
        return TETRA4
    elseif s == "TETRA10"
        return TETRA10
    elseif s == "HEX8" || s == "HEX"
        return HEX8
    else
        norma_abort("Unknown element type string: $s")
    end
end

default_num_int_pts(::Val{BAR2}) = 1
default_num_int_pts(::Val{TRI3}) = 3
# 7-point Dunavant degree-5 rule: all positive weights and 7 > 6 nodes so
# the element mass matrix Σ_p w_p N_p N_p^T is guaranteed full-rank for
# TRI6.  The former default of 4 produced a rank-deficient (rank ≤ 4) 6×6
# element mass matrix, causing a SingularException when the non-overlapping
# Schwarz Neumann projector H = ∫ N N^T dA was inverted with H \ I.
default_num_int_pts(::Val{TRI6}) = 7
default_num_int_pts(::Val{QUAD4}) = 4
default_num_int_pts(::Val{TETRA4}) = 4
default_num_int_pts(::Val{TETRA10}) = 5
default_num_int_pts(::Val{HEX8}) = 8
default_num_int_pts(et::ElementType) = default_num_int_pts(Val(et))

get_element_dim_nodes(::Val{BAR2}) = (1, 2)
get_element_dim_nodes(::Val{TRI3}) = (2, 3)
get_element_dim_nodes(::Val{TRI6}) = (2, 6)
get_element_dim_nodes(::Val{QUAD4}) = (2, 4)
get_element_dim_nodes(::Val{TETRA4}) = (3, 4)
get_element_dim_nodes(::Val{TETRA10}) = (3, 10)
get_element_dim_nodes(::Val{HEX8}) = (3, 8)
get_element_dim_nodes(et::ElementType) = get_element_dim_nodes(Val(et))

is_barycentric(::Val{BAR2}) = false
is_barycentric(::Val{TRI3}) = true
is_barycentric(::Val{TRI6}) = true
is_barycentric(::Val{QUAD4}) = false
is_barycentric(::Val{TETRA4}) = true
is_barycentric(::Val{TETRA10}) = true
is_barycentric(::Val{HEX8}) = false
is_barycentric(et::ElementType) = is_barycentric(Val(et))

function get_element_type(dim::Integer, num_nodes::Integer)
    if dim == 1 && num_nodes == 2
        return BAR2
    elseif dim == 2 && num_nodes == 3
        return TRI3
    elseif dim == 2 && num_nodes == 6
        return TRI6
    elseif dim == 2 && num_nodes == 4
        return QUAD4
    elseif dim == 3 && num_nodes == 4
        return TETRA4
    elseif dim == 3 && num_nodes == 10
        return TETRA10
    elseif dim == 3 && num_nodes == 8
        return HEX8
    else
        norma_abort("Invalid dim = $dim and num_nodes = $num_nodes")
    end
end

function isoparametric(element_type::ElementType, G::Integer)
    D, N = get_element_dim_nodes(element_type)
    return if is_barycentric(element_type)
        barycentric(Val(D), Val(N), Val(G))
    else
        lagrangian(Val(D), Val(N), Val(G))
    end
end

function interpolate(element_type::ElementType, ξ::AbstractVector{Float64})
    D, N = get_element_dim_nodes(element_type)
    ξs = SVector{D}(ξ)
    if is_barycentric(element_type)
        return barycentric_shape_functions(Val(D), Val(N), ξs)
    else
        return lagrangian_shape_functions(Val(D), Val(N), ξs)
    end
end

using Symbolics
@variables t, x, y, z

function get_side_set_nodal_forces(nodal_coord::Matrix{Float64}, traction_fun::Function, time::Float64)
    _, num_side_nodes = size(nodal_coord)
    element_type = get_element_type(2, num_side_nodes)
    num_int_points = default_num_int_pts(element_type)
    N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
    nodal_force_component = zeros(num_side_nodes)
    for point in 1:num_int_points
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * nodal_coord'
        j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
        wₚ = w[point]
        point_coord = nodal_coord * Nₚ
        txzy = (time, point_coord[1], point_coord[2], point_coord[3])
        traction_val = traction_fun(txzy)
        nodal_force_component += traction_val * Nₚ * j * wₚ
    end
    return nodal_force_component
end

function get_side_set_nodal_stiffness(nodal_coord::Matrix{Float64}, time::Float64)
    _, num_side_nodes = size(nodal_coord)
    element_type = get_element_type(2, num_side_nodes)
    num_int_points = default_num_int_pts(element_type)
    N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
    nodal_stiffness_component = zeros(num_side_nodes, num_side_nodes)
    for point in 1:num_int_points
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * nodal_coord'
        j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
        wₚ = w[point]
        nodal_stiffness_component += Nₚ * Nₚ' * j * wₚ
    end
    return nodal_stiffness_component
end

function get_side_set_nodal_pressure(nodal_coord::Matrix{Float64}, pressure_fun::Function, time::Float64)
    _, num_side_nodes = size(nodal_coord)
    element_type = get_element_type(2, num_side_nodes)
    num_int_points = default_num_int_pts(element_type)
    N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
    nodal_force_component = zeros(3, num_side_nodes)
    for point in 1:num_int_points
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * nodal_coord'
        perp_vector = cross(dXdξ[1, :], dXdξ[2, :])
        normal = LinearAlgebra.normalize(perp_vector)
        j = norm(perp_vector)
        wₚ = w[point]
        point_coord = nodal_coord * Nₚ
        txzy = (time, point_coord[1], point_coord[2], point_coord[3])
        pressure_val = pressure_fun(txzy)
        nodal_force_component += normal * (pressure_val * Nₚ * j * wₚ)'
    end
    return nodal_force_component
end

function map_to_parametric(element_type::ElementType, nodes::Matrix{Float64}, point::Vector{Float64})
    tol = 1.0e-08
    max_iters = 1024
    ξ = MVector{3,Float64}(0.0, 0.0, 0.0)
    hessian = MMatrix{3,3,Float64}(undef)
    for _ in 1:max_iters
        N, dN, _ = interpolate(element_type, ξ)
        trial_point = nodes * N
        residual = trial_point - point
        hessian = nodes * dN'
        δ = -hessian \ residual
        ξ = ξ + δ
        err = norm(δ)
        if err <= tol
            break
        end
    end
    return ξ
end

function is_inside_parametric(element_type::ElementType, ξ::AbstractVector{Float64}, tol::Float64=1e-6)
    factor = 1.0 + tol
    if element_type == BAR2
        return -factor ≤ ξ[1] ≤ factor
    elseif element_type == TRI3 || element_type == TRI6 || element_type == TETRA4 || element_type == TETRA10
        return sum(ξ) ≤ factor
    elseif element_type == QUAD4
        return all(-factor .≤ ξ[1:2] .≤ factor)
    elseif element_type == HEX8
        return all(-factor .≤ ξ[1:3] .≤ factor)
    else
        norma_abort("Unsupported element type: $element_type")
    end
end

function is_inside(element_type::ElementType, nodes::Matrix{Float64}, point::Vector{Float64}, tol::Float64=1.0e-06)
    ξ = @SVector [0.0, 0.0, 0.0]
    if in_bounding_box(nodes, point, 0.1) == false
        return ξ, false
    end
    ξ = map_to_parametric(element_type, nodes, point)
    return ξ, is_inside_parametric(element_type, ξ, tol)
end

function in_bounding_box(nodes::Matrix{Float64}, point::Vector{Float64}, tol::Float64=1.0e-06)::Bool
    for i in 1:3
        coord_min = minimum(nodes[i, :])
        coord_max = maximum(nodes[i, :])
        range_d = coord_max - coord_min
        lower_bound = coord_min - tol * range_d
        upper_bound = coord_max + tol * range_d
        if point[i] < lower_bound || point[i] > upper_bound
            return false
        end
    end
    return true
end

function closest_face_to_point(point::Vector{Float64}, model::SolidMechanics, side_set_id::Integer)
    mesh = model.mesh
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    ss_node_index = 1
    closest_face_nodes = Array{Float64}(undef, 0)
    closest_face_node_indices = Array{Int64}(undef, 0)
    minimum_nodal_distance = Inf
    for num_nodes_side in num_nodes_sides
        face_node_indices = side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)]
        face_nodes = model.reference[:, face_node_indices] + model.displacement[:, face_node_indices]
        nodal_distance = get_minimum_distance_to_nodes(face_nodes, point)
        if nodal_distance < minimum_nodal_distance
            minimum_nodal_distance = nodal_distance
            closest_face_nodes = face_nodes
            closest_face_node_indices = face_node_indices
        end
        ss_node_index += num_nodes_side
    end
    return closest_face_nodes, closest_face_node_indices, minimum_nodal_distance
end

# Find the minimum distance of a point to the nodes of each face on the side set
# and then project the point that closest face in the side set.
# This is done in place of a strict search because the contact surfaces may be deformed
# and not match each other exactly. We assume that we know the contact surfaces in advance
function project_point_to_side_set(point::Vector{Float64}, model::SolidMechanics, side_set_id::Integer)
    face_nodes, face_node_indices, _ = closest_face_to_point(point, model, side_set_id)
    new_point, ξ, surface_distance, normal = closest_point_projection(face_nodes, point)
    return new_point, ξ, face_nodes, face_node_indices, normal, surface_distance
end

function project_point_to_side_set(
    point::Vector{Float64}, face_nodes::Matrix{Float64}, face_node_indices::Vector{Int32}
)
    new_point, ξ, surface_distance, normal = closest_point_projection(face_nodes, point)
    return new_point, ξ, face_nodes, face_node_indices, normal, surface_distance
end

function get_minimum_distance_to_nodes(nodes::Matrix{Float64}, point::Vector{Float64})
    distances = norm.(eachcol(nodes) .- Ref(point))
    return minimum(distances)
end

# How far outside the parametric domain a projected point landed (0 = inside).
function parametric_overshoot(element_type::ElementType, ξ::AbstractVector{Float64})
    if element_type == BAR2
        return max(0.0, abs(ξ[1]) - 1.0)
    elseif element_type == TRI3 || element_type == TRI6
        return max(0.0, -ξ[1], -ξ[2], ξ[1] + ξ[2] - 1.0)
    elseif element_type == QUAD4
        return max(0.0, abs(ξ[1]) - 1.0, abs(ξ[2]) - 1.0)
    else
        norma_abort("Unsupported side element type: $element_type")
    end
end

# Static-array version of the closest point projection for a facet with N
# nodes: the same Newton iteration as the array version below, without heap
# allocation. Used by the containing-facet search, which calls it for every
# candidate facet of every destination integration point.
function closest_point_projection(nodes::SMatrix{3,N,Float64}, x::SVector{3,Float64}; strict::Bool=true) where {N}
    element_type = get_element_type(2, N)
    ξ = zero(SVector{2,Float64})
    y = zero(SVector{3,Float64})
    yx = zero(SVector{3,Float64})
    tol = 1.0e-10
    iteration = 1
    max_iterations = 64
    converged = false
    while true
        Nξ, dN, ddN = interpolate(element_type, ξ)
        y = nodes * Nξ
        dydξ = dN * nodes'
        yx = y - x
        residual = dydξ * yx
        ddyddξyx = SMatrix{2,2,Float64,4}(ntuple(m -> begin
            i = (m - 1) % 2 + 1
            j = (m - 1) ÷ 2 + 1
            total = 0.0
            @inbounds for k in 1:3
                row = 0.0
                for l in 1:N
                    row += ddN[i, j, l] * nodes[k, l]
                end
                total += row * yx[k]
            end
            total
        end, Val(4)))
        hessian = ddyddξyx + dydξ * dydξ'
        δ = -(hessian \ residual)
        ξ = ξ + δ
        if norm(δ) <= tol
            converged = true
            break
        end
        iteration += 1
        if iteration > max_iterations
            break
        end
    end
    if converged == false
        if strict == true
            norma_abort("Closest point projection failed to converge")
        end
        return y, ξ, Inf, zero(SVector{3,Float64})
    end
    _, dN, _ = interpolate(element_type, ξ)
    dxdξ = dN * nodes'
    perp_vec = cross(dxdξ[1, :], dxdξ[2, :])
    normal = perp_vec / norm(perp_vec)
    distance = -copysign(norm(yx), dot(yx, normal))
    return y, ξ, distance, normal
end

# Side set data gathered once for a containing-facet search over many points:
# the side set node list and current coordinates were read from the Exodus
# file and rebuilt for every point before, which made the interface projector
# setup quadratic with a large constant.
struct SideSetSearch
    node_coords::Vector{SVector{3,Float64}}      # one entry per side set node list entry
    facet_node_indices::Vector{Vector{Int64}}
    facet_coords::Vector{Matrix{Float64}}
    box_min::Vector{SVector{3,Float64}}
    box_max::Vector{SVector{3,Float64}}
    half_diagonal::Vector{Float64}
end

function build_side_set_search(model::SolidMechanics, side_set_id::Integer)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(model.mesh, side_set_id)
    coords = model.reference + model.displacement
    node_coords = [SVector{3,Float64}(view(coords, :, i)) for i in side_set_node_indices]
    num_facets = length(num_nodes_sides)
    facet_node_indices = Vector{Vector{Int64}}(undef, num_facets)
    facet_coords = Vector{Matrix{Float64}}(undef, num_facets)
    box_min = Vector{SVector{3,Float64}}(undef, num_facets)
    box_max = Vector{SVector{3,Float64}}(undef, num_facets)
    half_diagonal = Vector{Float64}(undef, num_facets)
    ss_node_index = 1
    for (facet, num_nodes_side) in enumerate(num_nodes_sides)
        indices = Int64.(side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)])
        ss_node_index += num_nodes_side
        nodes = coords[:, indices]
        lo = SVector{3,Float64}(minimum(nodes[i, :]) for i in 1:3)
        hi = SVector{3,Float64}(maximum(nodes[i, :]) for i in 1:3)
        facet_node_indices[facet] = indices
        facet_coords[facet] = nodes
        box_min[facet] = lo
        box_max[facet] = hi
        half_diagonal[facet] = 0.5 * norm(hi - lo)
    end
    return SideSetSearch(node_coords, facet_node_indices, facet_coords, box_min, box_max, half_diagonal)
end

# Function barrier on the facet node count for the static projection.
function project_onto_facet(nodes::Matrix{Float64}, point::SVector{3,Float64})
    n = size(nodes, 2)
    if n == 4
        return closest_point_projection(SMatrix{3,4,Float64,12}(nodes), point; strict=false)
    elseif n == 3
        return closest_point_projection(SMatrix{3,3,Float64,9}(nodes), point; strict=false)
    else
        y, ξ, distance, normal = closest_point_projection(nodes, Vector(point); strict=false)
        return SVector{3,Float64}(y), SVector{2,Float64}(ξ), distance, SVector{3,Float64}(normal)
    end
end

function project_point_to_containing_facet(point::SVector{3,Float64}, search::SideSetSearch)
    inside_tol = 1.0e-06
    best_key = (Inf, Inf)
    best_facet = 0
    best_point = zero(SVector{3,Float64})
    best_ξ = zero(SVector{2,Float64})
    best_normal = zero(SVector{3,Float64})
    best_distance = Inf
    # See the array version below for the rationale of the node gap margin.
    node_gap = Inf
    @inbounds for c in search.node_coords
        node_gap = min(node_gap, norm(point - c))
    end
    @inbounds for facet in eachindex(search.facet_coords)
        margin = node_gap + search.half_diagonal[facet]
        lo = search.box_min[facet]
        hi = search.box_max[facet]
        outside = false
        for i in 1:3
            if point[i] < lo[i] - margin || point[i] > hi[i] + margin
                outside = true
                break
            end
        end
        outside && continue
        nodes = search.facet_coords[facet]
        new_point, ξ, distance, normal = project_onto_facet(nodes, point)
        isfinite(distance) || continue
        element_type = get_element_type(2, size(nodes, 2))
        overshoot = parametric_overshoot(element_type, ξ)
        key = (overshoot ≤ inside_tol ? 0.0 : overshoot, abs(distance))
        if key < best_key
            best_key = key
            best_facet = facet
            best_point, best_ξ, best_normal, best_distance = new_point, ξ, normal, distance
        end
    end
    if best_facet == 0
        # Every facet's projection diverged (strongly deformed interface):
        # fall back to the nearest-node facet.
        nearest = 0
        nearest_distance = Inf
        for facet in eachindex(search.facet_coords)
            d = get_minimum_distance_to_nodes(search.facet_coords[facet], Vector(point))
            if d < nearest_distance
                nearest_distance = d
                nearest = facet
            end
        end
        nodes = search.facet_coords[nearest]
        y, ξ, distance, normal = closest_point_projection(nodes, Vector(point))
        return y, ξ, nodes, search.facet_node_indices[nearest], normal, distance
    end
    return Vector(best_point),
    Vector(best_ξ), search.facet_coords[best_facet], search.facet_node_indices[best_facet], Vector(best_normal),
    best_distance
end

# Containing-facet projection for interface transfer assembly. The
# nearest-node heuristic of closest_face_to_point can select a facet that
# does not contain the point whenever the two side sets are refined
# differently: near a shared corner the nearest node's facet is frequently
# the wrong one (with exact ties broken by traversal order), and the
# unclamped Newton projection then extrapolates that facet's shape functions
# (values outside [0, 1]), silently corrupting the transfer matrix. Constants
# and linear fields hide the corruption on flat interfaces because their
# extrapolation is still exact. Here every facet of the side set is scanned:
# facets whose parametric projection lies inside are preferred, with the true
# surface distance as tie-breaker; if no facet contains the point (gapped or
# strongly deformed interfaces), the facet with least parametric overshoot is
# used, which is the closest facet region.
function project_point_to_containing_facet(point::Vector{Float64}, model::SolidMechanics, side_set_id::Integer)
    mesh = model.mesh
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    coords = model.reference + model.displacement
    inside_tol = 1.0e-06
    best_key = (Inf, Inf)
    local best_point, best_ξ, best_nodes, best_indices, best_normal, best_distance
    found = false
    # Distance from the point to the nearest side-set node. The closest-point
    # projection onto the containing facet lands no farther than its nearest
    # node, so this upper-bounds the projection distance and lets the prefilter
    # below prune facets that cannot possibly contain the point without ever
    # discarding the one that does — even across a wide tied-contact gap, where
    # the destination node sits a full facet-width off the source surface.
    node_gap = minimum(norm(point - coords[:, i]) for i in side_set_node_indices)
    ss_node_index = 1
    for num_nodes_side in num_nodes_sides
        face_node_indices = side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)]
        ss_node_index += num_nodes_side
        face_nodes = coords[:, face_node_indices]
        # A facet can contain the point's projection only if the point lies
        # within (nearest-node distance + half the facet diagonal) of its
        # bounding box. The node-gap term admits legitimately distant tied and
        # gapped interfaces; the diagonal term keeps the flat dimension of a
        # planar facet from being pinched to zero by floating-point fuzz. On
        # curved side sets this also skips distant, badly oriented facets on
        # which the unclamped Newton projection can diverge.
        margin = node_gap + 0.5 * norm([maximum(face_nodes[i, :]) - minimum(face_nodes[i, :]) for i in 1:3])
        outside = any(
            point[i] < minimum(face_nodes[i, :]) - margin || point[i] > maximum(face_nodes[i, :]) + margin for i in 1:3
        )
        outside && continue
        # strict=false reports a diverging Newton iteration as an infinite
        # distance so the candidate is disqualified instead of aborting the run.
        new_point, ξ, distance, normal = closest_point_projection(face_nodes, point; strict=false)
        isfinite(distance) || continue
        element_type = get_element_type(2, Int64(num_nodes_side))
        overshoot = parametric_overshoot(element_type, ξ)
        key = (overshoot ≤ inside_tol ? 0.0 : overshoot, abs(distance))
        if key < best_key
            best_key = key
            best_point, best_ξ, best_nodes, best_indices, best_normal, best_distance = new_point,
            ξ, face_nodes, face_node_indices, normal,
            distance
            found = true
        end
    end
    if found == false
        # Every facet's projection diverged (strongly deformed interface):
        # fall back to the nearest-node facet search.
        face_nodes, face_node_indices, _ = closest_face_to_point(point, model, side_set_id)
        new_point, ξ, distance, normal = closest_point_projection(face_nodes, point)
        return new_point, ξ, face_nodes, face_node_indices, normal, distance
    end
    return best_point, best_ξ, best_nodes, best_indices, best_normal, best_distance
end

# The maps between the global indices of the nodes of a side set and their
# local indices, numbered in the order of first appearance in the node list
# of the side set: local from global as a dictionary, and global from local
# as a vector.  The node list is read once.
function get_side_set_node_maps(mesh::ExodusDatabase, side_set_id::Integer)
    global_from_local_map = Int64.(unique(Exodus.read_side_set_node_list(mesh, side_set_id)[2]))
    local_from_global_map = Dict{Int64,Int64}()
    for (i, node) in enumerate(global_from_local_map)
        local_from_global_map[node] = i
    end
    return local_from_global_map, global_from_local_map
end

function get_side_set_local_from_global_map(mesh::ExodusDatabase, side_set_id::Integer)
    return get_side_set_node_maps(mesh, side_set_id)[1]
end

function get_square_projection_matrix(model::SolidMechanics, bc::SolidMechanicsSchwarzBoundaryCondition)
    num_nodes_sides = bc.num_nodes_sides
    side_set_node_indices = bc.side_set_node_indices
    local_from_global_map = bc.local_from_global_map
    num_nodes = length(local_from_global_map)
    coords = model.reference
    square_projection_matrix = zeros(num_nodes, num_nodes)
    side_set_node_index = 1
    for num_nodes_side in num_nodes_sides
        side_nodes = side_set_node_indices[side_set_node_index:(side_set_node_index + num_nodes_side - 1)]
        side_coordinates = coords[:, side_nodes]
        side_matrix = get_side_set_nodal_stiffness(side_coordinates, 0.0)
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        square_projection_matrix[local_indices, local_indices] += side_matrix
        side_set_node_index += num_nodes_side
    end
    return square_projection_matrix
end

# Rectangular interface projection matrix L_ij = ∫ N_dst_i N_src_j dΓ,
# integrated over the SOURCE side's facets. The facet quadrature is exact
# only where the integrand is smooth, and the two integration directions
# protect different invariants of the transfer on non-conforming interfaces:
# - Integrating over the destination facets (get_rectangular_projection_matrix)
#   makes the Dirichlet projector W⁻¹L reproduce constants exactly for any
#   quadrature, because the source partition of unity Σ_j N_src_j = 1 holds
#   pointwise at every quadrature point. Use it for Dirichlet transfers.
# - Integrating over the source facets (this function) makes the Neumann
#   projector LH⁻¹ exactly conservative — the total of any transferred nodal
#   force is preserved, since column sums ∫ N_src_j are then integrated
#   exactly while Σ_i N_dst_i = 1 holds pointwise. Use it for Neumann
#   (traction/reaction) transfers.
# When the source trace mesh is nested inside the destination facets, the
# source-side integration here is also pointwise exact. The destination-side
# rule always keeps its partition-of-unity invariant, but its entries carry
# facet-quadrature error wherever the source basis kinks inside a destination
# facet (the integrand is only piecewise smooth there).
function get_src_integrated_rectangular_projection_matrix(
    dst_model::SolidMechanics,
    dst_bc::SolidMechanicsSchwarzBoundaryCondition,
    src_model::SolidMechanics,
    src_bc::SolidMechanicsSchwarzBoundaryCondition,
)
    L, W_src = get_rectangular_projection_matrix(src_model, src_bc, dst_model, dst_bc)
    return Matrix(transpose(L)), W_src
end

# Rectangular interface projection matrix L_ij = ∫ N_dst_i N_src_j dΓ over the
# facets of the destination side set, by segment-based quadrature: each
# destination facet is split into the pieces on which the source basis is
# smooth, the intersections of the facet with the source facets in the plane
# of the destination facet, and each piece is integrated with a seven-point
# triangle rule. The product of two bilinear shape functions on planar
# parallelogram facets is a polynomial of degree four, so L is exact there
# whatever the mesh ratio, both integration directions give the same matrix,
# and W⁻¹L reproduces linear fields. A facet that the pieces do not cover
# (gap, partial overlap, or a curved interface on which a projection fails)
# is integrated with its own Gauss rule and the containing-facet projection
# of each quadrature point (add_facet_point_quadrature!), which carries the
# quadrature error of the kinks of the source basis inside the facet.
#
# Returns L and the destination's facet mass matrix W_ij = ∫ N_dst_i N_dst_j dΓ
# integrated facet by facet with the same quadrature as L. Since the source
# basis sums to one at every quadrature point, the row sums of L equal those
# of this W exactly on every facet, so W⁻¹L reproduces constants and the
# transposed transfer conserves the total force for any quadrature error.
# On planar parallelogram facets this W equals the facet mass matrix of the
# facet's own Gauss rule (get_square_projection_matrix) to roundoff; on other
# facets, where the physical-coordinate rule is inexact, the two differ by
# the quadrature error (up to 1e-2 on the paved cap facets of the nested
# cylinders), and only this W is consistent with L.
function get_rectangular_projection_matrix(
    dst_model::SolidMechanics,
    dst_bc::SolidMechanicsSchwarzBoundaryCondition,
    src_model::SolidMechanics,
    src_bc::SolidMechanicsSchwarzBoundaryCondition,
)
    src_local_from_global_map = src_bc.local_from_global_map
    src_num_nodes = length(src_local_from_global_map)
    src_side_set_id = src_bc.side_set_id
    dst_local_from_global_map = dst_bc.local_from_global_map
    dst_num_nodes_sides = dst_bc.num_nodes_sides
    dst_side_set_node_indices = dst_bc.side_set_node_indices
    dst_num_nodes = length(dst_local_from_global_map)
    dst_coords = dst_model.reference
    dst_side_set_node_index = 1
    rectangular_projection_matrix = zeros(dst_num_nodes, src_num_nodes)
    square_projection_matrix = zeros(dst_num_nodes, dst_num_nodes)
    src_search = build_side_set_search(src_model, src_side_set_id)
    num_segmented = 0
    for dst_num_nodes_side in dst_num_nodes_sides
        dst_side_nodes = dst_side_set_node_indices[dst_side_set_node_index:(dst_side_set_node_index + dst_num_nodes_side - 1)]
        dst_local_indices = get.(Ref(dst_local_from_global_map), dst_side_nodes, 0)
        dst_side_coordinates = dst_coords[:, dst_side_nodes]
        segmented = add_facet_segmented_quadrature!(
            rectangular_projection_matrix,
            square_projection_matrix,
            dst_local_indices,
            dst_side_coordinates,
            src_search,
            src_local_from_global_map,
        )
        if segmented
            num_segmented += 1
        else
            add_facet_point_quadrature!(
                rectangular_projection_matrix,
                square_projection_matrix,
                dst_local_indices,
                dst_side_coordinates,
                src_search,
                src_local_from_global_map,
            )
        end
        dst_side_set_node_index += dst_num_nodes_side
    end
    num_facets = length(dst_num_nodes_sides)
    if num_segmented < num_facets
        norma_logf(
            0,
            :warn,
            "Interface transfer '%s'/'%s': %d of %d facets are not covered by the projected facets of the other " *
            "side (gap, partial overlap, or curved interface) and were integrated with the facet's Gauss rule, " *
            "which is inexact where the facets of the two sides are not nested.",
            dst_bc.name,
            src_bc.name,
            num_facets - num_segmented,
            num_facets,
        )
    end
    return rectangular_projection_matrix, square_projection_matrix
end

# The point rule: the facet's own Gauss points, each located on the source
# side set by the containing-facet projection.
function add_facet_point_quadrature!(
    L::Matrix{Float64},
    W::Matrix{Float64},
    dst_local_indices::Vector{Int64},
    dst_side_coordinates::Matrix{Float64},
    src_search::SideSetSearch,
    src_local_from_global_map::Dict{Int64,Int64},
)
    dst_num_nodes_side = size(dst_side_coordinates, 2)
    dst_element_type = get_element_type(2, dst_num_nodes_side)
    dst_num_int_points = default_num_int_pts(dst_element_type)
    dst_N, dst_dNdξ, dst_w, _ = isoparametric(dst_element_type, dst_num_int_points)
    for dst_point in 1:dst_num_int_points
        dst_Nₚ = dst_N[:, dst_point]
        dst_dNdξₚ = dst_dNdξ[:, :, dst_point]
        dst_dXdξ = dst_dNdξₚ * dst_side_coordinates'
        dst_j = norm(cross(dst_dXdξ[1, :], dst_dXdξ[2, :]))
        dst_wₚ = dst_w[dst_point]
        dst_int_point_coord = dst_side_coordinates * dst_Nₚ
        _, ξ, src_side_coordinates, src_side_nodes, _, _ = project_point_to_containing_facet(
            SVector{3,Float64}(dst_int_point_coord), src_search
        )
        src_side_element_type = get_element_type(2, size(src_side_coordinates)[2])
        src_Nₚ, _, _ = interpolate(src_side_element_type, ξ)
        src_local_indices = get.(Ref(src_local_from_global_map), src_side_nodes, 0)
        L[dst_local_indices, src_local_indices] += dst_Nₚ * src_Nₚ' * dst_j * dst_wₚ
        W[dst_local_indices, dst_local_indices] += dst_Nₚ * dst_Nₚ' * dst_j * dst_wₚ
    end
    return nothing
end

# Relative area below which a clipped piece is discarded, and relative
# coverage error above which a facet falls back to the point rule.
const SEGMENT_AREA_TOLERANCE = 1.0e-12
const SEGMENT_COVERAGE_TOLERANCE = 1.0e-8

# Segment-based quadrature of one destination facet with three or four nodes.
# Returns false, having added nothing, when the facet is of another type or
# the pieces do not cover it.
function add_facet_segmented_quadrature!(
    L::Matrix{Float64},
    W::Matrix{Float64},
    dst_local_indices::Vector{Int64},
    dst_side_coordinates::Matrix{Float64},
    src_search::SideSetSearch,
    src_local_from_global_map::Dict{Int64,Int64},
)
    dst_num_nodes_side = size(dst_side_coordinates, 2)
    dst_num_nodes_side == 3 || dst_num_nodes_side == 4 || return false
    dst_nodes = facet_static_nodes(dst_side_coordinates)
    center, e1, e2, plane_normal = facet_plane(dst_side_coordinates)
    dst_polygon = project_polygon(dst_side_coordinates, center, e1, e2)
    dst_area = polygon_area(dst_polygon)
    dst_area > 0.0 || return false
    dst_lo = SVector{3,Float64}(minimum(dst_side_coordinates[i, :]) for i in 1:3)
    dst_hi = SVector{3,Float64}(maximum(dst_side_coordinates[i, :]) for i in 1:3)
    margin = 0.5 * norm(dst_hi - dst_lo)
    ξ_rule, w_rule = barycentric_quadrature(Val(2), Val(6), Val(7))
    contributions = Vector{Tuple{Vector{Int64},Matrix{Float64}}}()
    square_block = zeros(dst_num_nodes_side, dst_num_nodes_side)
    covered_area = 0.0
    @inbounds for src_facet in eachindex(src_search.facet_coords)
        lo = src_search.box_min[src_facet]
        hi = src_search.box_max[src_facet]
        separated = false
        for i in 1:3
            if hi[i] < dst_lo[i] - margin || lo[i] > dst_hi[i] + margin
                separated = true
                break
            end
        end
        separated && continue
        src_side_coordinates = src_search.facet_coords[src_facet]
        src_num_nodes_side = size(src_side_coordinates, 2)
        src_num_nodes_side == 3 || src_num_nodes_side == 4 || return false
        src_polygon = project_polygon(src_side_coordinates, center, e1, e2)
        polygon_area(src_polygon) > SEGMENT_AREA_TOLERANCE * dst_area || continue
        piece = clip_convex_polygon(src_polygon, dst_polygon)
        length(piece) ≥ 3 || continue
        piece_area = polygon_area(piece)
        piece_area > SEGMENT_AREA_TOLERANCE * dst_area || continue
        src_nodes = facet_static_nodes(src_side_coordinates)
        src_element_type = get_element_type(2, src_num_nodes_side)
        block = zeros(dst_num_nodes_side, src_num_nodes_side)
        centroid = sum(piece) / length(piece)
        for k in eachindex(piece)
            a = piece[k]
            b = piece[k == length(piece) ? 1 : k + 1]
            triangle_area =
                0.5 * ((a[1] - centroid[1]) * (b[2] - centroid[2]) - (a[2] - centroid[2]) * (b[1] - centroid[1]))
            triangle_area > 0.0 || continue
            for q in eachindex(w_rule)
                λ1, λ2 = ξ_rule[1, q], ξ_rule[2, q]
                p = (1.0 - λ1 - λ2) * centroid + λ1 * a + λ2 * b
                x = center + p[1] * e1 + p[2] * e2
                _, ξ_dst, dst_distance, dst_normal = project_onto_facet_nodes(dst_nodes, x)
                _, ξ_src, src_distance, _ = project_onto_facet_nodes(src_nodes, x)
                (isfinite(dst_distance) && isfinite(src_distance)) || return false
                # Area of the facet per area of its projection onto the plane.
                slope = abs(dot(dst_normal, plane_normal))
                slope > 0.5 || return false
                dst_Nₚ = interpolate(get_element_type(2, dst_num_nodes_side), ξ_dst)[1]
                src_Nₚ = interpolate(src_element_type, ξ_src)[1]
                area = 2.0 * triangle_area * w_rule[q] / slope
                block += dst_Nₚ * src_Nₚ' * area
                square_block += dst_Nₚ * dst_Nₚ' * area
            end
        end
        covered_area += piece_area
        src_local_indices = get.(Ref(src_local_from_global_map), src_search.facet_node_indices[src_facet], 0)
        push!(contributions, (src_local_indices, block))
    end
    abs(covered_area - dst_area) ≤ SEGMENT_COVERAGE_TOLERANCE * dst_area || return false
    for (src_local_indices, block) in contributions
        L[dst_local_indices, src_local_indices] += block
    end
    W[dst_local_indices, dst_local_indices] += square_block
    return true
end

function facet_static_nodes(coordinates::Matrix{Float64})
    n = size(coordinates, 2)
    return n == 4 ? SMatrix{3,4,Float64,12}(coordinates) : SMatrix{3,3,Float64,9}(coordinates)
end

function project_onto_facet_nodes(nodes::SMatrix{3,N,Float64}, x::SVector{3,Float64}) where {N}
    return closest_point_projection(nodes, x; strict=false)
end

# Center, two in-plane unit vectors, and unit normal of the plane fitted to
# a facet: the normal of the diagonals of a quadrilateral or of the edges of
# a triangle.
function facet_plane(coordinates::Matrix{Float64})
    n = size(coordinates, 2)
    center = SVector{3,Float64}(sum(coordinates; dims=2) / n)
    x = [SVector{3,Float64}(coordinates[:, i]) for i in 1:n]
    d1 = n == 4 ? x[3] - x[1] : x[2] - x[1]
    d2 = n == 4 ? x[4] - x[2] : x[3] - x[1]
    normal = cross(d1, d2)
    normal = normal / norm(normal)
    e1 = d1 - dot(d1, normal) * normal
    e1 = e1 / norm(e1)
    e2 = cross(normal, e1)
    return center, e1, e2, normal
end

# Nodes of a facet projected onto the plane, as a counterclockwise polygon.
function project_polygon(
    coordinates::Matrix{Float64}, center::SVector{3,Float64}, e1::SVector{3,Float64}, e2::SVector{3,Float64}
)
    polygon = [
        SVector{2,Float64}(
            dot(SVector{3,Float64}(coordinates[:, i]) - center, e1),
            dot(SVector{3,Float64}(coordinates[:, i]) - center, e2),
        ) for i in 1:size(coordinates, 2)
    ]
    polygon_area(polygon) < 0.0 && reverse!(polygon)
    return polygon
end

function polygon_area(polygon::Vector{SVector{2,Float64}})
    area = 0.0
    n = length(polygon)
    for i in 1:n
        a = polygon[i]
        b = polygon[i == n ? 1 : i + 1]
        area += a[1] * b[2] - a[2] * b[1]
    end
    return 0.5 * area
end

# Sutherland-Hodgman clipping of a convex polygon by a convex counterclockwise
# clipper. A point within a relative 1e-12 of the clipper's diameter of an
# edge line counts as inside, so that subject edges lying on clipper edges,
# which conforming and nested facets produce, are kept whole rather than
# intersected with a parallel line.
function clip_convex_polygon(subject::Vector{SVector{2,Float64}}, clipper::Vector{SVector{2,Float64}})
    output = subject
    n = length(clipper)
    diameter = maximum(norm(clipper[i] - clipper[j]) for i in 1:n for j in 1:n)
    for i in 1:n
        isempty(output) && break
        a = clipper[i]
        b = clipper[i == n ? 1 : i + 1]
        edge = b - a
        tolerance = 1.0e-12 * diameter * norm(edge)
        inside(p) = edge[1] * (p[2] - a[2]) - edge[2] * (p[1] - a[1]) ≥ -tolerance
        input = output
        output = SVector{2,Float64}[]
        m = length(input)
        for j in 1:m
            current = input[j]
            previous = input[j == 1 ? m : j - 1]
            current_inside = inside(current)
            previous_inside = inside(previous)
            if current_inside
                if !previous_inside
                    push!(output, segment_intersection(previous, current, a, b))
                end
                push!(output, current)
            elseif previous_inside
                push!(output, segment_intersection(previous, current, a, b))
            end
        end
    end
    return output
end

# Intersection of the segment pq with the line ab.
function segment_intersection(
    p::SVector{2,Float64}, q::SVector{2,Float64}, a::SVector{2,Float64}, b::SVector{2,Float64}
)
    d = q - p
    e = b - a
    denominator = d[1] * e[2] - d[2] * e[1]
    # Parallel within roundoff: the crossing point lies on p within the
    # tolerance of the inside test.
    abs(denominator) ≤ 1.0e-14 * norm(d) * norm(e) && return p
    t = ((a[1] - p[1]) * e[2] - (a[2] - p[2]) * e[1]) / denominator
    return p + clamp(t, 0.0, 1.0) * d
end

function get_overlap_rectangular_projection_matrix(
    dst_model::SolidMechanics,
    dst_bc::SolidMechanicsOverlapSchwarzBoundaryCondition,
    src_model::SolidMechanics,
    coupled_block_name::String,
    tol::Float64,
)
    src_mesh = src_model.mesh
    src_block_id = block_id_from_name(coupled_block_name, src_mesh)
    src_element_type_string = Exodus.read_block_parameters(src_mesh, src_block_id)[1]
    src_element_type = element_type_from_string(src_element_type_string)
    src_num_nodes = size(src_model.reference, 2)
    dst_local_from_global_map = dst_bc.local_from_global_map
    dst_num_nodes_sides = dst_bc.num_nodes_sides
    dst_side_set_node_indices = dst_bc.side_set_node_indices
    dst_num_nodes = length(dst_local_from_global_map)
    dst_coords = dst_model.reference
    dst_side_set_node_index = 1
    rectangular_projection_matrix = zeros(dst_num_nodes, src_num_nodes)
    for dst_num_nodes_side in dst_num_nodes_sides
        dst_side_nodes = dst_side_set_node_indices[dst_side_set_node_index:(dst_side_set_node_index + dst_num_nodes_side - 1)]
        dst_local_indices = get.(Ref(dst_local_from_global_map), dst_side_nodes, 0)
        dst_side_coordinates = dst_coords[:, dst_side_nodes]
        dst_element_type = get_element_type(2, Int64(dst_num_nodes_side))
        dst_num_int_points = default_num_int_pts(dst_element_type)
        dst_N, dst_dNdξ, dst_w, _ = isoparametric(dst_element_type, dst_num_int_points)
        for dst_point in 1:dst_num_int_points
            dst_Nₚ = Vector(dst_N[:, dst_point])
            dst_dNdξₚ = Matrix(dst_dNdξ[:, :, dst_point])
            dst_wₚ = dst_w[dst_point]
            dst_dXdξ = dst_dNdξₚ * dst_side_coordinates'
            dst_j = norm(cross(dst_dXdξ[1, :], dst_dXdξ[2, :]))
            dst_int_point_coord = dst_side_coordinates * dst_Nₚ
            src_node_indices, ξ, found = find_point_in_mesh(dst_int_point_coord, src_model, src_block_id, tol)
            if found == false
                norma_abortf(
                    "Weak overlap projector: could not find point (%.4e, %.4e, %.4e) in source mesh",
                    dst_int_point_coord[1],
                    dst_int_point_coord[2],
                    dst_int_point_coord[3],
                )
            end
            src_Nₚ = interpolate(src_element_type, ξ)[1]
            rectangular_projection_matrix[dst_local_indices, src_node_indices] += dst_Nₚ * src_Nₚ' * dst_j * dst_wₚ
        end
        dst_side_set_node_index += dst_num_nodes_side
    end
    return rectangular_projection_matrix
end

function compute_normal(mesh::ExodusDatabase, side_set_id::Int64, model::SolidMechanics)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    # The normal of the configuration in which equilibrium is stated: the
    # current configuration for finite kinematics, the reference configuration
    # for infinitesimal kinematics.
    if model.kinematics == Finite
        coords = model.reference .+ model.displacement
    else
        coords = model.reference
    end
    num_nodes = length(local_from_global_map)
    space_dim, _ = size(coords)
    normals = zeros(space_dim, num_nodes)
    local_indices = Array{Int64}(undef, 0)
    side_set_node_index = 1
    for num_nodes_side in num_nodes_sides
        side_nodes = side_set_node_indices[side_set_node_index:(side_set_node_index + num_nodes_side - 1)]
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        coordinates = coords[:, side_nodes]
        point_A = coordinates[:, 1]
        point_B = coordinates[:, 2]
        point_C = coordinates[:, end]
        BA = point_B - point_A
        CA = point_C - point_A
        N = cross(BA, CA)
        normals[:, local_indices] .= N / norm(N)
        side_set_node_index += num_nodes_side
    end
    return normals
end

function interpolate(tᵃ::Float64, tᵇ::Float64, xᵃ::Vector{Float64}, xᵇ::Vector{Float64}, t::Float64)
    Δt = tᵇ - tᵃ
    if Δt == 0.0
        return 0.5 * (xᵃ + xᵇ)
    end
    p = (tᵇ - t) / Δt
    q = 1.0 - p
    return p * xᵃ + q * xᵇ
end

function interpolate(param_hist::Vector{Float64}, value_hist::Vector{Vector{Float64}}, param::Float64)
    if param < param_hist[1] || isapprox(param, param_hist[1]; rtol=1.0e-06, atol=1.0e-12)
        return value_hist[1]
    end
    if param > param_hist[end] || isapprox(param, param_hist[end]; rtol=1.0e-06, atol=1.0e-12)
        return value_hist[end]
    end
    # Largest index with param_hist[index] ≤ param; the guards above ensure
    # param lies strictly inside the history, so 1 ≤ index < length and the
    # query interpolates on its own segment [index, index + 1].
    index = searchsortedlast(param_hist, param)
    return interpolate(param_hist[index], param_hist[index + 1], value_hist[index], value_hist[index + 1], param)
end

using Einsum

function closest_point_projection(nodes::Matrix{Float64}, x::Vector{Float64}; strict::Bool=true)
    num_nodes = size(nodes, 2)
    element_type = get_element_type(2, num_nodes)
    ξ = @MVector zeros(Float64, 2)
    y = @MVector zeros(Float64, 3)
    residual = @MVector zeros(Float64, 2)
    hessian = @MMatrix zeros(Float64, 2, 2)
    yx = @MVector zeros(Float64, 3)
    ddyddξ = MArray{Tuple{2,2,3},Float64,3}(undef)
    ddyddξyx = MMatrix{2,2,Float64}(undef)
    tol = 1.0e-10
    iteration = 1
    max_iterations = 64
    converged = false
    while true
        N, dN, ddN = interpolate(element_type, SVector(ξ))
        y = nodes * N
        dydξ = dN * nodes'
        yx = y - x
        residual = dydξ * yx
        @einsum ddyddξ[i, j, k] := ddN[i, j, l] * nodes[k, l]
        @einsum ddyddξyx[i, j] := ddyddξ[i, j, k] * yx[k]
        hessian = ddyddξyx + dydξ * dydξ'
        δ = -hessian \ residual
        ξ = ξ + δ
        err = norm(δ)
        if err <= tol
            converged = true
            break
        end
        iteration += 1
        if iteration > max_iterations
            break
        end
    end
    if converged == false
        if strict == true
            norma_abort("Closest point projection failed to converge")
        end
        # Non-strict callers scan many candidate facets; a diverging Newton on
        # a distant or badly oriented facet just disqualifies that candidate.
        return y, ξ, Inf, zeros(3)
    end
    _, dN, _ = interpolate(element_type, ξ)
    dxdξ = dN * nodes'
    perp_vec = cross(dxdξ[1, :], dxdξ[2, :])
    normal = perp_vec / norm(perp_vec)
    distance = -copysign(norm(yx), dot(yx, normal))
    return y, ξ, distance, normal
end
