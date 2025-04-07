# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using StaticArrays

function barycentricD2N3(ξ::SVector{2,T}) where {T<:Real}
    N = @SVector [one(T) - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = @SMatrix [
        -one(T) one(T) zero(T)
        -one(T) zero(T) one(T)
    ]
    ddN = zeros(SArray{Tuple{2,2,3},T,3})
    return N, dN, ddN
end

function barycentricD2N3G1()
    w = @SVector [1 / 2]
    N = MMatrix{3,1,Float64}(undef)
    dN = MArray{Tuple{2,3,1},Float64}(undef)
    ξ = @SMatrix [1 / 3; 1 / 3]
    for p in 1:1
        Np, dNp, _ = barycentricD2N3(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:3
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:2
            @inbounds @simd for j in 1:3
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function barycentricD2N3G3()
    w = @SVector [1 / 6, 1 / 6, 1 / 6]
    N = MMatrix{3,3,Float64}(undef)
    dN = MArray{Tuple{2,3,3},Float64}(undef)
    ξ = @SMatrix [
        1/6 4/6 1/6
        1/6 1/6 4/6
    ]
    for p in 1:3
        Np, dNp, _ = barycentricD2N3(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:3
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:2
            @inbounds @simd for j in 1:3
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function barycentricD3N4(ξ::SVector{3,T}) where {T<:Real}
    N = @SVector [one(T) - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3]]
    dN = @SMatrix [
        -one(T) one(T) zero(T) zero(T)
        -one(T) zero(T) one(T) zero(T)
        -one(T) zero(T) zero(T) one(T)
    ]
    ddN = zeros(SArray{Tuple{3,3,4},T,3})
    return N, dN, ddN
end

function barycentricD3N10(ξ::SVector{3,T}) where {T<:Real}
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

function barycentricD3N4G1()
    w = @SVector [1 / 6]
    N = MMatrix{4,1,Float64}(undef)
    dN = MArray{Tuple{3,4,1},Float64}(undef)
    ξ = @SMatrix [0.25; 0.25; 0.25]
    for p in 1:1
        Np, dNp, _ = barycentricD3N4(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:4
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:3
            @inbounds @simd for j in 1:4
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function barycentricD3N4G4()
    w = @SVector [1 / 24, 1 / 24, 1 / 24, 1 / 24]
    N = MMatrix{4,4,Float64}(undef)
    dN = MArray{Tuple{3,4,4},Float64}(undef)
    s = sqrt(5)
    a = (5 + 3 * s) / 20
    b = (5 - s) / 20
    ξ = @SMatrix [
        b a b b
        b b a b
        b b b a
    ]
    for p in 1:4
        Np, dNp, _ = barycentricD3N4(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:4
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:3
            @inbounds @simd for j in 1:4
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function barycentricD3N10G4()
    w = @SVector [1 / 24, 1 / 24, 1 / 24, 1 / 24]
    N = MMatrix{10,4,Float64}(undef)
    dN = MArray{Tuple{3,10,4},Float64}(undef)
    s = sqrt(5)
    a = (5 + 3 * s) / 20
    b = (5 - s) / 20
    ξ = @SMatrix [
        b a b b
        b b a b
        b b b a
    ]
    for p in 1:4
        Np, dNp, _ = barycentricD3N10(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:10
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:3
            @inbounds @simd for j in 1:10
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function barycentricD3N10G5()
    a = -2 / 15
    b = 3 / 40
    w = @SVector [a, b, b, b, b]
    N = MMatrix{10,5,Float64}(0)
    dN = MArray{Tuple{3,10,5},Float64}(undef)
    ξ = @SMatrix [
        1/4 1/6 1/6 1/6 1/2
        1/4 1/6 1/6 1/2 1/6
        1/4 1/6 1/2 1/6 1/6
    ]
    for p in 1:5
        Np, dNp, _ = barycentricD3N10(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:10
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:3
            @inbounds @simd for j in 1:10
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function lagrangianD1N2(ξ::SVector{1,T}) where {T<:Real}
    N = @SVector [0.5 * (1.0 - ξ[1]), 0.5 * (1.0 + ξ[1])]
    dN = @SMatrix [-0.5; 0.5]
    ddN = zeros(SArray{Tuple{1,1,2},T,3})
    return N, dN, ddN
end

function lagrangianD1N2G1()
    w = @SVector [2.0]
    N = MMatrix{2,1,Float64}(undef)
    dN = MArray{Tuple{1,2,1},Float64}(undef)
    ξ = @SMatrix [0.0]
    for p in 1:1
        Np, dNp, _ = lagrangianD1N2(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:2
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:1
            @inbounds @simd for j in 1:2
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function lagrangianD1N2G2()
    w = @SVector [1.0, 1.0]
    N = MMatrix{2,2,Float64}(undef)
    dN = MArray{Tuple{1,2,2},Float64}(undef)
    g = 1 / sqrt(3)
    ξ = @SMatrix [-g; g]
    for p in 1:2
        Np, dNp, _ = lagrangianD1N2(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:2
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:1
            @inbounds @simd for j in 1:2
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function lagrangianD2N4(ξ::SVector{2,T}) where {T<:Real}
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
        ddN[1, 1, p] = zero(T)
        ddN[1, 2, p] = 0.25 * ra[p] * sa[p]
        ddN[2, 1, p] = 0.25 * ra[p] * sa[p]
        ddN[2, 2, p] = zero(T)
    end
    return SVector(N), SMatrix(dN), SArray(ddN)
end

function lagrangianD2N4G4()
    w = @SVector [1.0, 1.0, 1.0, 1.0]
    N = MMatrix{4,4,Float64}(undef)
    dN = MArray{Tuple{2,4,4},Float64}(undef)
    g = 1 / sqrt(3)
    ξ = @SMatrix [
        -g g g -g
        -g -g g g
    ]
    for p in 1:4
        Np, dNp, _ = lagrangianD2N4(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:4
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:2
            @inbounds @simd for j in 1:4
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function lagrangianD2N4G9()
    a = 25 / 81
    b = 40 / 81
    c = 64 / 81
    w = @SVector [a, a, a, a, b, b, b, b, c]
    N = MMatrix{4,9,Float64}(undef)
    dN = MArray{Tuple{2,4,9},Float64}(undef)
    g = sqrt(3 / 5)
    ξ = @SMatrix [
        -g g g -g 0 g 0 -g 0
        -g -g g g -g 0 g 0 0
    ]
    for p in 1:9
        Np, dNp, _ = lagrangianD2N4(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:4
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:2
            @inbounds @simd for j in 1:4
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function lagrangianD3N8(ξ::SVector{3,T}) where {T<:Real}
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

function lagrangianD3N8G8()
    w = @SVector [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    N = MMatrix{8,8,Float64}(undef)
    dN = MArray{Tuple{3,8,8},Float64}(undef)
    g = 1 / sqrt(3)
    ξ = @SMatrix [
        -g g g -g -g g g -g
        -g -g g g -g -g g g
        -g -g -g -g g g g g
    ]
    for p in 1:8
        Np, dNp, _ = lagrangianD3N8(SVector(ξ[:, p]))
        @inbounds @simd for i in 1:8
            N[i, p] = Np[i]
        end
        @inbounds for i in 1:3
            @inbounds @simd for j in 1:8
                dN[i, j, p] = dNp[i, j]
            end
        end
    end
    return SMatrix(N), SArray(dN), w, ξ
end

function default_num_int_pts(element_type::String)
    if element_type == "BAR2"
        return 1
    elseif element_type == "TRI3"
        return 3
    elseif element_type == "QUAD4"
        return 4
    elseif element_type == "TETRA4"
        return 4
    elseif element_type == "TETRA10"
        return 4
    elseif element_type == "HEX8"
        return 8
    else
        error("Invalid element type: ", element_type)
    end
end

function get_element_type(dim::Integer, num_nodes::Integer)
    if dim == 1 && num_nodes == 2
        return "BAR2"
    elseif dim == 2 && num_nodes == 3
        return "TRI3"
    elseif dim == 2 && num_nodes == 4
        return "QUAD4"
    elseif dim == 3 && num_nodes == 4
        return "TETRA4"
    elseif dim == 3 && num_nodes == 10
        return "TETRA10"
    elseif dim == 3 && num_nodes == 8
        return "HEX8"
    else
        error("Invalid dimension : ", dim, " and number of nodes : ", num_nodes)
    end
end

#
# Compute isoparametric interpolation functions, their parametric
# derivatives, integration weights, and integration point locations
#
function isoparametric(element_type::String, num_int::Integer)
    msg1 = "Invalid number of integration points: "
    msg2 = " for element type: "
    if element_type == "BAR2"
        if num_int == 1
            return lagrangianD1N2G1()
        elseif num_int == 2
            return lagrangianD1N2G2()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TRI3"
        if num_int == 1
            return barycentricD2N3G1()
        elseif num_int == 3
            return barycentricD2N3G3()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "QUAD4"
        if num_int == 4
            return lagrangianD2N4G4()
        elseif num_int == 9
            return lagrangianD2N4G9()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TETRA4"
        if num_int == 1
            return barycentricD3N4G1()
        elseif num_int == 4
            return barycentricD3N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TETRA10"
        if num_int == 4
            return barycentricD3N10G4()
        elseif num_int == 5
            return barycentricD3N10G5()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "HEX8"
        if num_int == 8
            return lagrangianD3N8G8()
        else
            error(msg1, num_int, msg2, element_type)
        end
    else
        error("Invalid element type: ", element_type)
    end
end

using Symbolics
@variables t, x, y, z

function get_side_set_nodal_forces(nodal_coord::Matrix{Float64}, traction_num::Num, time::Float64)
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
        values = Dict(t => time, x => point_coord[1], y => point_coord[2], z => point_coord[3])
        traction_sym = substitute(traction_num, values)
        traction_val = extract_value(traction_sym)
        nodal_force_component += traction_val * Nₚ * j * wₚ
    end
    return nodal_force_component
end

function map_to_parametric(element_type::String, nodes::Matrix{Float64}, point::Vector{Float64})
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

function interpolate(element_type::String, ξ::AbstractVector{Float64})
    if element_type == "BAR2"
        return lagrangianD1N2(SVector(ξ))
    elseif element_type == "TRI3"
        return barycentricD2N3(SVector(ξ))
    elseif element_type == "QUAD4"
        return lagrangianD2N4(SVector(ξ))
    elseif element_type == "TETRA4"
        return barycentricD3N4(SVector(ξ))
    elseif element_type == "TETRA10"
        return barycentricD3N10(SVector(ξ))
    elseif element_type == "HEX8"
        return lagrangianD3N8(SVector(ξ))
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside_parametric(element_type::String, ξ::AbstractVector{Float64}, tol::Float64=1.0e-06)
    factor = 1.0 + tol
    if element_type == "BAR2"
        return -factor ≤ ξ ≤ factor
    elseif element_type == "TRI3" || element_type == "TETRA4" || element_type == "TETRA10"
        return sum(ξ) ≤ factor
    elseif element_type == "QUAD4"
        return reduce(*, -factor * ones(2) .≤ ξ .≤ factor * ones(2))
    elseif element_type == "HEX8"
        return reduce(*, -factor * ones(3) .≤ ξ .≤ factor * ones(3))
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside(element_type::String, nodes::Matrix{Float64}, point::Vector{Float64}, tol::Float64=1.0e-06)
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
    num_nodes_per_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    ss_node_index = 1
    closest_face_nodes = Array{Float64}(undef, 0)
    closest_face_node_indices = Array{Int64}(undef, 0)
    minimum_nodal_distance = Inf
    for num_nodes_side in num_nodes_per_sides
        face_node_indices = side_set_node_indices[ss_node_index:(ss_node_index + num_nodes_side - 1)]
        face_nodes = model.current[:, face_node_indices]
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

function get_minimum_distance_to_nodes(nodes::Matrix{Float64}, point::Vector{Float64})
    distances = norm.(eachcol(nodes) .- Ref(point))
    return minimum(distances)
end

function get_side_set_local_from_global_map(mesh::ExodusDatabase, side_set_id::Integer)
    num_nodes_per_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    unique_node_indices = unique(side_set_node_indices)
    num_nodes = length(unique_node_indices)
    local_from_global_map = Dict{Int64,Int64}()
    for i in 1:num_nodes
        local_from_global_map[Int64(unique_node_indices[i])] = i
    end
    return local_from_global_map, num_nodes_per_sides, Int64.(side_set_node_indices)
end

function get_side_set_global_from_local_map(mesh::ExodusDatabase, side_set_id::Integer)
    side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)[2]
    unique_node_indices = unique(side_set_node_indices)
    num_nodes = length(unique_node_indices)
    global_from_local_map = zeros(Int64, num_nodes)
    for i in 1:num_nodes
        global_from_local_map[i] = Int64(unique_node_indices[i])
    end
    return global_from_local_map
end

function get_square_projection_matrix(model::SolidMechanics, side_set_id::Integer)
    mesh = model.mesh
    local_from_global_map, num_nodes_sides, side_set_node_indices = get_side_set_local_from_global_map(
        mesh, side_set_id
    )
    num_nodes = length(local_from_global_map)
    if model.kinematics == Finite
        coords = model.reference
    else
        coords = model.current
    end
    square_projection_matrix = zeros(num_nodes, num_nodes)
    side_set_node_index = 1
    for num_nodes_side in num_nodes_sides
        side_nodes = side_set_node_indices[side_set_node_index:(side_set_node_index + num_nodes_side - 1)]
        side_coordinates = coords[:, side_nodes]
        element_type = get_element_type(2, Int64(num_nodes_side))
        num_int_points = default_num_int_pts(element_type)
        N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
        side_matrix = zeros(num_nodes_side, num_nodes_side)
        for point in 1:num_int_points
            Nₚ = N[:, point]
            dNdξₚ = dNdξ[:, :, point]
            dXdξ = dNdξₚ * side_coordinates'
            j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
            wₚ = w[point]
            side_matrix += Nₚ * Nₚ' * j * wₚ
        end
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        square_projection_matrix[local_indices, local_indices] += side_matrix
        side_set_node_index += num_nodes_side
    end
    return square_projection_matrix
end

function get_rectangular_projection_matrix(
    src_model::SolidMechanics, src_side_set_id::Integer, dst_model::SolidMechanics, dst_side_set_id::Integer
)
    src_mesh = src_model.mesh
    src_local_from_global_map, _, _ = get_side_set_local_from_global_map(src_mesh, src_side_set_id)
    src_num_nodes = length(src_local_from_global_map)
    dst_mesh = dst_model.mesh
    dst_local_from_global_map, dst_num_nodes_sides, dst_side_set_node_indices = get_side_set_local_from_global_map(
        dst_mesh, dst_side_set_id
    )
    dst_num_nodes = length(dst_local_from_global_map)
    if dst_model.kinematics == Finite
        dst_coords = dst_model.reference
    else
        dst_coords = dst_model.current
    end
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
            dst_Nₚ = dst_N[:, dst_point]
            dst_dNdξₚ = dst_dNdξ[:, :, dst_point]
            dst_dXdξ = dst_dNdξₚ * dst_side_coordinates'
            dst_j = norm(cross(dst_dXdξ[1, :], dst_dXdξ[2, :]))
            dst_wₚ = dst_w[dst_point]
            dst_int_point_coord = dst_side_coordinates * dst_Nₚ
            _, ξ, src_side_coordinates, src_side_nodes, _, _ = project_point_to_side_set(
                dst_int_point_coord, src_model, src_side_set_id
            )
            src_side_element_type = get_element_type(2, size(src_side_coordinates)[2])
            src_Nₚ, _, _ = interpolate(src_side_element_type, ξ)
            src_local_indices = get.(Ref(src_local_from_global_map), src_side_nodes, 0)
            rectangular_projection_matrix[dst_local_indices, src_local_indices] += dst_Nₚ * src_Nₚ' * dst_j * dst_wₚ
        end
        dst_side_set_node_index += dst_num_nodes_side
    end
    return rectangular_projection_matrix
end

function compute_normal(mesh::ExodusDatabase, side_set_id::Int64, model::SolidMechanics)
    local_from_global_map, num_nodes_sides, side_set_node_indices = get_side_set_local_from_global_map(
        mesh, side_set_id
    )
    if model.kinematics == Finite
        coords = model.reference
    else
        coords = model.current
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
    if param < param_hist[1]
        param = param_hist[1]
    end
    if param > param_hist[end]
        param = param_hist[end]
    end
    index = 1
    size = length(param_hist)
    while param_hist[index] < param
        if index == size
            break
        end
        index += 1
    end
    if index == 1
        return value_hist[1]
    elseif index == size
        return value_hist[size]
    else
        return interpolate(param_hist[index], param_hist[index + 1], value_hist[index], value_hist[index + 1], param)
    end
end

using Einsum

function closest_point_projection(nodes::Matrix{Float64}, x::Vector{Float64})
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
            break
        end
        iteration += 1
        if iteration > max_iterations
            error("Closest point projection failed to converge")
        end
    end
    _, dN, _ = interpolate(element_type, ξ)
    dxdξ = dN * nodes'
    perp_vec = cross(dxdξ[1, :], dxdξ[2, :])
    normal = perp_vec / norm(perp_vec)
    distance = -copysign(norm(yx), dot(yx, normal))
    return y, ξ, distance, normal
end
