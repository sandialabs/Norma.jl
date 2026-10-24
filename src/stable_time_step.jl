# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Critical time step of central difference (issue #243).
#
# The step is CFL times an estimate of the critical step 2/ω_max, where ω_max²
# is the largest eigenvalue of K φ = ω² M φ with the lumped mass M. Two
# estimates are available, selected by `stable time step` in the integrator:
#
#   eigenvalue  the smallest element critical step 2/ω_max,e over the mesh,
#               with ω_max,e² the largest eigenvalue of the element stiffness
#               (material tangent at the current configuration) against the
#               element lumped mass. The largest element eigenvalue bounds the
#               global one from above, so this step is a lower bound on the
#               critical step: CFL ≤ 1 is stable for a linear problem.
#   length      (default) the critical step of an undistorted reference element
#               of the block's type and material, scaled by the ratio of the
#               smallest element altitude of the block to the altitude of the
#               reference element. Exact on undistorted elements; on distorted
#               ones it follows the smallest altitude, which detects flat and
#               sliver elements, but it is not a bound.
#
# Both are recomputed every `stable time step interval` evaluations (default 1
# for length, 100 for eigenvalue), on the current configuration, so that a
# mesh that deforms changes the step.
#
# The quadratic tetrahedron (TETRA10) has negative row-sum lumped masses at its
# corners, for which ω_max is undefined. Its estimate uses the diagonal of the
# consistent mass scaled to the element mass, which is positive; the integrator
# still uses the row sums, so for this element the estimate is not a bound.

# Element volume by the block's quadrature rule.
function element_volume(X::SMatrix{3,NN,T}, dN::SArray{Tuple{3,NN,NP},T}, weights::AbstractVector) where {NN,NP,T}
    volume = zero(T)
    @inbounds for point in 1:NP
        volume += det(dN[:, :, point] * X') * weights[point]
    end
    return volume
end

triangle_area(p, q, r) = 0.5 * norm(cross(q - p, r - p))
# Area of a quadrilateral face from its diagonals; exact for a planar face.
quadrilateral_area(p, q, r, s) = 0.5 * norm(cross(r - p, s - q))

# Smallest altitude of an element. A tetrahedron (TETRA4, TETRA10, and the
# composite TETRA10, from the four vertices): the smallest of 3V/A over its
# faces. A hexahedron: the smallest of V divided by the mean area of each pair
# of opposite faces, which is the distance between the faces of a
# parallelepiped.
function element_altitude(element_type::ElementType, X::SMatrix{3,NN,T}, volume::T) where {NN,T}
    if element_type == HEX8
        x = ntuple(a -> X[:, a], Val(8))
        pairs = (
            (quadrilateral_area(x[1], x[2], x[3], x[4]), quadrilateral_area(x[5], x[6], x[7], x[8])),
            (quadrilateral_area(x[1], x[2], x[6], x[5]), quadrilateral_area(x[4], x[3], x[7], x[8])),
            (quadrilateral_area(x[1], x[4], x[8], x[5]), quadrilateral_area(x[2], x[3], x[7], x[6])),
        )
        return minimum(volume / (0.5 * (A1 + A2)) for (A1, A2) in pairs)
    end
    a, b, c, d = X[:, 1], X[:, 2], X[:, 3], X[:, 4]
    tetra_volume = abs(dot(b - a, cross(c - a, d - a))) / 6
    faces = (triangle_area(b, c, d), triangle_area(a, c, d), triangle_area(a, b, d), triangle_area(a, b, c))
    return minimum(3 * tetra_volume / A for A in faces)
end

# Largest eigenvalue ω² of the element stiffness against the element mass used
# for the estimate, at the reference coordinates X and current coordinates x.
# The material state and the lattice orientation are those of `state_element`
# of the block; the state is not updated.
function element_omega_squared(
    model::SolidMechanics,
    block_index::Int,
    state_element::Int,
    material::Material,
    density::Float64,
    element_type::ElementType,
    X::SMatrix{3,NN,T},
    x::SMatrix{3,NN,T},
    N::SMatrix{NN,NP,T},
    dN::SArray{Tuple{3,NN,NP},T},
    weights::AbstractVector,
    mass_N::SMatrix,
    mass_dN::SArray,
    mass_weights::AbstractVector,
) where {NN,NP,T}
    K = create_element_matrix(T, Val(NN))
    fill!(K, zero(T))
    m = create_element_vector(T, Val(NN))
    fill!(m, zero(T))
    consistent = element_type == TETRA10
    Mc = consistent ? create_element_matrix(T, Val(NN)) : nothing
    consistent && fill!(Mc, zero(T))
    Q = element_lattice_orientation(model, material, block_index, state_element)
    for point in 1:NP
        dXdξ = dN[:, :, point] * X'
        dNdX = dXdξ \ dN[:, :, point]
        F = x * dNdX'
        dvol = det(dXdξ) * weights[point]
        state = model.state_old[block_index][state_element][point]
        _, _, AA, _ = constitutive_update(material, F, Q, state, 0.0; need_tangent=true)
        add_stiffness!(K, dNdX, AA, dvol)
        if consistent
            add_mass!(Mc, N[:, point], density, dvol)
        elseif element_type != COMPOSITE_TETRA10
            add_lumped_mass!(m, N[:, point], density, dvol)
        end
    end
    if element_type == COMPOSITE_TETRA10
        for point in eachindex(mass_weights)
            dvol = det(mass_dN[:, :, point] * X') * mass_weights[point]
            add_lumped_mass!(m, mass_N[:, point], density, dvol)
        end
    elseif consistent
        diagonal = [Mc[i, i] for i in 1:(3 * NN)]
        total = sum(Mc) / 3
        scale = 3 * total / sum(diagonal)
        for i in 1:(3 * NN)
            m[i] = scale * diagonal[i]
        end
    end
    s = 1 ./ sqrt.(Vector(m))
    A = Symmetric(s .* Matrix(K) .* s')
    return eigmax(A)
end

# Undistorted reference element of each type, with unit altitude scale: the
# unit cube and the regular tetrahedron of unit edge (TETRA10 and the
# composite TETRA10 with midside nodes at the edge midpoints, in Exodus order).
function reference_element_coordinates(element_type::ElementType)
    if element_type == HEX8
        return Float64[
            0 1 1 0 0 1 1 0
            0 0 1 1 0 0 1 1
            0 0 0 0 1 1 1 1
        ]
    end
    v = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, sqrt(3) / 2, 0.0], [0.5, sqrt(3) / 6, sqrt(2 / 3)])
    element_type == TETRA4 && return reduce(hcat, v)
    if element_type in (TETRA10, COMPOSITE_TETRA10)
        mid(i, j) = (v[i] + v[j]) / 2
        return reduce(hcat, (v..., mid(1, 2), mid(2, 3), mid(3, 1), mid(1, 4), mid(2, 4), mid(3, 4)))
    end
    return norma_abort(
        "The stable time step of central difference is not implemented for element type $(element_type)."
    )
end

# Critical step of the reference element of a block divided by its altitude.
function reference_step_per_altitude(
    model::SolidMechanics, block_index::Int, block::ElementBlockData, N::SMatrix{NN,NP,T}, dN::SArray{Tuple{3,NN,NP},T}
) where {NN,NP,T}
    material = model.materials[block_index]
    X = SMatrix{3,NN,T}(reference_element_coordinates(block.element_type))
    volume = element_volume(X, dN, block.weights)
    ω² = element_omega_squared(
        model,
        block_index,
        1,
        material,
        material.ρ,
        block.element_type,
        X,
        X,
        N,
        dN,
        block.weights,
        block.mass_N,
        block.mass_dN,
        block.mass_weights,
    )
    return (2 / sqrt(ω²)) / element_altitude(block.element_type, X, volume)
end

@inline function element_coordinates(model::SolidMechanics, block::ElementBlockData, element::Int, ::Val{NN}) where {NN}
    nodes = view(block.connectivity, :, element)
    X = SMatrix{3,NN,Float64}(view(model.reference, :, nodes))
    x = X + SMatrix{3,NN,Float64}(view(model.displacement, :, nodes))
    return X, x
end

# Smallest element altitude of a block on the current configuration. The
# shape function tables of a block are stored untyped, so they are passed as
# arguments: read inside the loop they would make every element a dynamic
# call (205 MB allocated and 27 ms per call on 110 592 hexahedra before).
function minimum_block_altitude(
    model::SolidMechanics,
    block::ElementBlockData,
    ::SMatrix{NN,NP,T},
    dN::SArray{Tuple{3,NN,NP},T},
    weights::AbstractVector{T},
) where {NN,NP,T}
    element_type = block.element_type
    minima = fill(Inf, maxthreadid())
    @threads for element in 1:block.num_elements
        _, x = element_coordinates(model, block, element, Val(NN))
        h = element_altitude(element_type, x, element_volume(x, dN, weights))
        t = threadid()
        @inbounds minima[t] = min(minima[t], h)
    end
    return minimum(minima)
end

# Smallest element critical step of a block on the current configuration.
function minimum_block_eigenvalue_step(
    model::SolidMechanics,
    block_index::Int,
    block::ElementBlockData,
    N::SMatrix{NN,NP,T},
    dN::SArray{Tuple{3,NN,NP},T},
    weights::AbstractVector{T},
    mass_N::SMatrix,
    mass_dN::SArray,
    mass_weights::AbstractVector,
) where {NN,NP,T}
    material = model.materials[block_index]
    element_type = block.element_type
    density = material.ρ
    minima = fill(Inf, maxthreadid())
    @threads for element in 1:block.num_elements
        X, x = element_coordinates(model, block, element, Val(NN))
        ω² = element_omega_squared(
            model,
            block_index,
            element,
            material,
            density,
            element_type,
            X,
            x,
            N,
            dN,
            weights,
            mass_N,
            mass_dN,
            mass_weights,
        )
        t = threadid()
        @inbounds minima[t] = min(minima[t], 2 / sqrt(ω²))
    end
    return minimum(minima)
end

# Estimate of the critical step of central difference on the current
# configuration, before the CFL factor.
function critical_time_step(integrator::CentralDifference, model::SolidMechanics)
    critical = Inf
    for (block_index, block) in enumerate(model.blocks)
        if integrator.stable_method === :eigenvalue
            block_step = minimum_block_eigenvalue_step(
                model,
                block_index,
                block,
                block.N,
                block.dN,
                block.weights,
                block.mass_N,
                block.mass_dN,
                block.mass_weights,
            )
        else
            if length(integrator.stable_reference) < block_index || !isfinite(integrator.stable_reference[block_index])
                resize_reference!(integrator, length(model.blocks))
                integrator.stable_reference[block_index] = reference_step_per_altitude(
                    model, block_index, block, block.N, block.dN
                )
            end
            block_step =
                integrator.stable_reference[block_index] *
                minimum_block_altitude(model, block, block.N, block.dN, block.weights)
        end
        critical = min(critical, block_step)
    end
    return critical
end

function resize_reference!(integrator::CentralDifference, num_blocks::Int)
    n = length(integrator.stable_reference)
    n < num_blocks && append!(integrator.stable_reference, fill(NaN, num_blocks - n))
    return nothing
end

# The largest step the central difference scheme takes on the current
# configuration: CFL times the estimate of the critical step, recomputed every
# `stable time step interval` evaluations. Implicit and reduced integrators
# have no such limit.
stable_time_step(::TimeIntegrator, ::Model) = Inf

function stable_time_step(integrator::CentralDifference, model::SolidMechanics)
    if integrator.stable_count % integrator.stable_interval == 0 || !isfinite(integrator.stable_critical)
        integrator.stable_critical = critical_time_step(integrator, model)
    end
    integrator.stable_count += 1
    return integrator.CFL * integrator.stable_critical
end
