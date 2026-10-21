# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Anisotropic elastic constants and lattice orientations.
#
# Elastic constants are given in the lattice frame in Voigt notation, with
# the order 11, 22, 33, 23, 13, 12 and engineering shear strains, so that C44
# relates σ23 to 2ε23.
#
# The lattice orientation Q holds the lattice axes, expressed in the sample
# frame, as its columns. A lattice vector v_L has sample components Q v_L, and
# a lattice-frame tensor has sample components
#
#     C_ijkl = Q_ip Q_jq Q_kr Q_ls C_L_pqrs.
#
# The elastic tensor and the slip systems of a crystal are both rotated with
# this same Q. Bunge Euler angles (φ1, Φ, φ2) give Q = Rz(φ1) Rx(Φ) Rz(φ2), the
# active rotation that carries the sample axes onto the lattice axes.

using LinearAlgebra
using Random
using StaticArrays

const VOIGT_PAIRS = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

# Voigt index of the index pair (i, j): (2, 3) -> 4, (1, 3) -> 5, (1, 2) -> 6.
@inline voigt_index(i::Integer, j::Integer) = i == j ? Int(i) : 9 - Int(i) - Int(j)

# Fourth-order tensor with minor and major symmetries from a Voigt matrix.
function voigt_to_tensor(C::SMatrix{6,6,Float64,36})::SArray{Tuple{3,3,3,3},Float64,4,81}
    CC = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        CC[i, j, k, l] = C[voigt_index(i, j), voigt_index(k, l)]
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(CC)
end

function tensor_to_voigt(CC::SArray{Tuple{3,3,3,3},Float64,4,81})::SMatrix{6,6,Float64,36}
    C = MMatrix{6,6,Float64,36}(undef)
    for (a, (i, j)) in enumerate(VOIGT_PAIRS), (b, (k, l)) in enumerate(VOIGT_PAIRS)
        C[a, b] = CC[i, j, k, l]
    end
    return SMatrix{6,6,Float64,36}(C)
end

# Sample components of a lattice-frame fourth-order tensor,
# C_ijkl = Q_ip Q_jq Q_kr Q_ls C_L_pqrs. Contracting one index at a time costs
# 4·3⁵ multiplies instead of 3⁸.
function rotate_tensor4(Q::SMatrix{3,3,Float64,9}, CC::SArray{Tuple{3,3,3,3},Float64,4,81})
    T1 = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, q in 1:3, r in 1:3, s in 1:3
        v = 0.0
        for p in 1:3
            v += Q[i, p] * CC[p, q, r, s]
        end
        T1[i, q, r, s] = v
    end
    T2 = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, r in 1:3, s in 1:3
        v = 0.0
        for q in 1:3
            v += Q[j, q] * T1[i, q, r, s]
        end
        T2[i, j, r, s] = v
    end
    T3 = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, k in 1:3, s in 1:3
        v = 0.0
        for r in 1:3
            v += Q[k, r] * T2[i, j, r, s]
        end
        T3[i, j, k, s] = v
    end
    T4 = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        v = 0.0
        for s in 1:3
            v += Q[l, s] * T3[i, j, k, s]
        end
        T4[i, j, k, l] = v
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(T4)
end

# Double contraction CC : A of a fourth-order tensor with a second-order one,
# generic in the element type of A so that ForwardDiff can pass through.
function double_contraction(CC::SArray{Tuple{3,3,3,3},Float64,4,81}, A::SMatrix{3,3,T,9}) where {T<:Number}
    return SMatrix{3,3,T,9}(sum(CC[i, j, k, l] * A[k, l] for k in 1:3, l in 1:3) for i in 1:3, j in 1:3)
end

# The nine constants of an orthotropic material and their Voigt positions.
const ORTHOTROPIC_CONSTANTS = (
    ("C11", (1, 1)),
    ("C22", (2, 2)),
    ("C33", (3, 3)),
    ("C12", (1, 2)),
    ("C13", (1, 3)),
    ("C23", (2, 3)),
    ("C44", (4, 4)),
    ("C55", (5, 5)),
    ("C66", (6, 6)),
)

# All 21 independent constants Cij, i ≤ j, of a triclinic material.
const TRICLINIC_CONSTANTS = ["C$i$j" for i in 1:6 for j in i:6]

# Lattice-frame Voigt stiffness of a material from its `symmetry:` and the
# independent constants of that symmetry class:
#   isotropic             two of the isotropic constants (as for linear elastic)
#   cubic                 C11, C12, C44
#   hexagonal             C11, C12, C13, C33, C44, with C66 = (C11 - C12) / 2;
#                         the c axis is lattice axis 3 ("transversely
#                         isotropic" is accepted as a synonym)
#   tetragonal            C11, C12, C13, C33, C44, C66
#   orthotropic           C11, C22, C33, C12, C13, C23, C44, C55, C66
#   triclinic             any of the 21 Cij, i ≤ j; the others are zero
function anisotropic_elastic_constants(params::Parameters)::SMatrix{6,6,Float64,36}
    haskey(params, "symmetry") || norma_abort(
        "An anisotropic elastic material requires a `symmetry:` (isotropic, cubic, hexagonal, " *
        "tetragonal, orthotropic, or triclinic).",
    )
    symmetry = lowercase(string(params["symmetry"]))
    symmetry == "transversely isotropic" && (symmetry = "hexagonal")
    function constant(name::String)
        haskey(params, name) || norma_abort("Elastic constant $name is required for $symmetry symmetry.")
        return Float64(params[name])
    end
    C = zeros(MMatrix{6,6,Float64,36})
    if symmetry == "isotropic"
        _, _, _, λ, μ = elastic_constants(params)
        C[1, 1] = C[2, 2] = C[3, 3] = λ + 2μ
        C[1, 2] = C[1, 3] = C[2, 3] = λ
        C[4, 4] = C[5, 5] = C[6, 6] = μ
    elseif symmetry == "cubic"
        C11, C12, C44 = constant("C11"), constant("C12"), constant("C44")
        C[1, 1] = C[2, 2] = C[3, 3] = C11
        C[1, 2] = C[1, 3] = C[2, 3] = C12
        C[4, 4] = C[5, 5] = C[6, 6] = C44
    elseif symmetry == "hexagonal" || symmetry == "tetragonal"
        C11, C12, C13, C33, C44 = constant("C11"), constant("C12"), constant("C13"), constant("C33"), constant("C44")
        C[1, 1] = C[2, 2] = C11
        C[3, 3] = C33
        C[1, 2] = C12
        C[1, 3] = C[2, 3] = C13
        C[4, 4] = C[5, 5] = C44
        C[6, 6] = symmetry == "hexagonal" ? (C11 - C12) / 2 : constant("C66")
    elseif symmetry == "orthotropic"
        for (name, (a, b)) in ORTHOTROPIC_CONSTANTS
            C[a, b] = constant(name)
        end
    elseif symmetry == "triclinic"
        for name in TRICLINIC_CONSTANTS
            haskey(params, name) && (C[parse(Int, name[2]), parse(Int, name[3])] = Float64(params[name]))
        end
    else
        norma_abort(
            "Unknown elastic symmetry '$symmetry'; expected isotropic, cubic, hexagonal, " *
            "tetragonal, orthotropic, or triclinic.",
        )
    end
    for a in 1:6, b in (a + 1):6
        C[b, a] = C[a, b]
    end
    eigenvalues = eigvals(Symmetric(Matrix(C)))
    if minimum(eigenvalues) <= 0.0
        norma_abortf(
            "The elastic constants for %s symmetry are not positive definite (smallest eigenvalue of the " *
            "Voigt stiffness %.3e).",
            symmetry,
            minimum(eigenvalues),
        )
    end
    return SMatrix{6,6,Float64,36}(C)
end

# Largest eigenvalue of the acoustic tensor A_ik = C_ijkl n_j n_l.
function acoustic_maximum(CC::SArray{Tuple{3,3,3,3},Float64,4,81}, n::SVector{3,Float64})
    A = SMatrix{3,3,Float64,9}(sum(CC[i, j, k, l] * n[j] * n[l] for j in 1:3, l in 1:3) for i in 1:3, k in 1:3)
    return maximum(eigvals(Symmetric(A)))
end

# Largest wave modulus ρ c² over all propagation directions and modes, the
# anisotropic counterpart of the P-wave modulus λ + 2μ that limits the stable
# step of explicit dynamics. It is a property of the tensor and does not
# depend on the lattice orientation. The acoustic tensor is sampled on a
# Fibonacci hemisphere (A(n) = A(-n)) and the best samples are refined by a
# pattern search on the sphere.
function maximum_wave_modulus(CC::SArray{Tuple{3,3,3,3},Float64,4,81})
    num_samples = 2048
    golden_angle = π * (3 - sqrt(5))
    samples = Vector{Tuple{Float64,SVector{3,Float64}}}(undef, num_samples)
    for k in 0:(num_samples - 1)
        z = 1 - (k + 0.5) / num_samples
        r = sqrt(1 - z^2)
        n = SVector(r * cos(golden_angle * k), r * sin(golden_angle * k), z)
        samples[k + 1] = (acoustic_maximum(CC, n), n)
    end
    sort!(samples; by=first, rev=true)
    directions = (
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(-1.0, 0.0, 0.0),
        SVector(0.0, -1.0, 0.0),
        SVector(0.0, 0.0, -1.0),
    )
    best = 0.0
    for (value, n) in samples[1:8]
        step = 0.05
        while step > 1.0e-8
            improved = false
            for d in directions
                trial = normalize(n + step * d)
                trial_value = acoustic_maximum(CC, trial)
                if trial_value > value
                    value, n, improved = trial_value, trial, true
                end
            end
            improved || (step /= 2)
        end
        best = max(best, value)
    end
    return best
end

# Q = Rz(φ1) Rx(Φ) Rz(φ2) from Bunge Euler angles in radians.
function orientation_from_bunge_angles(φ1::Real, Φ::Real, φ2::Real)::SMatrix{3,3,Float64,9}
    c1, s1 = cos(φ1), sin(φ1)
    c, s = cos(Φ), sin(Φ)
    c2, s2 = cos(φ2), sin(φ2)
    Rz1 = @SMatrix [c1 -s1 0.0; s1 c1 0.0; 0.0 0.0 1.0]
    Rx = @SMatrix [1.0 0.0 0.0; 0.0 c -s; 0.0 s c]
    Rz2 = @SMatrix [c2 -s2 0.0; s2 c2 0.0; 0.0 0.0 1.0]
    return Rz1 * Rx * Rz2
end

# Bunge Euler angles (φ1, Φ, φ2) in radians of a rotation Q, the inverse of
# orientation_from_bunge_angles. At Φ = 0 or π only φ1 ± φ2 is defined, and
# φ2 is taken as zero.
function bunge_angles(Q::AbstractMatrix{<:Real})
    Φ = acos(clamp(Q[3, 3], -1.0, 1.0))
    if sin(Φ) > 1.0e-10
        φ1 = atan(Q[1, 3], -Q[2, 3])
        φ2 = atan(Q[3, 1], Q[3, 2])
    else
        φ1 = atan(Q[2, 1], Q[1, 1])
        φ2 = 0.0
    end
    return mod2pi(φ1), Φ, mod2pi(φ2)
end

# A rotation drawn uniformly over SO(3), from the unit quaternion of three
# uniform numbers (Shoemake, Graphics Gems III, 1992): the lattice
# orientation of a grain of an untextured polycrystal.
function uniform_random_rotation(rng::AbstractRNG)::SMatrix{3,3,Float64,9}
    u1, u2, u3 = rand(rng), rand(rng), rand(rng)
    a, b = sqrt(1 - u1), sqrt(u1)
    w, x, y, z = a * sin(2π * u2), a * cos(2π * u2), b * sin(2π * u3), b * cos(2π * u3)
    return @SMatrix [
        1-2(y^2+z^2) 2(x*y-z*w) 2(x*z+y*w)
        2(x*y+z*w) 1-2(x^2+z^2) 2(y*z-x*w)
        2(x*z-y*w) 2(y*z+x*w) 1-2(x^2+y^2)
    ]
end

# Q from the three lattice axes expressed in the sample frame (LCM's `Basis
# Vector 1/2/3`). Each vector is normalized; the triad must be orthogonal and
# right-handed to within a tolerance, and the closest rotation is returned so
# that rounded input still gives an exact rotation.
function orientation_from_basis_vectors(vectors::AbstractVector)::SMatrix{3,3,Float64,9}
    length(vectors) == 3 && all(v -> length(v) == 3, vectors) ||
        norma_abort("`basis vectors` must be a list of three 3-component vectors.")
    B = hcat((normalize(SVector{3,Float64}(Float64.(v))) for v in vectors)...)
    B = SMatrix{3,3,Float64,9}(B)
    deviation = maximum(abs, B' * B - I)
    deviation <= 1.0e-6 || norma_abortf("`basis vectors` are not orthogonal (|bᵢ·bⱼ| up to %.2e).", deviation)
    det(B) > 0.0 || norma_abort("`basis vectors` must form a right-handed triad, with b3 = b1 × b2.")
    F = svd(B)
    return SMatrix{3,3,Float64,9}(F.U * F.Vt)
end

function angle_unit_factor(orientation::AbstractDict)
    unit = lowercase(string(get(orientation, "angle unit", "degrees")))
    unit in ("degree", "degrees") && return π / 180
    unit in ("radian", "radians") && return 1.0
    return norma_abort("`angle unit` must be degrees or radians (got '$unit').")
end

# Lattice orientation from an `orientation:` entry given in the input file:
# either `basis vectors` (the lattice axes in the sample frame) or `Bunge
# Euler angles` [φ1, Φ, φ2], in the `angle unit` (degrees by default).
function lattice_orientation(orientation::AbstractDict)::SMatrix{3,3,Float64,9}
    has_basis = haskey(orientation, "basis vectors")
    has_angles = haskey(orientation, "Bunge Euler angles")
    has_basis + has_angles == 1 ||
        norma_abort("An `orientation:` takes exactly one of `basis vectors` or `Bunge Euler angles`.")
    has_basis && return orientation_from_basis_vectors(orientation["basis vectors"])
    angles = orientation["Bunge Euler angles"]
    length(angles) == 3 || norma_abort("`Bunge Euler angles` must be [φ1, Φ, φ2].")
    factor = angle_unit_factor(orientation)
    return orientation_from_bunge_angles((factor * Float64(a) for a in angles)...)
end

# Lattice orientation of every element of a block. Without an `orientation:`
# the lattice is aligned with the sample axes; with `mesh variables` the
# orientation of each element is read from element variables of the input
# mesh: three names give Bunge Euler angles (in the `angle unit`), nine names
# give Q row-major, Q11, Q12, ..., Q33.
function block_lattice_orientations(
    material_params::AbstractDict, mesh::ExodusDatabase, block_id::Integer, num_elements::Integer
)::Vector{SMatrix{3,3,Float64,9}}
    orientation = get(material_params, "orientation", nothing)
    orientation === nothing && return fill(SMatrix{3,3,Float64,9}(I), num_elements)
    haskey(orientation, "mesh variables") || return fill(lattice_orientation(orientation), num_elements)
    # The mesh would silently win over an orientation also given on the
    # material, as basis vectors and angles together already abort.
    (haskey(orientation, "basis vectors") || haskey(orientation, "Bunge Euler angles")) && norma_abort(
        "An `orientation:` takes exactly one of `basis vectors`, `Bunge Euler angles` and `mesh variables`."
    )
    names = String.(orientation["mesh variables"])
    length(names) in (3, 9) ||
        norma_abort("`mesh variables` must name three element variables (Bunge Euler angles) or nine (Q row-major).")
    available = Exodus.read_names(mesh, ElementVariable)
    for name in names
        name in available || norma_abort("Element variable '$name' for the lattice orientation is not in the mesh.")
    end
    time_index = Int(get(orientation, "time index", 1))
    num_times = length(Exodus.read_times(mesh))
    1 <= time_index <= num_times || norma_abortf(
        "The `time index` %d of the lattice orientation is not a step of the mesh, which has %d.",
        time_index,
        num_times,
    )
    columns = map(names) do name
        try
            return Float64.(Exodus.read_values(mesh, ElementVariable, time_index, block_id, name))
        catch error
            error isa Exodus.ExodusError || rethrow()
            norma_abortf(
                "Element variable '%s' for the lattice orientation is not defined on block %d of the mesh.",
                name,
                block_id,
            )
        end
    end
    orientations = Vector{SMatrix{3,3,Float64,9}}(undef, num_elements)
    if length(names) == 3
        factor = angle_unit_factor(orientation)
        for e in 1:num_elements
            orientations[e] = orientation_from_bunge_angles(
                factor * columns[1][e], factor * columns[2][e], factor * columns[3][e]
            )
        end
    else
        for e in 1:num_elements
            Q = SMatrix{3,3,Float64,9}(columns[3 * (i - 1) + j][e] for i in 1:3, j in 1:3)
            orientations[e] = orientation_from_basis_vectors([Q[:, 1], Q[:, 2], Q[:, 3]])
        end
    end
    return orientations
end
