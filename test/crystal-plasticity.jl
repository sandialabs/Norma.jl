# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "lcm-reference", "lcm_crystal_plasticity.jl"))

# The LCM cases a material point can replay: every one whose gold file
# records each step. SingleSlipHard_Explicit and
# MinisolverStep_ConjugateGradient write only every tenth and fifth step, so
# their steps are not known; the first is checked with the FEM instead.
const LCM_REPLAYABLE_CASES = [
    "CP-uniaxial",
    "ElasticShear2d",
    "PlasticShear2d",
    "ElasticShear3d",
    "ElasticShearRotated3d",
    "SingleSlip_Explicit",
    "SingleSlip_Implicit",
    "SingleSlipHard_Implicit",
    "SingleSlipSaturation",
    "MultiSlipPlaneHard_Implicit",
    "MultiSlipPlaneHard_Implicit_Active_Sets",
    "MinisolverStep_Newton",
    "MinisolverStep_NewtonLineSearch",
    "MinisolverStep_TrustRegion",
    "MinisolverStep_LineSearchRegularized",
    "ThermallyActivatedSlip",
    "MultiFamily",
    "DislocationDensityHardening",
    "DislocationDensityTwoFamilies",
    "QuadSlipDislocationDensityTraction",
]

# An FCC crystal with LCM's copper-like constants and its 12 {111}<110>
# systems. With a `drag` coefficient the flow rule is the power law with
# drag.
function fcc_crystal(;
    scheme="implicit",
    residual="slip",
    rate_exponent=20.0,
    reference_rate=1.0,
    drag=nothing,
    hardening=Dict(
        "type" => "linear minus recovery",
        "hardening modulus" => 355.0,
        "recovery modulus" => 2.9,
        "initial hardening state" => 122.0,
    ),
    systems="fcc",
)
    flow = Dict{String,Any}(
        "type" => "power law", "reference slip rate" => reference_rate, "rate exponent" => rate_exponent
    )
    if drag !== nothing
        flow["type"] = "power law with drag"
        flow["drag coefficient"] = drag
    end
    return Norma.Parameters(
        "model" => "crystal plasticity",
        "symmetry" => "cubic",
        "C11" => 204600.0,
        "C12" => 137700.0,
        "C44" => 126200.0,
        "integration scheme" => scheme,
        "residual" => residual,
        "slip families" => [Dict("slip systems" => systems, "flow rule" => flow, "hardening law" => hardening)],
    )
end

const HCP_FAMILIES = (
    "hcp basal", "hcp prismatic", "hcp pyramidal a", "hcp pyramidal c+a first order", "hcp pyramidal c+a second order"
)

# An HCP crystal with elastic constants and c/a close to magnesium's, and
# all five HCP families, in the order of HCP_FAMILIES, each with a constant
# slip resistance and a steep power law.
function hcp_crystal(; ca=1.624, resistances=(10.0, 40.0, 60.0, 80.0, 80.0), rate_exponent=100.0)
    return Norma.Parameters(
        "model" => "crystal plasticity",
        "symmetry" => "hexagonal",
        "C11" => 59400.0,
        "C12" => 25610.0,
        "C13" => 21440.0,
        "C33" => 61600.0,
        "C44" => 16400.0,
        "slip families" => [
            Dict(
                "slip systems" => name,
                "c/a ratio" => ca,
                "flow rule" =>
                    Dict("type" => "power law", "reference slip rate" => 1.0, "rate exponent" => rate_exponent),
                "hardening law" => Dict("type" => "none", "initial hardening state" => g),
            ) for (name, g) in zip(HCP_FAMILIES, resistances)
        ],
    )
end

function rotation_about(axis, θ)
    n = normalize(SVector{3,Float64}(axis))
    K = @SMatrix [0.0 -n[3] n[2]; n[3] 0.0 -n[1]; -n[2] n[1] 0.0]
    return SMatrix{3,3,Float64,9}(I + sin(θ) * K + (1 - cos(θ)) * K * K)
end

# The rotation that turns lattice direction d onto sample x.
function orientation_along_x(d)
    n = normalize(SVector{3,Float64}(d))
    axis = cross(n, SVector(1.0, 0.0, 0.0))
    norm(axis) < 1.0e-12 && return SMatrix{3,3,Float64,9}(I)
    return rotation_about(axis, acos(clamp(n[1], -1.0, 1.0)))
end

# Uniaxial stress along sample x for the material point driver.
const UNIAXIAL_FREE = [false true true; false true true; false false true]

@testset "Crystal Plasticity" begin
    identity3 = SMatrix{3,3,Float64,9}(I)

    @testset "LCM Gold At The Material Point: $case" for case in LCM_REPLAYABLE_CASES
        errors = lcm_crystal_errors(lcm_crystal_replay(case))
        @test maximum(values(errors)) < 1.0e-7
    end

    @testset "Slip System Tables" begin
        counts = Dict(
            "fcc" => 12,
            "bcc 110" => 12,
            "bcc 112" => 12,
            "bcc 123" => 24,
            "hcp basal" => 3,
            "hcp prismatic" => 3,
            "hcp pyramidal a" => 6,
            "hcp pyramidal c+a first order" => 12,
            "hcp pyramidal c+a second order" => 6,
        )
        for (name, count) in counts
            material = Norma.create_material(fcc_crystal(; systems=name))
            @test Norma.number_slip_systems(material) == count
            @test all(abs(dot(s, n)) < 1.0e-12 for (s, n) in zip(material.s, material.n))
            @test all(norm(s) ≈ 1.0 && norm(n) ≈ 1.0 for (s, n) in zip(material.s, material.n))
            schmid = [s * n' for (s, n) in zip(material.s, material.n)]
            @test all(!isapprox(schmid[i], schmid[j]) for i in 1:count for j in (i + 1):count)
        end
    end

    # The slip systems of a crystallographic family: every pair of a plane
    # normal and a direction in that plane, from all the index variants of
    # the family, as unit vectors and once up to the sense of slip.
    same_system((s1, n1), (s2, n2)) = abs(dot(s1, s2)) ≈ 1.0 && abs(dot(n1, n2)) ≈ 1.0
    function family_systems(normals, directions)
        systems = NTuple{2,SVector{3,Float64}}[]
        for n in normals, s in directions
            abs(dot(s, n)) < 1.0e-12 * norm(s) * norm(n) || continue
            system = (normalize(s), normalize(n))
            any(same_system(system, other) for other in systems) || push!(systems, system)
        end
        return systems
    end
    # The systems of a table and of a family correspond one to one.
    function same_systems(table, family)
        length(table) == length(family) || return false
        return all(count(f -> same_system(t, f), family) == 1 for t in table) &&
               all(count(t -> same_system(t, f), table) == 1 for f in family)
    end
    table_systems(material, systems=1:Norma.number_slip_systems(material)) =
        [(material.s[α], material.n[α]) for α in systems]
    permutations = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))

    @testset "Cubic Slip System Crystallography" begin
        # Each cubic table is its whole family: the permutations and sign
        # changes of the plane and direction indices.
        cubic_variants(v) = unique([
            SVector{3,Float64}(σ[1] * v[p[1]], σ[2] * v[p[2]], σ[3] * v[p[3]]) for p in permutations for
            σ in Iterators.product((1, -1), (1, -1), (1, -1))
        ])
        families = (
            "fcc" => ((1, 1, 1), (1, 1, 0)),
            "bcc 110" => ((1, 1, 0), (1, 1, 1)),
            "bcc 112" => ((1, 1, 2), (1, 1, 1)),
            "bcc 123" => ((1, 2, 3), (1, 1, 1)),
        )
        for (name, (plane, direction)) in families
            material = Norma.create_material(fcc_crystal(; systems=name))
            family = family_systems(cubic_variants(plane), cubic_variants(direction))
            @test same_systems(table_systems(material), family)
        end
    end

    @testset "Cubic Schmid Factors" begin
        # The largest Schmid factor cos φ cos λ in uniaxial loading along
        # [001], [011] and [111], and the number of systems that reach it.
        # {111}<110>, and {110}<111> with the roles of planes and directions
        # exchanged: 1/√6 ≈ 0.408 on 8 and 4 systems, and 2/(3√6) ≈ 0.272 on 6.
        # {112}<111>: along [001] the planes (±1 ±1 2), at cos φ = 2/√6, with
        # their <111> at cos λ = 1/√3, give √2/3 ≈ 0.471; along [111],
        # (1 1 2)[1 1 -1] gives (4/√18)(1/3) = 2√2/9 ≈ 0.314.
        # {123}<111>: along [001] the planes with the 3 along the load, at
        # cos φ = 3/√14, give 3/√42 ≈ 0.463; along [111], (1 2 3)[1 1 -1]
        # gives (6/√42)(1/3) = 2/√42 ≈ 0.309.
        loads = ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0])
        expected = (
            "fcc" => ((1 / sqrt(6.0), 8), (1 / sqrt(6.0), 4), (2 / (3 * sqrt(6.0)), 6)),
            "bcc 110" => ((1 / sqrt(6.0), 8), (1 / sqrt(6.0), 4), (2 / (3 * sqrt(6.0)), 6)),
            "bcc 112" => ((sqrt(2.0) / 3, 4), (sqrt(2.0) / 3, 2), (2 * sqrt(2.0) / 9, 3)),
            "bcc 123" => ((3 / sqrt(42.0), 8), (3 / sqrt(42.0), 4), (2 / sqrt(42.0), 6)),
        )
        for (name, values) in expected
            material = Norma.create_material(fcc_crystal(; systems=name))
            for (load, (m_max, count_max)) in zip(loads, values)
                l = normalize(SVector{3,Float64}(load))
                m = [abs(dot(l, s) * dot(l, n)) for (s, n) in table_systems(material)]
                @test maximum(m) ≈ m_max rtol = 1.0e-12
                @test count(x -> isapprox(x, m_max; rtol=1.0e-12), m) == count_max
            end
        end
    end

    @testset "Hexagonal Slip System Crystallography" begin
        # Miller-Bravais indices in the lattice frame of the tables, with
        # a = 1: the direction [u v t w] is u a1 + v a2 + t a3 + w c, and the
        # plane (h k i l) has the normal h b1 + k b2 + l b3, with b1, b2, b3
        # the reciprocal basis of a1, a2, c. A family takes the permutations
        # of the first three indices, both signs of those three, and both
        # signs of the fourth.
        a = (SVector(1.0, 0.0, 0.0), SVector(-0.5, sqrt(3.0) / 2, 0.0), SVector(-0.5, -sqrt(3.0) / 2, 0.0))
        direction(v, ca) = v[1] * a[1] + v[2] * a[2] + v[3] * a[3] + SVector(0.0, 0.0, v[4] * ca)
        normal(v, ca) = SVector(v[1], (v[1] + 2 * v[2]) / sqrt(3.0), v[4] / ca)
        bravais_variants(v) = unique([
            (σ * v[p[1]], σ * v[p[2]], σ * v[p[3]], σ4 * v[4]) for p in permutations for σ in (1, -1) for σ4 in (1, -1)
        ])
        families = (
            ((0, 0, 0, 1), (1, 1, -2, 0)), # basal {0001}<11-20>
            ((1, 0, -1, 0), (1, 1, -2, 0)), # prismatic {10-10}<11-20>
            ((1, 0, -1, 1), (1, 1, -2, 0)), # pyramidal {10-11}<11-20>
            ((1, 0, -1, 1), (1, 1, -2, 3)), # first-order pyramidal {10-11}<11-23>
            ((1, 1, -2, 2), (1, 1, -2, 3)), # second-order pyramidal {11-22}<11-23>
        )
        # Unit vectors along the a axes, <11-20>, and along <10-10>, at 30°
        # from them, in the basal plane.
        a_axes = [SVector(cos(k * π / 3), sin(k * π / 3)) for k in 0:5]
        m_axes = [SVector(cos((k + 0.5) * π / 3), sin((k + 0.5) * π / 3)) for k in 0:5]
        along(v, axes) = any(isapprox(normalize(SVector(v[1], v[2])), u; atol=1.0e-12) for u in axes)
        # Ti, ideal and Zn.
        for ca in (1.587, sqrt(8 / 3), 1.856)
            material = Norma.create_material(hcp_crystal(; ca=ca))
            for (family, (plane, slip)) in zip(material.families, families)
                expected = family_systems(
                    [normal(v, ca) for v in bravais_variants(plane)], [direction(v, ca) for v in bravais_variants(slip)]
                )
                @test same_systems(table_systems(material, family.systems), expected)
            end
            tables = [Norma.hexagonal_slip_systems(name, ca) for name in HCP_FAMILIES]
            # <a> slip: the directions are in the basal plane, along the a axes.
            @test all(s[3] == 0.0 && along(s, a_axes) for table in tables[1:3] for (s, _) in table)
            # <c+a> slip: the direction is a + c, with a along an a axis.
            # Scaled so that its basal part has the length a, its component
            # along the c axis is c and its length √(a² + c²).
            for table in tables[4:5], (s, _) in table
                @test along(s, a_axes)
                @test norm(s) / hypot(s[1], s[2]) ≈ sqrt(1 + ca^2) rtol = 1.0e-14
                @test s[3] / hypot(s[1], s[2]) ≈ ca rtol = 1.0e-14
            end
            # Normals: the basal plane is normal to c and the prism planes
            # contain it. The first-order pyramidal planes {10-11} are
            # inclined at atan(2c/(√3 a)) to the basal plane, with normals
            # over <10-10>; the second-order ones {11-22} at atan(c/a), with
            # normals over <11-20>.
            @test all(n[1] == 0.0 && n[2] == 0.0 for (_, n) in tables[1])
            @test all(n[3] == 0.0 && along(n, m_axes) for (_, n) in tables[2])
            for table in tables[[3, 4]], (_, n) in table
                @test along(n, m_axes)
                @test hypot(n[1], n[2]) / abs(n[3]) ≈ 2 * ca / sqrt(3.0) rtol = 1.0e-14
            end
            for (_, n) in tables[5]
                @test along(n, a_axes)
                @test hypot(n[1], n[2]) / abs(n[3]) ≈ ca rtol = 1.0e-14
            end
        end
    end

    @testset "LCM's First-Order Pyramidal <c+a> Direction" begin
        # LCM's create_slip_systems.py gives the sixth system of this family
        # (its HCP system 18) the direction (1, 0, 2c/a) = 1/3[2 -1 -1 6] on
        # the plane (-1 0 1 1), with normal (-√3 c/a, -c/a, √3). That
        # direction is not a <c+a> vector and is not in the plane; the
        # <c+a> direction a1 + c = 1/3[2 -1 -1 3] = (1, 0, c/a) is.
        ca = 1.624
        s, n = SVector{3,Float64}.(Norma.hexagonal_slip_systems("hcp pyramidal c+a first order", ca)[6])
        @test normalize(n) ≈ normalize(SVector(-sqrt(3.0) * ca, -ca, sqrt(3.0))) rtol = 1.0e-14
        @test normalize(s) ≈ normalize(SVector(1.0, 0.0, ca)) rtol = 1.0e-14
        @test dot(s, n) ≈ 0.0 atol = 1.0e-14
        @test dot(SVector(1.0, 0.0, 2 * ca), n) ≈ sqrt(3.0) * ca rtol = 1.0e-14
    end

    @testset "Hexagonal Schmid Factors" begin
        # Uniaxial loading along c: the <a> systems have no resolved shear
        # stress; every first-order <c+a> system has cos λ = (c/a)/√(1 + (c/a)²)
        # and cos φ = √3/√(3 + 4 (c/a)²), every second-order one
        # cos φ = 1/√(1 + (c/a)²): about 0.40 and 0.45 for magnesium.
        # Along a1: no basal slip; the prismatic system that slips along a1
        # has its plane parallel to the load, and the two others have
        # cos φ = √3/2 and cos λ = 1/2. At 45° from c towards a1, the basal
        # system along a1 has the largest possible factor, 1/2.
        ca = 1.624
        material = Norma.create_material(hcp_crystal(; ca=ca))
        function schmid_factors(load, family)
            l = normalize(SVector{3,Float64}(load))
            return sort([abs(dot(l, material.s[α]) * dot(l, material.n[α])) for α in material.families[family].systems])
        end
        c_axis, a1, inclined = [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]
        @test all(maximum(schmid_factors(c_axis, f)) < 1.0e-15 for f in 1:3)
        @test schmid_factors(c_axis, 4) ≈ fill(sqrt(3.0) * ca / sqrt((1 + ca^2) * (3 + 4 * ca^2)), 12) rtol = 1.0e-14
        @test schmid_factors(c_axis, 5) ≈ fill(ca / (1 + ca^2), 6) rtol = 1.0e-14
        @test maximum(schmid_factors(a1, 1)) < 1.0e-15
        @test schmid_factors(a1, 2) ≈ [0.0, sqrt(3.0) / 4, sqrt(3.0) / 4] atol = 1.0e-14
        @test schmid_factors(inclined, 1) ≈ [0.25, 0.25, 0.5] rtol = 1.0e-14
    end

    @testset "Elastic Limit" begin
        # Without slip systems the model is anisotropic St. Venant-Kirchhoff.
        params = fcc_crystal()
        params["slip families"] = []
        crystal = Norma.create_material(params)
        svk = Norma.create_material(
            Norma.Parameters(
                "model" => "anisotropic Saint-Venant Kirchhoff",
                "symmetry" => "cubic",
                "C11" => 204600.0,
                "C12" => 137700.0,
                "C44" => 126200.0,
            ),
        )
        Q = rotation_about([1.0, 2.0, 3.0], 0.8)
        F = @SMatrix [1.003 0.002 -0.001; 0.001 0.998 0.002; -0.002 0.001 1.001]
        W1, P1, A1, _ = Norma.constitutive_update(crystal, F, Q, Norma.initial_state(crystal), 0.1)
        W2, P2, A2, _ = Norma.constitutive_update(svk, F, Q, Float64[], 0.1)
        @test W1 ≈ W2 rtol = 1.0e-12
        @test P1 ≈ P2 rtol = 1.0e-12
        @test A1 ≈ A2 rtol = 1.0e-8
    end

    # One plastic step of an FCC crystal at a general orientation, from the
    # virgin state, with a soft flow rule so that the trial stress stays
    # within the range of the power law.
    Q = rotation_about([0.3, -0.5, 1.0], 0.7)
    F = SMatrix{3,3,Float64,9}(I + 1.5e-3 * [1.0 0.4 -0.3; 0.1 -0.4 0.25; -0.3 0.2 -0.5])
    soft = Dict(
        "type" => "linear minus recovery",
        "hardening modulus" => 355.0,
        "recovery modulus" => 2.9,
        "initial hardening state" => 60.0,
    )
    soft_crystal(; kwargs...) = fcc_crystal(; rate_exponent=10.0, reference_rate=1.0e-3, hardening=soft, kwargs...)

    @testset "Consistent Tangent, $flow, $scheme, $residual" for flow in ("power law", "power law with drag"),
        (scheme, residual) in (("implicit", "slip"), ("implicit", "slip and hardness"), ("explicit", "slip"))

        # The drag coefficient puts the crossover of the power law and the
        # drag, d |τ/g|^(m-1) = 1, within the range of τ/g of the step.
        drag = flow == "power law" ? nothing : 0.05
        material = Norma.create_material(soft_crystal(; scheme=scheme, residual=residual, drag=drag))
        state = Norma.initial_state(material)
        Δt = 0.01
        _, P, AA, state_new = Norma.constitutive_update(material, F, Q, state, Δt)
        N = Norma.number_slip_systems(material)
        @test maximum(abs, state_new[10:(9 + N)]) > 1.0e-6
        if drag !== nothing
            # The drag slows the slip down.
            plain = Norma.create_material(soft_crystal(; scheme=scheme, residual=residual))
            _, _, _, state_plain = Norma.constitutive_update(plain, F, Q, state, Δt; need_tangent=false)
            @test maximum(abs, state_new[10:(9 + N)]) < 0.8 * maximum(abs, state_plain[10:(9 + N)])
        end
        @test det(SMatrix{3,3,Float64,9}(state_new[1:9])) ≈ 1.0 atol = 1.0e-14
        h = 1.0e-8
        scale = maximum(abs, AA)
        for k in 1:3, l in 1:3
            dF = SMatrix{3,3,Float64,9}((i, j) == (k, l) ? h : 0.0 for i in 1:3, j in 1:3)
            _, Pp, _, _ = Norma.constitutive_update(material, F + dF, Q, state, Δt; need_tangent=false)
            _, Pm, _, _ = Norma.constitutive_update(material, F - dF, Q, state, Δt; need_tangent=false)
            @test maximum(abs, (Pp - Pm) / 2h - AA[:, :, k, l]) < 1.0e-5 * scale
        end
    end

    @testset "Frame Invariance" begin
        # Rotating the sample frame rotates the lattice, the deformation and
        # the stress together; the slips do not change. A model that rotated
        # its elasticity and its slip systems differently would fail this.
        material = Norma.create_material(soft_crystal())
        state = Norma.initial_state(material)
        R = rotation_about([-1.0, 0.4, 0.9], 1.1)
        _, P, _, s = Norma.constitutive_update(material, F, Q, state, 0.01)
        _, P_R, _, s_R = Norma.constitutive_update(material, R * F * R', R * Q, state, 0.01)
        σ = P * F' / det(F)
        σ_R = P_R * (R * F * R')' / det(F)
        @test σ_R ≈ R * σ * R' rtol = 1.0e-10
        N = Norma.number_slip_systems(material)
        @test s_R[10:(9 + N)] ≈ s[10:(9 + N)] rtol = 1.0e-10
        @test SMatrix{3,3,Float64,9}(s_R[1:9]) ≈ R * SMatrix{3,3,Float64,9}(s[1:9]) * R' rtol = 1.0e-10
    end

    @testset "Objectivity" begin
        material = Norma.create_material(soft_crystal())
        state = Norma.initial_state(material)
        R = rotation_about([2.0, -1.0, 0.5], 0.9)
        _, P, _, s = Norma.constitutive_update(material, F, Q, state, 0.01)
        _, P_R, _, s_R = Norma.constitutive_update(material, R * F, Q, state, 0.01)
        @test P_R ≈ R * P rtol = 1.0e-10
        @test s_R ≈ s rtol = 1.0e-10
    end

    @testset "Equivalent Orientations" begin
        # A quarter turn about a cube axis maps the FCC lattice and its slip
        # systems onto themselves.
        material = Norma.create_material(soft_crystal())
        state = Norma.initial_state(material)
        G = rotation_about([0.0, 0.0, 1.0], π / 2)
        W1, P1, _, s1 = Norma.constitutive_update(material, F, Q, state, 0.01)
        W2, P2, _, s2 = Norma.constitutive_update(material, F, Q * G, state, 0.01)
        @test W1 ≈ W2 rtol = 1.0e-10
        @test P1 ≈ P2 rtol = 1.0e-8
        N = Norma.number_slip_systems(material)
        @test sort(abs.(s1[10:(9 + N)])) ≈ sort(abs.(s2[10:(9 + N)])) rtol = 1.0e-8
        @test s1[end] ≈ s2[end] rtol = 1.0e-8
    end

    # A single slip system along x on the plane normal to y, sheared along it.
    single_slip(; hardening) = Norma.create_material(
        fcc_crystal(; systems=[Dict("direction" => [1.0, 0.0, 0.0], "normal" => [0.0, 1.0, 0.0])], hardening=hardening),
    )
    shear(rate) = t -> SMatrix{3,3,Float64,9}([1.0 rate*t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

    @testset "Rate Sensitivity" begin
        # Once the elastic strain saturates, the slip rate equals the applied
        # shear rate and τ = g (γ̇/γ̇0)^(1/m).
        g = 122.0
        material = single_slip(; hardening=Dict("type" => "none", "initial hardening state" => g))
        names = Norma.internal_variable_names(material)
        for rate in (1.0e-3, 1.0e-1)
            times = collect(range(0.0, 0.02 / rate; length=201))
            history = Norma.material_point_history(material, shear(rate), times)
            τ = history.state[end][findfirst(==("tau_1"), names)]
            @test τ ≈ g * rate^(1 / 20) rtol = 1.0e-3
        end
    end

    @testset "Power Law With Drag" begin
        # The power law and the linear drag γ̇0 x/d, x = τ/g, in series:
        # 1/γ̇ = 1/γ̇_power + 1/γ̇_drag. The slower mechanism controls the rate:
        # the power law below the crossover d |x|^(m-1) = 1, at x ≈ 1.04 here,
        # and the drag above it.
        γ̇0, m, d, g = 1.0e-3, 20.0, 0.5, 122.0
        material = Norma.create_material(fcc_crystal(; rate_exponent=m, reference_rate=γ̇0, drag=d))
        rule = material.families[1].flow
        @test rule isa Norma.PowerLawDragFlow
        @test (rule.γ̇0, rule.m, rule.drag) == (γ̇0, m, d)
        power = Norma.PowerLawFlow(γ̇0, m, Norma.power_law_ratio_range(m)...)
        rate(flow, x) = Norma.slip_rate(flow, x * g, g)
        drag_rate(x) = γ̇0 * x / d
        for x in (-6.0, -1.1, -0.7, 0.3, 0.9, 1.0, 1.04, 1.1, 1.5, 3.0, 6.0)
            @test rate(rule, x) ≈ 1 / (1 / rate(power, x) + 1 / drag_rate(x)) rtol = 1.0e-14
            @test rate(rule, -x) == -rate(rule, x)
        end
        @test issorted([rate(rule, x) for x in range(-5.0, 5.0; length=101)])
        # Where the drag is negligible, d |x|^(m-1) ≪ 1, the rule is the power
        # law: the relative difference is d |x|^(m-1)/(1 + d |x|^(m-1)), and
        # none once d |x|^(m-1) is below the round-off.
        for d_small in (1.0e-3, 1.0e-6, 1.0e-9)
            small = Norma.PowerLawDragFlow(γ̇0, m, d_small, Norma.power_law_ratio_range(m)...)
            ratio = d_small * 1.2^(m - 1)
            @test 1 - rate(small, 1.2) / rate(power, 1.2) ≈ ratio / (1 + ratio) rtol = 1.0e-6
        end
        negligible = Norma.PowerLawDragFlow(γ̇0, m, 1.0e-30, rule.min_ratio, rule.max_ratio)
        @test rate(negligible, 1.2) ≈ rate(power, 1.2) rtol = 1.0e-15
        @test rate(rule, 0.1) ≈ rate(power, 0.1) rtol = 1.0e-15
        # At large x the rule is the drag, linear in τ.
        for x in (20.0, 40.0)
            @test rate(rule, x) ≈ drag_rate(x) rtol = 1.0e-14
            @test rate(rule, 2 * x) ≈ 2 * rate(rule, x) rtol = 1.0e-14
        end
        # Beyond the range of the power law it stays the drag, where the power
        # law fails; below the range both are zero.
        @test_throws DomainError rate(power, 2 * rule.max_ratio)
        @test rate(rule, 2 * rule.max_ratio) ≈ drag_rate(2 * rule.max_ratio) rtol = 1.0e-14
        @test rate(rule, 0.5 * rule.min_ratio) == 0.0
    end

    @testset "Power Law With Drag At A Material Point" begin
        # Single slip sheared at the rate γ̇: once the elastic strain saturates
        # the slip rate is γ̇, and x = τ/g solves γ̇0 x^m/(1 + d x^(m-1)) = γ̇.
        # With d = 1 the crossover is at x = 1, reached at γ̇ = γ̇0/2. Well
        # below it τ follows the power law, g (γ̇/γ̇0)^(1/m); well above it the
        # drag, g d γ̇/γ̇0.
        g, m, d = 122.0, 20.0, 1.0
        material = Norma.create_material(
            fcc_crystal(;
                systems=[Dict("direction" => [1.0, 0.0, 0.0], "normal" => [0.0, 1.0, 0.0])],
                drag=d,
                hardening=Dict("type" => "none", "initial hardening state" => g),
            ),
        )
        names = Norma.internal_variable_names(material)
        # The steady x by bisection; x^m/(1 + d x^(m-1)) increases with x and
        # is at least x/(1 + d) for x ≥ 1.
        function steady_ratio(rate)
            lower, upper = 0.0, max(1.0, rate * (1 + d))
            for _ in 1:200
                middle = (lower + upper) / 2
                middle^m / (1 + d * middle^(m - 1)) < rate ? (lower = middle) : (upper = middle)
            end
            return (lower + upper) / 2
        end
        for (rate, limit) in ((1.0e-4, g * 1.0e-4^(1 / m)), (0.5, g), (10.0, g * d * 10.0))
            times = collect(range(0.0, 0.2 / rate; length=201))
            history = Norma.material_point_history(material, shear(rate), times)
            τ = history.state[end][findfirst(==("tau_1"), names)]
            γ̇ = history.state[end][findfirst(==("gamma_dot_1"), names)]
            @test γ̇ ≈ rate rtol = 1.0e-6
            @test τ ≈ g * steady_ratio(rate) rtol = 1.0e-6
            @test τ ≈ limit rtol = 1.0e-4
        end
    end

    @testset "Linear Minus Recovery In Closed Form" begin
        H, R, g0 = 355.0, 2.9, 122.0
        hardening = Dict(
            "type" => "linear minus recovery",
            "hardening modulus" => H,
            "recovery modulus" => R,
            "initial hardening state" => g0,
        )
        material = single_slip(; hardening=hardening)
        names = Norma.internal_variable_names(material)
        history = Norma.material_point_history(material, shear(1.0e-2), collect(range(0.0, 5.0; length=101)))
        γ = [s[findfirst(==("gamma_1"), names)] for s in history.state]
        g = [s[findfirst(==("tau_hard_1"), names)] for s in history.state]
        @test γ[end] > 0.02
        @test g ≈ g0 .+ H / R .* (1 .- exp.(-R .* abs.(γ))) rtol = 1.0e-12
    end

    # Uniaxial stress along sample x at the strain rate 1e-3.
    stretch = t -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + 1.0e-3 * t, 1.0, 1.0]))

    @testset "Schmid Factors, $systems" for (systems, N_001, m_001, N_111, m_111) in (
        ("fcc", 8, 1 / sqrt(6.0), 6, 2 / (3 * sqrt(6.0))),
        ("bcc 112", 4, sqrt(2.0) / 3, 3, 2 * sqrt(2.0) / 9),
        ("bcc 123", 8, 3 / sqrt(42.0), 6, 2 / sqrt(42.0)),
    )
        # In uniaxial stress along d at a strain rate ε̇, the N systems with the
        # largest Schmid factor m_s carry the plastic flow, each at the rate
        # ε̇/(N m_s), so the steady flow stress is (g/m_s) (ε̇/(N m_s γ̇0))^(1/m).
        # N and m_s are those of "Cubic Schmid Factors"; in all three sets the
        # Schmid factors along [001] and [111] are in the ratio 3/2, so the
        # flow stresses are close to it.
        material = Norma.create_material(
            fcc_crystal(;
                systems=systems,
                rate_exponent=100.0,
                hardening=Dict("type" => "none", "initial hardening state" => 100.0),
            ),
        )
        times = collect(range(0.0, 10.0; length=101))
        flow_stress(d) = Norma.material_point_history(
            material, stretch, times; orientation=orientation_along_x(d), free=UNIAXIAL_FREE
        ).σ[end][
            1, 1
        ]
        σ_001 = flow_stress([0.0, 0.0, 1.0])
        σ_111 = flow_stress([1.0, 1.0, 1.0])
        flow(N, m_s) = 100.0 / m_s * (1.0e-3 / (N * m_s))^(1 / 100)
        @test σ_001 ≈ flow(N_001, m_001) rtol = 2.0e-3
        @test σ_111 ≈ flow(N_111, m_111) rtol = 2.0e-3
        @test σ_111 / σ_001 ≈ 1.5 rtol = 0.02
    end

    @testset "HCP Crystal In Uniaxial Stress, $load" for (load, d, family, m_s, N) in (
        ("45° from c", [1.0, 0.0, 1.0], 1, 0.5, 1),
        ("along a1", [1.0, 0.0, 0.0], 2, sqrt(3.0) / 4, 2),
        ("along c", [0.0, 0.0, 1.0], 5, 1.624 / (1 + 1.624^2), 6),
    )
        # The family that yields first has the smallest g/m_s, from the
        # Schmid factors of "Hexagonal Schmid Factors" and the resistances
        # of hcp_crystal, 10, 40, 60, 80 and 80 MPa:
        #   45° from c: basal, at 20 MPa; every other family above 180 MPa.
        #   Along a1: no basal slip; prismatic at 92 MPa, pyramidal <a> at
        #     157 MPa, the <c+a> families above 179 MPa.
        #   Along c: no <a> slip; the second-order pyramidal <c+a> family at
        #     179 MPa, the first-order one at 200 MPa.
        # The resistances do not harden, so the stress settles at the flow
        # stress of that family and no other family yields. The steps are
        # short enough for the soft basal family: with m = 100 and the default
        # maximum slip rate, ln(max Float64), the local update fails beyond
        # τ/g = 709^(1/100) ≈ 1.07, and a strain step of 2.5e-5 raises the
        # basal τ/g by about 0.05.
        resistances = (10.0, 40.0, 60.0, 80.0, 80.0)
        material = Norma.create_material(hcp_crystal(; ca=1.624, resistances=resistances))
        times = collect(range(0.0, 10.0; length=401))
        history = Norma.material_point_history(
            material, stretch, times; orientation=orientation_along_x(d), free=UNIAXIAL_FREE
        )
        family_slips(state) = [sum(abs, state[9 .+ f.systems]) for f in material.families]
        # At the first yield and at the end, the family carries the slip.
        onset = findfirst(state -> sum(family_slips(state)) > 1.0e-4, history.state)
        @test onset !== nothing
        for state in (history.state[onset], history.state[end])
            slips = family_slips(state)
            @test slips[family] > 0.999 * sum(slips)
        end
        # The N systems of the family with the largest Schmid factor share the
        # slip, and the stress settles at their flow stress. The flow rule
        # sees the Mandel stress, which for small elastic strains is the
        # Kirchhoff stress J σ; the Cauchy stress is lower by the elastic
        # volume change, up to 1.6e-3 along c.
        γ = abs.(history.state[end][9 .+ material.families[family].systems])
        @test count(>(0.5 * maximum(γ)), γ) == N
        flow = resistances[family] / m_s * (1.0e-3 / (N * m_s))^(1 / 100)
        @test det(history.F[end]) * history.σ[end][1, 1] ≈ flow rtol = 1.0e-3
    end

    @testset "Explicit And Implicit Converge Together" begin
        # Both schemes are first order; the explicit one is only stable for
        # Δt (∂γ̇/∂τ) μ < 2, hence the soft rate exponent.
        hardening = Dict(
            "type" => "linear minus recovery",
            "hardening modulus" => 355.0,
            "recovery modulus" => 2.9,
            "initial hardening state" => 122.0,
        )
        function final_stress(scheme, steps)
            material = Norma.create_material(
                fcc_crystal(;
                    scheme=scheme,
                    rate_exponent=5.0,
                    systems=[Dict("direction" => [1.0, 0.0, 0.0], "normal" => [0.0, 1.0, 0.0])],
                    hardening=hardening,
                ),
            )
            return Norma.material_point_history(material, shear(1.0e-2), collect(range(0.0, 2.0; length=steps + 1))).σ[end][
                1, 2
            ]
        end
        gap_coarse = abs(final_stress("explicit", 400) - final_stress("implicit", 400))
        gap_fine = abs(final_stress("explicit", 1600) - final_stress("implicit", 1600))
        @test 3 < gap_coarse / gap_fine < 5
    end

    @testset "Local Failure Signals A Step Cut" begin
        # A slip rate beyond the family's maximum throws a DomainError, which
        # the model turns into a cut of the global step.
        material = Norma.create_material(fcc_crystal(; scheme="explicit"))
        F_large = @SMatrix [1.0 0.5 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test_throws DomainError Norma.constitutive_update(
            material, F_large, identity3, Norma.initial_state(material), 0.01
        )
    end

    @testset "Input Errors" begin
        Norma.norma_log(0, :info, "Testing bad input. Mock code abort ...")
        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException Norma.create_material(
                fcc_crystal(; systems=[Dict("direction" => [1.0, 1.0, 0.0], "normal" => [1.0, 0.0, 0.0])])
            )
            @test_throws Norma.NormaAbortException Norma.create_material(fcc_crystal(; systems="diamond"))
            bad_flow = fcc_crystal()
            bad_flow["slip families"][1]["flow rule"]["type"] = "creep"
            @test_throws Norma.NormaAbortException Norma.create_material(bad_flow)
            bad_scheme = fcc_crystal()
            bad_scheme["integration scheme"] = "midpoint"
            @test_throws Norma.NormaAbortException Norma.create_material(bad_scheme)
        finally
            Norma.NORMA_TEST_MODE[] = false
        end
    end
end

# Inputs that would silently give another model are rejected, with a
# message, at construction.
@testset "Crystal Plasticity Input Checks" begin
    shear(t) = @SMatrix [1.0 t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    times = collect(0.0:2.5e-4:5.0e-3)
    hardness(material) = Norma.material_point_history(material, shear, times).state[end][(10 + 24):(9 + 36)]
    Norma.NORMA_TEST_MODE[] = true
    try
        # The explicit update hardens with the slip rates of the previous
        # step, which `slip predictor: none` sets to zero: it would never
        # harden. The implicit scheme only starts its Newton iteration there,
        # and converges to the same state as with `rate`.
        explicit_none = fcc_crystal(; scheme="explicit")
        explicit_none["slip predictor"] = "none"
        @test_throws Norma.NormaAbortException Norma.create_material(explicit_none)
        implicit_none = fcc_crystal()
        implicit_none["slip predictor"] = "none"
        h_none = hardness(Norma.create_material(implicit_none))
        h_rate = hardness(Norma.create_material(fcc_crystal()))
        @test maximum(h_rate) > 123.0
        @test maximum(abs.(h_none .- h_rate)) < 1.0e-8 * maximum(h_rate)
        # Every family needs a flow rule, and a flow rule or a hardening law
        # its type: without them the crystal was silently elastic or did not
        # harden. LCM's capitalized `Type` is not a type.
        no_flow = fcc_crystal()
        delete!(no_flow["slip families"][1], "flow rule")
        @test_throws Norma.NormaAbortException Norma.create_material(no_flow)
        untyped_flow = fcc_crystal()
        delete!(untyped_flow["slip families"][1]["flow rule"], "type")
        @test_throws Norma.NormaAbortException Norma.create_material(untyped_flow)
        capitalized = fcc_crystal()
        flow = capitalized["slip families"][1]["flow rule"]
        flow["Type"] = pop!(flow, "type")
        @test_throws Norma.NormaAbortException Norma.create_material(capitalized)
        untyped_hardening = fcc_crystal()
        delete!(untyped_hardening["slip families"][1]["hardening law"], "type")
        @test_throws Norma.NormaAbortException Norma.create_material(untyped_hardening)
        # The input validation reports misspelled keys inside the slip
        # families, the flow rules, the hardening laws, the slip systems and
        # the local solver.
        typos = fcc_crystal()
        family = typos["slip families"][1]
        family["flow rules"] = family["flow rule"]
        family["hardening law"]["Recovery Modulus"] = 2.9
        family["flow rule"]["rate exponant"] = 20.0
        family["slip systems"] = [Dict("direction" => [-1.0, 1.0, 0.0], "normals" => [1.0, 1.0, 1.0])]
        typos["local solver"] = Dict("maximum iteration" => 5)
        messages = Norma.validate_material!(
            String[], Dict("blocks" => Dict("grain" => "copper"), "copper" => typos), "test"
        )
        for key in ("flow rules", "Recovery Modulus", "rate exponant", "normals", "maximum iteration")
            @test any(occursin("\"$key\"", message) for message in messages)
        end
        @test isempty(
            Norma.validate_material!(
                String[], Dict("blocks" => Dict("grain" => "copper"), "copper" => fcc_crystal()), "test"
            ),
        )
        # Parameters out of their range, and families or systems that are
        # not lists, abort at construction with a message that names them,
        # instead of failing at the first step, being ignored, or raising a
        # Julia error.
        family(params) = params["slip families"][1]
        changes = (
            params -> (family(params)["flow rule"]["rate exponent"] = -2.0),
            params -> (family(params)["flow rule"]["rate exponent"] = 0.0),
            params -> (family(params)["flow rule"]["reference slip rate"] = 0.0),
            params -> (family(params)["hardening law"]["initial hardening state"] = 0.0),
            params -> (family(params)["hardening law"]["recovery modulus"] = -1.0),
            params -> (family(params)["maximum slip rate"] = 0.0),
            params -> (params["slip families"] = Dict("copper" => family(params))),
            params ->
                (family(params)["slip systems"] = Dict("direction" => [-1.0, 1.0, 0.0], "normal" => [1.0, 1.0, 1.0])),
        )
        for change in changes
            params = fcc_crystal()
            change(params)
            @test_throws Norma.NormaAbortException Norma.create_material(params)
        end
        @test_throws Norma.NormaAbortException Norma.create_material(fcc_crystal(; drag=0.0))
        @test_throws Norma.NormaAbortException Norma.create_material(hcp_crystal(; ca=0.0))
        saturation = Dict(
            "type" => "saturation",
            "hardening rate" => 100.0,
            "initial saturation stress" => 202.0,
            "saturation exponent" => -0.1,
            "initial hardening state" => 90.0,
        )
        @test_throws Norma.NormaAbortException Norma.create_material(fcc_crystal(; hardening=saturation))
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
end

# A NaN slip resistance gave a zero slip rate, since every comparison of the
# flow rules is false for it, so the point turned elastic without an error.
# Saturation hardening with its saturation stress at the initial hardening
# state gave 0/0 at the first plastic step.
@testset "Crystal Plasticity Degenerate Slip Resistance" begin
    shear(t) = @SMatrix [1.0 t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    times = collect(0.0:2.5e-4:1.0e-2)
    saturated = fcc_crystal(;
        hardening=Dict(
            "type" => "saturation",
            "hardening rate" => 100.0,
            "initial saturation stress" => 60.0,
            "initial hardening state" => 60.0,
        ),
    )
    material = Norma.create_material(saturated)
    history = Norma.material_point_history(material, shear, times)
    h = history.state[end][34:45]
    @test all(isfinite, h) && all(h .== 60.0)
    # Perfect plasticity at 60 MPa: the shear stress stays near the flow
    # stress of the most stressed systems instead of growing elastically.
    @test maximum(abs(σ[1, 2]) for σ in history.σ) < 1.5 * 60.0 / 0.4
    # A NaN slip resistance is a math error, which cuts the step.
    state = Norma.initial_state(material)
    state[34:45] .= NaN
    F = shear(1.0e-3)
    identity3 = SMatrix{3,3,Float64,9}(I)
    @test_throws DomainError Norma.constitutive_update(material, F, identity3, state, 1.0e-3)
    # A saturation stress below the initial hardening state hardens without
    # bound, and is rejected.
    Norma.NORMA_TEST_MODE[] = true
    try
        below = deepcopy(saturated)
        below["slip families"][1]["hardening law"]["initial saturation stress"] = 50.0
        @test_throws Norma.NormaAbortException Norma.create_material(below)
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
end

# Norma evaluates materials with Δt = 0 in the initial equilibrium solve, the
# Schwarz blended energy and at the first time of the drivers. The slips
# cannot change then, and the response is the elastic one on Fp_n, however
# large the resolved shear: the slip-rate cap used to reject it.
@testset "Crystal Plasticity At A Zero Time Step" begin
    F = @SMatrix [1.0 0.01 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    identity3 = SMatrix{3,3,Float64,9}(I)
    elastic = Norma.create_material(
        Norma.Parameters(
            "model" => "anisotropic Saint-Venant Kirchhoff",
            "symmetry" => "cubic",
            "C11" => 204600.0,
            "C12" => 137700.0,
            "C44" => 126200.0,
        ),
    )
    _, P_elastic, _, _ = Norma.constitutive_update(elastic, F, identity3, Float64[], 0.0)
    for scheme in ("implicit", "explicit")
        material = Norma.create_material(fcc_crystal(; scheme=scheme))
        _, P, _, state = Norma.constitutive_update(material, F, identity3, Norma.initial_state(material), 0.0)
        @test maximum(abs.(P - P_elastic)) < 1.0e-10 * maximum(abs.(P_elastic))
        @test state[10:21] == zeros(12)
    end
end

# The explicit update had no bound on the slip increment: past the range of
# the exponential map it gave NaN, which every guard of the stress update
# let through, instead of the error that cuts the step.
@testset "Crystal Plasticity Explicit Overflow Cuts The Step" begin
    F = SMatrix{3,3,Float64,9}(I + 0.0005 * [0.0 1.0 0.0; 0.0 0.2 0.0; 0.0 0.0 -0.2])
    Q = Norma.orientation_from_bunge_angles(0.3, 0.7, 1.1)
    material = Norma.create_material(fcc_crystal(; scheme="explicit", rate_exponent=1.0))
    # Fp of order 1e177, and slips of order 1e3, where the map overflows.
    for Δt in (1.0e3, 1.0e4)
        @test_throws DomainError Norma.constitutive_update(material, F, Q, Norma.initial_state(material), Δt)
    end
end

# Saturation hardening with the slip residual (and the explicit scheme)
# integrates the hardness by forward Euler, which overshoots the saturation
# stress once a step slips more than about (gs - g0)/(2 h0), and diverges
# beyond twice that. The overshoot now cuts the step; the slip and hardness
# residual, which integrates it implicitly, takes the same step. The
# slip-rate cap is lifted, as in the LCM comparisons, so that the local
# Newton iteration can take these steps. From the rate predictor, it needs
# more than its default 50 iterations for so long a step: with 50, Julia 1.13
# on macOS stops short of convergence where 1.12 converges.
@testset "Crystal Plasticity Saturation Hardening Does Not Overshoot" begin
    saturation = Dict(
        "type" => "saturation",
        "hardening rate" => 20000.0,
        "initial saturation stress" => 70.0,
        "initial hardening state" => 60.0,
    )
    crystal(residual) = begin
        params = fcc_crystal(; hardening=saturation, residual=residual)
        params["slip families"][1]["maximum slip rate"] = floatmax(Float64)
        params["local solver"] = Dict{String,Any}("maximum iterations" => 200)
        return Norma.create_material(params)
    end
    shear(t) = @SMatrix [1.0 t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    identity3 = SMatrix{3,3,Float64,9}(I)
    hardness(state) = state[34:45]
    # One step of 2e-3 in shear slips 2.5e-3 in all, which forward Euler
    # takes to a hardness of about 85, beyond the saturation stress of 70.
    slip = crystal("slip")
    @test_throws DomainError Norma.constitutive_update(slip, shear(2.0e-3), identity3, Norma.initial_state(slip), 0.01)
    implicit = crystal("slip and hardness")
    _, _, _, state = Norma.constitutive_update(implicit, shear(2.0e-3), identity3, Norma.initial_state(implicit), 0.01)
    @test all(60.0 .<= hardness(state) .<= 70.0)
    @test maximum(hardness(state)) > 62.0
    # Short steps stay below the saturation stress with either residual.
    times = collect(0.0:2.5e-4:2.0e-2)
    for material in (slip, implicit)
        history = Norma.material_point_history(material, shear, times)
        @test all(all(60.0 .<= hardness(state) .<= 70.0) for state in history.state)
        @test maximum(hardness(history.state[end])) > 61.0
    end
end

# The thermal activation rule keeps the rate at γ̇0 beyond ξ = 1, where the
# barrier vanishes; LCM returns NaN there for a fractional q exponent.
@testset "Crystal Plasticity Thermal Activation Beyond The Barrier" begin
    for q in (1.5, 2.0)
        rule = Norma.create_flow_rule(
            Dict(
                "type" => "thermal activation",
                "reference slip rate" => 2.0,
                "activation energy" => 3.0,
                "thermal resistance" => 10.0,
                "p exponent" => 0.5,
                "q exponent" => q,
            ),
        )
        g = 50.0
        @test Norma.slip_rate(rule, 40.0, g) == 0.0
        @test Norma.slip_rate(rule, g, g) == 0.0
        @test 0.0 < Norma.slip_rate(rule, 55.0, g) < 2.0
        @test Norma.slip_rate(rule, 60.0 - 1.0e-9, g) ≈ 2.0 rtol = 1.0e-6
        @test Norma.slip_rate(rule, 60.0, g) == 2.0
        @test Norma.slip_rate(rule, 75.0, g) == 2.0
        @test Norma.slip_rate(rule, -75.0, g) == -2.0
        @test Norma.slip_rate(rule, -55.0, g) == -Norma.slip_rate(rule, 55.0, g)
    end
end

# The explicit scheme takes the slip resistance from the hardening that the
# previous step's rates give (the rate predictor). That matters for the laws
# whose rate depends on the hardness, saturation and dislocation density, and
# no explicit LCM case uses one.
@testset "Crystal Plasticity Explicit Scheme Predicts The Hardening" begin
    saturation = Dict(
        "type" => "saturation",
        "hardening rate" => 100.0,
        "initial saturation stress" => 70.0,
        "initial hardening state" => 60.0,
    )
    material = Norma.create_material(fcc_crystal(; scheme="explicit", hardening=saturation))
    N = 12
    state = Norma.initial_state(material)
    γ̇_n = [1.0e-3 * k for k in 1:N]
    state[(10 + N):(9 + 2N)] .= γ̇_n
    h_n = state[(10 + 2N):(9 + 3N)]
    Δt = 1.0e-3
    F = @SMatrix [1.0 1.0e-3 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    identity3 = SMatrix{3,3,Float64,9}(I)
    _, _, _, state_new = Norma.constitutive_update(material, F, identity3, state, Δt)
    # The definition: g = h_pred, the hardening of the previous rates, and
    # γ̇ = γ̇(τ at Fp_n, g).
    h_pred, _ = Norma.harden(material, Δt, γ̇_n, h_n, h_n)
    @test maximum(abs, h_pred .- h_n) > 1.0e-3
    P = [SMatrix{3,3,Float64,9}(material.s[α] * material.n[α]') for α in 1:N]
    Fp_n = SMatrix{3,3,Float64,9}(state[1:9])
    _, _, τ_trial, _ = Norma.crystal_stress(material.CC, F, Fp_n, P)
    @test state_new[(10 + N):(9 + 2N)] ≈ Norma.slip_rates(material, τ_trial, h_pred) rtol = 1.0e-12
end

# LCM's dislocation density law does not harden when annihilation exceeds
# generation (the driver is max(k1 √ρF - k2 ρ, 0)), so a crystal that starts
# above saturation keeps its densities while it slips.
@testset "Crystal Plasticity Dislocation Density Above Saturation" begin
    density = Dict(
        "type" => "dislocation density",
        "geometric factor" => 0.4,
        "initial hardening state" => 80.0,
        "generation factor" => 5000.0,
        "annihilation factor" => 50000.0,
        "shear modulus" => 70700.0,
        "burgers vector magnitude" => 2.87e-4,
    )
    material = Norma.create_material(fcc_crystal(; hardening=density))
    shear(t) = @SMatrix [1.0 t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    history = Norma.material_point_history(material, shear, collect(0.0:5.0e-4:1.0e-2))
    @test maximum(abs, history.state[end][10:21]) > 1.0e-3
    @test all(state[34:45] == fill(80.0, 12) for state in history.state)
end

# The `local solver` options take effect: with tolerances that the first
# iterate meets, `minimum iterations` still runs Newton's method to the
# solution, and without it the update keeps the predictor.
@testset "Crystal Plasticity Local Solver Minimum Iterations" begin
    F = @SMatrix [1.0 2.0e-3 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    identity3 = SMatrix{3,3,Float64,9}(I)
    slips(params) = begin
        material = Norma.create_material(params)
        Norma.constitutive_update(material, F, identity3, Norma.initial_state(material), 1.0e-2)[4][10:21]
    end
    crystal() = begin
        params = fcc_crystal()
        params["slip families"][1]["maximum slip rate"] = floatmax(Float64)
        return params
    end
    reference = slips(crystal())
    @test maximum(abs, reference) > 1.0e-5
    loose(minimum_iterations) = begin
        params = crystal()
        params["local solver"] = Dict(
            "absolute tolerance" => 1.0e100,
            "relative tolerance" => 1.0e100,
            "minimum iterations" => minimum_iterations,
            "maximum iterations" => 20,
        )
        return params
    end
    @test maximum(abs, slips(loose(0))) == 0.0
    @test slips(loose(10)) ≈ reference rtol = 1.0e-8
end

# The update infers a concrete return type, as the assembly's guard expects:
# the Jacobians of the local solve are typed, and the closures of the update
# capture variables assigned once.
@testset "Crystal Plasticity Update Is Type-Stable" begin
    F = @SMatrix [1.0 2.0e-3 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    identity3 = SMatrix{3,3,Float64,9}(I)
    for scheme in ("implicit", "explicit")
        material = Norma.create_material(fcc_crystal(; scheme=scheme))
        result = @inferred Norma.constitutive_update(material, F, identity3, Norma.initial_state(material), 1.0e-2)
        @test result[4] isa Vector{Float64}
    end
end

# The translation of LCM materials (norma_crystal_parameters): the cases it
# maps, the inputs it rejects instead of turning them into another material,
# and its branch for a left-handed LCM basis, which no LCM file has.
@testset "LCM Translation" begin
    basis = [[1.0, 1.0, 0.0] / sqrt(2), [-1.0, 1.0, 0.0] / sqrt(2), [0.0, 0.0, -1.0]]
    lcm_material(elasticity=Dict{String,Any}()) = Dict{String,Any}(
        "Crystal Elasticity" => merge(
            Dict{String,Any}(
                "C11" => 300.0,
                "C12" => 100.0,
                "C44" => 100.0,
                "Basis Vector 1" => basis[1],
                "Basis Vector 2" => basis[2],
                "Basis Vector 3" => basis[3],
            ),
            elasticity,
        ),
        "Integration Scheme" => "Implicit",
        "Number of Slip Systems" => 2,
        "Slip System 1" => Dict{String,Any}("Slip Direction" => [1.0, -1.0, 0.0], "Slip Normal" => [1.0, 1.0, 1.0]),
        "Slip System 2" => Dict{String,Any}("Slip Direction" => [0.0, 1.0, -1.0], "Slip Normal" => [1.0, 1.0, 1.0]),
        "Slip System Family 0" => Dict{String,Any}(
            "Flow Rule" =>
                Dict{String,Any}("Type" => "Power Law", "Reference Slip Rate" => 1.0, "Rate Exponent" => 10.0),
            "Hardening Law" => Dict{String,Any}(
                "Type" => "Linear Minus Recovery",
                "Hardening Modulus" => 10.0,
                "Recovery Modulus" => 0.0,
                "Initial Hardening State" => 50.0,
            ),
        ),
    )
    # A left-handed basis: Norma's orientation is proper, and the sample-frame
    # slip vectors are LCM's, B s and B n with the improper B; this isotropic
    # cubic elasticity is unchanged by the reflection.
    params, _ = norma_crystal_parameters(lcm_material())
    material = Norma.create_material(params)
    Q = Norma.orientation_from_basis_vectors(params["orientation"]["basis vectors"])
    @test det(Q) ≈ 1.0
    B = hcat(basis...)
    for (α, system) in enumerate(("Slip System 1", "Slip System 2"))
        lcm = lcm_material()[system]
        @test Q * material.s[α] ≈ normalize(B * lcm["Slip Direction"])
        @test Q * material.n[α] ≈ normalize(B * lcm["Slip Normal"])
    end
    # LCM's tensor with C33 = C11 and C13 != C12 is tetragonal, with C66 = C44.
    params, _ = norma_crystal_parameters(lcm_material(Dict{String,Any}("C13" => 80.0)))
    @test params["symmetry"] == "tetragonal" && params["C13"] == 80.0 && params["C66"] == 100.0
    C = Norma.anisotropic_elastic_constants(params)
    @test C[1, 1] == C[2, 2] == C[3, 3] == 300.0 && C[1, 2] == 100.0 && C[1, 3] == C[2, 3] == 80.0
    @test C[4, 4] == C[5, 5] == C[6, 6] == 100.0
    # LCM's slip predictor is translated; the unsupported options stop it.
    none = lcm_material()
    none["Slip Predictor"] = "None"
    @test norma_crystal_parameters(none)[1]["slip predictor"] == "none"
    unsupported = (
        material -> (material["Crystal Elasticity"]["M11"] = 0.1),
        material -> (material["Crystal Elasticity"]["Reference Temperature"] = 300.0),
        material -> (material["Residual Type"] = "Constrained Slip Hardness"),
        material -> (material["Slip Predictor"] = "Solve"),
    )
    for change in unsupported
        material = lcm_material()
        change(material)
        @test_throws ErrorException norma_crystal_parameters(material)
    end
end

# A restart snapshot holds only the nodal displacement and velocity, so a
# restarted crystal would start again from Fp = I, zero slips and the initial
# hardness on the restored displacement. A restart of any material with
# internal state variables is rejected, in a single domain and in a Schwarz
# restart; materials without state still restart.
@testset "Crystal Plasticity Restart Is Rejected" begin
    crystal = fcc_crystal()
    crystal["density"] = 1000.0
    j2 = Norma.Parameters(
        "model" => "j2 plasticity",
        "elastic modulus" => 1.0e9,
        "Poisson's ratio" => 0.25,
        "density" => 1000.0,
        "yield stress" => 1.0e7,
        "hardening modulus" => 1.0e8,
    )
    elastic = Norma.Parameters(
        "model" => "neohookean", "elastic modulus" => 1.0e9, "Poisson's ratio" => 0.25, "density" => 1000.0
    )
    for schwarz in (false, true)
        message = restart_abort_message(crystal; schwarz=schwarz)
        @test message !== nothing && occursin("`crystal plasticity`", message)
        message = restart_abort_message(j2; schwarz=schwarz)
        @test message !== nothing && occursin("`j2 plasticity`", message)
        @test restart_abort_message(elastic; schwarz=schwarz) === nothing
    end
end
