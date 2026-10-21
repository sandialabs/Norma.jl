# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays

const ID3 = SMatrix{3,3,Float64,9}(I)

# Copper constants of LCM's SingleElement tests.
function copper_params(model)
    return Norma.Parameters(
        "model" => model, "symmetry" => "cubic", "C11" => 168.4, "C12" => 121.4, "C44" => 75.4, "density" => 1.0
    )
end

# Rotation by angle θ about a unit axis (Rodrigues).
function axis_rotation(axis::AbstractVector, θ::Real)
    n = normalize(SVector{3,Float64}(axis))
    K = @SMatrix [0.0 -n[3] n[2]; n[3] 0.0 -n[1]; -n[2] n[1] 0.0]
    return SMatrix{3,3,Float64,9}(I + sin(θ) * K + (1 - cos(θ)) * K * K)
end

# Uniaxial stress along the sample x axis: F11 prescribed, F21, F31 and F32
# fixed to remove the rigid rotation, the other components free.
const UNIAXIAL_STRESS_FREE = [false true true; false true true; false false true]

@testset "Anisotropic Elasticity" begin
    @testset "Voigt Conversion" begin
        C6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters(
                "symmetry" => "triclinic", ["C$i$j" => (i == j ? 10.0 * i : 0.1 * (i + j)) for i in 1:6 for j in i:6]...
            ),
        )
        CC = Norma.voigt_to_tensor(C6)
        @test Norma.tensor_to_voigt(CC) == C6
        @test all(
            CC[i, j, k, l] == CC[j, i, k, l] == CC[i, j, l, k] == CC[k, l, i, j] for
            i in 1:3, j in 1:3, k in 1:3, l in 1:3
        )
        @test CC[2, 3, 2, 3] == C6[4, 4]
        @test CC[1, 3, 1, 2] == C6[5, 6]
    end

    @testset "Symmetry Classes" begin
        C6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters(
                "symmetry" => "hexagonal", "C11" => 162.4, "C12" => 92.0, "C13" => 69.0, "C33" => 180.7, "C44" => 46.7
            ),
        )
        @test C6[1, 1] == C6[2, 2] == 162.4 && C6[3, 3] == 180.7
        @test C6[1, 3] == C6[2, 3] == 69.0 && C6[4, 4] == C6[5, 5] == 46.7
        @test C6[6, 6] ≈ (162.4 - 92.0) / 2
        T6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters(
                "symmetry" => "transversely isotropic",
                "C11" => 162.4,
                "C12" => 92.0,
                "C13" => 69.0,
                "C33" => 180.7,
                "C44" => 46.7,
            ),
        )
        @test T6 == C6
        G6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters(
                "symmetry" => "tetragonal",
                "C11" => 72.0,
                "C12" => 58.0,
                "C13" => 36.0,
                "C33" => 88.0,
                "C44" => 22.0,
                "C66" => 24.0,
            ),
        )
        @test G6[6, 6] == 24.0 && G6[1, 2] == 58.0 && G6[3, 3] == 88.0
        O6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters(
                "symmetry" => "orthotropic",
                "C11" => 10.0,
                "C22" => 11.0,
                "C33" => 12.0,
                "C12" => 1.0,
                "C13" => 2.0,
                "C23" => 3.0,
                "C44" => 4.0,
                "C55" => 5.0,
                "C66" => 6.0,
            ),
        )
        @test diag(O6) == [10.0, 11.0, 12.0, 4.0, 5.0, 6.0]
        @test O6[2, 3] == O6[3, 2] == 3.0 && O6[1, 4] == 0.0
        E, ν = 200.0, 0.3
        λ, μ = E * ν / (1 + ν) / (1 - 2ν), E / 2 / (1 + ν)
        S6 = Norma.anisotropic_elastic_constants(
            Norma.Parameters("symmetry" => "isotropic", "elastic modulus" => E, "Poisson's ratio" => ν)
        )
        @test S6[1, 1] ≈ λ + 2μ && S6[1, 2] ≈ λ && S6[4, 4] ≈ μ
    end

    @testset "Rotation" begin
        CC = Norma.create_material(copper_params("anisotropic linear elastic")).CC
        R1 = axis_rotation([1.0, 2.0, -0.5], 0.7)
        R2 = axis_rotation([-0.3, 0.4, 1.0], -1.1)
        @test Norma.rotate_tensor4(ID3, CC) == CC
        @test Norma.rotate_tensor4(R1, Norma.rotate_tensor4(R2, CC)) ≈ Norma.rotate_tensor4(R1 * R2, CC)
        # Cubic symmetry: quarter turns about the cube axes leave the tensor unchanged.
        for axis in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])
            @test Norma.rotate_tensor4(axis_rotation(axis, π / 2), CC) ≈ CC atol = 1.0e-12
        end
        @test !isapprox(Norma.rotate_tensor4(axis_rotation([0, 0, 1.0], π / 4), CC), CC; atol=1.0)
        # Transverse isotropy: any rotation about the c axis.
        H = Norma.voigt_to_tensor(
            Norma.anisotropic_elastic_constants(
                Norma.Parameters(
                    "symmetry" => "hexagonal",
                    "C11" => 162.4,
                    "C12" => 92.0,
                    "C13" => 69.0,
                    "C33" => 180.7,
                    "C44" => 46.7,
                ),
            ),
        )
        @test Norma.rotate_tensor4(axis_rotation([0, 0, 1.0], 0.37), H) ≈ H atol = 1.0e-12
        # Rotating the tensor and the strain together rotates the stress.
        ϵ = @SMatrix [1.0e-3 2.0e-4 -1.0e-4; 2.0e-4 -5.0e-4 3.0e-4; -1.0e-4 3.0e-4 2.0e-4]
        @test Norma.double_contraction(Norma.rotate_tensor4(R1, CC), R1 * ϵ * R1') ≈
            R1 * Norma.double_contraction(CC, ϵ) * R1'
    end

    @testset "Lattice Orientation Input" begin
        # LCM's SingleElement3d basis vectors are the columns of the Bunge
        # rotation for (0.12, 0.03, 0.28) rad, rounded to eight digits.
        lcm_basis = [
            [0.92107588, 0.38929489, 0.00828943],
            [-3.89366574e-01, 0.92063166, 0.02882734],
            [0.00359083, -2.97797916e-02, 0.99955003],
        ]
        Q_basis = Norma.lattice_orientation(Dict("basis vectors" => lcm_basis))
        Q_bunge = Norma.lattice_orientation(Dict("Bunge Euler angles" => [0.12, 0.03, 0.28], "angle unit" => "radians"))
        @test maximum(abs.(Q_basis - Q_bunge)) < 1.0e-8
        @test Q_basis' * Q_basis ≈ ID3 atol = 1.0e-14
        Q_degrees = Norma.lattice_orientation(Dict("Bunge Euler angles" => rad2deg.([0.12, 0.03, 0.28])))
        @test Q_degrees ≈ Q_bunge
        # φ1 = 90° alone turns lattice axis 1 onto sample y.
        Q_90 = Norma.lattice_orientation(Dict("Bunge Euler angles" => [90.0, 0.0, 0.0]))
        @test Q_90 * [1.0, 0.0, 0.0] ≈ [0.0, 1.0, 0.0] atol = 1.0e-15
        # Unnormalized input is accepted.
        @test Norma.lattice_orientation(Dict("basis vectors" => [[2.0, 2.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 3.0]])) ≈
            axis_rotation([0, 0, 1.0], π / 4)
    end

    @testset "Input Errors" begin
        Norma.norma_log(0, :info, "Testing bad input. Mock code abort ...")
        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException Norma.anisotropic_elastic_constants(Norma.Parameters("C11" => 1.0))
            @test_throws Norma.NormaAbortException Norma.anisotropic_elastic_constants(
                Norma.Parameters("symmetry" => "rhombic", "C11" => 1.0)
            )
            @test_throws Norma.NormaAbortException Norma.anisotropic_elastic_constants(
                Norma.Parameters("symmetry" => "cubic", "C11" => 1.0, "C12" => 1.0)
            )
            # C12 > C11 is not positive definite.
            @test_throws Norma.NormaAbortException Norma.anisotropic_elastic_constants(
                Norma.Parameters("symmetry" => "cubic", "C11" => 1.0, "C12" => 2.0, "C44" => 1.0)
            )
            @test_throws Norma.NormaAbortException Norma.lattice_orientation(
                Dict("basis vectors" => [[1.0, 0.0, 0.0], [0.1, 1.0, 0.0], [0.0, 0.0, 1.0]])
            )
            @test_throws Norma.NormaAbortException Norma.lattice_orientation(
                Dict("basis vectors" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
            )
            @test_throws Norma.NormaAbortException Norma.lattice_orientation(
                Dict(
                    "basis vectors" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    "Bunge Euler angles" => [0, 0, 0],
                ),
            )
            @test_throws Norma.NormaAbortException Norma.lattice_orientation(
                Dict("Bunge Euler angles" => [0, 0, 0], "angle unit" => "gradians")
            )
        finally
            Norma.NORMA_TEST_MODE[] = false
        end
    end

    @testset "Maximum Wave Modulus" begin
        copper = Norma.create_material(copper_params("anisotropic linear elastic"))
        @test Norma.get_p_wave_modulus(copper) ≈ (168.4 + 2 * 121.4 + 4 * 75.4) / 3 rtol = 1.0e-10
        # Zener ratio below one: the stiffest direction is [100].
        soft = Norma.create_material(
            Norma.Parameters(
                "model" => "anisotropic linear elastic",
                "symmetry" => "cubic",
                "C11" => 200.0,
                "C12" => 50.0,
                "C44" => 50.0,
            ),
        )
        @test Norma.get_p_wave_modulus(soft) ≈ 200.0 rtol = 1.0e-10
        E, ν = 200.0, 0.3
        λ, μ = E * ν / (1 + ν) / (1 - 2ν), E / 2 / (1 + ν)
        isotropic = Norma.create_material(
            Norma.Parameters(
                "model" => "anisotropic linear elastic",
                "symmetry" => "isotropic",
                "elastic modulus" => E,
                "Poisson's ratio" => ν,
            ),
        )
        @test Norma.get_p_wave_modulus(isotropic) ≈ λ + 2μ rtol = 1.0e-10
    end

    @testset "Isotropic Limit" begin
        # Cubic constants with C11 - C12 = 2 C44 describe an isotropic material,
        # whatever the lattice orientation.
        E, ν = 200.0, 0.3
        λ, μ = E * ν / (1 + ν) / (1 - 2ν), E / 2 / (1 + ν)
        Q = axis_rotation([1.0, -2.0, 0.5], 0.9)
        F = @SMatrix [1.05 0.02 -0.03; 0.01 0.97 0.04; 0.02 -0.01 1.02]
        for (isotropic_model, anisotropic_model) in (
            ("Saint-Venant Kirchhoff", "anisotropic Saint-Venant Kirchhoff"),
            ("linear elastic", "anisotropic linear elastic"),
        )
            isotropic = Norma.create_material(
                Norma.Parameters("model" => isotropic_model, "elastic modulus" => E, "Poisson's ratio" => ν)
            )
            cubic = Norma.create_material(
                Norma.Parameters(
                    "model" => anisotropic_model, "symmetry" => "cubic", "C11" => λ + 2μ, "C12" => λ, "C44" => μ
                ),
            )
            W1, P1, A1 = Norma.constitutive(isotropic, F)
            W2, P2, A2, _ = Norma.constitutive_update(cubic, F, Q, Float64[], 0.0)
            @test W2 ≈ W1 rtol = 1.0e-12
            @test P2 ≈ P1 rtol = 1.0e-12
            @test A2 ≈ A1 rtol = 1.0e-12
        end
    end

    @testset "Saint-Venant Kirchhoff Response" begin
        copper = Norma.create_material(copper_params("anisotropic Saint-Venant Kirchhoff"))
        Q = axis_rotation([0.3, 1.0, -0.2], 0.6)
        F = @SMatrix [1.04 0.03 -0.02; -0.01 0.98 0.05; 0.02 0.01 1.01]
        W, P, AA, _ = Norma.constitutive_update(copper, F, Q, Float64[], 0.0)
        # Objectivity: a superposed rotation rotates P and leaves W unchanged.
        R = axis_rotation([1.0, 1.0, 1.0], 1.3)
        W_R, P_R, _, _ = Norma.constitutive_update(copper, R * F, Q, Float64[], 0.0)
        @test W_R ≈ W rtol = 1.0e-12
        @test P_R ≈ R * P rtol = 1.0e-12
        # The stress is the derivative of the energy, and the tangent the
        # derivative of the stress (central differences).
        h = 1.0e-6
        for k in 1:3, l in 1:3
            dF = SMatrix{3,3,Float64,9}((i, j) == (k, l) ? h : 0.0 for i in 1:3, j in 1:3)
            Wp, Pp, _, _ = Norma.constitutive_update(copper, F + dF, Q, Float64[], 0.0; need_tangent=false)
            Wm, Pm, _, _ = Norma.constitutive_update(copper, F - dF, Q, Float64[], 0.0; need_tangent=false)
            @test (Wp - Wm) / 2h ≈ P[k, l] rtol = 1.0e-7
            @test (Pp - Pm) / 2h ≈ AA[:, :, k, l] rtol = 1.0e-6 atol = 1.0e-6
        end
        # The closed-form stress and tangent agree with the automatic
        # differentiation of the lattice-frame energy.
        F_L = Q' * F * Q
        W_ad, P_ad = Norma.stress_ad(copper, F_L)
        _, P_L, AA_L = Norma.constitutive(copper, F_L)
        @test P_L ≈ P_ad rtol = 1.0e-12
        @test AA_L ≈ Norma.tangent_ad(copper, F_L) rtol = 1.0e-10
        # At small strain the finite and infinitesimal models agree.
        linear = Norma.create_material(copper_params("anisotropic linear elastic"))
        Fs = ID3 + 1.0e-6 * (F - ID3)
        _, P_s, _, _ = Norma.constitutive_update(copper, Fs, Q, Float64[], 0.0)
        _, σ_s, _, _ = Norma.constitutive_update(linear, Fs, Q, Float64[], 0.0)
        @test P_s ≈ σ_s rtol = 1.0e-5
    end

    @testset "Directional Young's Modulus" begin
        # Uniaxial stress along a lattice direction d of a cubic crystal:
        # 1/E(d) = S11 - 2 (S11 - S12 - S44/2) (d1² d2² + d2² d3² + d3² d1²).
        copper = Norma.create_material(copper_params("anisotropic linear elastic"))
        S = inv(Matrix(Norma.tensor_to_voigt(copper.CC)))
        for d in ([1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
            n = normalize(d)
            expected =
                1 / (
                    S[1, 1] -
                    2 * (S[1, 1] - S[1, 2] - S[4, 4] / 2) * (n[1]^2 * n[2]^2 + n[2]^2 * n[3]^2 + n[3]^2 * n[1]^2)
                )
            # The orientation that turns lattice direction d onto sample x.
            axis = cross(n, [1.0, 0.0, 0.0])
            Q = norm(axis) < 1.0e-12 ? ID3 : axis_rotation(axis, acos(clamp(n[1], -1.0, 1.0)))
            @test Q * n ≈ [1.0, 0.0, 0.0] atol = 1.0e-14
            history = Norma.material_point_history(
                copper,
                s -> @SMatrix([1.0 + 1.0e-4 * s 0 0; 0 1 0; 0 0 1]),
                [0.0, 1.0];
                orientation=Q,
                free=UNIAXIAL_STRESS_FREE,
            )
            σ = history.σ[end]
            @test σ[1, 1] / 1.0e-4 ≈ expected rtol = 1.0e-10
            @test maximum(abs, σ - Diagonal([σ[1, 1], 0.0, 0.0])) < 1.0e-10 * abs(σ[1, 1])
        end
    end

    @testset "Material Point Uniaxial Stress" begin
        # For isotropic St. Venant-Kirchhoff in uniaxial stress, S11 = E E11
        # and the lateral Green strain is -ν E11.
        E, ν = 200.0, 0.3
        svk = Norma.create_material(
            Norma.Parameters("model" => "Saint-Venant Kirchhoff", "elastic modulus" => E, "Poisson's ratio" => ν)
        )
        history = Norma.material_point_history(
            svk, s -> @SMatrix([1.0 + s 0 0; 0 1 0; 0 0 1]), 0.0:0.02:0.1; free=UNIAXIAL_STRESS_FREE
        )
        for (s, F, P) in zip(0.0:0.02:0.1, history.F, history.P)
            E11 = ((1 + s)^2 - 1) / 2
            @test P[1, 1] ≈ (1 + s) * E * E11 atol = 1.0e-12
            @test F[2, 2] ≈ sqrt(1 - 2ν * E11) atol = 1.0e-12
            @test F[3, 3] ≈ F[2, 2] atol = 1.0e-12
        end
    end

    @testset "LCM Simple Shear" begin
        # LCM's SingleElement tests shear a copper crystal to F12 = 0.1. The
        # rotated case has the Bunge angles (0.12, 0.03, 0.28) rad, which LCM
        # receives as basis vectors rounded to eight digits.
        copper = Norma.create_material(copper_params("anisotropic Saint-Venant Kirchhoff"))
        Q = Norma.orientation_from_bunge_angles(0.12, 0.03, 0.28)
        for (case, orientation, tolerance) in
            (("ElasticShear3d", ID3, 1.0e-10), ("ElasticShearRotated3d", Q, 1.0e-7), ("ElasticShear2d", ID3, 1.0e-10))
            rows = lcm_reference(case)
            history = Norma.material_point_history(
                copper, s -> @SMatrix([1.0 s 0; 0 1 0; 0 0 1]), [row["time"] for row in rows]; orientation=orientation
            )
            σ_lcm = [lcm_tensor(row, "Cauchy_Stress") for row in rows]
            # Compare the components LCM writes (the in-plane ones in 2D).
            known = [isfinite(v) for v in σ_lcm[end]]
            scale = maximum(abs, σ_lcm[end][known])
            @test maximum(maximum(abs, (σ - σ_ref)[known]) for (σ, σ_ref) in zip(history.σ, σ_lcm)) < tolerance * scale
        end
        # The consistent rotation, checked against an index-by-index evaluation.
        F = @SMatrix [1.0 0.1 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        history = Norma.material_point_history(copper, s -> F, [0.0]; orientation=Q)
        @test history.σ[end] ≈ cubic_svk_cauchy(168.4, 121.4, 75.4, Q, F) rtol = 1.0e-12
    end
end

# Orientation inputs that were resolved silently or failed with an opaque
# Exodus error: an orientation from the mesh given together with one on the
# material (the mesh won), orientation variables that the truth table leaves
# undefined on a block, and a `time index` that is not a step of the mesh.
@testset "Lattice Orientation Input Checks" begin
    mesh_file = "orientation-checks.g"
    source = Exodus.ExodusDatabase(joinpath("lcm-reference", "meshes", "RubiksCube.g"), "r")
    Exodus.copy(source, mesh_file)
    close(source)
    mesh = Exodus.ExodusDatabase(mesh_file, "rw")
    names = ["phi1", "Phi", "phi2"]
    Exodus.write_number_of_variables(mesh, ElementVariable, 3)
    Exodus.write_names(mesh, ElementVariable, names)
    num_blocks = length(Exodus.read_ids(mesh, Block))
    truth_table = zeros(Cint, 3 * num_blocks)
    truth_table[1:3] .= 1  # the variables are defined on the first block only
    Exodus.LibExodus.ex_put_truth_table(
        Exodus.get_file_id(mesh), Exodus.LibExodus.EX_ELEM_BLOCK, num_blocks, 3, truth_table
    )
    Exodus.write_time(mesh, 1, 0.0)
    first_block, second_block = Exodus.read_ids(mesh, Block)[1:2]
    for (name, value) in zip(names, (30.0, 45.0, 0.0))
        Exodus.write_values(mesh, ElementVariable, 1, first_block, name, [value])
    end
    close(mesh)
    mesh = Exodus.ExodusDatabase(mesh_file, "r")
    from_mesh(extra) = Dict{String,Any}("orientation" => merge(Dict{String,Any}("mesh variables" => names), extra))
    Norma.NORMA_TEST_MODE[] = true
    try
        Q = Norma.block_lattice_orientations(from_mesh(Dict()), mesh, first_block, 1)
        @test Q[1] ≈ Norma.orientation_from_bunge_angles(deg2rad(30.0), deg2rad(45.0), 0.0)
        @test_throws Norma.NormaAbortException Norma.block_lattice_orientations(
            from_mesh(Dict("Bunge Euler angles" => [90.0, 0.0, 0.0])), mesh, first_block, 1
        )
        @test_throws Norma.NormaAbortException Norma.block_lattice_orientations(
            from_mesh(Dict()), mesh, second_block, 1
        )
        for time_index in (0, 2)
            @test_throws Norma.NormaAbortException Norma.block_lattice_orientations(
                from_mesh(Dict("time index" => time_index)), mesh, first_block, 1
            )
        end
    finally
        Norma.NORMA_TEST_MODE[] = false
        close(mesh)
        rm(mesh_file; force=true)
    end
end

# Each symmetry class reads only its own constants: a measured C66 given with
# `symmetry: hexagonal`, or C22 with cubic, was dropped without a word.
@testset "Unused Elastic Constants" begin
    hexagonal = Dict{String,Any}(
        "model" => "anisotropic linear elastic",
        "symmetry" => "hexagonal",
        "C11" => 59400.0,
        "C12" => 25610.0,
        "C13" => 21440.0,
        "C33" => 61600.0,
        "C44" => 16400.0,
    )
    validate(props) = Norma.validate_material!(String[], Dict("blocks" => Dict("b" => "m"), "m" => props), "test")
    @test isempty(validate(hexagonal))
    messages = validate(merge(hexagonal, Dict("C66" => 20000.0)))
    @test length(messages) == 1 && occursin("\"C66\"", messages[1]) && occursin("hexagonal", messages[1])
    cubic = Dict{String,Any}(
        "model" => "crystal plasticity", "symmetry" => "cubic", "C11" => 1.0, "C12" => 0.5, "C44" => 0.4, "C22" => 1.0
    )
    messages = validate(cubic)
    @test any(m -> occursin("\"C22\"", m) && occursin("cubic", m), messages)
    triclinic = Dict{String,Any}(
        "model" => "anisotropic linear elastic", "symmetry" => "triclinic", "C16" => 1.0, "C11" => 2.0
    )
    @test isempty(validate(triclinic))
end

# The adaptivity loop of mesh smoothing accepts or rejects a topology change
# by the energies of the elements. They ignored the lattice orientation of an
# anisotropic material, whose oriented energy the smoothing solve minimizes,
# and a material with internal state gave a MethodError.
@testset "Adaptivity Energy Of Oriented And Stateful Materials" begin
    cp("../examples/ems/cube/cube.g", "adaptivity-cube.g"; force=true)
    base = YAML.load_file("../examples/ems/cube/cube.yaml"; dicttype=Norma.Parameters)
    base["name"] = "adaptivity-cube.yaml"
    base["input mesh file"] = "adaptivity-cube.g"
    base["output mesh file"] = "adaptivity-cube.e"
    function energies(material)
        params = deepcopy(base)
        params["model"]["material"] = Norma.Parameters(
            "blocks" => Norma.Parameters("cube" => "solid"), "solid" => material
        )
        sim = Norma.create_simulation(params)
        try
            model = sim.model
            return Norma.element_energies(model, 1, model.blocks[1].connectivity, model.reference)
        finally
            Exodus.close(sim.params["input_mesh"])
            Exodus.close(sim.params["output_mesh"])
        end
    end
    cubic = Norma.Parameters(
        "model" => "anisotropic Saint-Venant Kirchhoff",
        "symmetry" => "cubic",
        "C11" => 2000.0,
        "C12" => 1000.0,
        "C44" => 1500.0,
        "density" => 1000.0,
    )
    oriented = merge(cubic, Norma.Parameters("orientation" => Dict("Bunge Euler angles" => [30.0, 45.0, 10.0])))
    # The same crystal with its constants rotated into the sample frame.
    Q = Norma.orientation_from_bunge_angles(deg2rad(30.0), deg2rad(45.0), deg2rad(10.0))
    C = Norma.tensor_to_voigt(
        Norma.rotate_tensor4(Q, Norma.voigt_to_tensor(Norma.anisotropic_elastic_constants(cubic)))
    )
    triclinic = Norma.Parameters(
        "model" => "anisotropic Saint-Venant Kirchhoff", "symmetry" => "triclinic", "density" => 1000.0
    )
    for i in 1:6, j in i:6
        triclinic["C$i$j"] = C[i, j]
    end
    try
        rotated = energies(triclinic)
        @test maximum(abs.(energies(oriented) .- rotated)) < 1.0e-10 * maximum(abs, rotated)
        @test maximum(abs.(energies(cubic) .- rotated)) > 1.0e-3 * maximum(abs, rotated)
        j2 = Norma.Parameters(
            "model" => "j2 plasticity",
            "elastic modulus" => 1.0e3,
            "Poisson's ratio" => 0.25,
            "yield stress" => 10.0,
            "hardening modulus" => 1.0,
            "density" => 1000.0,
        )
        Norma.NORMA_TEST_MODE[] = true
        @test_throws Norma.NormaAbortException energies(j2)
    finally
        Norma.NORMA_TEST_MODE[] = false
        rm("adaptivity-cube.g"; force=true)
        rm("adaptivity-cube.e"; force=true)
    end
end
