# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays

# The composite of LCM's AHD material point tests: half matrix (E = 2e5,
# ν = 0.3, damage saturating at 10 %) and half fibers along x (k = 1e4,
# q = 100, complete damage), plus a second fiber family of zero volume
# fraction.
function lcm_ahd_params()
    fiber(φ) = Dict(
        "direction" => [1.0, 0.0, 0.0],
        "k" => 1.0e4,
        "q" => 100.0,
        "volume fraction" => φ,
        "maximum damage" => 1.0,
        "damage saturation" => 1.0e4,
    )
    return Norma.Parameters(
        "model" => "anisotropic hyperelastic damage",
        "elastic modulus" => 2.0e5,
        "Poisson's ratio" => 0.3,
        "matrix volume fraction" => 0.5,
        "matrix maximum damage" => 0.1,
        "matrix damage saturation" => 1.0e4,
        "fibers" => [fiber(0.5), fiber(0.0)],
    )
end

@testset "Anisotropic Hyperelastic Damage" begin
    identity3 = SMatrix{3,3,Float64,9}(I)

    @testset "LCM Material Point $case" for case in ("AHD-uniaxial", "AHD-shear", "AHD-hydrostatic")
        material = Norma.create_material(lcm_ahd_params())
        rows = lcm_reference(case)
        # The LCM Material Point Simulator writes the step index as the time;
        # the model is rate independent, so the time increment does not matter.
        F_history = [lcm_tensor(row, "F") for row in rows]
        history = Norma.material_point_history(material, s -> F_history[Int(s) + 1], [row["time"] for row in rows])
        σ_lcm = [lcm_tensor(row, "Cauchy_Stress") for row in rows]
        scale = maximum(maximum(abs, σ) for σ in σ_lcm)
        @test maximum(maximum(abs, σ - σ_ref) for (σ, σ_ref) in zip(history.σ, σ_lcm)) < 1.0e-10 * scale
        # Under monotonic loading the energy history is the current energy.
        names = Norma.internal_variable_names(material)
        for (name, column) in (
            ("matrix_energy", "Matrix_Energy"),
            ("fiber_1_energy", "F1_Energy"),
            ("fiber_2_energy", "F2_Energy"),
            ("matrix_damage", "Matrix_Damage"),
            ("fiber_1_damage", "F1_Damage"),
            ("fiber_2_damage", "F2_Damage"),
        )
            k = findfirst(==(name), names)
            values = [state[k] for state in history.state]
            reference = [row[column] for row in rows]
            @test values ≈ reference rtol = 1.0e-10 atol = 1.0e-12 * max(1.0, maximum(abs, reference))
        end
    end

    @testset "Matrix Alone Is Neohookean" begin
        E, ν = 2.0e5, 0.3
        matrix = Norma.create_material(
            Norma.Parameters(
                "model" => "anisotropic hyperelastic damage", "elastic modulus" => E, "Poisson's ratio" => ν
            ),
        )
        neohookean = Norma.create_material(
            Norma.Parameters("model" => "neohookean", "elastic modulus" => E, "Poisson's ratio" => ν)
        )
        F = @SMatrix [1.08 0.05 -0.02; -0.03 0.95 0.04; 0.02 0.01 1.03]
        W, P, AA, _ = Norma.constitutive_update(matrix, F, identity3, Norma.initial_state(matrix), 0.0)
        W_n, P_n, AA_n = Norma.constitutive(neohookean, F)
        @test W ≈ W_n rtol = 1.0e-12
        @test P ≈ P_n rtol = 1.0e-12
        @test AA ≈ AA_n rtol = 1.0e-10
    end

    @testset "Fiber Under Uniaxial Stretch" begin
        k, q = 1.0e4, 100.0
        fiber = Norma.create_material(
            Norma.Parameters(
                "model" => "anisotropic hyperelastic damage",
                "elastic modulus" => 1.0,
                "Poisson's ratio" => 0.0,
                "matrix volume fraction" => 0.0,
                "fibers" => [Dict("direction" => [1.0, 0.0, 0.0], "k" => k, "q" => q, "volume fraction" => 1.0)],
            ),
        )
        for λ in (0.95, 1.0, 1.05, 1.1)
            F = SMatrix{3,3,Float64,9}(Diagonal([λ, 1.0, 1.0]))
            W, P, _, _ = Norma.constitutive_update(fiber, F, identity3, Norma.initial_state(fiber), 0.0)
            e = exp(q * (λ^2 - 1)^2)
            @test W ≈ k / q * (e - 1) atol = 1.0e-12
            @test P[1, 1] ≈ λ * 4k * (λ^2 - 1) * e atol = 1.0e-9
            @test maximum(abs, P - Diagonal([P[1, 1], 0.0, 0.0])) < 1.0e-12
        end
    end

    @testset "Damage Does Not Recover" begin
        # Loading to 10 % and back: every phase keeps the largest energy it
        # reached and the damage it caused, the matrix as well as the fibers.
        # This is the deliberate difference from LCM, whose damage follows
        # the energy of the previous step and recovers on unloading.
        material = Norma.create_material(lcm_ahd_params())
        stretch = s -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + 0.1 * (1 - abs(1 - s)), 1.0, 1.0]))
        history = Norma.material_point_history(material, stretch, 0.0:0.1:2.0)
        names = Norma.internal_variable_names(material)
        @test history.state[11][findfirst(==("fiber_1_damage"), names)] > 0.5
        @test history.state[11][findfirst(==("matrix_damage"), names)] > 0.0
        for (k, name) in enumerate(names)
            values = [state[k] for state in history.state]
            @test issorted(values)
            @test all(values[11:end] .== values[11])
        end
    end

    @testset "Wave Modulus" begin
        # The modulus of the explicit stable step, φ_m (λ + 2μ) + Σ 8 φ_f k_f,
        # is the longitudinal acoustic modulus along the fibers at F = I.
        material = Norma.create_material(lcm_ahd_params())
        modulus = material.φ_m * (material.λ + 2.0 * material.μ) + sum(8.0 * f.φ * f.k for f in material.fibers)
        @test Norma.get_p_wave_modulus(material) ≈ modulus rtol = 1.0e-14
        _, _, AA, _ = Norma.constitutive_update(material, identity3, identity3, Norma.initial_state(material), 0.0)
        @test AA[1, 1, 1, 1] ≈ modulus rtol = 1.0e-10
    end

    @testset "Orientation Rotates The Fibers" begin
        params = lcm_ahd_params()
        material = Norma.create_material(params)
        n = normalize(SVector(1.0, 2.0, -0.5))
        axis = cross(SVector(1.0, 0.0, 0.0), n)
        θ = acos(n[1])
        K = (@SMatrix [0.0 -axis[3] axis[2]; axis[3] 0.0 -axis[1]; -axis[2] axis[1] 0.0]) / norm(axis)
        Q = SMatrix{3,3,Float64,9}(I + sin(θ) * K + (1 - cos(θ)) * K * K)
        @test Q * SVector(1.0, 0.0, 0.0) ≈ n
        rotated = deepcopy(params)
        for fiber in rotated["fibers"]
            fiber["direction"] = collect(n)
        end
        direct = Norma.create_material(rotated)
        F = @SMatrix [1.06 0.03 -0.02; -0.01 0.97 0.05; 0.02 0.01 1.02]
        state = Norma.initial_state(material)
        W1, P1, AA1, s1 = Norma.constitutive_update(material, F, Q, state, 0.0)
        W2, P2, AA2, s2 = Norma.constitutive_update(direct, F, identity3, state, 0.0)
        @test W1 ≈ W2 rtol = 1.0e-12
        @test P1 ≈ P2 rtol = 1.0e-12
        @test AA1 ≈ AA2 rtol = 1.0e-10
        @test s1 ≈ s2 rtol = 1.0e-12
    end

    @testset "Consistent Tangent" begin
        # Central differences of the stress at a fixed start-of-step state,
        # while every phase is damaging.
        material = Norma.create_material(lcm_ahd_params())
        F = @SMatrix [1.06 0.03 -0.02; -0.01 0.97 0.05; 0.02 0.01 1.02]
        state = Norma.initial_state(material)
        _, P, AA, state_new = Norma.constitutive_update(material, F, identity3, state, 0.0)
        @test all(state_new[4:5] .> 0.0)
        h = 1.0e-7
        for k in 1:3, l in 1:3
            dF = SMatrix{3,3,Float64,9}((i, j) == (k, l) ? h : 0.0 for i in 1:3, j in 1:3)
            _, Pp, _, _ = Norma.constitutive_update(material, F + dF, identity3, state, 0.0; need_tangent=false)
            _, Pm, _, _ = Norma.constitutive_update(material, F - dF, identity3, state, 0.0; need_tangent=false)
            @test (Pp - Pm) / 2h ≈ AA[:, :, k, l] rtol = 1.0e-6 atol = 1.0e-6 * maximum(abs, AA)
        end
    end
end

# A restart snapshot holds only the nodal displacement and velocity, so a
# restarted AHD material would lose its damage. The restart is rejected, in a
# single domain and in a Schwarz restart.
# A HEX8 bar along x, nx × ny × nz elements on L × B × B, with node sets on
# its faces x = 0 (nsx-), x = L (nsx+), y = 0 (nsy-) and z = 0 (nsz-).
function write_ahd_bar(path; nx=20, ny=2, nz=2, L=1.0, B=0.1)
    node(i, j, k) = 1 + i + (nx + 1) * (j + (ny + 1) * k)
    coords = zeros(3, (nx + 1) * (ny + 1) * (nz + 1))
    for k in 0:nz, j in 0:ny, i in 0:nx
        coords[:, node(i, j, k)] = [L * i / nx, B * j / ny, B * k / nz]
    end
    connectivity = Int32[]
    for k in 0:(nz - 1), j in 0:(ny - 1), i in 0:(nx - 1)
        append!(
            connectivity,
            [
                node(i, j, k),
                node(i + 1, j, k),
                node(i + 1, j + 1, k),
                node(i, j + 1, k),
                node(i, j, k + 1),
                node(i + 1, j, k + 1),
                node(i + 1, j + 1, k + 1),
                node(i, j + 1, k + 1),
            ],
        )
    end
    num_elements = nx * ny * nz
    init = Exodus.Initialization{Int32}(
        Int32(3), Int32(size(coords, 2)), Int32(num_elements), Int32(1), Int32(4), Int32(0)
    )
    rm(path; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(path, "w", init)
    Exodus.write_coordinates(exo, coords)
    Exodus.write_block(exo, 1, "HEX8", reshape(connectivity, 8, num_elements))
    Exodus.write_name(exo, Block, 1, "bar")
    for (k, (name, axis, value)) in enumerate((("nsx-", 1, 0.0), ("nsx+", 1, L), ("nsy-", 2, 0.0), ("nsz-", 3, 0.0)))
        nodes = Int32.(findall(n -> abs(coords[axis, n] - value) < 1.0e-12, axes(coords, 2)))
        node_set = Exodus.NodeSet(Int32(k), nodes)
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, name)
    end
    Exodus.close(exo)
    return nothing
end

# exp(q (I4 - 1)^2) overflows at a fiber stretch that a Newton trial iterate
# easily reaches. The energy, stress and tangent were then Inf or NaN with no
# error, so the step was not cut, and with the direct solver the run aborted
# on a singular factorization.
@testset "AHD Overflowing Fiber Energy Cuts The Step" begin
    identity3 = SMatrix{3,3,Float64,9}(I)
    material = Norma.create_material(lcm_ahd_params())
    state = Norma.initial_state(material)
    stretched(λ) = SMatrix{3,3,Float64,9}(Diagonal([λ, 1.0, 1.0]))
    W, P, AA, _ = Norma.constitutive_update(material, stretched(1.2), identity3, state, 0.0)
    @test isfinite(W) && all(isfinite, P) && all(isfinite, AA)
    @test_throws DomainError Norma.constitutive_update(material, stretched(7.0), identity3, state, 0.0)
    # A bar pulled by 20 % in one step: the first iterate stretches the last
    # layer of elements 5 times, where the energy overflows, and the
    # adaptive step halves until it converges.
    write_ahd_bar("ahd-bar.g")
    fiber = Dict{String,Any}(
        "direction" => [1.0, 0.0, 0.0],
        "k" => 10.0,
        "q" => 2.0,
        "volume fraction" => 0.5,
        "maximum damage" => 0.5,
        "damage saturation" => 10.0,
    )
    params = Norma.Parameters(
        "type" => "single",
        "name" => "ahd-bar.yaml",
        "input mesh file" => "ahd-bar.g",
        "output mesh file" => "ahd-bar.e",
        "Exodus output interval" => 0.0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "solid mechanics",
            "material" => Dict{String,Any}(
                "blocks" => Dict{String,Any}("bar" => "tissue"),
                "tissue" => Dict{String,Any}(
                    "model" => "anisotropic hyperelastic damage",
                    "elastic modulus" => 100.0,
                    "Poisson's ratio" => 0.3,
                    "density" => 1.0,
                    "matrix volume fraction" => 0.5,
                    "fibers" => [fiber],
                ),
            ),
        ),
        "time integrator" => Dict{String,Any}(
            "type" => "quasi static",
            "initial time" => 0.0,
            "final time" => 1.0,
            "time step" => 1.0,
            "minimum time step" => 1.0e-3,
            "maximum time step" => 1.0,
            "decrease factor" => 0.5,
            "increase factor" => 1.0,
        ),
        "boundary conditions" => Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => "nsx-", "component" => "x", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsy-", "component" => "y", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsz-", "component" => "z", "function" => "0.0"),
                Dict{String,Any}("node set" => "nsx+", "component" => "x", "function" => "0.2 * t"),
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => "direct",
            "minimum iterations" => 1,
            "maximum iterations" => 25,
            "relative tolerance" => 1.0e-10,
            "absolute tolerance" => 1.0e-8,
        ),
    )
    sim = Norma.run(params)
    model = sim.model
    @test sim.integrator.time ≈ 1.0
    @test all(isfinite, model.displacement)
    @test maximum(model.displacement[1, :]) ≈ 0.2
    @test all(all(isfinite, point) for element in model.state[1] for point in element)
    rm("ahd-bar.g"; force=true)
    rm("ahd-bar.e"; force=true)
end

@testset "AHD Restart Is Rejected" begin
    material = lcm_ahd_params()
    material["density"] = 1000.0
    for schwarz in (false, true)
        message = restart_abort_message(material; schwarz=schwarz)
        @test message !== nothing && occursin("`anisotropic hyperelastic damage`", message)
    end
end

# Two fiber families with their own directions and volume fractions. Every
# other test gives the second family a zero volume fraction or the first one's
# direction, so its stress, and any use of a direction or fraction per family,
# went unchecked. The phases do not interact: the stress of the two-family
# material is the sum of the one-family ones less the matrix, and the stress
# of a family is φ (1 - d) F S with S = 4 k (I4 - 1) exp(q (I4 - 1)²) M ⊗ M
# and d = d_max (1 - exp(-W/s)) from its energy W = k/q (exp(q (I4 - 1)²) - 1).
@testset "AHD Two Fiber Families" begin
    identity3 = SMatrix{3,3,Float64,9}(I)
    family(direction, k, φ) = Dict(
        "direction" => direction,
        "k" => k,
        "q" => 50.0,
        "volume fraction" => φ,
        "maximum damage" => 0.8,
        "damage saturation" => 20.0,
    )
    x, yz = [1.0, 0.0, 0.0], [0.0, 1.0, 1.0] / sqrt(2)
    function material(φ_x, φ_yz)
        params = lcm_ahd_params()
        params["fibers"] = [family(x, 1.0e4, φ_x), family(yz, 3.0e4, φ_yz)]
        return Norma.create_material(params)
    end
    F = @SMatrix [1.05 0.02 0.0; 0.01 1.03 0.04; 0.0 0.02 0.98]
    P(m) = Norma.constitutive_update(m, F, identity3, Norma.initial_state(m), 1.0)[2]
    P_both, P_x, P_yz, P_matrix = P(material(0.3, 0.2)),
    P(material(0.3, 0.0)), P(material(0.0, 0.2)),
    P(material(0.0, 0.0))
    @test maximum(abs.(P_both - (P_x + P_yz - P_matrix))) < 1.0e-10 * maximum(abs.(P_both))
    C = F' * F
    M = yz
    I4 = dot(M, C * M)
    W = 3.0e4 / 50.0 * (exp(50.0 * (I4 - 1)^2) - 1)
    d = 0.8 * (1 - exp(-W / 20.0))
    S = 4 * 3.0e4 * (I4 - 1) * exp(50.0 * (I4 - 1)^2) * (M * M')
    @test maximum(abs.((P_yz - P_matrix) - 0.2 * (1 - d) * F * S)) < 1.0e-10 * maximum(abs.(F * S))
    @test 0.1 < d < 0.8
end

# The example of examples/materials/ahd: pulled along its fibers on rollers,
# the cube deforms homogeneously, so every point has the state of the
# material point in uniaxial stress at the same stretch.
@testset "AHD Example Cube" begin
    for file in ("cube.yaml", "cube.g")
        cp(joinpath("../examples/materials/ahd", file), file; force=true)
    end
    sim = Norma.run("cube.yaml")
    for file in ("cube.yaml", "cube.g", "cube.e")
        rm(file; force=true)
    end
    @test sim.integrator.time ≈ 1.0
    material = sim.model.materials[1]
    point = Norma.material_point_history(
        material,
        τ -> SMatrix{3,3,Float64,9}(Diagonal([1.0, 1.0, 1.0 + τ])),
        collect(0.0:0.01:0.08);
        free=[true true true; false true true; false false false],
    )
    expected = point.state[end]
    names = Norma.internal_variable_names(material)
    @test 0.1 < expected[findfirst(==("fiber_1_damage"), names)] < 0.2
    for block in sim.model.state_old, element in block, state in element
        @test state ≈ expected rtol = 1.0e-6
    end
end
