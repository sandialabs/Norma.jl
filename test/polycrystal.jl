# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using Random
using StaticArrays
using Statistics

@testset "Polycrystal" begin
    @testset "Uniform Random Rotations" begin
        rng = Xoshiro(7)
        rotations = [Norma.uniform_random_rotation(rng) for _ in 1:20000]
        @test all(Q' * Q ≈ I && det(Q) ≈ 1 for Q in rotations)
        # Over SO(3) every entry of Q has mean 0 and mean square 1/3.
        @test maximum(abs, mean(rotations)) < 0.02
        @test maximum(abs, mean(Q .^ 2 for Q in rotations) .- 1 / 3) < 0.01
        for Q in rotations[1:100]
            @test Norma.orientation_from_bunge_angles(Norma.bunge_angles(Q)...) ≈ Q atol = 1.0e-12
        end
        @test Norma.orientation_from_bunge_angles(Norma.bunge_angles(SMatrix{3,3,Float64,9}(I))...) ≈ I
    end

    # A displacement-based FEM solution is bounded above by the Voigt average,
    # the energy of the affine field, which the FEM space contains, and below
    # by the Reuss average, the complementary energy of the uniform stress:
    # the exact solution is stiffer than Reuss, and the FEM than the exact
    # solution.
    @testset "Voigt And Reuss Bound The Elastic Polycrystal: $mesh" for (folder, mesh) in (
        ("polycrystal", "polycrystal-tet4.g"), ("polycrystal", "polycrystal-hex8.g"), ("rubiks-cube", "rubiks-cube-2.g")
    )
        cp(joinpath(POLYCRYSTAL_EXAMPLES, folder, mesh), mesh; force=true)
        elastic = merge(copper_elastic(), Dict("orientation" => mesh_orientation()))
        materials = Dict(b => elastic for b in mesh_block_names(mesh))
        ε = 1.0e-3
        # Uniaxial strain on the whole boundary.
        sim = run_polycrystal(mesh, materials, affine_uniaxial_strain(ε))
        C_voigt, S_reuss, volume = voigt_reuss(sim.model, sim.model.materials[1].CC)
        @test volume ≈ 1.0 rtol = 1.0e-10
        C_reuss = inv(S_reuss)
        C11 = 2 * sim.model.strain_energy / (ε^2 * volume)
        @test C_reuss[1, 1] < C11 < C_voigt[1, 1]
        # Uniaxial stress on rollers: the Voigt bound relaxes the two lateral
        # normal strains only, the deformations the rollers admit.
        sim = run_polycrystal(mesh, materials, uniaxial_rollers("$ε * t"))
        E = axial_reaction(sim.model, "nsx+") / ε
        lateral = SVector(2, 3)
        E_voigt = C_voigt[1, 1] - C_voigt[1, lateral]' * (C_voigt[lateral, lateral] \ C_voigt[lateral, 1])
        E_reuss = 1 / S_reuss[1, 1]
        @test E_reuss < E < E_voigt
        rm(mesh; force=true)
    end

    # Equal grains of RubiksCube in uniaxial stress on rollers. The Taylor
    # (iso-strain) estimate stretches every grain alike, by an affine field
    # that the rollers and the FEM space admit, and the Sachs (iso-stress)
    # estimate loads every grain with the same uniaxial stress, which is
    # statically admissible. In the elastic range the two are the Voigt and
    # Reuss bounds; with plasticity they are estimates, which a correct FEM
    # solution is expected to lie between, but no theorem orders them.
    @testset "Sachs And Taylor Bracket The Plastic Polycrystal" begin
        mesh = "rubiks-cube-1.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        crystal = merge(rubiks_copper(), Dict("orientation" => mesh_orientation()))
        final_time, time_step = 5.0e-3, 2.5e-4
        sim = run_polycrystal(
            mesh,
            Dict(b => crystal for b in mesh_block_names(mesh)),
            uniaxial_rollers("1.0 * t");
            final_time=final_time,
            time_step=time_step,
            linear_solver="direct",
            adaptive=true,
        )
        @test sim.integrator.time ≈ final_time
        σ_fem = axial_reaction(sim.model, "nsx+")
        eqps = [state[end] for block in sim.model.state_old for element in block for state in element]
        @test maximum(eqps) > final_time / 2
        materials = sim.model.materials
        orientations = [only(o) for o in sim.model.lattice_orientation]
        weights = ones(length(materials))
        times = collect(range(0.0, final_time; step=time_step))
        taylor = Norma.taylor_history(
            materials,
            orientations,
            weights,
            t -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + t, 1.0, 1.0])),
            times;
            free=[false false false; false true false; false false true],
        )
        sachs = Norma.sachs_history(materials, orientations, weights, t -> t, times)
        @test all(abs(P[2, 2]) + abs(P[3, 3]) < 1.0e-8 * abs(P[1, 1]) for P in taylor.P[2:end])
        @test sachs.F[end][1, 1] - 1 ≈ final_time rtol = 1.0e-8
        C_voigt, S_reuss, _ = voigt_reuss(sim.model, materials[1].CC)
        lateral = SVector(2, 3)
        E_voigt = C_voigt[1, 1] - C_voigt[1, lateral]' * (C_voigt[lateral, lateral] \ C_voigt[lateral, 1])
        E_reuss = 1 / S_reuss[1, 1]
        @test taylor.P[2][1, 1] / time_step ≈ E_voigt rtol = 5.0e-3
        @test sachs.P[2][1, 1] / time_step ≈ E_reuss rtol = 5.0e-3
        σ_taylor = taylor.P[end][1, 1]
        σ_sachs = sachs.P[end][1, 1]
        @test σ_sachs < σ_fem < σ_taylor
        # Plastic flow: both estimates end well below their elastic response.
        @test σ_taylor < 0.5 * E_voigt * final_time
        @test σ_sachs < 0.5 * E_reuss * final_time
        rm(mesh; force=true)
    end

    @testset "Orientations By Block And By Element Agree" begin
        mesh = "rubiks-cube-1.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        database = Exodus.ExodusDatabase(mesh, "r")
        block_names = Exodus.read_names(database, Block)
        block_ids = Exodus.read_ids(database, Block)
        by_block = Dict{String,Any}()
        for (name, id) in zip(block_names, block_ids)
            angles = [
                Float64(only(Exodus.read_values(database, ElementVariable, 1, id, v))) for v in ("phi1", "Phi", "phi2")
            ]
            by_block[name] = merge(rubiks_copper(), Dict("orientation" => Dict("Bunge Euler angles" => angles)))
        end
        close(database)
        crystal = merge(rubiks_copper(), Dict("orientation" => mesh_orientation()))
        dirichlet = uniaxial_rollers("1.0 * t")
        run_two_steps(materials) =
            run_polycrystal(mesh, materials, dirichlet; final_time=5.0e-4, time_step=2.5e-4, linear_solver="direct")
        sim_element = run_two_steps(Dict(b => crystal for b in block_names))
        sim_block = run_two_steps(by_block)
        @test all(
            sim_element.model.lattice_orientation[b] ≈ sim_block.model.lattice_orientation[b] for
            b in eachindex(block_names)
        )
        @test sim_element.model.stress ≈ sim_block.model.stress rtol = 1.0e-10
        @test sim_element.model.state_old ≈ sim_block.model.state_old rtol = 1.0e-10
        rm(mesh; force=true)
    end
end

# A math error of the material in the first evaluation of a driver step
# escaped as a bare DomainError, without the time of the step, while those of
# the Newton corrections aborted with the driver's message.
@testset "Driver Failures Name The Step" begin
    database = Exodus.ExodusDatabase(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", "rubiks-cube-1.g"), "r")
    orientations = map(Exodus.read_ids(database, Block)) do id
        angles = [
            Float64(only(Exodus.read_values(database, ElementVariable, 1, id, v))) for v in ("phi1", "Phi", "phi2")
        ]
        return Norma.orientation_from_bunge_angles(deg2rad.(angles)...)
    end
    close(database)
    materials = [Norma.create_material(Norma.Parameters(rubiks_copper())) for _ in orientations]
    weights = ones(length(materials))
    stretch(t) = SMatrix{3,3,Float64,9}(Diagonal([1.0 + t, 1.0, 1.0]))
    lateral = [false false false; false true false; false false true]
    failure(f) =
        try
            f()
            nothing
        catch error
            error
        end
    Norma.NORMA_TEST_MODE[] = true
    try
        # Steps of 1e-3 take the grains beyond their slip-rate cap at the
        # first evaluation of the step where they start to flow.
        error = failure(
            () -> Norma.taylor_history(
                materials, orientations, weights, stretch, collect(0.0:1.0e-3:5.0e-3); free=lateral
            ),
        )
        @test error isa Norma.NormaAbortException &&
            occursin("Taylor aggregate", error.msg) &&
            occursin("t = ", error.msg)
        shear(t) = @SMatrix [1.0 t 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        error = failure(() -> Norma.material_point_history(materials[1], shear, [0.0, 0.5]))
        @test error isa Norma.NormaAbortException &&
            occursin("Material point", error.msg) &&
            occursin("t = ", error.msg)
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
end

# The tools of examples/crystal-plasticity that no other test runs: the voxel
# path of the MicroStructPy converter with one block per grain and a seed, and
# compare_estimates.jl, which compares a polycrystal run with its Sachs and
# Taylor estimates and checks that the run is loaded as they assume.
@testset "Polycrystal Tools" begin
    converter = Module()
    Base.include(converter, abspath(POLYCRYSTAL_EXAMPLES, "polycrystal", "microstructpy_to_exodus.jl"))
    voxels(seed) = (
        mesh="microstructpy-voxels-$seed.g";
        converter.main([joinpath("microstructpy", "periodic-rastermesh.txt"), mesh, "--blocks", "--seed", "$seed"]);
        mesh
    )
    angles = Dict{Int,Vector{Vector{Float64}}}()
    for seed in (1, 5)
        mesh = voxels(seed)
        database = Exodus.ExodusDatabase(mesh, "r")
        ids = Exodus.read_ids(database, Block)
        names = Exodus.read_names(database, Block)
        @test all(startswith("grain_"), names)
        num_elements = 0
        angles[seed] = Vector{Float64}[]
        for id in ids
            element_type, n = Exodus.read_block_parameters(database, id)[1:2]
            @test element_type == "HEX8"
            num_elements += n
            grain = [Exodus.read_values(database, ElementVariable, 1, id, v) for v in ("phi1", "Phi", "phi2")]
            # One orientation per grain, so per block.
            @test all(all(values .== values[1]) for values in grain)
            push!(angles[seed], [values[1] for values in grain])
        end
        @test num_elements == 64
        # The i-th node of an upper face set is the image of the i-th node of
        # the lower one.
        X = Exodus.read_coordinates(database)
        sets = Dict(zip(Exodus.read_names(database, NodeSet), Exodus.read_sets(database, NodeSet)))
        for (axis, lower, upper) in ((1, "nsx-", "nsx+"), (2, "nsy-", "nsy+"), (3, "nsz-", "nsz+"))
            shift = zeros(3)
            shift[axis] = 1.0
            @test maximum(abs.(X[:, sets[upper].nodes] .- X[:, sets[lower].nodes] .- shift)) < 1.0e-12
        end
        close(database)
        rm(mesh; force=true)
    end
    @test angles[1] != angles[5]

    estimates = Module()
    Base.include(estimates, abspath(POLYCRYSTAL_EXAMPLES, "compare_estimates.jl"))
    cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", "rubiks-cube-1.g"), "rubiks-cube-1.g"; force=true)
    params = YAML.load_file(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", "rubiks-cube.yaml"))
    params["input mesh file"] = "rubiks-cube-1.g"
    params["output mesh file"] = "rubiks-cube-estimates.e"
    params["time integrator"]["final time"] = 2.0e-3
    YAML.write_file("rubiks-cube-estimates.yaml", params)
    try
        estimates.main(["rubiks-cube-estimates.yaml", "rubiks-cube-estimates.csv"])
        table, _ = readdlm("rubiks-cube-estimates.csv", ',', Float64, '\n'; header=true)
        @test size(table, 1) == 5
        # Sachs below the finite element response, below Taylor, at every stop
        # after the first, with the strain at the stop.
        @test all(table[2:end, 2] .< table[2:end, 3] .< table[2:end, 4])
        @test table[end, 1] ≈ 2.0e-3 rtol = 1.0e-12
        # A load that is not linear in time from zero is not what the
        # estimates assume.
        params["boundary conditions"]["Dirichlet"][end]["function"] = "500.0 * t^2"
        YAML.write_file("rubiks-cube-estimates.yaml", params)
        @test_throws ErrorException estimates.main(["rubiks-cube-estimates.yaml", "rubiks-cube-estimates.csv"])
    finally
        for file in
            ("rubiks-cube-1.g", "rubiks-cube-estimates.yaml", "rubiks-cube-estimates.e", "rubiks-cube-estimates.csv")
            rm(file; force=true)
        end
    end
end
