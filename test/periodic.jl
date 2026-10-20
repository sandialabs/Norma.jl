# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays
using YAML

const PERIODIC_PAIRS = [["nsx-", "nsx+"], ["nsy-", "nsy+"], ["nsz-", "nsz+"]]

# A periodic condition on the faces of the cube with the components of F̄ in
# `deformation` and the stress-controlled components in `stress`.
function periodic_conditions(deformation::AbstractDict, stress::AbstractDict=Dict{String,Any}())
    return Dict(
        "Periodic" => [
            Dict{String,Any}(
                "node set pairs" => PERIODIC_PAIRS, "deformation gradient" => deformation, "stress" => stress
            ),
        ],
    )
end

# Uniaxial stress along x at the strain rate `rate`: F̄11 prescribed and the
# components in `free` with zero stress. The lower shears stay prescribed, so
# the rigid rotation is fixed.
function periodic_uniaxial_stress(rate; free=("yy", "zz", "xy", "xz", "yz"))
    return periodic_conditions(Dict{String,Any}("xx" => "1.0 + $rate * t"), Dict{String,Any}(c => "0.0" for c in free))
end

# A periodic condition on the pairs `pairs` of PERIODIC_PAIRS, with Dirichlet
# conditions (node set, component, function) in `dirichlet`.
function periodic_cell(pairs, deformation::AbstractDict, stress::AbstractDict=Dict{String,Any}(); dirichlet=())
    conditions = Dict{String,Any}(
        "Periodic" => [
            Dict{String,Any}(
                "node set pairs" => PERIODIC_PAIRS[pairs], "deformation gradient" => deformation, "stress" => stress
            ),
        ],
    )
    if !isempty(dirichlet)
        conditions["Dirichlet"] = [
            Dict{String,Any}("node set" => n, "component" => c, "function" => f) for (n, c, f) in dirichlet
        ]
    end
    return conditions
end

# The largest violation of u(X⁺) - u(X⁻) = (F̄ - I)(X⁺ - X⁻) over the nodes of
# the faces, each paired here with the node of the opposite face nearest to
# its image, relative to the largest displacement.
function periodic_violation(model, mesh)
    H = Norma.periodic_deformation_gradient(model) - I
    X, u = model.reference, model.displacement
    database = Exodus.ExodusDatabase(mesh, "r")
    ids = Exodus.read_ids(database, NodeSet)
    sets = Dict(
        n => Int.(Exodus.read_node_set_nodes(database, id)) for
        (n, id) in zip(Exodus.read_names(database, NodeSet), ids)
    )
    close(database)
    violation = 0.0
    for (axis, (lower, upper)) in enumerate(PERIODIC_PAIRS)
        shift = [j == axis ? maximum(X[axis, :]) - minimum(X[axis, :]) : 0.0 for j in 1:3]
        for a in sets[lower]
            b = argmin(n -> norm(X[:, n] - X[:, a] - shift), sets[upper])
            violation = max(violation, norm(u[:, b] - u[:, a] - H * (X[:, b] - X[:, a])))
        end
    end
    return violation / maximum(abs, u)
end

# The displacement of every node minus that of the origin and the affine field
# (F̄ - I)(X - X_O), relative to the largest displacement.
function fluctuation(model)
    H = Norma.periodic_deformation_gradient(model) - I
    origin = Norma.periodic_boundary_condition(model).origin
    X, u = model.reference, model.displacement
    return maximum(abs, u .- u[:, origin] - H * (X .- X[:, origin])) / maximum(abs, u)
end

# The fluctuation of the first Newton iterate of every solve of `run()`.
function predicted_fluctuations(run)
    predicted = Float64[]
    Norma.SOLVE_ITERATION_OBSERVER[] = function (model, iteration)
        iteration == 0 && push!(predicted, fluctuation(model))
        return nothing
    end
    sim = try
        run()
    finally
        Norma.SOLVE_ITERATION_OBSERVER[] = nothing
    end
    return sim, predicted
end

# A copy of the mesh `source` with its coordinates X replaced by map(X).
function mapped_mesh(source, target, map)
    rm(target; force=true)
    database = Exodus.ExodusDatabase(source, "r")
    Exodus.copy(database, target)
    close(database)
    database = Exodus.ExodusDatabase(target, "rw")
    Exodus.write_coordinates(database, map(Exodus.read_coordinates(database)))
    close(database)
    return target
end

function isotropic_steel()
    return Dict{String,Any}(
        "model" => "linear elastic", "elastic modulus" => 200.0, "Poisson's ratio" => 0.25, "density" => 1.0
    )
end

@testset "Periodic Boundary Conditions" begin
    # A homogeneous cell deforms affinely, and the predictor applies the
    # increment of F̄ to the whole cell: the first iterate is the solution.
    @testset "Homogeneous Cell Deforms Affinely" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        E, ν = 200.0, 0.25
        λ, μ = E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
        ε = 1.0e-3
        deformation = Dict{String,Any}("xx" => "1.0 + $ε * t", "xy" => "$(ε / 2) * t", "zx" => "$(-ε / 3) * t")
        sim = run_polycrystal(
            mesh,
            materials,
            periodic_conditions(deformation);
            time_step=0.5,
            linear_solver="direct",
            exodus_interval=0.5,
        )
        F = Norma.periodic_deformation_gradient(sim.model)
        P = Norma.periodic_stress(sim.model)
        @test F ≈ [1+ε ε/2 0; 0 1 0; -ε/3 0 1] rtol = 1.0e-14
        @test fluctuation(sim.model) < 1.0e-10
        @test periodic_violation(sim.model, mesh) < 1.0e-12
        strain = (F - I + (F - I)') / 2
        σ = λ * tr(strain) * I + 2μ * strain
        @test P ≈ σ rtol = 1.0e-10
        # F̄ and P̄ are global variables of the output, at t = 0, 0.5 and 1.
        output = Exodus.ExodusDatabase("rubiks-cube-2.e", "r")
        @test Exodus.read_names(output, GlobalVariable) == Norma.PERIODIC_GLOBAL_NAMES
        @test Exodus.read_values(output, GlobalVariable, 3) ≈ vcat(vec(permutedims(F)), vec(permutedims(P)))
        close(output)
        rm("rubiks-cube-2.e"; force=true)
        # Uniaxial stress: the free components relax to zero stress.
        sim = run_polycrystal(mesh, materials, periodic_uniaxial_stress(ε); time_step=0.5, linear_solver="direct")
        F = Norma.periodic_deformation_gradient(sim.model)
        P = Norma.periodic_stress(sim.model)
        @test F ≈ Diagonal([1 + ε, 1 - ν * ε, 1 - ν * ε]) rtol = 1.0e-12
        @test P[1, 1] ≈ E * ε rtol = 1.0e-10
        @test maximum(abs, P - Diagonal([P[1, 1], 0, 0])) < 1.0e-10 * P[1, 1]
        @test fluctuation(sim.model) < 1.0e-10
        # Periodic along x and y only, with free faces normal to z: the cell
        # stretches along x with F̄_yy held, and contracts freely along z. The
        # column of F̄ along z, which the condition does not control, is the
        # volume average of F, in the output too.
        conditions = Dict(
            "Periodic" => [
                Dict{String,Any}(
                    "node set pairs" => PERIODIC_PAIRS[1:2],
                    "deformation gradient" => Dict{String,Any}("xx" => "1.0 + $ε * t"),
                ),
            ],
        )
        sim = run_polycrystal(mesh, materials, conditions; linear_solver="direct", exodus_interval=1.0)
        ε_zz = -λ * ε / (λ + 2μ)
        X, u = sim.model.reference, sim.model.displacement
        @test u[3, :] ≈ ε_zz * (X[3, :] .- minimum(X[3, :])) atol = 1.0e-12
        @test Norma.periodic_stress(sim.model) ≈ Diagonal([(λ + 2μ) * ε + λ * ε_zz, λ * (ε + ε_zz), 0]) atol = 1.0e-12
        F = Norma.periodic_deformation_gradient(sim.model)
        @test F ≈ Diagonal([1 + ε, 1, 1 + ε_zz]) atol = 1.0e-12
        output = Exodus.ExodusDatabase("rubiks-cube-2.e", "r")
        @test Exodus.read_values(output, GlobalVariable, 2)[1:9] ≈ vec(permutedims(F)) atol = 1.0e-12
        close(output)
        rm("rubiks-cube-2.e"; force=true)
        rm(mesh; force=true)
    end

    # The lengths L_k and the volume V scale the corner displacements
    # (F̄ - I) L_k e_k, the corner forces V P̄ᵢⱼ / L_k and P̄ = (1/V) Σ f ⊗ X,
    # which the other cells, unit cubes with zero stress targets, do not
    # show. The cell of rubiks-cube-2.g scaled to 2 × 1.5 × 0.5 and moved off
    # the origin, with F̄_xx and F̄_xy prescribed, P̄_yy = 0.1 t, and the other
    # components of the upper triangle under zero stress: the homogeneous cell
    # deforms affinely, every first iterate included, with F̄ and P̄ of the
    # material point.
    @testset "Scaled And Translated Cell With A Stress Target" begin
        mesh = mapped_mesh(
            joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", "rubiks-cube-2.g"),
            "rubiks-cube-2-mapped.g",
            X -> Diagonal([2.0, 1.5, 0.5]) * X .+ [1.0, -2.0, 3.0],
        )
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        E, ν = 200.0, 0.25
        λ, μ = E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
        ε, s = 1.0e-3, 0.1
        deformation = Dict{String,Any}("xx" => "1.0 + $ε * t", "xy" => "$(ε / 2) * t")
        stress = Dict{String,Any}("yy" => "$s * t", "zz" => "0.0", "xz" => "0.0", "yz" => "0.0")
        sim, predicted = predicted_fluctuations(
            () -> run_polycrystal(
                mesh, materials, periodic_conditions(deformation, stress); time_step=0.5, linear_solver="direct"
            ),
        )
        bc = Norma.periodic_boundary_condition(sim.model)
        @test bc.lengths ≈ [2.0, 1.5, 0.5] rtol = 1.0e-14
        @test bc.volume ≈ 1.5 rtol = 1.0e-14
        @test sim.model.reference[:, bc.origin] ≈ [0.0, -2.75, 2.75] atol = 1.0e-14
        H_yy, H_zz = [λ+2μ λ; λ λ+2μ] \ [s - λ * ε, -λ * ε]
        H = [ε ε/2 0; 0 H_yy 0; 0 0 H_zz]
        strain = (H + H') / 2
        @test maximum(abs, Norma.periodic_deformation_gradient(sim.model) - I - H) < 1.0e-12
        @test maximum(abs, Norma.periodic_stress(sim.model) - (λ * tr(strain) * I + 2μ * strain)) < 1.0e-10 * s
        @test length(predicted) == 2
        @test maximum(predicted) < 1.0e-10
        @test fluctuation(sim.model) < 1.0e-10
        @test periodic_violation(sim.model, mesh) < 1.0e-12
        rm(mesh; force=true)
    end

    # The increment of F̄ applied to the whole cell is the prediction of
    # `constant`, and of a step the other predictors do not predict, such as
    # the first step of `linear`, which has no increment to extrapolate. A
    # step that `linear` or `tangent` predicts does not add the increment of
    # F̄ again. On a homogeneous cell every predicted iterate is the affine
    # solution, also when F̄ at the initial time is not I: without the initial
    # equilibrium, the cell starts from the affine field of F̄(t₀). From the
    # corners alone, `tangent` left F̄(t₀) - I on the elements next to the
    # upper faces at the first step, and `linear` extrapolated it into the
    # second.
    @testset "Predictors Of The Time Integrator" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        deformation = Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t", "yz" => "5.0e-4 * t")
        initial = Dict{String,Any}("xx" => "1.002 + 1.0e-3 * t", "yz" => "5.0e-4 * t")
        for gradient in (deformation, initial), predictor in ("constant", "linear", "tangent")
            sim, predicted = predicted_fluctuations(
                () -> run_polycrystal(
                    mesh,
                    materials,
                    periodic_conditions(gradient);
                    time_step=0.25,
                    linear_solver="direct",
                    predictor=predictor,
                ),
            )
            @test length(predicted) == 4
            @test maximum(predicted) < 1.0e-10
            @test fluctuation(sim.model) < 1.0e-10
        end
        params = polycrystal_parameters(mesh, materials, periodic_conditions(initial); time_step=0.25)
        sim = Norma.create_simulation(params)
        try
            Norma.sync_control_time(sim)
            Norma.initialize(sim)
            @test Norma.periodic_deformation_gradient(sim.model)[1, 1] ≈ 1.002 rtol = 1.0e-14
            @test fluctuation(sim.model) < 1.0e-12
            @test periodic_violation(sim.model, mesh) < 1.0e-12
        finally
            Norma.finalize_writing(sim)
        end
        # Before its first step, the tangent predictor evaluates the start of
        # the step: the state with the prescribed degrees of freedom set back
        # to their start values. The evaluation must leave that state as it is,
        # and not move the corners to the new F̄.
        params = polycrystal_parameters(
            mesh, materials, periodic_conditions(deformation); time_step=0.25, predictor="tangent"
        )
        sim = Norma.create_simulation(params)
        try
            Norma.sync_control_time(sim)
            Norma.initialize(sim)
            u = sim.integrator.displacement
            start = copy(u)
            sim.model.time = 0.25
            Norma.apply_bcs(sim.model)
            fixed = .!sim.model.free_dofs
            @test u[fixed] != start[fixed]
            u[fixed] = start[fixed]
            Norma.evaluate(sim.integrator, sim.solver, sim.model)
            @test u == start
        finally
            Norma.finalize_writing(sim)
        end
        rm("rubiks-cube-2.e"; force=true)
        rm(mesh; force=true)
    end

    # The relative tolerance is measured from the residual of the iterate
    # without the increment of F̄, as for a predictor. Measured from the
    # residual of the affine iterate, which is round-off for a homogeneous cell
    # (2e-5 in SI units), a relative tolerance of 1e-10 with an absolute one of
    # 1e-12 could not be met, and the run aborted.
    @testset "Relative Tolerance Of A Homogeneous Cell In SI Units" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        steel = Dict{String,Any}(
            "model" => "linear elastic", "elastic modulus" => 2.0e11, "Poisson's ratio" => 0.25, "density" => 7800.0
        )
        materials = Dict(b => steel for b in mesh_block_names(mesh))
        conditions = periodic_conditions(Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t"))
        for predictor in ("constant", "linear", "tangent")
            params = polycrystal_parameters(
                mesh, materials, conditions; time_step=0.5, linear_solver="direct", predictor=predictor
            )
            merge!(params["solver"], Dict("absolute tolerance" => 1.0e-12, "relative tolerance" => 1.0e-10))
            sim = Norma.run(params)
            rm("rubiks-cube-2.e"; force=true)
            @test sim.integrator.time ≈ 1.0
            @test Norma.UNCONVERGED_SOLVE_COUNT[] == 0
            # The first iterate of each step is the solution, which one Newton
            # iteration confirms.
            @test Norma.NEWTON_ITERATION_COUNT[] == 2
            @test fluctuation(sim.model) < 1.0e-10
        end
        rm(mesh; force=true)
    end

    # A swap hands the converged state to a new simulation, whose periodic
    # condition measures the next increment of F̄ from it. Applied first to
    # the new model at zero displacement, before the state was copied in, the
    # condition counted F̄ - I as an increment, and the first step after the
    # swap predicted it a second time. A time swap of a homogeneous cell to the
    # same input: every first iterate is the affine solution.
    @testset "Prediction After A Swap" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        conditions = periodic_conditions(Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t", "xy" => "5.0e-4 * t"))
        for predictor in ("constant", "linear", "tangent")
            options = (time_step=0.25, linear_solver="direct", predictor=predictor)
            replacement = polycrystal_parameters(mesh, materials, conditions; options...)
            delete!(replacement, "name")
            replacement["output mesh file"] = "rubiks-cube-2-swapped.e"
            YAML.write_file("rubiks-cube-2-swapped.yaml", replacement)
            params = polycrystal_parameters(mesh, materials, conditions; options...)
            params["swaps"] = [
                Dict{String,Any}(
                    "replacement" => "rubiks-cube-2-swapped.yaml",
                    "criterion" => Dict{String,Any}("type" => "time", "t_swap" => 0.45),
                ),
            ]
            sim, predicted = predicted_fluctuations(() -> Norma.run(params))
            rm("rubiks-cube-2.e"; force=true)
            rm("rubiks-cube-2-swapped.e"; force=true)
            @test sim.name == "rubiks-cube-2-swapped"
            @test sim.integrator.time ≈ 1.0
            @test length(predicted) == 4
            @test maximum(predicted) < 1.0e-10
            @test fluctuation(sim.model) < 1.0e-10
        end
        rm("rubiks-cube-2-swapped.yaml"; force=true)
        rm(mesh; force=true)
    end

    # A periodic cell of one crystal deforms homogeneously, so it is a
    # material point: in uniaxial stress, F̄ and P̄ follow those of the
    # material point driver under the same mixed conditions.
    @testset "Single Crystal Cell Is A Material Point" begin
        mesh = "rubiks-cube-1.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        crystal = merge(rubiks_copper(), Dict("orientation" => Dict("Bunge Euler angles" => [20.0, 50.0, 70.0])))
        final_time, time_step = 5.0e-3, 2.5e-4
        sim = run_polycrystal(
            mesh,
            Dict(b => crystal for b in mesh_block_names(mesh)),
            periodic_uniaxial_stress(1.0);
            final_time=final_time,
            time_step=time_step,
            linear_solver="direct",
        )
        @test sim.integrator.time ≈ final_time
        F = Norma.periodic_deformation_gradient(sim.model)
        P = Norma.periodic_stress(sim.model)
        history = Norma.material_point_history(
            sim.model.materials[1],
            τ -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + τ, 1.0, 1.0])),
            collect(range(0.0, final_time; step=time_step));
            orientation=sim.model.lattice_orientation[1][1],
            free=[false true true; false true true; false false true],
        )
        eqps = [state[end] for block in sim.model.state_old for element in block for state in element]
        @test minimum(eqps) > final_time / 2
        @test F ≈ history.F[end] rtol = 1.0e-10
        @test P ≈ history.P[end] atol = 1.0e-7 * abs(P[1, 1])
        @test fluctuation(sim.model) < 1.0e-8
        rm(mesh; force=true)
    end

    # A Dirichlet condition may hold the origin. The corners and the dependent
    # nodes follow it, and since nothing else holds the cell, a moving origin
    # translates the solution with the origin fixed: F̄ and P̄ are the same. On
    # the elastic voxel polycrystal, with F̄_xx and F̄_xy prescribed and F̄_yy
    # under stress.
    @testset "Dirichlet Condition On The Origin" begin
        mesh = "polycrystal-hex8.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "polycrystal", mesh), mesh; force=true)
        elastic = merge(copper_elastic(), Dict("orientation" => mesh_orientation()))
        materials = Dict(b => elastic for b in mesh_block_names(mesh))
        deformation = Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t", "xy" => "5.0e-4 * t")
        stress = Dict{String,Any}("yy" => "0.0")
        fixed = run_polycrystal(mesh, materials, periodic_conditions(deformation, stress))
        moving = periodic_conditions(deformation, stress)
        moving["Dirichlet"] = [
            Dict{String,Any}("node set" => "nsx-y-z-", "component" => c, "function" => "0.01 * t") for
            c in ("x", "y", "z")
        ]
        sim = run_polycrystal(mesh, materials, moving)
        origin = Norma.periodic_boundary_condition(sim.model).origin
        @test sim.model.displacement[:, origin] == [0.01, 0.01, 0.01]
        @test maximum(abs, sim.model.displacement .- 0.01 - fixed.model.displacement) < 1.0e-12
        @test periodic_violation(sim.model, mesh) < 1.0e-12
        @test Norma.periodic_deformation_gradient(sim.model) ≈ Norma.periodic_deformation_gradient(fixed.model)
        @test Norma.periodic_stress(sim.model) ≈ Norma.periodic_stress(fixed.model) rtol = 1.0e-8
        rm(mesh; force=true)
    end

    # The prediction moves the whole cell with an origin that a Dirichlet
    # condition moves, u ← u + Δu_O + ΔF̄ (X - X_O), so on a homogeneous cell
    # the first iterate of every step, a translation and a homogeneous
    # deformation, is the solution. Without Δu_O, only the origin and its
    # images moved with it, and the elements next to them took the jump.
    @testset "Prediction With A Moving Origin" begin
        mesh = "polycrystal-hex8.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "polycrystal", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        conditions = periodic_conditions(Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t", "xy" => "5.0e-4 * t"))
        conditions["Dirichlet"] = [
            Dict{String,Any}("node set" => "nsx-y-z-", "component" => c, "function" => "1.0e-3 * t") for
            c in ("x", "y", "z")
        ]
        for predictor in ("constant", "linear", "tangent")
            sim, predicted = predicted_fluctuations(
                () -> run_polycrystal(
                    mesh, materials, conditions; time_step=0.25, linear_solver="direct", predictor=predictor
                ),
            )
            origin = Norma.periodic_boundary_condition(sim.model).origin
            @test sim.model.displacement[:, origin] ≈ fill(1.0e-3, 3) rtol = 1.0e-14
            @test length(predicted) == 4
            @test maximum(predicted) < 1.0e-10
            @test fluctuation(sim.model) < 1.0e-10
        end
        rm(mesh; force=true)
    end

    # The apparent stiffness of a periodic cell lies between those of the same
    # cell under uniform displacement (KUBC) and uniform traction (SUBC)
    # conditions. In uniaxial strain, the periodic fields contain the affine
    # field of KUBC, so the periodic C11 is below KUBC; in uniaxial stress,
    # the traction conditions admit every periodic field, so the periodic E is
    # above SUBC. Both lie within the Voigt and Reuss bounds.
    @testset "Elastic Periodic Polycrystal Between KUBC And SUBC: $mesh" for mesh in (
        "polycrystal-tet4.g", "polycrystal-hex8.g"
    )
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "polycrystal", mesh), mesh; force=true)
        elastic = merge(copper_elastic(), Dict("orientation" => mesh_orientation()))
        materials = Dict(b => elastic for b in mesh_block_names(mesh))
        ε = 1.0e-3
        # Uniaxial strain: KUBC and periodic.
        sim = run_polycrystal(mesh, materials, affine_uniaxial_strain(ε))
        C_voigt, S_reuss, volume = voigt_reuss(sim.model, sim.model.materials[1].CC)
        C11_kubc = 2 * sim.model.strain_energy / (ε^2 * volume)
        sim = run_polycrystal(mesh, materials, periodic_conditions(Dict{String,Any}("xx" => "1.0 + $ε * t")))
        P = Norma.periodic_stress(sim.model)
        C11_periodic = P[1, 1] / ε
        # Hill-Mandel: the average stress power is the stored energy.
        @test C11_periodic ≈ 2 * sim.model.strain_energy / (ε^2 * volume) rtol = 1.0e-8
        @test periodic_violation(sim.model, mesh) < 1.0e-12
        @test inv(S_reuss)[1, 1] < C11_periodic < C11_kubc < C_voigt[1, 1]
        # The corners carry the volume average of F of the periodic field.
        H = Norma.integrated_displacement_gradient(sim.model) / volume
        @test maximum(abs, H - (Norma.periodic_deformation_gradient(sim.model) - I)) < 1.0e-12 * ε
        # Uniaxial stress: periodic and SUBC. The traction conditions hold the
        # corners against the rigid motions only.
        sim = run_polycrystal(mesh, materials, periodic_uniaxial_stress(ε))
        P = Norma.periodic_stress(sim.model)
        E_periodic = P[1, 1] / ε
        @test maximum(abs, P - Diagonal([P[1, 1], 0, 0])) < 1.0e-8 * P[1, 1]
        s = 0.1
        traction = Dict(
            "Neumann" => [
                Dict{String,Any}("side set" => "ssx+", "component" => "x", "function" => "$s * t"),
                Dict{String,Any}("side set" => "ssx-", "component" => "x", "function" => "$(-s) * t"),
            ],
            "Dirichlet" => [
                Dict{String,Any}("node set" => n, "component" => c, "function" => "0.0") for (n, c) in (
                    ("nsx-y-z-", "x"),
                    ("nsx-y-z-", "y"),
                    ("nsx-y-z-", "z"),
                    ("nsx+y-z-", "y"),
                    ("nsx+y-z-", "z"),
                    ("nsx-y+z-", "z"),
                )
            ],
        )
        sim = run_polycrystal(mesh, materials, traction)
        E_subc = s^2 * volume / (2 * sim.model.strain_energy)
        @test 1 / S_reuss[1, 1] < E_subc < E_periodic < 1 / inv(C_voigt)[1, 1]
        rm(mesh; force=true)
    end

    # The MicroStructPy converter of examples/crystal-plasticity/polycrystal,
    # with quadratic tetrahedra, on a coarse periodic mesh of three grains
    # written by MicroStructPy (test/microstructpy/periodic-polycrystal.xml).
    # For every axis, the i-th node of the upper face set is the image of the
    # i-th node of the lower one, mid-edge nodes included, and a node set holds
    # each corner of the cube. On the TETRA10 cell, a homogeneous material
    # deforms affinely, with the exact P̄, under every predictor, and the
    # elastic polycrystal satisfies Hill-Mandel and lies between its bounds.
    @testset "Periodic TETRA10 Polycrystal From MicroStructPy" begin
        converter = Module()
        Base.include(converter, abspath(POLYCRYSTAL_EXAMPLES, "polycrystal", "microstructpy_to_exodus.jl"))
        mesh = "microstructpy-tet10.g"
        converter.main([joinpath("microstructpy", "periodic-trimesh.txt"), mesh, "--tet10"])
        database = Exodus.ExodusDatabase(mesh, "r")
        X = Exodus.read_coordinates(database)
        block = only(Exodus.read_ids(database, Block))
        element_type, num_elements, num_element_nodes = Exodus.read_block_parameters(database, block)[1:3]
        connectivity = reshape(
            Exodus.read_block_connectivity(database, block, num_elements * num_element_nodes),
            num_element_nodes,
            num_elements,
        )
        sets = Dict(
            n => Int.(Exodus.read_node_set_nodes(database, id)) for
            (n, id) in zip(Exodus.read_names(database, NodeSet), Exodus.read_ids(database, NodeSet))
        )
        close(database)
        @test element_type == "TETRA10"
        midnodes = Set(connectivity[5:10, :])
        for (axis, (lower, upper)) in enumerate(PERIODIC_PAIRS)
            shift = [j == axis ? 1.0 : 0.0 for j in 1:3]
            @test length(sets[lower]) == length(sets[upper])
            @test maximum(norm(X[:, b] - X[:, a] - shift) for (a, b) in zip(sets[lower], sets[upper])) < 1.0e-12
            @test any(in(midnodes), sets[lower])
        end
        for corner in Iterators.product((0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
            name = "ns" * join(("x", "y", "z") .* map(c -> c == 0.0 ? "-" : "+", corner))
            @test X[:, only(sets[name])] == collect(corner)
        end
        E, ν = 200.0, 0.25
        λ, μ = E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
        ε = 1.0e-3
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        deformation = Dict{String,Any}("xx" => "1.0 + $ε * t", "xy" => "$(ε / 2) * t", "zx" => "$(-ε / 3) * t")
        for predictor in ("constant", "linear", "tangent")
            sim = run_polycrystal(
                mesh,
                materials,
                periodic_conditions(deformation);
                time_step=0.5,
                linear_solver="direct",
                predictor=predictor,
            )
            F = Norma.periodic_deformation_gradient(sim.model)
            strain = (F - I + (F - I)') / 2
            @test F ≈ [1+ε ε/2 0; 0 1 0; -ε/3 0 1] rtol = 1.0e-14
            @test fluctuation(sim.model) < 1.0e-10
            @test periodic_violation(sim.model, mesh) < 1.0e-12
            @test Norma.periodic_stress(sim.model) ≈ λ * tr(strain) * I + 2μ * strain rtol = 1.0e-10
        end
        elastic = merge(copper_elastic(), Dict("orientation" => mesh_orientation()))
        materials = Dict(b => elastic for b in mesh_block_names(mesh))
        sim = run_polycrystal(mesh, materials, affine_uniaxial_strain(ε); linear_solver="direct")
        C_voigt, S_reuss, volume = voigt_reuss(sim.model, sim.model.materials[1].CC)
        C11_kubc = 2 * sim.model.strain_energy / (ε^2 * volume)
        stretch = periodic_conditions(Dict{String,Any}("xx" => "1.0 + $ε * t"))
        sim = run_polycrystal(mesh, materials, stretch; linear_solver="direct")
        C11_periodic = Norma.periodic_stress(sim.model)[1, 1] / ε
        @test C11_periodic ≈ 2 * sim.model.strain_energy / (ε^2 * volume) rtol = 1.0e-8
        @test inv(S_reuss)[1, 1] < C11_periodic < C11_kubc < C_voigt[1, 1]
        rm(mesh; force=true)
    end

    # With `linear` and `tangent`, the predictor moves the corners of the
    # stress-controlled components with the rest of the field, and a failed
    # step is predicted again from its start. The RubiksCube of
    # rubiks-cube-2.g in uniaxial stress along x, with its five other
    # components of P̄ zero, at LCM's steps of 5e-4: every predictor halves the
    # step at the onset of flow, [1e-3, 1.5e-3], and reaches the same stress,
    # `linear` and `tangent` in fewer Newton iterations than `constant`. None
    # halves the first step, which `linear` cannot predict and starts with the
    # increment of F̄; without it, the whole increment fell on the elements
    # next to the upper faces and the step failed.
    @testset "Predictors With Stress Control And A Halved Step" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        crystal = merge(rubiks_copper(), Dict("orientation" => mesh_orientation()))
        materials = Dict(b => crystal for b in mesh_block_names(mesh))
        stress = Dict{String,Matrix{Float64}}()
        iterations = Dict{String,Int}()
        for predictor in ("constant", "linear", "tangent")
            times = Float64[]
            Norma.SOLVE_ITERATION_OBSERVER[] = function (model, iteration)
                iteration == 0 && push!(times, model.time)
                return nothing
            end
            sim = try
                run_polycrystal(
                    mesh,
                    materials,
                    periodic_uniaxial_stress(1.0);
                    final_time=2.0e-3,
                    time_step=5.0e-4,
                    linear_solver="direct",
                    adaptive=true,
                    predictor=predictor,
                )
            finally
                Norma.SOLVE_ITERATION_OBSERVER[] = nothing
            end
            @test sim.integrator.time ≈ 2.0e-3
            @test any(t -> isapprox(t, 1.25e-3; rtol=1.0e-9), times)
            @test !any(t -> isapprox(t, 2.5e-4; rtol=1.0e-9), times)
            stress[predictor] = Norma.periodic_stress(sim.model)
            iterations[predictor] = Norma.NEWTON_ITERATION_COUNT[]
            P = stress[predictor]
            @test maximum(abs, P - Diagonal([P[1, 1], 0, 0])) < 1.0e-8 * P[1, 1]
        end
        for predictor in ("linear", "tangent")
            @test stress[predictor] ≈ stress["constant"] rtol = 1.0e-6
            @test iterations[predictor] < iterations["constant"]
        end
        rm(mesh; force=true)
    end

    # The grains of RubiksCube in a periodic cell, in uniaxial stress with
    # the shears held, lie between the Sachs and Taylor estimates, as
    # expected though not guaranteed with plasticity (see polycrystal.jl);
    # the Taylor estimate has the same free components.
    @testset "Sachs And Taylor Bracket The Periodic Plastic Polycrystal" begin
        mesh = "rubiks-cube-1.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        crystal = merge(rubiks_copper(), Dict("orientation" => mesh_orientation()))
        final_time, time_step = 5.0e-3, 2.5e-4
        sim = run_polycrystal(
            mesh,
            Dict(b => crystal for b in mesh_block_names(mesh)),
            periodic_uniaxial_stress(1.0; free=("yy", "zz"));
            final_time=final_time,
            time_step=time_step,
            linear_solver="direct",
            adaptive=true,
        )
        @test sim.integrator.time ≈ final_time
        P = Norma.periodic_stress(sim.model)
        @test abs(P[2, 2]) + abs(P[3, 3]) < 1.0e-8 * P[1, 1]
        @test periodic_violation(sim.model, mesh) < 1.0e-12
        materials = sim.model.materials
        orientations = [only(o) for o in sim.model.lattice_orientation]
        weights = ones(length(materials))
        times = collect(range(0.0, final_time; step=time_step))
        taylor = Norma.taylor_history(
            materials,
            orientations,
            weights,
            τ -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + τ, 1.0, 1.0])),
            times;
            free=[false false false; false true false; false false true],
        )
        sachs = Norma.sachs_history(materials, orientations, weights, τ -> τ, times)
        @test sachs.P[end][1, 1] < P[1, 1] < taylor.P[end][1, 1]
        rm(mesh; force=true)
    end

    # A prescribed component of F̄ in a plane holds the rigid rotation in that
    # plane. A cell with one or two pairs has planes with one such component,
    # or none, so the rotation is free unless another condition holds it, and
    # the tangent is singular. A cell periodic along x only, and one periodic
    # along x and y with F̄_zx under stress, are rejected, and are accepted
    # when a Dirichlet condition on a node away from the origin holds the
    # rotation. The homogeneous cell then deforms without rotation, in
    # uniaxial stress along x, or in plane strain along y with F̄_zx zero.
    @testset "Rigid Rotations Of Cells With One Or Two Pairs" begin
        mesh = "polycrystal-hex8.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "polycrystal", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        E, ν = 200.0, 0.25
        λ, μ = E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
        ε = 1.0e-3
        stretch = Dict{String,Any}("xx" => "1.0 + $ε * t")
        free_zx = Dict{String,Any}("zx" => "0.0")
        Norma.NORMA_TEST_MODE[] = true
        try
            @test_throws Norma.NormaAbortException run_polycrystal(mesh, materials, periodic_cell(1:1, stretch))
            @test_throws Norma.NormaAbortException run_polycrystal(
                mesh, materials, periodic_cell(1:2, stretch, free_zx)
            )
        finally
            Norma.NORMA_TEST_MODE[] = false
        end
        held_x = periodic_cell(1:1, stretch; dirichlet=[("nsx-y+z-", "z", "0.0")])
        sim = run_polycrystal(mesh, materials, held_x; linear_solver="direct")
        X = sim.model.reference .- sim.model.reference[:, Norma.periodic_boundary_condition(sim.model).origin]
        @test maximum(abs, sim.model.displacement - Diagonal([ε, -ν * ε, -ν * ε]) * X) < 1.0e-12
        @test maximum(abs, Norma.periodic_stress(sim.model) - Diagonal([E * ε, 0, 0])) < 1.0e-10
        F = Norma.periodic_deformation_gradient(sim.model)
        @test maximum(abs, F - Diagonal([1 + ε, 1 - ν * ε, 1 - ν * ε])) < 1.0e-12
        held_xz = periodic_cell(1:2, stretch, free_zx; dirichlet=[("nsx-y-z+", "x", "0.0")])
        sim = run_polycrystal(mesh, materials, held_xz; linear_solver="direct")
        ε_zz = -λ * ε / (λ + 2μ)
        @test maximum(abs, sim.model.displacement - Diagonal([ε, 0, ε_zz]) * X) < 1.0e-12
        F = Norma.periodic_deformation_gradient(sim.model)
        @test maximum(abs, F - Diagonal([1 + ε, 1, 1 + ε_zz])) < 1.0e-12
        rm(mesh; force=true)
    end

    # The origin of a cell with one or two pairs is the lowest node, by x,
    # then y, then z, of the face or edge common to the first sets, with
    # coordinates within the matching tolerance equal. Compared exactly, the
    # coordinates along the pair axes, which the common nodes share up to
    # roundoff, picked the origin: x = -1e-16 at a node inside the face x-, or
    # a signed zero at the end of the edge (x-, y-), made that node the
    # origin, and the documented holds acted on the wrong nodes. On
    # polycrystal-hex8.g with those two coordinates, the origin is the corner
    # nsx-y-z-, and the cells of the last test, and a cell periodic along x
    # and y with nsx-y-z- held, take the solutions of the exact mesh.
    @testset "Origin Of Cells With One Or Two Pairs Under Roundoff" begin
        node(X, point) = findfirst(n -> norm(X[:, n] - point) < 1.0e-12, axes(X, 2))
        mesh = mapped_mesh(
            joinpath(POLYCRYSTAL_EXAMPLES, "polycrystal", "polycrystal-hex8.g"),
            "polycrystal-hex8-roundoff.g",
            function (X)
                X = copy(X)
                X[1, node(X, [0.0, 0.5, 0.5])] = -1.0e-16
                X[1, node(X, [0.0, 0.0, 1.0])] = -0.0
                return X
            end,
        )
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        E, ν = 200.0, 0.25
        λ, μ = E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
        ε = 1.0e-3
        ε_zz = -λ * ε / (λ + 2μ)
        stretch = Dict{String,Any}("xx" => "1.0 + $ε * t")
        held_x = periodic_cell(1:1, stretch; dirichlet=[("nsx-y+z-", "z", "0.0")])
        free_zx = Dict{String,Any}("zx" => "0.0")
        held_xz = periodic_cell(1:2, stretch, free_zx; dirichlet=[("nsx-y-z+", "x", "0.0")])
        held_origin = periodic_cell(1:2, stretch; dirichlet=[("nsx-y-z-", c, "0.0") for c in ("x", "y", "z")])
        cells = (
            held_x => Diagonal([ε, -ν * ε, -ν * ε]),
            held_xz => Diagonal([ε, 0, ε_zz]),
            held_origin => Diagonal([ε, 0, ε_zz]),
        )
        for (conditions, H) in cells
            sim = run_polycrystal(mesh, materials, conditions; linear_solver="direct")
            X = sim.model.reference
            @test X[:, Norma.periodic_boundary_condition(sim.model).origin] == zeros(3)
            @test maximum(abs, sim.model.displacement - H * X) < 1.0e-12
            @test maximum(abs, Norma.periodic_deformation_gradient(sim.model) - I - H) < 1.0e-12
            P = Norma.periodic_stress(sim.model)
            @test maximum(abs, P - Diagonal([P[1, 1], P[2, 2], 0])) < 1.0e-10
        end
        rm(mesh; force=true)
    end

    # P̄ is the volume average of the stress of the material, which does not
    # include the forces of the Robin and penalty Surface conditions that
    # evaluate adds to the internal force. A layer periodic along x and y,
    # loaded on its upper face by a traction of 0.1 along -z and held on its
    # lower face by springs, or by a penalty Surface: the material carries
    # P̄_zz = -0.1. With the forces of the springs, P̄_zz was -0.052, and it
    # changed with the position of the mesh along z.
    @testset "Average Stress Without The Forces Of Robin And Surface Conditions" begin
        mesh = "rubiks-cube-2.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        traction = Dict{String,Any}("side set" => "ssz+", "component" => "z", "function" => "-0.1 * t")
        supports = (
            "Robin" => Dict{String,Any}(
                "side set" => "ssz-", "component" => "z", "function" => "0.0", "robin parameter" => 1000.0
            ),
            "Surface" => Dict{String,Any}(
                "side set" => "ssz-", "function" => "z + 0.5", "enforcement" => "penalty", "penalty" => 1.0e4
            ),
        )
        for (type, support) in supports
            conditions = periodic_cell(1:2, Dict{String,Any}())
            conditions["Neumann"] = [traction]
            conditions[type] = [support]
            sim = run_polycrystal(mesh, materials, conditions; linear_solver="direct")
            @test Norma.periodic_stress(sim.model)[3, 3] ≈ -0.1 rtol = 1.0e-8
        end
        rm(mesh; force=true)
    end

    @testset "Input Errors" begin
        mesh = "rubiks-cube-1.g"
        cp(joinpath(POLYCRYSTAL_EXAMPLES, "rubiks-cube", mesh), mesh; force=true)
        materials = Dict(b => isotropic_steel() for b in mesh_block_names(mesh))
        stretch = Dict{String,Any}("xx" => "1.0 + 1.0e-3 * t")
        rejected = [
            # A component both prescribed and stress controlled.
            periodic_conditions(stretch, Dict{String,Any}("xx" => "0.0")),
            # Both shears of a plane free: the rigid rotation is free.
            periodic_conditions(stretch, Dict{String,Any}("xy" => "0.0", "yx" => "0.0")),
            # An expression of position.
            periodic_conditions(Dict{String,Any}("xx" => "1.0 + 1.0e-3 * x")),
            # Node sets that are not images of each other along an axis.
            Dict("Periodic" => [Dict{String,Any}("node set pairs" => [["nsx-", "nsy+"]])]),
        ]
        # A Dirichlet condition on nodes that the periodic condition moves,
        # and what the message says it holds. On a second node set, its nodes
        # and its corner. On a whole first node set, its edges on the second
        # sets of the other pairs and their corners.
        dirichlet = [
            "nsx+" => ["nodes of nsx+", "[nsx-, nsx+]"],
            "nsx-" => ["nodes of nsy+", "nodes of nsz+", "[nsy-, nsy+]", "[nsz-, nsz+]"],
        ]
        Norma.NORMA_TEST_MODE[] = true
        try
            for conditions in rejected
                @test_throws Norma.NormaAbortException run_polycrystal(mesh, materials, conditions)
            end
            for (node_set, held) in dirichlet
                conditions = periodic_cell(1:3, stretch; dirichlet=[(node_set, "y", "0.0")])
                message = try
                    run_polycrystal(mesh, materials, conditions)
                    ""
                catch error
                    error isa Norma.NormaAbortException || rethrow()
                    error.msg
                end
                @test all(text -> occursin(text, message), held)
                @test !occursin("nodes of $node_set", message) || node_set == "nsx+"
            end
            # Matrix output, whose stiffness would not be the condensed one
            # the solve uses.
            params = polycrystal_parameters(mesh, materials, periodic_conditions(stretch))
            params["matrix output"] = Dict{String,Any}("matrices" => ["stiffness"])
            message = try
                Norma.run(params)
                ""
            catch error
                error isa Norma.NormaAbortException || rethrow()
                error.msg
            end
            @test occursin("`matrix output` does not support periodic boundary conditions", message)
        finally
            Norma.NORMA_TEST_MODE[] = false
        end
        rm(mesh; force=true)
        rm("rubiks-cube-1.e"; force=true)
    end
end
