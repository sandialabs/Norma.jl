# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# The composite tet 10 against LCM's (`Use Composite Tet 10: true`), on the
# reference data of lcm-reference/composite-tet (see its README section):
# LCM's ProjectTet10 gold with its nodal projection, its ClampedSDBC1DCt test
# (implicit with the consistent analytic mass, explicit with the lumped one)
# and a Neohookean cylinder of curved elements stretched by 30 %.
using Exodus
using LinearAlgebra
using StaticArrays

const COMPOSITE_TET_REFERENCE = joinpath(@__DIR__, "lcm-reference", "composite-tet")

composite_tet_reference(case) = lcm_reference(joinpath("composite-tet", case))

# Position key of a point, to pair LCM's nodes with Norma's.
composite_tet_key(x) = Tuple(round.(Int, x .* 1.0e9))

function composite_tet_parameters(case, mesh, material, integrator, dirichlet; solver=nothing, ics=nothing)
    cp(joinpath(@__DIR__, "lcm-reference", "meshes", mesh), mesh; force=true)
    params = Norma.Parameters(
        "name" => case * ".yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => case * ".e",
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "solid mechanics",
            "composite tet 10" => true,
            "material" => Dict{String,Any}("blocks" => Dict(material[1] => "solid"), "solid" => material[2]),
        ),
        "time integrator" => integrator,
        "boundary conditions" => Dict(
            "Dirichlet" =>
                [Dict{String,Any}("node set" => n, "component" => c, "function" => f) for (n, c, f) in dirichlet],
        ),
        "solver" => solver === nothing ? Dict{String,Any}("type" => "explicit solver", "step" => "explicit") : solver,
    )
    ics === nothing || (params["initial conditions"] = ics)
    return params
end

function composite_tet_cleanup(case, mesh)
    rm(mesh; force=true)
    rm(case * ".e"; force=true)
    rm(case * ".log"; force=true)
    return nothing
end

# Largest deviation of Norma's nodal displacements u (3 × nodes, Norma's
# numbering) from the LCM rows of one step, paired by position.
function composite_tet_nodal_error(u, reference, rows)
    index = Dict(composite_tet_key(reference[:, n]) => n for n in axes(reference, 2))
    error = 0.0
    for row in rows
        n = index[composite_tet_key([row["X"], row["Y"], row["Z"]])]
        error = max(error, abs(u[1, n] - row["ux"]), abs(u[2, n] - row["uy"]), abs(u[3, n] - row["uz"]))
    end
    return error
end

rows_at_step(rows, step) = filter(row -> row["step"] == step, rows)

const COMPOSITE_TET_DIRECT_SOLVER = Dict{String,Any}(
    "type" => "Hessian minimizer",
    "step" => "full Newton",
    "linear solver" => "direct",
    "minimum iterations" => 1,
    "maximum iterations" => 20,
)

const CLAMPED_CT_DIRICHLET = [
    ("NodeSet0", "x", "0.0"),
    ("NodeSet1", "x", "0.0"),
    ("NodeSet2", "y", "0.0"),
    ("NodeSet3", "y", "0.0"),
    ("NodeSet4", "z", "0.0"),
    ("NodeSet5", "z", "0.0"),
]

# LCM's `Gaussian Z` initial condition with data [0.01, 0.5, 0.02], and its
# linear elastic material; LCM's default cubature degree 3 is the 5-point rule.
const CLAMPED_CT_ICS = Dict(
    "displacement" => [
        Dict{String,Any}(
            "node set" => "NodeSetAll",
            "component" => "z",
            "function" => "a=0.01; b=0.5; c=0.02; a*exp(-(z-b)^2/c^2/2)",
        ),
    ],
)
const CLAMPED_CT_MATERIAL = Dict{String,Any}(
    "model" => "linear elastic", "elastic modulus" => 1.0e9, "Poisson's ratio" => 0.0, "density" => 1.0e3
)

@testset "Composite Tet 10 LCM" begin
    @testset "ProjectTet10" begin
        # LCM's gold: a cube of 6 composite tets stretched by 1 %, 4 points,
        # and the projection of point values onto the nodes that the test is
        # of (`Mass Matrix Type: Full`).
        case, mesh = "composite-tet-project", "cube-single-tet10.g"
        material = Dict{String,Any}(
            "model" => "linear elastic", "elastic modulus" => 1.0, "Poisson's ratio" => 0.25, "density" => 1.0
        )
        integrator = Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0
        )
        dirichlet = [
            ("nodelist_1", "x", "0.0"),
            ("nodelist_3", "y", "0.0"),
            ("nodelist_5", "z", "0.0"),
            ("nodelist_2", "x", "0.01"),
        ]
        solver = merge(
            COMPOSITE_TET_DIRECT_SOLVER, Dict("absolute tolerance" => 1.0e-14, "relative tolerance" => 1.0e-12)
        )
        params = composite_tet_parameters(case, mesh, ("fine", material), integrator, dirichlet; solver=solver)
        params["model"]["num integration points"] = Dict("fine" => 4)
        params["model"]["nodal recovery"] = Dict{String,Any}("method" => "consistent", "stress" => true)
        sim = Norma.create_simulation(params)
        Norma.evolve(sim)
        model = sim.model
        @test model.blocks[1].element_type == Norma.COMPOSITE_TETRA10
        # The solution is affine (uniaxial stress of 0.01 along x on the cube
        # [-0.5, 0.5]³), which Norma reproduces to round-off. LCM's gold is
        # converged to its residual tolerance of 1e-10 only: its stress is
        # 9.1e-11 off the exact one, and so is Norma's from the gold.
        exact = similar(model.reference)
        exact[1, :] = 0.01 * (model.reference[1, :] .+ 0.5)
        exact[2, :] = -0.0025 * (model.reference[2, :] .+ 0.5)
        exact[3, :] = -0.0025 * (model.reference[3, :] .+ 0.5)
        @test maximum(abs.(model.displacement - exact)) < 1.0e-15
        @test all(maximum(abs.(σ - [0.01, 0, 0, 0, 0, 0])) < 1.0e-15 for element in model.stress[1] for σ in element)
        nodal = composite_tet_reference("ProjectTet10-nodal")
        @test composite_tet_nodal_error(model.displacement, model.reference, nodal) < 1.0e-10
        error = 0.0
        for row in composite_tet_reference("ProjectTet10")
            σ = model.stress[1][Int(row["element"])][Int(row["point"])]
            σ_lcm = lcm_tensor(row, "Cauchy_Stress")
            voigt = (σ_lcm[1, 1], σ_lcm[2, 2], σ_lcm[3, 3], σ_lcm[2, 3], σ_lcm[1, 3], σ_lcm[1, 2])
            error = max(error, maximum(abs.(σ .- voigt)))
        end
        @test error < 1.0e-10
        # The consistent recovery gives the uniform stress at every node. The
        # gold's projected stress is not compared: unlike its point stress, it
        # is not uniform.
        @test maximum(abs.(model.recovered_stress .- [0.01, 0, 0, 0, 0, 0])) < 1.0e-14
        # The linear field of LCM's test, 1.5x - 0.8y + 1.1z at the points,
        # comes back at the nodes, as from LCM's projection
        # (`proj_nodal_linear`, which its test checks to 1e-10).
        linear(x) = 1.5 * x[1] - 0.8 * x[2] + 1.1 * x[3]
        block = model.blocks[1]
        N, _, _, _ = Norma.isoparametric(block.element_type, block.num_points)
        for e in 1:(block.num_elements), q in 1:(block.num_points)
            point = model.reference[:, block.connectivity[:, e]] * N[:, q]
            model.stress[1][e][q] .= (linear(point), 0.0, 0.0, 0.0, 0.0, 0.0)
        end
        Norma.recover_stress!(model)
        recovered = model.recovered_stress[1, :]
        @test maximum(abs(recovered[n] - linear(model.reference[:, n])) for n in axes(model.reference, 2)) < 1.0e-14
        index = Dict(composite_tet_key(model.reference[:, n]) => n for n in axes(model.reference, 2))
        lcm = [(index[composite_tet_key([row["X"], row["Y"], row["Z"]])], row["proj_nodal_linear"]) for row in nodal]
        @test maximum(abs(recovered[n] - value) for (n, value) in lcm) < 1.0e-14
        # Close the meshes before removing them (an open file cannot be
        # removed on Windows).
        Norma.finalize_writing(sim)
        composite_tet_cleanup(case, mesh)
    end

    @testset "ClampedSDBC1DCt explicit" begin
        # Central difference with LCM's lumped analytic composite tet mass.
        case, mesh = "composite-tet-clamped-explicit", "clamped-ct.g"
        integrator = Dict{String,Any}(
            "type" => "central difference",
            "initial time" => 0.0,
            "final time" => 1.5e-5,
            "time step" => 1.0e-7,
            "CFL" => 1.0,
            "γ" => 0.5,
        )
        params = composite_tet_parameters(
            case, mesh, ("Block0", CLAMPED_CT_MATERIAL), integrator, CLAMPED_CT_DIRICHLET; ics=CLAMPED_CT_ICS
        )
        params["model"]["num integration points"] = Dict("Block0" => 5)
        params["Exodus output interval"] = 5.0e-6
        sim = Norma.run(params)
        exo = ExodusDatabase(case * ".e", "r")
        times = Exodus.read_times(exo)
        rows = composite_tet_reference("ClampedCt_Explicit-nodal")
        for step in unique(row["step"] for row in rows)
            step_rows = rows_at_step(rows, step)
            k = findfirst(t -> isapprox(t, step_rows[1]["time"]; atol=1.0e-12), times)
            u = vcat((Exodus.read_values(exo, NodalVariable, k, 1, "disp_" * c)' for c in ("x", "y", "z"))...)
            # The amplitude of the initial displacement is 0.01.
            @test composite_tet_nodal_error(u, sim.model.reference, step_rows) < 1.0e-14
        end
        close(exo)
        composite_tet_cleanup(case, mesh)
    end

    @testset "ClampedSDBC1DCt implicit" begin
        # Newmark (β = 1/4) with LCM's consistent analytic composite tet mass,
        # over the 15 steps of LCM's test, which checks the last. LCM's Tempus
        # Newmark d-form starts from a zero acceleration instead of solving
        # M a = f, so this run starts from zero acceleration too.
        case, mesh = "composite-tet-clamped-implicit", "clamped-ct.g"
        integrator = Dict{String,Any}(
            "type" => "Newmark",
            "initial time" => 0.0,
            "final time" => 1.5e-7,
            "time step" => 1.0e-8,
            "β" => 0.25,
            "γ" => 0.5,
        )
        solver = merge(
            COMPOSITE_TET_DIRECT_SOLVER, Dict("absolute tolerance" => 1.0e-6, "relative tolerance" => 1.0e-9)
        )
        params = composite_tet_parameters(
            case,
            mesh,
            ("Block0", CLAMPED_CT_MATERIAL),
            integrator,
            CLAMPED_CT_DIRICHLET;
            solver=solver,
            ics=CLAMPED_CT_ICS,
        )
        params["model"]["num integration points"] = Dict("Block0" => 5)
        params["Exodus output interval"] = 5.0e-8
        sim = Norma.create_simulation(params)
        Norma.sync_control_time(sim)
        Norma.initialize(sim)
        sim.integrator.acceleration .= 0.0
        history = [copy(reshape(sim.integrator.displacement, 3, :))]
        while true
            Norma.advance_control_time(sim)
            Norma.sync_control_time(sim)
            Norma.advance_control(sim)
            push!(history, copy(reshape(sim.integrator.displacement, 3, :)))
            Norma.stop_evolve(sim) && break
        end
        # Close the meshes before removing them (an open file cannot be
        # removed on Windows).
        Norma.finalize_writing(sim)
        composite_tet_cleanup(case, mesh)
        rows = composite_tet_reference("ClampedCt_Implicit-nodal")
        for step in unique(row["step"] for row in rows)
            error = composite_tet_nodal_error(history[Int(step)], sim.model.reference, rows_at_step(rows, step))
            @test error < 1.0e-14
        end
    end

    @testset "Neohookean cylinder" begin
        # Curved composite tets, clamped ends, 30 % stretch in 10 steps, 4 points.
        case, mesh = "composite-tet-cylinder", "cylinder-tet10.g"
        material = Dict{String,Any}(
            "model" => "neohookean", "elastic modulus" => 1000.0, "Poisson's ratio" => 0.3, "density" => 1.0
        )
        integrator = Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 0.1
        )
        dirichlet = [
            ("nodelist_1", "x", "0.0"),
            ("nodelist_1", "y", "0.0"),
            ("nodelist_1", "z", "0.0"),
            ("nodelist_2", "x", "0.0"),
            ("nodelist_2", "y", "0.0"),
            ("nodelist_2", "z", "0.6 * t"),
        ]
        solver = merge(
            COMPOSITE_TET_DIRECT_SOLVER, Dict("absolute tolerance" => 1.0e-11, "relative tolerance" => 1.0e-14)
        )
        params = composite_tet_parameters(case, mesh, ("block_1", material), integrator, dirichlet; solver=solver)
        params["Exodus output interval"] = 0.5
        sim = Norma.run(params)
        model = sim.model
        exo = ExodusDatabase(case * ".e", "r")
        times = Exodus.read_times(exo)
        rows = composite_tet_reference("Cylinder-nodal")
        for step in unique(row["step"] for row in rows)
            step_rows = rows_at_step(rows, step)
            k = findfirst(t -> isapprox(t, step_rows[1]["time"]; atol=1.0e-12), times)
            u = vcat((Exodus.read_values(exo, NodalVariable, k, 1, "disp_" * c)' for c in ("x", "y", "z"))...)
            # The largest displacement is 0.6.
            @test composite_tet_nodal_error(u, model.reference, step_rows) < 1.0e-14
        end
        close(exo)
        composite_tet_cleanup(case, mesh)
        # The deformation gradient and the Cauchy stress at every point of the
        # last step. LCM keeps the element order of this mesh.
        block = model.blocks[1]
        N, dN, _ = Norma.isoparametric(block.element_type, block.num_points)
        σ_error = 0.0
        F_error = 0.0
        σ_max = 0.0
        for row in composite_tet_reference("Cylinder")
            e, q = Int(row["element"]), Int(row["point"])
            nodes = block.connectivity[:, e]
            X = model.reference[:, nodes]
            x = X + model.displacement[:, nodes]
            dNdX = (dN[:, :, q] * X') \ dN[:, :, q]
            F_error = max(F_error, maximum(abs.(x * dNdX' - lcm_tensor(row, "F"))))
            σ_lcm = lcm_tensor(row, "Cauchy_Stress")
            voigt = (σ_lcm[1, 1], σ_lcm[2, 2], σ_lcm[3, 3], σ_lcm[2, 3], σ_lcm[1, 3], σ_lcm[1, 2])
            σ_error = max(σ_error, maximum(abs.(model.stress[1][e][q] .- voigt)))
            σ_max = max(σ_max, maximum(abs, voigt))
        end
        @test F_error < 1.0e-13
        @test σ_error < 1.0e-12 * σ_max
    end
end
