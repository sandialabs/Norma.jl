# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays
using YAML

include(joinpath(@__DIR__, "lcm-reference", "lcm_crystal_plasticity.jl"))

# A Norma run of an LCM crystal plasticity FEM case, on the case's mesh with
# its materials, its Dirichlet conditions and its fixed step. Returns the
# simulation and, for each block id, the position in Norma's numbering of
# each LCM slip system.
function run_lcm_fem_case(case::String, mesh::String, dirichlet, time_step::Float64, final_time::Float64)
    lcm = YAML.load_file(lcm_material_file(case))["LCM"]
    material = Dict{String,Any}("blocks" => Dict{String,Any}())
    positions = Dict{Int,Vector{Int}}()
    for (block_name, block) in lcm["ElementBlocks"]
        label = block["material"]
        params, position = norma_crystal_parameters(lcm["Materials"][label])
        material["blocks"][block_name] = label
        material[label] = params
        positions[parse(Int, split(block_name, "_")[end])] = position
    end
    cp(joinpath(@__DIR__, "lcm-reference", "meshes", mesh), mesh; force=true)
    params = Norma.Parameters(
        "name" => case * ".yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => case * ".e",
        "Exodus output interval" => final_time,
        "CSV output interval" => 0,
        "model" => Dict("type" => "solid mechanics", "material" => material),
        "time integrator" =>
            Dict("type" => "quasi static", "initial time" => 0.0, "final time" => final_time, "time step" => time_step),
        "boundary conditions" => Dict(
            "Dirichlet" =>
                [Dict{String,Any}("node set" => n, "component" => c, "function" => f) for (n, c, f) in dirichlet],
        ),
        "solver" => Dict(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => "direct",
            "minimum iterations" => 1,
            "maximum iterations" => 20,
            "absolute tolerance" => 1.0e-10,
            "relative tolerance" => 1.0e-14,
        ),
    )
    sim = Norma.run(params)
    rm(mesh; force=true)
    rm(case * ".e"; force=true)
    return sim, positions
end

# Largest deviations at the last step, over every integration point, of the
# Cauchy stress (relative to the largest LCM stress) and of the slips
# (relative to the largest LCM slip, or the elastic strain when larger).
function lcm_fem_deviations(sim, positions, case::String)
    model = sim.model
    rows = lcm_final_rows(case)
    σ_scale = maximum(maximum(abs, lcm_tensor(row, "Cauchy_Stress")) for row in rows)
    σ_deviation = 0.0
    γ_deviation = 0.0
    γ_scale = 0.0
    for row in rows
        block = Int(row["block"])
        block_index = findfirst(b -> b.id == block, model.blocks)
        element = Int(row["element"])
        material = model.materials[block_index]
        names = Norma.internal_variable_names(material)
        σ_lcm = lcm_tensor(row, "Cauchy_Stress")
        position = positions[block]
        for (point, voigt) in enumerate(model.stress[block_index][element])
            σ = @SMatrix [voigt[1] voigt[6] voigt[5]; voigt[6] voigt[2] voigt[4]; voigt[5] voigt[4] voigt[3]]
            σ_deviation = max(σ_deviation, maximum(abs, σ - σ_lcm))
            state = model.state_old[block_index][element][point]
            for i in eachindex(position)
                haskey(row, "gamma_$i") || continue
                γ = state[findfirst(==("gamma_$(position[i])"), names)]
                γ_deviation = max(γ_deviation, abs(γ - row["gamma_$i"]))
                γ_scale = max(γ_scale, abs(row["gamma_$i"]))
            end
        end
    end
    stiffness = maximum(maximum(abs, Norma.tensor_to_voigt(m.CC)) for m in model.materials)
    return σ_deviation / σ_scale, γ_deviation / max(γ_scale, σ_scale / stiffness)
end

# The two single-element blocks of LCM's SingleSlip mesh, pulled along x at
# 0.3 per unit time: x = 0 faces fixed along x, and one edge of each fixed
# along y and one along z.
const SINGLE_SLIP_UNIAXIAL = [
    ("nodelist_1", "x", "0.0"),
    ("nodelist_2", "x", "0.3 * t"),
    ("nodelist_3", "y", "0.0"),
    ("nodelist_4", "z", "0.0"),
    ("nodelist_11", "x", "0.0"),
    ("nodelist_12", "x", "0.3 * t"),
    ("nodelist_13", "y", "0.0"),
    ("nodelist_14", "z", "0.0"),
]

# The same blocks with both x faces fixed laterally: uniaxial strain.
const SINGLE_SLIP_LATERALLY_FIXED = vcat(
    [(n, c, "0.0") for n in ("nodelist_1", "nodelist_11") for c in ("x", "y", "z")],
    [(n, c, "0.0") for n in ("nodelist_2", "nodelist_12") for c in ("y", "z")],
    [("nodelist_2", "x", "1.0 * t"), ("nodelist_12", "x", "1.0 * t")],
)

# LCM's MultiSlipPlaneHard element, pulled along x at 0.5 per unit time.
const MULTI_SLIP_UNIAXIAL = [
    ("nodelist_11", "x", "0.0"),
    ("nodelist_12", "x", "0.5 * t"),
    ("nodelist_13", "y", "0.0"),
    ("nodelist_14", "z", "0.0"),
]

@testset "Single Crystal LCM Cases" begin
    @testset "$case" for (case, mesh, dirichlet, time_step, final_time, tolerance) in (
        # LCM's first step of this case is not in equilibrium: only F11 moves,
        # and the free lateral faces carry stresses of a third of σ11. The
        # difference decays with the viscous relaxation, to 8.2e-6 at the end.
        ("SingleSlip_Explicit", "SingleSlip.g", SINGLE_SLIP_UNIAXIAL, 1.0e-3, 0.05, 5.0e-5),
        ("SingleSlip_Implicit", "SingleSlip.g", SINGLE_SLIP_LATERALLY_FIXED, 2.0e-4, 0.01, 1.0e-8),
        # LCM writes only every tenth step of this case, which the material
        # point replay cannot follow; the FEM takes LCM's own steps.
        ("SingleSlipHard_Explicit", "SingleSlip.g", SINGLE_SLIP_UNIAXIAL, 2.5e-4, 0.05, 1.0e-8),
        ("MultiSlipPlaneHard_Implicit_Active_Sets", "MultiSlipPlaneHard.g", MULTI_SLIP_UNIAXIAL, 1.0e-3, 0.03, 1.0e-8),
    )
        sim, positions = run_lcm_fem_case(case, mesh, dirichlet, time_step, final_time)
        @test sim.integrator.time ≈ final_time
        σ_error, γ_error = lcm_fem_deviations(sim, positions, case)
        @test σ_error < tolerance
        @test γ_error < tolerance
    end
end
