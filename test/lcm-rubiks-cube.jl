# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays
using YAML

include(joinpath(@__DIR__, "lcm-reference", "lcm_crystal_plasticity.jl"))

# LCM's RubiksCube (tests/LCM/CrystalPlasticity/RubiksCube, implicit): 27
# FCC grains of one HEX8 element each, pulled along x to 0.5 % on rollers.
const RUBIKS_CUBE_DIRICHLET = [
    ("nodelist_1", "x", "0.0"), ("nodelist_2", "x", "1.0 * t"), ("nodelist_3", "y", "0.0"), ("nodelist_5", "z", "0.0")
]

# The time of every LCM step, in order: from <case>-times.csv when the case
# keeps only some of its steps (extract.jl), else from the case itself.
function lcm_step_times(case::AbstractString)
    path = joinpath(@__DIR__, "lcm-reference", case * "-times.csv")
    if isfile(path)
        data, _ = readdlm(path, ',', Float64, '\n'; header=true, comments=true)
        return data[sortperm(data[:, 1]), 2]
    end
    times = Dict{Int,Float64}()
    for row in lcm_reference(case)
        times[Int(row["step"])] = row["time"]
    end
    return [times[step] for step in sort!(collect(keys(times)))]
end

# A Norma run of an LCM crystal plasticity FEM case at LCM's own times: one
# step from each output time to the next, with no cut allowed. Records, at
# every time, the Cauchy stress and the state of every integration point and
# the axial reaction on `loaded`. The slip-rate cap of each family is lifted
# unless the LCM material sets one (Max Slip Increment). LCM has the same cap
# with the same default, log(floatmax), and both codes only use it to give up
# on a local iterate and cut the step; it plays no part in a converged state.
# The two codes' Newton iterates differ, so a step that LCM had to cut is not
# always one that Norma would, and vice versa; without the cap, Norma takes
# every step that LCM took.
function run_lcm_fem_stops(case::String, mesh::String, dirichlet, times::Vector{Float64}, loaded::String)
    lcm = YAML.load_file(lcm_material_file(case))["LCM"]
    material = Dict{String,Any}("blocks" => Dict{String,Any}())
    positions = Dict{Int,Vector{Int}}()
    for (block_name, block) in lcm["ElementBlocks"]
        label = block["material"]
        params, position = norma_crystal_parameters(lcm["Materials"][label])
        for family in params["slip families"]
            haskey(family, "maximum slip rate") || (family["maximum slip rate"] = floatmax(Float64))
        end
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
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}("type" => "solid mechanics", "material" => material),
        "time integrator" => Dict(
            "type" => "quasi static",
            "initial time" => times[1],
            "final time" => times[end],
            "time step" => times[2] - times[1],
            "predictor" => "tangent",
        ),
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
            "relative tolerance" => 1.0e-13,
        ),
    )
    sim = Norma.create_simulation(params)
    model = sim.model
    stresses = [deepcopy(model.stress)]
    states = [deepcopy(model.state_old)]
    reactions = [0.0]
    controller, integrator = sim.controller, sim.integrator
    Norma.sync_control_time(sim)
    Norma.initialize(sim)
    for k in 2:length(times)
        # One stop per LCM step, taken as a single step of its own length.
        controller.prev_time = controller.time
        controller.time = times[k]
        controller.stop = k - 1
        Δt = times[k] - times[k - 1]
        integrator.time_step = integrator.minimum_time_step = integrator.maximum_time_step = Δt
        Norma.sync_control_time(sim)
        Norma.advance_control(sim)
        push!(stresses, deepcopy(model.stress))
        push!(states, deepcopy(model.state_old))
        push!(reactions, axial_reaction(model, loaded))
    end
    # Close the meshes before removing them, as on Windows an open file
    # cannot be removed; the model keeps everything the comparisons need.
    Norma.finalize_writing(sim)
    rm(mesh; force=true)
    rm(case * ".e"; force=true)
    return (; sim, positions, stresses, states, reactions)
end

# Nodal displacements of an LCM case at every written step, keyed by step,
# 3 × num_nodes each, in
# the node order of `reference` (the coordinates of the Norma model). LCM
# writes its output with its own node numbering, so the nodes are matched by
# their coordinates.
function lcm_nodal_displacements(case::AbstractString, reference::AbstractMatrix)
    path = joinpath(@__DIR__, "lcm-reference", case * "-nodal.csv")
    data, header = readdlm(path, ',', Float64, '\n'; header=true, comments=true)
    names = String.(vec(header))
    column(name) = data[:, findfirst(==(name), names)]
    steps = Int.(column("step"))
    first_step = steps .== minimum(steps)
    X = hcat(column("x"), column("y"), column("z"))[first_step, :]
    lcm_node = map(axes(reference, 2)) do n
        return only(findall(i -> maximum(abs, X[i, :] - reference[:, n]) < 1.0e-12, axes(X, 1)))
    end
    u = hcat(column("displacement_x"), column("displacement_y"), column("displacement_z"))
    return Dict(s => permutedims(u[steps .== s, :][lcm_node, :]) for s in unique(steps))
end

# The Norma integration point of each LCM point (Intrepid2's order), found by
# matching the deformation gradient LCM writes at every point with Norma's
# at its points for the same nodal displacements. Checked over every element
# and step; LCM must not average J for the match to be exact.
function lcm_point_order(case::AbstractString, model)
    rows = lcm_reference(case)
    displacements = lcm_nodal_displacements(case, model.reference)
    order = nothing
    worst = 0.0
    for row in rows
        Int(row["step"]) == 1 && continue
        block_index = findfirst(b -> b.id == Int(row["block"]), model.blocks)
        block = model.blocks[block_index]
        nodes = view(block.connectivity, :, Int(row["element"]))
        X = Norma.gather_nodal(model.reference, nodes, block.N)
        u = displacements[Int(row["step"])][:, nodes]
        F_lcm = lcm_tensor(row, "F")
        distances = map(1:block.num_points) do p
            dN = block.dN[:, :, p]
            F = I + u * ((dN * X') \ dN)'
            return maximum(abs, F - F_lcm)
        end
        p = argmin(distances)
        worst = max(worst, distances[p])
        order === nothing && (order = zeros(Int, block.num_points))
        q = Int(row["point"])
        order[q] == 0 || order[q] == p || error("LCM point $q is Norma point $(order[q]) and $p")
        order[q] = p
    end
    return order, worst
end

# Largest normalized deviations from LCM over every step, element and point,
# with where they occur: the Cauchy stress over the largest LCM stress; the
# slips over the largest LCM slip, or the elastic strain (stress over the
# stiffest elastic constant) when larger; Fp - I over the same scale.
function lcm_fem_history_deviations(run, case::AbstractString, order::Vector{Int})
    model = run.sim.model
    rows = lcm_reference(case)
    σ_scale = maximum(maximum(abs, lcm_tensor(row, "Cauchy_Stress")) for row in rows)
    γ_magnitude = maximum(abs(row["gamma_$i"]) for row in rows for i in 1:12 if haskey(row, "gamma_$i"); init=0.0)
    stiffness = maximum(maximum(abs, Norma.tensor_to_voigt(m.CC)) for m in model.materials)
    strain_scale = max(γ_magnitude, σ_scale / stiffness)
    worst = Dict(f => (0.0, 0, 0, 0) for f in ("Cauchy_Stress", "gamma", "Fp"))
    function record!(field, deviation, row)
        deviation > worst[field][1] || return nothing
        worst[field] = (deviation, Int(row["step"]), Int(row["block"]), Int(row["point"]))
        return nothing
    end
    for row in rows
        step, block = Int(row["step"]), Int(row["block"])
        block_index = findfirst(b -> b.id == block, model.blocks)
        element, point = Int(row["element"]), order[Int(row["point"])]
        voigt = run.stresses[step][block_index][element][point]
        σ = @SMatrix [voigt[1] voigt[6] voigt[5]; voigt[6] voigt[2] voigt[4]; voigt[5] voigt[4] voigt[3]]
        record!("Cauchy_Stress", maximum(abs, σ - lcm_tensor(row, "Cauchy_Stress")) / σ_scale, row)
        state = run.states[step][block_index][element][point]
        names = Norma.internal_variable_names(model.materials[block_index])
        position = run.positions[block]
        γ = [state[findfirst(==("gamma_$(position[i])"), names)] - row["gamma_$i"] for i in eachindex(position)]
        record!("gamma", maximum(abs, γ) / strain_scale, row)
        if haskey(row, "Fp_1")
            Fp = SMatrix{3,3,Float64,9}(state[1:9])
            record!("Fp", maximum(abs, Fp - lcm_tensor(row, "Fp")) / strain_scale, row)
        end
    end
    return worst
end

# Axial reaction of LCM on the face x = `position` at every step, from LCM's
# own integration point fields: the x component of the internal force
# f_a = Σ_q J σ F⁻ᵀ ∇N_a w_q summed over the nodes of the face, with the F and
# σ that LCM writes (F̄ when it averages J); NaN at the steps the case does
# not keep. LCM's nodal "force" is the residual, zero on Dirichlet rows, so it
# does not carry the reaction.
function lcm_axial_reaction(case::AbstractString, model, order::Vector{Int}, position::Float64)
    rows = lcm_reference(case)
    num_steps = length(lcm_step_times(case))
    reaction = fill(NaN, num_steps)
    for step in unique(Int(row["step"]) for row in rows)
        reaction[step] = 0.0
    end
    for row in rows
        block_index = findfirst(b -> b.id == Int(row["block"]), model.blocks)
        block = model.blocks[block_index]
        nodes = view(block.connectivity, :, Int(row["element"]))
        X = Norma.gather_nodal(model.reference, nodes, block.N)
        on_face = [abs(model.reference[1, n] - position) < 1.0e-12 for n in nodes]
        any(on_face) || continue
        p = order[Int(row["point"])]
        dN = block.dN[:, :, p]
        J = dN * X'
        F = lcm_tensor(row, "F")
        P = det(F) * lcm_tensor(row, "Cauchy_Stress") * inv(F)'
        f = P * (J \ dN) * (det(J) * block.weights[p])
        reaction[Int(row["step"])] += sum(f[1, on_face])
    end
    return reaction
end

@testset "LCM RubiksCube" begin
    # Without LCM's volume averaging the two codes solve the same discrete
    # problem. Norma takes every LCM step, and LCM's integration point order
    # is recovered from its deformation gradients. Every point is compared
    # at every fourth LCM step and the last, which the case keeps; they
    # include the steps of the largest deviations over all of them. The
    # slips and Fp differ by 1.8e-8 and 2.2e-8 at the end, in the center
    # grain, and the stress by 1.2e-8, at 0.16 % strain; neither tightening
    # LCM's global and local tolerances (to 1e-12 and 1e-14) nor Norma's
    # changes this. The slips start to differ at the onset of flow, where
    # the rate exponent 20 makes them 20 times as sensitive to the resolved
    # shear as the stress. The axial reaction agrees to 2.7e-9, at 0.19 %.
    case = "RubiksCube_NoAverage"
    times = lcm_step_times(case)
    run = run_lcm_fem_stops(case, "RubiksCube.g", RUBIKS_CUBE_DIRICHLET, times, "nodelist_2")
    @test run.sim.integrator.time ≈ times[end]
    order, F_mismatch = lcm_point_order(case, run.sim.model)
    @test sort(order) == 1:8
    @test F_mismatch < 1.0e-12
    worst = lcm_fem_history_deviations(run, case, order)
    @test worst["Cauchy_Stress"][1] < 3.0e-8
    @test worst["gamma"][1] < 5.0e-8
    @test worst["Fp"][1] < 5.0e-8
    lcm_reaction = lcm_axial_reaction(case, run.sim.model, order, 0.5)
    written = findall(isfinite, lcm_reaction)
    @test length(written) == 9
    @test maximum(abs, run.reactions[written] - lcm_reaction[written]) < 1.0e-8 * maximum(abs, lcm_reaction[written])

    # RubiksCube as LCM ships it averages J and the pressure over every
    # element (Weighted Volume Average J, Volume Average Pressure), which
    # Norma does not. Its F̄-based response is softer: at 0.5 % the axial
    # reaction is 4.0 % below Norma's, which matches the unaveraged LCM run.
    # The case keeps only what this reads: the last step of the elements on
    # the loaded face.
    case = "RubiksCube"
    times = lcm_step_times(case)
    averaged = run_lcm_fem_stops(case, "RubiksCube.g", RUBIKS_CUBE_DIRICHLET, times, "nodelist_2")
    @test averaged.sim.integrator.time ≈ times[end]
    lcm_averaged = lcm_axial_reaction(case, averaged.sim.model, order, 0.5)
    @test averaged.reactions[end] ≈ lcm_reaction[end] rtol = 1.0e-4
    @test 0.035 < averaged.reactions[end] / lcm_averaged[end] - 1 < 0.045
end
