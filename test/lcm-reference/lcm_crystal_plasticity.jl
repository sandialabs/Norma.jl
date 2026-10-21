# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Replays LCM's crystal plasticity gold cases at the material point: each
# LCM material is translated to Norma's input, and each element (or
# integration point) is driven along the deformation gradient LCM recorded,
# at LCM's own times. Needs helpers.jl.

using LinearAlgebra
using StaticArrays
using YAML

lcm_material_file(case::AbstractString) = joinpath(@__DIR__, "materials", case * ".yaml")

# Exodus block id => block name, from the comment line of the case's CSV.
function lcm_block_names(case::AbstractString)
    for line in eachline(joinpath(@__DIR__, case * ".csv"))
        startswith(line, "#") || break
        startswith(line, "# block names: ") || continue
        pairs = split(line[(length("# block names: ") + 1):end], ", ")
        return Dict(parse(Int, strip(id)) => String(strip(name)) for (id, name) in (split(p, " => ") for p in pairs))
    end
    return Dict{Int,String}()
end

function lowercase_entries(d::AbstractDict)
    return Dict{String,Any}(lowercase(k) => (v isa AbstractString ? lowercase(v) : v) for (k, v) in d)
end

# Norma parameters of an LCM CrystalPlasticity material, and the position in
# Norma's numbering of each LCM slip system (Norma numbers the systems family
# by family). A left-handed LCM basis is the rotation with its third axis
# flipped composed with the reflection of lattice axis 3; the reflection
# leaves cubic, tetragonal and hexagonal elasticity unchanged and is applied
# to the slip vectors instead. LCM inputs that Norma does not support stop
# the translation instead of becoming another material: temperature-dependent
# elasticity, the constrained slip and hardness residual, and the solved slip
# predictor.
function norma_crystal_parameters(lcm::AbstractDict)
    elasticity = lcm["Crystal Elasticity"]
    for key in ("M11", "M12", "M13", "M33", "M44", "Reference Temperature")
        haskey(elasticity, key) &&
            error("norma_crystal_parameters: LCM's temperature-dependent elasticity (`$key`) is not supported.")
    end
    C11, C12, C44 = Float64(elasticity["C11"]), Float64(elasticity["C12"]), Float64(elasticity["C44"])
    C13 = Float64(get(elasticity, "C13", C12))
    C33 = Float64(get(elasticity, "C33", C11))
    params = Norma.Parameters("model" => "crystal plasticity", "density" => Float64(get(lcm, "Density", 1.0)))
    # LCM's tensor (CrystalPlasticityCore_Def.hpp): C11 = C22 = C33 when C33 =
    # C11, with C13 = C23 and C66 = C44, which is cubic only when C13 = C12;
    # otherwise C66 = (C11 - C12)/2, hexagonal.
    if C33 == C11 && C13 == C12
        merge!(params, Dict("symmetry" => "cubic", "C11" => C11, "C12" => C12, "C44" => C44))
    elseif C33 == C11
        merge!(
            params,
            Dict(
                "symmetry" => "tetragonal",
                "C11" => C11,
                "C12" => C12,
                "C13" => C13,
                "C33" => C33,
                "C44" => C44,
                "C66" => C44,
            ),
        )
    else
        merge!(
            params,
            Dict("symmetry" => "hexagonal", "C11" => C11, "C12" => C12, "C13" => C13, "C33" => C33, "C44" => C44),
        )
    end
    B = hcat((normalize(Float64.(elasticity["Basis Vector $k"])) for k in 1:3)...)
    reflect = det(B) < 0.0
    reflect && (B[:, 3] .*= -1.0)
    D = reflect ? Diagonal([1.0, 1.0, -1.0]) : Diagonal(ones(3))
    params["orientation"] = Dict("basis vectors" => [B[:, k] for k in 1:3])
    params["integration scheme"] = lowercase(get(lcm, "Integration Scheme", "Explicit"))
    residual = lowercase(get(lcm, "Residual Type", "Slip"))
    residual in ("slip", "slip hardness") ||
        error("norma_crystal_parameters: LCM's `Residual Type: $(lcm["Residual Type"])` is not supported.")
    params["residual"] = residual == "slip hardness" ? "slip and hardness" : "slip"
    predictor = lowercase(get(lcm, "Slip Predictor", "Rate"))
    predictor in ("rate", "none") ||
        error("norma_crystal_parameters: LCM's `Slip Predictor: $(lcm["Slip Predictor"])` is not supported.")
    params["slip predictor"] = predictor
    num_families = Int(get(lcm, "Number of Slip Families", 1))
    num_systems = Int(get(lcm, "Number of Slip Systems", 0))
    members = [Int[] for _ in 1:num_families]
    for i in 1:num_systems
        push!(members[Int(get(lcm["Slip System $i"], "Slip Family", 0)) + 1], i)
    end
    position = zeros(Int, num_systems)
    families = Dict{String,Any}[]
    for k in 1:num_families
        isempty(members[k]) && continue
        lcm_family = lcm["Slip System Family $(k - 1)"]
        systems = map(members[k]) do i
            system = lcm["Slip System $i"]
            entry = Dict{String,Any}(
                "direction" => D * Float64.(system["Slip Direction"]), "normal" => D * Float64.(system["Slip Normal"])
            )
            haskey(system, "Initial Hardening State") &&
                (entry["initial hardening state"] = system["Initial Hardening State"])
            return entry
        end
        family = Dict{String,Any}("slip systems" => systems)
        family["flow rule"] =
            haskey(lcm_family, "Flow Rule") ? lowercase_entries(lcm_family["Flow Rule"]) : Dict("type" => "none")
        family["hardening law"] = lowercase_entries(lcm_family["Hardening Law"])
        haskey(lcm_family, "Max Slip Increment") && (family["maximum slip rate"] = lcm_family["Max Slip Increment"])
        for i in members[k]
            position[i] = count(>(0), position) + 1
        end
        push!(families, family)
    end
    params["slip families"] = families
    return params, position
end

# Material point histories along LCM's recorded deformation, one per
# (block, element, point) of the case's CSV, with LCM's rows for comparison.
function lcm_crystal_replay(case::AbstractString)
    rows = lcm_reference(case)
    lcm = YAML.load_file(lcm_material_file(case))["LCM"]
    block_names = lcm_block_names(case)
    materials = lcm["Materials"]
    element_blocks = lcm["ElementBlocks"]
    replays = []
    for key in unique((Int(row["block"]), Int(row["element"]), Int(row["point"])) for row in rows)
        block, element, point = key
        point_rows = filter(row -> (Int(row["block"]), Int(row["element"]), Int(row["point"])) == key, rows)
        block_name = haskey(block_names, block) ? block_names[block] : only(keys(element_blocks))
        lcm_material = materials[element_blocks[block_name]["material"]]
        params, position = norma_crystal_parameters(lcm_material)
        material = Norma.create_material(params)
        Q = Norma.lattice_orientation(params["orientation"])
        # The Material Point Simulator writes the step index as the time.
        simulator = get(lcm_material, "Material Point Simulator", nothing)
        scale = simulator === nothing ? 1.0 : Float64(simulator["Step Size"])
        times = [scale * row["time"] for row in point_rows]
        field = haskey(first(point_rows), "F_1") ? "F" : "F_nodal"
        haskey(first(point_rows), "$(field)_1") || error("LCM case $case records no deformation gradient.")
        F_history = [lcm_tensor(row, field) for row in point_rows]
        # 2D cases: the out-of-plane components of F are those of plane strain.
        F_history = [
            SMatrix{3,3,Float64,9}(isfinite(F[i, j]) ? F[i, j] : (i == j ? 1.0 : 0.0) for i in 1:3, j in 1:3) for
            F in F_history
        ]
        history = Norma.material_point_history(
            material, t -> F_history[searchsortedfirst(times, t)], times; orientation=Q
        )
        push!(replays, (; block, element, point, rows=point_rows, history, material, position))
    end
    return replays
end

# Errors of a case against LCM, field by field, over all its points and
# steps. Each error is the largest deviation over a scale that keeps fields
# near zero meaningful: stresses over the largest stress of the case, the
# hardening state over its own largest value, and slips, Fp - I and eqps over
# their largest value or at least the elastic strain of the case (stress over
# the stiffest elastic constant); slip rates over theirs or that strain over
# the duration of the case.
function lcm_crystal_errors(replays)
    deviation = Dict{String,Float64}()
    magnitude = Dict{String,Float64}()
    function record!(field, norma, lcm)
        deviation[field] = max(get(deviation, field, 0.0), maximum(abs, norma - lcm))
        return magnitude[field] = max(get(magnitude, field, 0.0), maximum(abs, lcm))
    end
    row_col(k) = ((k - 1) ÷ 3 + 1, (k - 1) % 3 + 1)
    column_major(k) = (k - 1) % 3 * 3 + (k - 1) ÷ 3 + 1
    δ(k) = k in (1, 5, 9) ? 1.0 : 0.0
    duration = 0.0
    stiffness = 0.0
    for replay in replays
        rows = replay.rows
        columns = keys(first(rows))
        names = Norma.internal_variable_names(replay.material)
        index(name) = findfirst(==(name), names)
        written(field) = [k for k in 1:9 if "$(field)_$k" in columns && isfinite(first(rows)["$(field)_$k"])]
        duration = max(duration, replay.history.time[end] - replay.history.time[1])
        stiffness = max(stiffness, maximum(abs, Norma.tensor_to_voigt(replay.material.CC)))
        for (step, row) in enumerate(rows)
            state = replay.history.state[step]
            ks = written("Cauchy_Stress")
            isempty(ks) || record!(
                "Cauchy_Stress",
                [replay.history.σ[step][row_col(k)...] for k in ks],
                [row["Cauchy_Stress_$k"] for k in ks],
            )
            ks = written("Fp")
            isempty(ks) || record!("Fp", [state[column_major(k)] - δ(k) for k in ks], [row["Fp_$k"] - δ(k) for k in ks])
            for field in ("gamma", "gamma_dot", "tau_hard", "tau")
                systems = [i for i in eachindex(replay.position) if "$(field)_$i" in columns]
                isempty(systems) && continue
                record!(
                    field,
                    [state[index("$(field)_$(replay.position[i])")] for i in systems],
                    [row["$(field)_$i"] for i in systems],
                )
            end
            "eqps" in columns && record!("eqps", [state[index("eqps")]], [row["eqps"]])
        end
    end
    stress = get(magnitude, "Cauchy_Stress", 0.0)
    strain = stress / stiffness
    scale = Dict(
        "Cauchy_Stress" => stress,
        "tau" => stress,
        "tau_hard" => get(magnitude, "tau_hard", 0.0),
        "gamma" => max(get(magnitude, "gamma", 0.0), strain),
        "Fp" => max(get(magnitude, "Fp", 0.0), strain),
        "eqps" => max(get(magnitude, "eqps", 0.0), strain),
        "gamma_dot" => max(get(magnitude, "gamma_dot", 0.0), strain / duration),
    )
    return Dict(field => deviation[field] / scale[field] for field in keys(deviation))
end
