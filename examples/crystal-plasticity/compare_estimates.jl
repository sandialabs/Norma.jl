# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Runs a polycrystal pulled along x (rubiks-cube.yaml, polycrystal-tension.yaml,
# polycrystal-periodic.yaml) and compares its axial nominal stress with the
# Sachs (iso-stress) and Taylor (iso-strain) estimates of the same grains:
#
#   julia --project=<Norma.jl> compare_estimates.jl <input.yaml> [output.csv]
#
# The input either prescribes the x displacement of node set nsx+ and holds
# nsx-, nsy- and nsz- on rollers, or has a periodic boundary condition that
# prescribes F̄_xx. In the periodic case the strain and stress are F̄_xx - 1
# and P̄_xx, and the Taylor estimate has the free components of the periodic
# condition; on rollers it has F_yy and F_zz free. Run it from the folder of
# the input. The CSV (default: <input>-estimates.csv) has a row per stop: the
# axial strain and the Sachs, finite element and Taylor stresses.
#
# The estimates assume the loading of the examples: an axial strain linear in
# time from zero at t = 0, and, for a periodic cell, zero stress on its free
# components and its other prescribed components of F̄ at the identity. The
# script checks these against the run and stops otherwise.

using LinearAlgebra
using Norma
using Printf
using StaticArrays

function element_volumes(model, b)
    block = model.blocks[b]
    volumes = zeros(block.num_elements)
    for e in 1:block.num_elements
        X = Norma.gather_nodal(model.reference, view(block.connectivity, :, e), block.N)
        for p in 1:block.num_points
            volumes[e] += det(block.dN[:, :, p] * X') * block.weights[p]
        end
    end
    return volumes
end

# One grain per block and lattice orientation, weighted by its volume.
function grains(model)
    weights = Dict{Tuple{Int,SMatrix{3,3,Float64,9}},Float64}()
    for b in eachindex(model.blocks), (e, v) in enumerate(element_volumes(model, b))
        key = (b, model.lattice_orientation[b][e])
        weights[key] = get(weights, key, 0.0) + v
    end
    keys_sorted = collect(keys(weights))
    materials = [model.materials[b] for (b, _) in keys_sorted]
    return materials, [Q for (_, Q) in keys_sorted], [weights[k] for k in keys_sorted]
end

function main(args)
    isempty(args) && error("usage: compare_estimates.jl <input.yaml> [output.csv]")
    input = args[1]
    output = length(args) > 1 ? args[2] : splitext(input)[1] * "-estimates.csv"
    sim = Norma.create_simulation(Norma.load_input(input))
    model = sim.model
    periodic = Norma.periodic_boundary_condition(model)
    free = [false false false; false true false; false false true]
    if periodic === nothing
        loaded = only(
            bc for bc in model.boundary_conditions if
            bc isa Norma.SolidMechanicsDirichletBoundaryCondition && bc.name == "nsx+" && bc.offset == 1
        )
        lo, hi = vec(minimum(model.reference; dims=2)), vec(maximum(model.reference; dims=2))
        length_x = hi[1] - lo[1]
        area = (hi[2] - lo[2]) * (hi[3] - lo[3])
        axial_strain = () -> model.displacement[1, loaded.node_set_node_indices[1]] / length_x
        axial_stress = () -> sum(model.internal_force[3 * (n - 1) + 1] for n in loaded.node_set_node_indices) / area
    else
        free = falses(3, 3)
        free[periodic.components[periodic.stress_controlled]] .= true
        axial_strain = () -> Norma.periodic_deformation_gradient(model)[1, 1] - 1
        axial_stress = () -> Norma.periodic_stress(model)[1, 1]
    end
    times, strains, stresses = [0.0], [0.0], [0.0]
    # Norma.evolve, recording the reaction at every stop.
    Norma.sync_control_time(sim)
    Norma.initialize(sim)
    Norma.initialize_writing(sim)
    Norma.write_stop(sim)
    try
        while true
            Norma.advance_control_time(sim)
            Norma.maybe_apply_swaps!(sim)
            Norma.sync_control_time(sim)
            Norma.advance_control(sim)
            Norma.write_stop(sim)
            push!(times, sim.integrator.time)
            push!(strains, axial_strain())
            push!(stresses, axial_stress())
            Norma.stop_evolve(sim) && break
        end
    finally
        Norma.finalize_writing(sim)
    end
    materials, orientations, weights = grains(model)
    @printf("%d grains\n", length(materials))
    # The estimates assume the loading of the examples (see the header).
    times[1] == 0.0 || error("compare_estimates.jl: the estimates start at t = 0; this input starts at $(times[1]).")
    rate = strains[end] / times[end]
    maximum(abs.(strains .- rate .* times)) <= 1.0e-9 * abs(strains[end]) ||
        error("compare_estimates.jl: the estimates assume an axial strain linear in time from zero; this input is not.")
    if periodic !== nothing
        P̄ = Norma.periodic_stress(model)
        F̄ = Norma.periodic_deformation_gradient(model)
        for (c, stress_controlled) in zip(periodic.components, periodic.stress_controlled)
            c == CartesianIndex(1, 1) && continue
            if stress_controlled
                abs(P̄[c]) <= 1.0e-9 * abs(P̄[1, 1]) || error(
                    "compare_estimates.jl: the estimates assume zero stress on the free components of F̄; P̄$(Tuple(c)) = $(P̄[c]).",
                )
            else
                abs(F̄[c] - (c[1] == c[2])) <= 1.0e-12 || error(
                    "compare_estimates.jl: the estimates assume the prescribed components of F̄ other than xx at the identity; F̄$(Tuple(c)) = $(F̄[c]).",
                )
            end
        end
    end
    # The estimates integrate each grain in steps of at most 5e-5 strain
    # between the stops of the run: the local integration of the flow rule
    # starts from the slip rates of the previous step, and at the onset of flow
    # a step of the run's size can take it beyond the range of the flow rule.
    estimate_times = [0.0]
    stop_index = [1]
    for n in 2:length(times)
        substeps = max(1, ceil(Int, rate * (times[n] - times[n - 1]) / 5.0e-5))
        for k in 1:(substeps - 1)
            push!(estimate_times, times[n - 1] + k * (times[n] - times[n - 1]) / substeps)
        end
        push!(estimate_times, times[n])
        push!(stop_index, length(estimate_times))
    end
    taylor = Norma.taylor_history(
        materials,
        orientations,
        weights,
        t -> SMatrix{3,3,Float64,9}(Diagonal([1.0 + rate * t, 1.0, 1.0])),
        estimate_times;
        free=free,
    )
    sachs = Norma.sachs_history(materials, orientations, weights, t -> rate * t, estimate_times)
    σ_sachs = [sachs.P[i][1, 1] for i in stop_index]
    σ_taylor = [taylor.P[i][1, 1] for i in stop_index]
    open(output, "w") do io
        println(io, "strain,sachs,fem,taylor")
        for n in eachindex(times)
            @printf(io, "%.6e,%.6e,%.6e,%.6e\n", strains[n], σ_sachs[n], stresses[n], σ_taylor[n])
        end
    end
    @printf("%12s %12s %12s %12s\n", "strain", "Sachs", "FEM", "Taylor")
    for n in eachindex(times)
        @printf("%12.4e %12.4e %12.4e %12.4e\n", strains[n], σ_sachs[n], stresses[n], σ_taylor[n])
    end
    return nothing
end

# Run as a script; the tests include the file and call main.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
