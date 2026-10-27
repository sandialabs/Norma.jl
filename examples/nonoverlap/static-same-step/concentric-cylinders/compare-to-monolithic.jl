# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Runs the quasi-static nonoverlapping Schwarz cases of the quarter of the three
# nested cylinders, quarter/nonoverlap/cylinders-N.yaml for N = 2, 3 or 4
# subdomains, and the monolithic run quarter/monolithic/cylinders.yaml, and
# compares the mean z displacement of the core at t = 1, the measure of test
# 154. For each case it prints whether the run converged, the Schwarz
# iterations of every stop, and the deviation from the monolithic run:
#
#     julia --project=/path/to/Norma.jl compare-to-monolithic.jl      # N = 2, 3, 4
#     julia --project=/path/to/Norma.jl compare-to-monolithic.jl 4    # N = 4 only
#
# The runs are made in a temporary copy of quarter/, so no output is written
# into this directory. A case that aborts is reported and the next one runs.

using Norma
using Printf

function mean_core_displacement_z(model)
    nodes = unique(vec(model.blocks[1].connectivity))
    return sum(model.displacement[3, nodes]) / length(nodes)
end

function run_in(dir::String, input::String)
    return cd(() -> Norma.run(input), dir)
end

function main(args::Vector{String})
    counts = isempty(args) ? [2, 3, 4] : parse.(Int, args)
    work = mktempdir()
    cp(joinpath(@__DIR__, "quarter"), joinpath(work, "quarter"))
    mono = run_in(joinpath(work, "quarter", "monolithic"), "cylinders.yaml")
    reference = mean_core_displacement_z(mono.model)
    results = String[]
    for n in counts
        input = "cylinders-$n.yaml"
        sim = try
            run_in(joinpath(work, "quarter", "nonoverlap"), input)
        catch e
            push!(results, @sprintf("%-17s aborted: %s", input, sprint(showerror, e)))
            continue
        end
        value = mean_core_displacement_z(sim.subsims[1].model)
        deviation = 100 * abs(value - reference) / abs(reference)
        iters = join(sim.controller.schwarz_iters, ", ")
        converged = sim.controller.converged
        line = @sprintf("%-17s converged = %s, Schwarz iterations per stop = [%s], ", input, converged, iters)
        push!(results, line * @sprintf("mean core u_z = %.4e (%.2f%% from monolithic)", value, deviation))
    end
    @printf("cylinders.yaml    monolithic, mean core u_z = %.4e\n", reference)
    foreach(println, results)
    return nothing
end

main(ARGS)
