# Summarize the runs of the dynamic stability study.
#
#   julia --project=../.. collect.jl [--runs DIR] [--output FILE]
#
# Reads every case directory under the runs directory (default: runs) and
# writes one row per case to summary.csv (or FILE), then prints a table.
# A case contributes what it has: a run that stopped early still gives its
# energy up to the time it reached.
#
# The energy is the blended total energy of run-energy.csv, which counts
# an overlap region once, normalized by its value at the initial time.  The
# columns are:
#   status        completed, failed, or not run
#   reached       last time in the energy history
#   E/E0@t        energy ratio at each checkpoint of the problem (blank if
#                 not reached)
#   E/E0 max, min extremes of the energy ratio over the run
#   t(1.1)        first time the ratio exceeds 1.1 (blank if never)
#   t(2)          first time the ratio exceeds 2 (blank if never)
#   iterations    mean and maximum Schwarz iterations per stop, and the
#                 number of stops that reached the iteration limit
#   outcome       conserving, dissipating, growing, or failed (see classify)
#   wall time     seconds

include(joinpath(@__DIR__, "matrix.jl"))
using DelimitedFiles
using Printf

# Tolerances of the outcome classes, on the energy ratio.
const GROWTH_THRESHOLD = 1.05
const LOSS_THRESHOLD = 0.95

function read_key_values(file)
    values = Dict{String,String}()
    isfile(file) || return values
    for line in eachline(file)
        key, value = strip.(split(line, ':'; limit=2))
        values[key] = value
    end
    return values
end

function energy_history(dir)
    file = joinpath(dir, "run-energy.csv")
    isfile(file) || return nothing
    data, header = readdlm(file, ','; header=true)
    size(data, 1) == 0 && return nothing
    time = Float64.(data[:, findfirst(==("time"), vec(header))])
    total = Float64.(data[:, findfirst(==("total_energy"), vec(header))])
    return time, total ./ total[1]
end

function schwarz_iterations(dir, limit)
    file = joinpath(dir, "run.log")
    isfile(file) || return (mean=NaN, max=0, at_limit=0)
    counts = Int[]
    for line in eachline(file)
        m = match(r"Performed (\d+) Schwarz Iterations?", line)
        m === nothing || push!(counts, parse(Int, m[1]))
    end
    isempty(counts) && return (mean=NaN, max=0, at_limit=0)
    return (mean=sum(counts) / length(counts), max=maximum(counts), at_limit=count(>=(limit), counts))
end

function first_time_above(time, ratio, threshold)
    i = findfirst(>(threshold), ratio)
    return i === nothing ? NaN : time[i]
end

function ratio_at(time, ratio, t)
    t > time[end] * (1 + 1.0e-9) && return NaN
    return ratio[argmin(abs.(time .- t))]
end

# A run that ends early is failed whatever its energy did before; a
# completed run is growing if its energy ratio ever exceeded the growth
# threshold, dissipating if it ended below the loss threshold, and
# conserving otherwise.
function classify(status, maximum_ratio, final_ratio)
    status != "completed" && return status == "not run" ? "not run" : "failed"
    maximum_ratio > GROWTH_THRESHOLD && return "growing"
    final_ratio < LOSS_THRESHOLD && return "dissipating"
    return "conserving"
end

format(x) = isnan(x) ? "" : @sprintf("%.4g", x)

function summarize(dir)
    info = read_key_values(joinpath(dir, "case.txt"))
    isempty(info) && return nothing
    problem = info["problem"]
    p = problem_parameters(problem)
    status_info = read_key_values(joinpath(dir, "status.txt"))
    status = get(status_info, "status", "not run")
    history = energy_history(dir)
    checkpoints = p.checkpoints
    row = Dict{String,Any}(
        "case" => info["name"], "problem" => problem, "coupling" => info["coupling"], "pair" => info["pair"],
        "level" => info["level"], "ratio" => info["ratio"] == "NaN" ? "" : info["ratio"], "status" => status,
        "wall time" => get(status_info, "wall time", ""),
    )
    if history === nothing
        time, ratio = [0.0], [NaN]
    else
        time, ratio = history
    end
    row["reached"] = format(time[end])
    for (k, t) in enumerate(checkpoints)
        row["E/E0 $k"] = format(ratio_at(time, ratio, t))
    end
    row["E/E0 max"] = format(maximum(ratio))
    row["E/E0 min"] = format(minimum(ratio))
    row["t(1.1)"] = format(first_time_above(time, ratio, 1.1))
    row["t(2)"] = format(first_time_above(time, ratio, 2.0))
    iterations = schwarz_iterations(dir, p.maximum_iterations)
    row["iterations mean"] = format(iterations.mean)
    row["iterations max"] = string(iterations.max)
    row["stops at limit"] = string(iterations.at_limit)
    row["outcome"] = classify(status, maximum(ratio), ratio[end])
    return row
end

const COLUMNS = [
    "case", "problem", "coupling", "pair", "level", "ratio", "status", "outcome", "reached",
    "E/E0 1", "E/E0 2", "E/E0 3", "E/E0 4", "E/E0 max", "E/E0 min", "t(1.1)", "t(2)",
    "iterations mean", "iterations max", "stops at limit", "wall time",
]

function main(args)
    runs = joinpath(@__DIR__, "runs")
    output = joinpath(@__DIR__, "summary.csv")
    i = 1
    while i <= length(args)
        args[i] == "--runs" && (runs = abspath(args[i + 1]))
        args[i] == "--output" && (output = abspath(args[i + 1]))
        i += 2
    end
    rows = filter(!isnothing, [summarize(joinpath(runs, d)) for d in sort(readdir(runs)) if isdir(joinpath(runs, d))])
    open(output, "w") do io
        println(io, join(COLUMNS, ","))
        for row in rows
            println(io, join((row[c] for c in COLUMNS), ","))
        end
    end
    # The checkpoints differ by problem, so each problem gets its own table.
    for problem in ("beam", "cyl")
        subset = filter(r -> r["problem"] == problem, rows)
        isempty(subset) && continue
        checkpoints = problem_parameters(problem).checkpoints
        times = join((@sprintf("%8s", @sprintf("%.3g ms", 1.0e3 * t)) for t in checkpoints), " ")
        println("\n$problem: E/E0 at $times")
        for row in subset
            ratios = join((@sprintf("%8s", row["E/E0 $k"]) for k in 1:4), " ")
            @printf("  %-28s %-11s %s  max %-7s  it %-5s %s\n", row["case"], row["outcome"], ratios, row["E/E0 max"],
                row["iterations mean"], row["status"] == "failed" ? "at " * row["reached"] : "")
        end
    end
    println("\n$(length(rows)) cases; summary written to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
