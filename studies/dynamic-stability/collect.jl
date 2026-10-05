# Summarize the runs of the dynamic stability study.
#
#   julia --project=../.. collect.jl [--runs DIR] [--output FILE]
#
# Reads every case directory under the runs directory (default: runs) and
# writes one row per case to summary.csv (or FILE), then prints a table.
# A case contributes what it has: a run that stopped early still gives its
# energy up to the time it reached.
#
# The energy E is the blended total energy of run-energy.csv (kinetic plus
# strain; for central difference the kinetic energy is the staggered form
# that the scheme conserves), normalized by its value at the initial time.
# The pseudo-energy is the energy functional applied to the velocity and the
# acceleration, the quantity that the velocity constraint conserves.  The
# columns are:
#   status            completed, failed, or not run
#   outcome           conserving to roundoff, conserving, dissipating,
#                     growing, or failed (see classify)
#   reached           last time in the energy history (s)
#   E/E0@t            energy ratio at each checkpoint (blank if not reached)
#   E/E0 max, min     extremes of the energy ratio over the run
#   dE                max |E/E0 - 1| over the run
#   dpseudo           max |pseudo/pseudo0 - 1| over the run (blank for a
#                     monolithic run without the column)
#   jump u, jump v    root mean square interface displacement and velocity
#                     jump at the end of the run (m, m/s; blank if none)
#   impulse           largest relative interface impulse residual (subcycled
#                     constrained cases only)
#   iterations        mean and maximum Schwarz iterations per stop, and the
#                     number of stops that reached the iteration limit
#   wall time         seconds

include(joinpath(@__DIR__, "matrix.jl"))
using DelimitedFiles
using Printf

# Tolerances of the outcome classes, on the energy ratio.
const GROWTH_THRESHOLD = 1.05
const LOSS_THRESHOLD = 0.95
# A completed run whose ratio stays within this distance of 1 over the whole
# run conserves to roundoff and Schwarz tolerance.
const ROUNDOFF_THRESHOLD = 1.0e-8

function read_key_values(file)
    values = Dict{String,String}()
    isfile(file) || return values
    for line in eachline(file)
        occursin(':', line) || continue
        key, value = strip.(split(line, ':'; limit=2))
        values[key] = value
    end
    return values
end

# The energy history and the other columns of run-energy.csv, as named
# vectors; nothing without a history.
function read_history(dir)
    file = joinpath(dir, "run-energy.csv")
    isfile(file) || return nothing
    data, header = readdlm(file, ','; header=true)
    size(data, 1) == 0 && return nothing
    names = String.(vec(header))
    columns = Dict{String,Vector{Float64}}()
    for (j, name) in enumerate(names)
        columns[name] = [x isa Number ? Float64(x) : NaN for x in data[:, j]]
    end
    return columns
end

# The first column whose name starts with a prefix, or nothing.
function column(columns, prefix)
    for (name, values) in columns
        startswith(name, prefix) && return values
    end
    return nothing
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

function ratio_at(time, ratio, t)
    t > time[end] * (1 + 1.0e-9) && return NaN
    return ratio[argmin(abs.(time .- t))]
end

# A run that ends early is failed whatever its energy did before; a
# completed run conserves to roundoff if its energy ratio stayed within
# ROUNDOFF_THRESHOLD of 1, is growing if the ratio ever exceeded the growth
# threshold, dissipating if it ended below the loss threshold, and
# conserving otherwise.
function classify(status, maximum_ratio, minimum_ratio, final_ratio)
    status != "completed" && return status == "not run" ? "not run" : "failed"
    max(maximum_ratio - 1.0, 1.0 - minimum_ratio) < ROUNDOFF_THRESHOLD && return "conserving to roundoff"
    maximum_ratio > GROWTH_THRESHOLD && return "growing"
    final_ratio < LOSS_THRESHOLD && return "dissipating"
    return "conserving"
end

format(x) = (x === nothing || isnan(x)) ? "" : @sprintf("%.4g", x)
finite(v) = filter(!isnan, v)
max_deviation(v) = (f = finite(v); isempty(f) ? NaN : maximum(abs.(f ./ f[1] .- 1)))
last_finite(v) = (f = finite(v); isempty(f) ? NaN : f[end])
max_finite(v) = (f = finite(v); isempty(f) ? NaN : maximum(f))

function summarize(dir)
    info = read_key_values(joinpath(dir, "case.txt"))
    isempty(info) && return nothing
    status_info = read_key_values(joinpath(dir, "status.txt"))
    status = get(status_info, "status", "not run")
    columns = read_history(dir)
    row = Dict{String,Any}("status" => status, "wall time" => get(status_info, "wall time", ""))
    for key in ("name", "coupling", "pair", "level", "ratio", "tolerance", "steps", "role", "solver", "dirichlet side")
        row[key] = get(info, key, "")
    end
    row["case"] = row["name"]
    row["ratio"] = row["ratio"] == "NaN" ? "" : row["ratio"]
    if columns === nothing
        time, ratio = [0.0], [NaN]
        row["dpseudo"] = row["jump u"] = row["jump v"] = row["impulse"] = ""
    else
        time = columns["time"]
        total = columns["total_energy"]
        ratio = total ./ total[1]
        pseudo = get(columns, "pseudo_energy_total", nothing)
        row["dpseudo"] = format(pseudo === nothing ? nothing : max_deviation(pseudo))
        jump_u = column(columns, "displacement_jump_rms")
        jump_v = column(columns, "velocity_jump_rms")
        impulse = column(columns, "impulse_residual_relative")
        row["jump u"] = format(jump_u === nothing ? nothing : last_finite(jump_u))
        row["jump v"] = format(jump_v === nothing ? nothing : last_finite(jump_v))
        row["impulse"] = format(impulse === nothing ? nothing : max_finite(impulse))
    end
    row["reached"] = format(time[end])
    for (k, t) in enumerate(BEAM.checkpoints)
        row["E/E0 $k"] = format(ratio_at(time, ratio, t))
    end
    row["E/E0 max"] = format(max_finite(ratio))
    row["E/E0 min"] = format(-max_finite(-ratio))
    row["dE"] = format(max_deviation(ratio))
    iterations = schwarz_iterations(dir, BEAM.maximum_iterations)
    row["iterations mean"] = format(iterations.mean)
    row["iterations max"] = string(iterations.max)
    row["stops at limit"] = string(iterations.at_limit)
    row["outcome"] = classify(status, max_finite(ratio), -max_finite(-ratio), last_finite(ratio))
    return row
end

const COLUMNS = [
    "case", "coupling", "pair", "level", "ratio", "tolerance", "steps", "role", "solver", "dirichlet side",
    "status", "outcome", "reached", "E/E0 1", "E/E0 2", "E/E0 3", "E/E0 4", "E/E0 max", "E/E0 min", "dE",
    "dpseudo", "jump u", "jump v", "impulse", "iterations mean", "iterations max", "stops at limit", "wall time",
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
    times = join((@sprintf("%8s", @sprintf("%.3g ms", 1.0e3 * t)) for t in BEAM.checkpoints), " ")
    println("\nE/E0 at $times      dE       dpseudo  it")
    for row in rows
        ratios = join((@sprintf("%8s", row["E/E0 $k"]) for k in 1:4), " ")
        @printf("  %-34s %-22s %s  %-8s %-8s %-5s %s\n", row["case"], row["outcome"], ratios, row["dE"], row["dpseudo"],
            row["iterations mean"], row["status"] == "failed" ? "ended at " * row["reached"] * " s" : "")
    end
    println("\n$(length(rows)) cases; summary written to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
