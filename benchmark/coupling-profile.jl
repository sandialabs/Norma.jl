# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Where the time of a coupled run goes:
#   julia -t N --project=.. coupling-profile.jl case.yaml
# The case runs once to compile and once under the sampling profiler. The
# samples of thread 1, which drives the Schwarz iteration and the solves, are
# attributed inclusively to the source files of Norma (a sample counts for
# every file on its stack) and to the functions of schwarz.jl and
# simulation.jl, and exclusively to the innermost Norma function on the stack;
# the chunked loops appear under the closure of `chunk_results`, which also
# holds the wait for the other threads.
using Norma, Profile
input = ARGS[1]
Norma.run(input)
Profile.init(; n=10^8, delay=0.001)
Profile.clear()
Profile.@profile Norma.run(input)
data = Profile.fetch(; include_meta=true)
lidict = Profile.getdict(unique(data))
nmeta = Profile.nmeta
samples = Vector{Vector{UInt}}()
sleeping = 0
let i = 1, start = 1
    while i <= length(data)
        if Profile.is_block_end(data, i)
            if data[i - Profile.META_OFFSET_THREADID] == 1
                push!(samples, data[start:(i - nmeta - 1)])
                sleeping += data[i - Profile.META_OFFSET_SLEEPSTATE] == 2
            end
            start = i + 1
        end
        i += 1
    end
end
total = length(samples)
println("SAMPLES thread1=", total, " sleeping=", sleeping, " threads=", Threads.nthreads())
frames_of(ip) = [f for f in get(lidict, ip, Profile.StackFrame[]) if !f.from_c]
file_inclusive = Dict{String,Int}()
func_inclusive = Dict{String,Int}()
func_exclusive = Dict{String,Int}()
for ips in samples
    files = Set{String}()
    funcs = Set{String}()
    innermost = nothing
    for ip in ips
        for f in frames_of(ip)
            path = string(f.file)
            occursin("Norma", path) || continue
            file = basename(path)
            name = string(f.func) * " (" * file * ")"
            innermost === nothing && (innermost = name)
            push!(files, file)
            push!(funcs, name)
        end
    end
    for file in files
        file_inclusive[file] = get(file_inclusive, file, 0) + 1
    end
    for name in funcs
        func_inclusive[name] = get(func_inclusive, name, 0) + 1
    end
    innermost === nothing || (func_exclusive[innermost] = get(func_exclusive, innermost, 0) + 1)
end
percent(count) = round(100 * count / total; digits=1)
println("==== INCLUSIVE BY FILE, percent of the samples of thread 1")
for (name, count) in sort(collect(file_inclusive); by=last, rev=true)
    count / total < 0.005 && break
    println(rpad(name, 32), lpad(percent(count), 6))
end
println("==== INCLUSIVE, functions of schwarz.jl and simulation.jl")
for (name, count) in sort(collect(func_inclusive); by=last, rev=true)
    (occursin("schwarz.jl", name) || occursin("simulation.jl", name)) || continue
    count / total < 0.003 && continue
    println(rpad(name, 64), lpad(percent(count), 6))
end
println("==== EXCLUSIVE, innermost Norma function, top 30")
for (name, count) in first(sort(collect(func_exclusive); by=last, rev=true), 30)
    println(rpad(name, 64), lpad(percent(count), 6))
end
