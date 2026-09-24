# Run cases of the dynamic stability study, several at a time, each in its
# own Julia process.
#
#   julia --project=../.. run.jl [options] TIER|CASE|all ...
#
# Every argument is a tier letter, a case name, or `all` (every case
# directory under the runs directory).  The cases must have been generated
# first.  A case whose status.txt says completed is skipped unless --force
# is given, so an interrupted tier can be resumed by running the same
# command again.  The output of each run goes to run.log in its directory.
#
# Options:
#   --runs DIR     directory of the case directories (default: runs)
#   --jobs N       cases run at the same time (default 2)
#   --threads T    Julia threads per case (default 4)
#   --force        rerun cases that already completed
#
# Choose jobs times threads at most the number of cores.  The level 2
# cylinder cases need several gigabytes of memory each.

include(joinpath(@__DIR__, "matrix.jl"))

function parse_arguments(args)
    options = Dict{String,Any}("runs" => joinpath(@__DIR__, "runs"), "jobs" => 2, "threads" => 4, "force" => false)
    selections = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--runs"
            options["runs"] = abspath(args[i + 1]); i += 2
        elseif arg == "--jobs"
            options["jobs"] = parse(Int, args[i + 1]); i += 2
        elseif arg == "--threads"
            options["threads"] = parse(Int, args[i + 1]); i += 2
        elseif arg == "--force"
            options["force"] = true; i += 1
        else
            push!(selections, arg); i += 1
        end
    end
    isempty(selections) && error("Give one or more tiers ($(join(first.(TIERS), ", "))), case names, or all")
    return selections, options
end

function case_directories(selections, runs)
    names = String[]
    for s in selections
        if s == "all"
            append!(names, filter(n -> isfile(joinpath(runs, n, "run.yaml")), readdir(runs)))
        elseif find_case(s) !== nothing
            push!(names, s)
        else
            append!(names, case_name.(tier_cases(s)))
        end
    end
    dirs = [joinpath(runs, n) for n in unique(names)]
    missing_dirs = filter(d -> !isfile(joinpath(d, "run.yaml")), dirs)
    isempty(missing_dirs) || error("Not generated yet: $(join(basename.(missing_dirs), ", ")); run generate.jl first")
    return dirs
end

function completed(dir)
    file = joinpath(dir, "status.txt")
    return isfile(file) && occursin("status: completed", read(file, String))
end

function main(args)
    selections, options = parse_arguments(args)
    dirs = case_directories(selections, options["runs"])
    options["force"] || (dirs = filter(!completed, dirs))
    project = normpath(joinpath(@__DIR__, "..", ".."))
    script = joinpath(@__DIR__, "run_case.jl")
    julia = Base.julia_cmd()
    println("$(length(dirs)) cases to run, $(options["jobs"]) at a time with $(options["threads"]) threads each")
    finished = Ref(0)
    asyncmap(dirs; ntasks=options["jobs"]) do dir
        rm(joinpath(dir, "status.txt"); force=true)
        command = `$julia --project=$project -t $(options["threads"]) $script $dir`
        command = addenv(command, "OPENBLAS_NUM_THREADS" => "1")
        start = time()
        run(pipeline(ignorestatus(command); stdout=joinpath(dir, "run.log"), stderr=joinpath(dir, "run.log")))
        status = isfile(joinpath(dir, "status.txt")) ? split(readline(joinpath(dir, "status.txt")))[end] : "crashed"
        finished[] += 1
        println("[$(finished[])/$(length(dirs))] $(basename(dir)): $status in $(round(Int, time() - start)) s")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
