# Run one case of the dynamic stability study in its directory and record
# the outcome in status.txt: completed or failed, the wall time, and the
# error message of a failed run.  run.jl calls this once per case; it can
# also be run by hand:
#
#   julia --project=../.. -t 4 run_case.jl runs/<case>

using Norma

function run_case(dir::AbstractString)
    cd(dir)
    status = "completed"
    message = ""
    start = time()
    try
        Norma.run("run.yaml")
    catch error
        status = "failed"
        message = replace(sprint(showerror, error), '\n' => ' ')
    end
    wall_time = time() - start
    write("status.txt", "status: $status\nwall time: $(round(wall_time; digits=1))\nmessage: $(first(message, 400))\n")
    return status
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_case(ARGS[1])
end
