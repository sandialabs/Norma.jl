# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

module Norma

using Logging

# Enable debugging for a specific module in the environment variable JULIA_DEBUG (all for all modules)
function configure_logger()
    debug_level = get(ENV, "JULIA_DEBUG", "")

    if debug_level != ""
        # Enable debug logging
        global_logger(ConsoleLogger(stderr, Logging.Debug))
        @info "Debugging enabled for: $debug_level"
    else
        # Default to Info level
        global_logger(ConsoleLogger(stderr, Logging.Info))
    end
end

configure_logger()

include("minitensor.jl")
include("simulation.jl")
include("evolve.jl")

function format_time(seconds::Float64)::String
    days = floor(Int, seconds / 86400)  # 86400 seconds in a day
    seconds %= 86400

    hours = floor(Int, seconds / 3600)  # 3600 seconds in an hour
    seconds %= 3600

    minutes = floor(Int, seconds / 60)  # 60 seconds in a minute
    seconds %= 60

    # Build the formatted string based on non-zero values
    time_str = []
    if days > 0
        push!(time_str, @sprintf("%dd", days))
    end
    if hours > 0 || days > 0  # Show hours if days > 0 for completeness
        push!(time_str, @sprintf("%dh", hours))
    end
    if minutes > 0 || hours > 0 || days > 0  # Show minutes if higher units are non-zero
        push!(time_str, @sprintf("%dm", minutes))
    end
    push!(time_str, @sprintf("%.1fs", seconds))

    return join(time_str, " ")
end

function run(input_file::String)
    start_time = time()
    sim = create_simulation(input_file)
    evolve(sim)
    elapsed_time = time() - start_time
    println("Simulation completed in ", format_time(elapsed_time))
    return sim
end

function run(params::Parameters, name::String)
    start_time = time()
    sim = create_simulation(params, name)
    evolve(sim)
    elapsed_time = time() - start_time
    println("Simulation completed in ", format_time(elapsed_time))
    return sim
end

for input_file in ARGS
    run(input_file)
end

end
