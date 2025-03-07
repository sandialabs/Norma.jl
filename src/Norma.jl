# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
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
        @info "Debugging disabled"
    end
end

configure_logger()

include("minitensor.jl")
include("simulation.jl")
include("evolve.jl")

function run(input_file::String)
    sim = create_simulation(input_file)
    evolve(sim)
    return sim
end

function run(params::Parameters, name::String)
    sim = create_simulation(params, name)
    evolve(sim)
    return sim
end

for input_file in ARGS
    run(input_file)
end

end
