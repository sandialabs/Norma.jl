# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

module Norma

include("utils.jl")
include("simulation.jl")
include("evolve.jl")

function run(input_file::String)
    run(create_simulation(input_file))
end

function run(params::Parameters)
    run(create_simulation(params))
end

function run(sim::Simulation)
    start_time = time()
    println("📐 Norma.jl")
    if get(sim.params, "enable FPE", false) == true
        enable_fpe_traps()
    end
    evolve(sim)
    elapsed_time = time() - start_time
    println("⏹️  Simulation Complete")
    println("⌚️ Total Time = ", format_time(elapsed_time))
    return sim
end

configure_logger()

if abspath(PROGRAM_FILE) == @__FILE__
    input_file = parse_args()
    run(input_file)
end

end
