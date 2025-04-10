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
    start_time = time()
    sim = create_simulation(input_file)
    if get(sim.params, "enable FPE", false) == true
        enable_fpe_traps()
    end
    evolve(sim)
    elapsed_time = time() - start_time
    println("Simulation completed in ", format_time(elapsed_time))
    return sim
end

function run(params::Parameters, name::String)
    start_time = time()
    sim = create_simulation(params, name)
    if get(sim.params, "enable FPE", false) == true
        enable_fpe_traps()
    end
    evolve(sim)
    elapsed_time = time() - start_time
    println("Simulation completed in ", format_time(elapsed_time))
    return sim
end

configure_logger()

if abspath(PROGRAM_FILE) == @__FILE__
    input_file = parse_args()
    run(input_file)
end

end
