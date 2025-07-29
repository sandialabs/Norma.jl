# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

module Norma

include("utils.jl")
include("simulation.jl")

function run(input_file::String)
    norma_log(0, :norma, "BEGIN SIMULATION")
    return run(create_simulation(input_file))
end

function run(params::Parameters)
    norma_log(0, :norma, "BEGIN SIMULATION")
    return run(create_simulation(params))
end

function run(sim::Simulation)
    start_time = time()
    if get(sim.params, "enable FPE", false) == true
        enable_fpe_traps()
    end
    evolve(sim)
    elapsed_time = time() - start_time
    norma_log(0, :done, "Simulation Complete")
    norma_log(0, :time, "Run Time = " * format_time(elapsed_time))
    norma_log(0, :norma, "END SIMULATION")
    return sim
end

configure_logger()

if abspath(PROGRAM_FILE) == @__FILE__
    input_file = parse_args()
    run(input_file)
end

end
