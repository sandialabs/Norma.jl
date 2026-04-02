# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

module Norma

# --- Core utilities ---
include("utils.jl")
include("minitensor.jl")

# --- Types (order matters due to cross-dependencies) ---
include("base_types.jl")                # Simulation, Parameters
include("interpolation_types.jl")       # @enum ElementType
include("constitutive_types.jl")
include("boundary_conditions_types.jl") # needs Simulation
include("model_types.jl")              # needs BoundaryCondition
include("time_integrator_types.jl")
include("solver_types.jl")
include("simulation_types.jl")         # needs TimeIntegrator, Solver, Model

# --- Methods ---
include("constitutive.jl")
include("interpolation.jl")
include("opinf/opinf_model.jl")
include("model.jl")
include("boundary_conditions.jl")
include("schwarz.jl")
include("opinf/opinf_ics_bcs.jl")
include("time_integrator.jl")
include("solver.jl")
include("io.jl")
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
