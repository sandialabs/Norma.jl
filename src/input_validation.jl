# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Input keyword validation, in the fashion of Carina's input_parsing.jl.
#
# Norma reads its YAML input with plain Dict lookups, so an unrecognized key
# is silently ignored: `adjoint paring: false` leaves adjoint pairing on,
# `relaxation type: classical` sets nothing, and the run proceeds with
# defaults the user believes were overridden. Every loaded input file is
# therefore walked against the known-key sets below, and any key that no
# parser reads produces a warning with a "did you mean" suggestion based on
# Levenshtein distance. Warnings, not aborts: a stale key in an old deck
# should not stop a run that has always worked, and the parsers keep their
# own aborts for missing required keys and invalid values.
#
# The known-key sets mirror what the parsers actually read (see
# constitutive.jl, model.jl, solver.jl, time_integrator.jl,
# boundary_conditions.jl, simulation.jl, swap.jl, io.jl, opinf/, kroms/).
# When adding an input key to a parser, add it to the matching set here;
# the docs tables under docs/src/reference are the user-facing statement of
# the same schema.

# Levenshtein distance for "did you mean" suggestions. Indexing is over
# characters, not bytes: input keys contain non-ASCII (β, γ, HHT α, Lamé's
# first constant), where byte indexing would throw.
function levenshtein_distance(a::AbstractString, b::AbstractString)
    ca = collect(a)
    cb = collect(b)
    la = length(ca)
    lb = length(cb)
    la == 0 && return lb
    lb == 0 && return la
    prev = collect(0:lb)
    curr = similar(prev)
    for i in 1:la
        curr[1] = i
        for j in 1:lb
            cost = ca[i] == cb[j] ? 0 : 1
            curr[j + 1] = min(prev[j + 1] + 1, curr[j] + 1, prev[j] + cost)
        end
        prev, curr = curr, prev
    end
    return prev[lb + 1]
end

function suggest_key(key::AbstractString, known::AbstractSet{String}; max_distance::Int=5)
    best_distance = max_distance + 1
    best_key = ""
    for k in known
        d = levenshtein_distance(lowercase(key), lowercase(k))
        if d < best_distance
            best_distance = d
            best_key = k
        end
    end
    return best_distance <= max_distance ? best_key : ""
end

# --- Known keys per input-file section ---

# Keys written into params at runtime ("name") are tolerated so that
# validation order relative to the bookkeeping writes does not matter.
const TOP_LEVEL_SINGLE_KEYS = Set([
    "type",
    "name",
    "input mesh file",
    "output mesh file",
    "time integrator",
    "solver",
    "model",
    "boundary conditions",
    "initial conditions",
    "restart",
    "swaps",
    "Exodus output interval",
    "CSV output interval",
    "CSV write sidesets",
    "mesh smoothing",
])

const TOP_LEVEL_MULTI_KEYS = Set([
    "type",
    "name",
    "domains",
    "minimum iterations",
    "maximum iterations",
    "absolute tolerance",
    "relative tolerance",
    "initial time",
    "final time",
    "time step",
    "relaxation",
    "relaxation parameter",
    "aitken N0 parameter",
    "naive stabilized",
    "interface predictor",
    "stalled interface jump action",
    "unconverged step action",
    "blended energy output",
    "Exodus output interval",
    "CSV output interval",
    "CSV write sidesets",
    "restart",
    "swaps",
])

const RESTART_KEYS = Set(["index"])

const TIME_INTEGRATOR_COMMON_KEYS = Set([
    "type",
    "time step",
    "initial time",
    "final time",
    "minimum time step",
    "maximum time step",
    "decrease factor",
    "increase factor",
])

const TIME_INTEGRATOR_TYPE_KEYS = Dict(
    "quasi static" => Set(["initial equilibrium"]),
    "Newmark" => Set(["β", "γ", "HHT α", "HHT alpha"]),
    "central difference" => Set(["CFL", "γ"]),
)

const SOLVER_KEYS = Set([
    "type",
    "step",
    "minimum iterations",
    "maximum iterations",
    "absolute tolerance",
    "relative tolerance",
    "linear solver absolute tolerance",
    "linear solver relative tolerance",
    "use line search",
    "line search backtrack factor",
    "line search decrease factor",
    "line search maximum iterations",
    "step length",
    "memory",
    "energy stagnation window",
    "energy stagnation tolerance",
    "rom linear solver",
])

# `stress recovery` and `recover internal variables` are legacy keys whose
# presence the model parser rejects with its own targeted message; they are
# known here so that the parser's abort is not preceded by a spurious
# unknown-key warning.
const MODEL_KEYS = Set([
    "type",
    "material",
    "num integration points",
    "smooth reference",
    "size field",
    "nodal recovery",
    "stress recovery",
    "recover internal variables",
    "model-file",
    "model-directory",
    "ensemble-size",
    "rbf-type",
])

const NODAL_RECOVERY_KEYS = Set(["method", "stress", "von mises stress", "internal variables", "deformation gradient"])

const MATERIAL_PROPS_COMMON_KEYS = Set([
    "model", "elastic modulus", "Poisson's ratio", "bulk modulus", "Lamé's first constant", "shear modulus", "density"
])

const MATERIAL_PROPS_MODEL_KEYS = Dict(
    "seth-hill" => Set(["m", "n"]), "j2 plasticity" => Set(["yield stress", "hardening modulus"])
)

# Material model names accepted by create_material (constitutive.jl). A
# recognized model validates against the common keys plus its own extras; an
# unrecognized model (the material parser's abort to report) validates
# against the union so only genuinely foreign keys warn.
const MATERIAL_MODEL_NAMES = Set([
    "linear elastic", "Saint-Venant Kirchhoff", "neohookean", "r-neohookean", "seth-hill", "hencky", "j2 plasticity"
])

const BC_ENTRY_KEYS = Dict(
    "Dirichlet" => Set(["node set", "side set", "function", "component"]),
    "OpInf Dirichlet" => Set(["node set", "function", "component", "model-directory", "ensemble-size"]),
    "Neumann" => Set(["side set", "function", "component"]),
    "Neumann pressure" => Set(["side set", "function"]),
    "Surface" => Set(["side set", "function", "enforcement", "penalty"]),
    "Robin" => Set(["side set", "function", "component", "robin parameter"]),
    "Schwarz contact" => Set(["source", "side set", "source side set", "friction type", "swap BC types"]),
    "Schwarz overlap" => Set([
        "source",
        "side set",
        "source block",
        "source side set",
        "search tolerance",
        "weak",
        "compute overlap L2 relative error",
    ]),
    "Schwarz DN nonoverlap" => Set([
        "source",
        "side set",
        "source side set",
        "source block",
        "search tolerance",
        "default BC type",
        "swap BC types",
    ]),
    "Schwarz impedance nonoverlap" => Set([
        "source",
        "side set",
        "source side set",
        "source block",
        "search tolerance",
        "robin parameter",
        "impedance scale",
        "adjoint pairing",
    ]),
    "Schwarz impedance overlap" => Set([
        "source",
        "side set",
        "source block",
        "source side set",
        "search tolerance",
        "robin parameter",
        "impedance scale",
        "partner traction",
        "transfer",
        "transfer quadrature subdivisions",
        "content aware absorption",
        "representable dashpot",
    ]),
    "OpInf Schwarz overlap" => Set([
        "source",
        "side set",
        "source block",
        "search tolerance",
        "model-directory",
        "ensemble-size",
        "compute overlap L2 relative error",
    ]),
)
# The classical Robin-Robin condition accepts a subset of the impedance keys:
# `impedance scale` is rejected by its parser with a hard abort that points at
# `Schwarz impedance nonoverlap`. Share the impedance key set here so that
# abort is not preceded by a spurious unknown-key warning.
BC_ENTRY_KEYS["Schwarz RR nonoverlap"] = BC_ENTRY_KEYS["Schwarz impedance nonoverlap"]

const BC_TYPE_KEYS = Set(keys(BC_ENTRY_KEYS))

const IC_ENTRY_KEYS = Set(["node set", "function", "component"])

const IC_TYPE_KEYS = Set(["displacement", "velocity"])

const SWAP_PLAN_KEYS = Set(["subsim", "criterion", "replacement"])

const SWAP_CRITERION_KEYS = Set(["type", "t_swap", "tolerance", "direction"])

# --- Validation walk ---

"""
Warn for every key of `dict` that is not in `known`, with a "did you mean"
suggestion. `section` names the location inside `file` for the message.
Returns the messages so callers (and tests) can inspect what was reported.
"""
function warn_unknown_keys!(
    messages::Vector{String}, dict::AbstractDict, known::AbstractSet{String}, section::String, file::String
)
    for key in keys(dict)
        key isa String || continue
        key in known && continue
        message = "Input file '$file': unknown key \"$key\" in $section."
        suggestion = suggest_key(key, known)
        if !isempty(suggestion)
            message *= " Did you mean \"$suggestion\"?"
        end
        norma_log(0, :warning, message)
        push!(messages, message)
    end
    return messages
end

function validate_time_integrator!(messages::Vector{String}, params::Parameters, file::String)
    integrator_params = get(params, "time integrator", nothing)
    integrator_params isa AbstractDict || return messages
    integrator_type = get(integrator_params, "type", "")
    known = if haskey(TIME_INTEGRATOR_TYPE_KEYS, integrator_type)
        union(TIME_INTEGRATOR_COMMON_KEYS, TIME_INTEGRATOR_TYPE_KEYS[integrator_type])
    else
        # Unknown or absent type: the integrator parser aborts with its own
        # message; validate against the union so only genuinely foreign keys
        # warn.
        union(TIME_INTEGRATOR_COMMON_KEYS, values(TIME_INTEGRATOR_TYPE_KEYS)...)
    end
    return warn_unknown_keys!(messages, integrator_params, known, "time integrator", file)
end

function validate_solver!(messages::Vector{String}, params::Parameters, file::String)
    solver_params = get(params, "solver", nothing)
    solver_params isa AbstractDict || return messages
    return warn_unknown_keys!(messages, solver_params, SOLVER_KEYS, "solver", file)
end

function validate_material!(messages::Vector{String}, material_params::AbstractDict, file::String)
    blocks = get(material_params, "blocks", nothing)
    labels = Set{String}()
    if blocks isa AbstractDict
        for value in values(blocks)
            value isa AbstractString && push!(labels, String(value))
        end
    end
    # Beyond `blocks`, every key of the material section should be a property
    # dict that some block references; anything else is typically a renamed
    # block whose old material was left behind.
    warn_unknown_keys!(messages, material_params, union(Set(["blocks"]), labels), "model material", file)
    for label in labels
        props = get(material_params, label, nothing)
        props isa AbstractDict || continue
        material_model = get(props, "model", "")
        known = if material_model in MATERIAL_MODEL_NAMES
            union(MATERIAL_PROPS_COMMON_KEYS, get(MATERIAL_PROPS_MODEL_KEYS, material_model, Set{String}()))
        else
            union(MATERIAL_PROPS_COMMON_KEYS, values(MATERIAL_PROPS_MODEL_KEYS)...)
        end
        warn_unknown_keys!(messages, props, known, "material \"$label\"", file)
    end
    return messages
end

function validate_model!(messages::Vector{String}, params::Parameters, file::String)
    model_params = get(params, "model", nothing)
    model_params isa AbstractDict || return messages
    warn_unknown_keys!(messages, model_params, MODEL_KEYS, "model", file)
    material_params = get(model_params, "material", nothing)
    if material_params isa AbstractDict
        validate_material!(messages, material_params, file)
    end
    recovery_params = get(model_params, "nodal recovery", nothing)
    if recovery_params isa AbstractDict
        warn_unknown_keys!(messages, recovery_params, NODAL_RECOVERY_KEYS, "nodal recovery", file)
    end
    return messages
end

function validate_boundary_conditions!(messages::Vector{String}, params::Parameters, file::String)
    bc_params = get(params, "boundary conditions", nothing)
    bc_params isa AbstractDict || return messages
    warn_unknown_keys!(messages, bc_params, BC_TYPE_KEYS, "boundary conditions", file)
    for (bc_type, entries) in bc_params
        haskey(BC_ENTRY_KEYS, bc_type) || continue
        entries isa AbstractVector || continue
        for (index, entry) in enumerate(entries)
            entry isa AbstractDict || continue
            warn_unknown_keys!(messages, entry, BC_ENTRY_KEYS[bc_type], "$bc_type boundary condition $index", file)
        end
    end
    return messages
end

function validate_initial_conditions!(messages::Vector{String}, params::Parameters, file::String)
    ic_params = get(params, "initial conditions", nothing)
    ic_params isa AbstractDict || return messages
    warn_unknown_keys!(messages, ic_params, IC_TYPE_KEYS, "initial conditions", file)
    for (ic_type, entries) in ic_params
        ic_type in IC_TYPE_KEYS || continue
        entries isa AbstractVector || continue
        for (index, entry) in enumerate(entries)
            entry isa AbstractDict || continue
            warn_unknown_keys!(messages, entry, IC_ENTRY_KEYS, "$ic_type initial condition $index", file)
        end
    end
    return messages
end

function validate_restart!(messages::Vector{String}, params::Parameters, file::String)
    restart_params = get(params, "restart", nothing)
    restart_params isa AbstractDict || return messages
    return warn_unknown_keys!(messages, restart_params, RESTART_KEYS, "restart", file)
end

function validate_swaps!(messages::Vector{String}, params::Parameters, file::String)
    swap_plans = get(params, "swaps", nothing)
    swap_plans isa AbstractVector || return messages
    for (index, plan) in enumerate(swap_plans)
        plan isa AbstractDict || continue
        warn_unknown_keys!(messages, plan, SWAP_PLAN_KEYS, "swap plan $index", file)
        criterion = get(plan, "criterion", nothing)
        criterion isa AbstractDict || continue
        warn_unknown_keys!(messages, criterion, SWAP_CRITERION_KEYS, "swap plan $index criterion", file)
    end
    return messages
end

"""
Validate the keys of one loaded input file against the sets the parsers
actually read. `file` names the file in the warnings. Unknown keys warn (with
a "did you mean" suggestion) and never abort; the parsers retain their own
aborts for missing required keys and invalid values. Returns the warning
messages, empty when every key is recognized.
"""
function validate_input_parameters(params::Parameters, file::AbstractString)
    messages = String[]
    file_name = String(file)
    simulation_type = get(params, "type", "")
    if simulation_type == "multi"
        warn_unknown_keys!(messages, params, TOP_LEVEL_MULTI_KEYS, "the top level", file_name)
        validate_restart!(messages, params, file_name)
        validate_swaps!(messages, params, file_name)
    else
        # `type: single` files, subdomain files, and swap replacement files
        # all share the single-domain format. An absent or unrecognized
        # `type` is the simulation parser's abort to report.
        warn_unknown_keys!(messages, params, TOP_LEVEL_SINGLE_KEYS, "the top level", file_name)
        validate_time_integrator!(messages, params, file_name)
        validate_solver!(messages, params, file_name)
        validate_model!(messages, params, file_name)
        validate_boundary_conditions!(messages, params, file_name)
        validate_initial_conditions!(messages, params, file_name)
        validate_restart!(messages, params, file_name)
        validate_swaps!(messages, params, file_name)
    end
    return messages
end
