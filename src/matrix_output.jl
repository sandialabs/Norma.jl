# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Matrix output: writes the assembled global operators of a single-domain
# solid mechanics simulation so the discrete problem can be reproduced
# outside Norma. Requested with an input block such as
#
#   matrix output:
#     matrices: [stiffness, mass, external force]   # default: all three
#     steps: [0]                                    # default; or all, [0, 10], {every: 5}
#
# Files, with the simulation name as prefix and the stop index padded to 4:
#
#   <name>-stiffness-NNNN.mm        tangent stiffness K (Matrix Market)
#   <name>-mass.mm                  consistent mass M, written once
#   <name>-lumped-mass.mm           lumped mass as a diagonal matrix, written once
#   <name>-external-force-NNNN.csv  body + boundary (Neumann) force, nodal layout
#   <name>-free-dofs-NNNN.csv       one line per DOF: true = free, false = constrained
#
# DOF numbering is node-major: DOF 3(n-1)+c holds component c (1=x, 2=y, 3=z)
# of node n (1-based, the Exodus node order). Each line of a nodal CSV holds
# the x,y,z values of one node. Matrices include every DOF, constrained or
# not; the free-dofs file tells the two apart.
#
# Only monolithic (type: single) runs are supported; Schwarz subdomains abort.

const MATRIX_OUTPUT_NAMES = ("stiffness", "mass", "external force")

struct MatrixOutputOptions
    matrices::Vector{String}
    all_steps::Bool
    step_list::Vector{Int64}
    every::Int64
end

_is_positive_integer(value) = value isa Integer && value > 0

function parse_matrix_output(params::Parameters)
    block = params["matrix output"]
    block isa AbstractDict || norma_abort("`matrix output` must be a block with `matrices` and/or `steps`")
    matrices_param = get(block, "matrices", collect(MATRIX_OUTPUT_NAMES))
    matrices_param isa AbstractVector ||
        norma_abort("`matrix output: matrices` must be a list, e.g. [stiffness, mass, external force]")
    matrices = String[]
    for name in matrices_param
        name_string = string(name)
        if !(name_string in MATRIX_OUTPUT_NAMES)
            norma_abort(
                "Unknown `matrix output` entry '$name_string'. Valid entries: " * join(MATRIX_OUTPUT_NAMES, ", ")
            )
        end
        name_string in matrices || push!(matrices, name_string)
    end
    steps = get(block, "steps", [0])
    all_steps = false
    step_list = Int64[]
    every = 0
    if steps == "all"
        all_steps = true
    elseif steps isa Integer
        push!(step_list, steps)
    elseif steps isa AbstractVector && all(s -> s isa Integer, steps)
        append!(step_list, steps)
    elseif steps isa AbstractDict && length(steps) == 1 && _is_positive_integer(get(steps, "every", 0))
        every = steps["every"]
    else
        norma_abort(
            "`matrix output: steps` must be `all`, a stop index, a list of stop indices, or `{every: n}` with n > 0"
        )
    end
    any(<(0), step_list) && norma_abort("`matrix output: steps` must not contain negative stop indices")
    return MatrixOutputOptions(matrices, all_steps, step_list, every)
end

function is_matrix_output_stop(options::MatrixOutputOptions, stop::Integer)
    options.all_steps && return true
    options.every > 0 && return stop % options.every == 0
    return stop in options.step_list
end

function reject_matrix_output_in_schwarz(params::Parameters, file::AbstractString)
    if haskey(params, "matrix output")
        norma_abort(
            "`matrix output` is only supported for monolithic (type: single) simulations, " *
            "but it appears in the Schwarz (type: multi) input '$file'. Remove it from that file.",
        )
    end
    return nothing
end

# Matrix Market coordinate format, general storage, 1-based indices.
# Every stored entry of the sparse matrix is written, with 17 significant
# digits so that the values round-trip exactly.
function write_matrix_market(filename::String, A::SparseMatrixCSC{Float64,Int64}, description::String)
    num_rows, num_cols = size(A)
    rows = rowvals(A)
    values = nonzeros(A)
    open(filename, "w") do io
        println(io, "%%MatrixMarket matrix coordinate real general")
        println(io, "% Norma ", description)
        println(io, "% DOF 3*(node-1)+component, node-major, 1-based; components 1=x 2=y 3=z")
        @printf(io, "%d %d %d\n", num_rows, num_cols, nnz(A))
        for col in 1:num_cols
            for k in nzrange(A, col)
                @printf(io, "%d %d %.16e\n", rows[k], col, values[k])
            end
        end
    end
    return nothing
end

function write_matrix_output(sim::SingleDomainSimulation)
    haskey(sim.params, "matrix output") || return nothing
    is_coupled(sim) && reject_matrix_output_in_schwarz(sim.params, sim.name)
    return write_matrix_output(sim, sim.model)
end

function write_matrix_output(sim::SingleDomainSimulation, ::Model)
    norma_abort("`matrix output` is only supported for the solid mechanics full-order model")
    return nothing
end

function write_matrix_output(sim::SingleDomainSimulation, model::SolidMechanics)
    params = sim.params
    options = parse_matrix_output(params)
    stop = sim.controller.stop
    is_matrix_output_stop(options, stop) || return nothing
    norma_log(0, :output, "Matrix Market operators for $(sim.name) [MM]")
    prefix = sim.name * "-"
    index_string = "-" * string(stop; pad=4)
    num_dofs = length(model.free_dofs)
    if "stiffness" in options.matrices
        if size(model.stiffness) != (num_dofs, num_dofs)
            norma_abort(
                "`matrix output: stiffness` needs an assembled stiffness matrix, which requires an " *
                "implicit time integrator (quasi-static or Newmark) with the Hessian minimizer solver",
            )
        end
        write_matrix_market(
            prefix * "stiffness" * index_string * ".mm",
            model.stiffness,
            "stiffness matrix K at stop $stop, time $(sim.controller.time)",
        )
    end
    if "mass" in options.matrices && get(params, "_matrix output mass written", false) == false
        if size(model.mass) == (num_dofs, num_dofs)
            write_matrix_market(prefix * "mass.mm", model.mass, "consistent mass matrix M")
        elseif length(model.lumped_mass) == num_dofs
            write_matrix_market(
                prefix * "lumped-mass.mm", spdiagm(0 => model.lumped_mass), "lumped mass matrix M (diagonal)"
            )
        else
            norma_abort("`matrix output: mass` needs a dynamic time integrator; this simulation has no mass matrix")
        end
        params["_matrix output mass written"] = true
    end
    if "external force" in options.matrices
        external_force = copy(model.boundary_force)
        length(model.body_force) == num_dofs && (external_force .+= model.body_force)
        writedlm_nodal_array(prefix * "external-force" * index_string * ".csv", reshape(external_force, 3, :))
    end
    # The free/constrained map is always written: the operators cannot be
    # used without it, and CSV output may be off.
    writedlm(prefix * "free-dofs" * index_string * ".csv", model.free_dofs)
    return nothing
end
