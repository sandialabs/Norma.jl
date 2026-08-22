# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Printf
using YAML


function create_simulation(input_file::String)
    norma_log(0, :setup, "Reading from " * input_file)
    params = YAML.load_file(input_file; dicttype=Parameters)
    basename = stripped_name(input_file)
    params["name"] = basename
    return create_simulation(params)
end

# Whether a `model: type:` string supports restart is now resolved via
# supports_restart(model_type_for(model_type)) (model.jl / model_types.jl)
# rather than a hand-kept list here -- see process_restart!() below. "solid
# mechanics" is the FOM model; every ROM model type (see model_type_for() in
# model.jl for the canonical list of recognized `model: type:` strings) also
# supports it. Each ROM type wraps an internal FOM model
# (`fom_model::SolidMechanics`) that is restored from the restart snapshot
# the same way a standalone FOM model is; see process_restart!() and
# apply_ics(::Parameters, ::RomModel, ...) for how the restored FOM state is
# then projected onto the reduced basis.
#
# This applies uniformly to single-domain and multi-domain (Schwarz-coupled)
# simulations: any combination of FOM and ROM subdomains (including mixed
# FOM-ROM pairings) may restart. See process_multidomain_restart!() below for
# the multi-domain-specific checks.

# Restart and mid-run mesh swapping (`swaps:`, see swap.jl) are mutually
# incompatible in the current implementation and must never be combined.
# Restart seeds the initial displacement/velocity fields with a *positional*
# read from the input mesh file (see process_restart!() below): node i of the
# restart snapshot is assigned directly to node i of the model, with no
# node-ID cross-reference, coordinate check, or interpolation. Swap plans can
# replace that model outright with one on a different mesh. The combination
# has not been validated end-to-end and gives no correctness guarantee, so it
# is rejected outright here rather than risking silently wrong results.
# Checked wherever a `params` (or subdomain `subparams`) dict could carry
# both a `restart:` and a `swaps:` block: SingleDomainSimulation() (covers
# both standalone single-domain runs and individual Schwarz subdomains,
# whose `restart:` may have been injected by process_multidomain_restart!())
# and MultiDomainSimulation() (covers a shared top-level `restart:` combined
# with top-level `swaps:` plans that target subdomains by `subsim:` name).
function _reject_restart_with_swaps!(params::Parameters, context::String)
    if haskey(params, "restart") && haskey(params, "swaps")
        norma_abort(
            "`restart:` and `swaps:` cannot be used together ($context). Restart seeds " *
            "the initial state with a positional (node-for-node) read from the input " *
            "mesh file and has no cross-mesh mapping, while `swaps:` may replace the " *
            "model with one on a different mesh mid-run; combining the two is not " *
            "supported. Remove one or the other.",
        )
    end
    return nothing
end

# Validate and normalize a `restart: index:` value. Integer-valued numbers
# (an `Integer`, or an integer-valued `AbstractFloat` such as YAML's parse of
# `2.0`) are accepted and converted to Int64. Anything else (e.g. `2.5`)
# previously reached `Int64(...)` directly, which throws a raw, uninformative
# `InexactError` instead of a clean, actionable abort message.
function _validate_restart_index(index_value)
    if !(index_value isa Integer || (index_value isa AbstractFloat && isinteger(index_value)))
        norma_abortf(
            "`restart: index:` must be an integer, got %s (a %s).",
            string(index_value),
            string(typeof(index_value)),
        )
    end
    return Int64(index_value)
end

# Read the restart snapshot (time + nodal displacement/velocity fields) for a
# single-domain simulation from its (already-opened) input mesh, validate the
# request, and stash the result on params under "restart_info" for SolidMechanics
# to consume. Also overwrites params["time integrator"]["initial time"] so that
# the controller and integrator are built with the restart time from the outset
# (num_stops depends on it). No-op (params["restart_info"] = nothing) when the
# input file has no `restart:` block.
function process_restart!(params::Parameters, input_mesh, basename::String)
    if haskey(params, "restart") == false
        params["restart_info"] = nothing
        return nothing
    end
    if haskey(params, "initial conditions")
        norma_abort(
            "Cannot specify both `restart:` and `initial conditions:` in the same " *
            "input file ($basename). The restart snapshot already supplies the " *
            "initial displacement and velocity fields; remove one or the other.",
        )
    end
    model_params = get(params, "model", Parameters())
    model_type = get(model_params, "type", "")
    if model_type == "mesh smoothing"
        norma_abort(
            "Restart is not supported for `model: type: mesh smoothing`: mesh smoothing " *
            "is not a stateful dynamic simulation, so there is no displacement/velocity " *
            "state for a restart snapshot to resume.",
        )
    end
    julia_model_type = model_type_for(model_type)
    if julia_model_type === nothing || !supports_restart(julia_model_type)
        norma_abortf(
            "Restart is only supported for `model: type: solid mechanics` (FOM) models " *
            "and ROM models (linear/quadratic/cubic OpInf or kernel ROM, neural network " *
            "OpInf ROM, RBF kernel ROM). Got `model: type: %s`.",
            model_type,
        )
    end
    # ROM restart is built on top of the FOM restart machinery: every ROM model
    # type constructs an internal `fom_model::SolidMechanics` (see opinf_model.jl
    # / krom_model.jl), which is restored from the restart snapshot exactly like
    # a standalone FOM model. apply_ics() then projects the restored FOM
    # displacement/velocity onto the reduced basis (see opinf_ics_bcs.jl). ROM
    # restart is therefore only as good as FOM restart — in particular it is
    # still subject to the underlying J2-plasticity internal-variable
    # limitation enforced in SolidMechanics().
    #
    # Multi-domain (Schwarz-coupled) restart supports both FOM
    # (`solid mechanics`) and ROM subdomains — including mixed FOM-ROM
    # pairings — since ROM restart is layered on the same FOM machinery used
    # here; no extra check is needed in this per-subdomain function. See
    # process_multidomain_restart!() for the multi-domain-specific checks
    # (shared restart index/time, no per-subdomain `restart:` blocks).
    restart_params = params["restart"]
    haskey(restart_params, "index") || norma_abort("`restart:` block must specify an `index`.")
    requested_restart_index = _validate_restart_index(restart_params["index"])
    num_steps = Exodus.read_number_of_time_steps(input_mesh)
    if num_steps < 1
        norma_abort("Restart mesh '$(input_mesh.file_name)' contains no time steps.")
    end
    # Negative indices count back from the last snapshot, Python-style:
    # -1 is the final snapshot, -2 the one before it, and so on.
    restart_index =
        requested_restart_index < 0 ? num_steps + requested_restart_index + 1 : requested_restart_index
    if restart_index < 1 || restart_index > num_steps
        norma_abortf(
            "Restart index %d is out of range for mesh '%s', which has %d time step(s) " *
            "(valid range is 1:%d, or -1:-%d counting back from the last snapshot).",
            requested_restart_index,
            input_mesh.file_name,
            num_steps,
            num_steps,
            num_steps,
        )
    end
    restart_time = Exodus.read_time(input_mesh, restart_index)
    integrator_params = params["time integrator"]
    # A restart at or past `final time` leaves nothing to integrate. Without
    # this check, num_stops = max(round((final_time - restart_time) /
    # time_step) + 1, 2) is forced to at least 2 regardless, and evolve()
    # always advances one control step before checking stop_evolve(), so the
    # run would silently take one step past final_time and write a bogus
    # extra stop there. Abort instead.
    final_time = get(integrator_params, "final time", nothing)
    if final_time !== nothing
        final_time = Float64(final_time)
        if restart_time > final_time || isapprox(restart_time, final_time; rtol=1e-9, atol=1e-12)
            norma_abortf(
                "Restart time %.6e (index %d of '%s') is at or past `time integrator: " *
                "final time` %.6e; there is nothing left to integrate. Choose an earlier " *
                "restart index, or increase `final time`.",
                restart_time,
                restart_index,
                input_mesh.file_name,
                final_time,
            )
        end
    end
    disp_x = Exodus.read_values(input_mesh, NodalVariable, restart_index, "disp_x")
    disp_y = Exodus.read_values(input_mesh, NodalVariable, restart_index, "disp_y")
    disp_z = Exodus.read_values(input_mesh, NodalVariable, restart_index, "disp_z")
    num_nodes = Exodus.num_nodes(input_mesh.init)
    displacement = Matrix{Float64}(undef, 3, num_nodes)
    displacement[1, :] = disp_x
    displacement[2, :] = disp_y
    displacement[3, :] = disp_z
    # Only a dynamic-integrator run ever writes velo_x/y/z to its Exodus
    # output (see is_dynamic(integrator) in initialize_writing(), io.jl); a
    # checkpoint written by a quasistatic (StaticTimeIntegrator) run has no
    # velocity fields at all, since velocity is not a meaningful quasistatic
    # quantity. Reading them unconditionally would throw a raw, uninformative
    # Exodus "variable not found" error. Use zero initial velocity instead
    # when the checkpoint has no velocity fields -- this is the correct
    # restart behavior for a quasistatic checkpoint, not just a fallback: a
    # quasistatic run can never have had a nonzero velocity to lose (nor can
    # `initial conditions:` have set one, since it is rejected above whenever
    # `restart:` is present).
    velocity = Matrix{Float64}(undef, 3, num_nodes)
    available_node_vars = Exodus.read_names(input_mesh, NodalVariable)
    has_velocity = all(v -> v in available_node_vars, ("velo_x", "velo_y", "velo_z"))
    if has_velocity
        velo_x = Exodus.read_values(input_mesh, NodalVariable, restart_index, "velo_x")
        velo_y = Exodus.read_values(input_mesh, NodalVariable, restart_index, "velo_y")
        velo_z = Exodus.read_values(input_mesh, NodalVariable, restart_index, "velo_z")
        velocity[1, :] = velo_x
        velocity[2, :] = velo_y
        velocity[3, :] = velo_z
    else
        velocity .= 0.0
        norma_log(
            0,
            :restart,
            "Restart checkpoint has no velocity fields (quasistatic); using zero initial velocity.",
        )
    end
    # Override the initial time unconditionally: the restart snapshot's time is
    # the only initial time consistent with the restart displacement/velocity
    # fields, regardless of what `time integrator: initial time` says in the
    # input file.
    requested_initial_time = get(integrator_params, "initial time", nothing)
    integrator_params["initial time"] = restart_time
    norma_log(0, :restart, "Restarting simulation from snapshot data.")
    norma_logf(0, :restart, "Restart index:  %d (of %d)", restart_index, num_steps)
    norma_logf(0, :restart, "Restart time:   %.6e", restart_time)
    if requested_initial_time !== nothing && !isapprox(Float64(requested_initial_time), restart_time; atol=1e-12)
        norma_logf(
            0,
            :restart,
            "Overriding `time integrator: initial time` (%.6e) with restart time (%.6e).",
            Float64(requested_initial_time),
            restart_time,
        )
    end
    params["restart_info"] = (
        index=restart_index,
        time=restart_time,
        displacement=displacement,
        velocity=velocity,
    )
    return nothing
end

# The list of `material: <name>: model:` strings actually referenced by a
# subdomain's `material: blocks:` mapping (block_name => material_name). Used
# by process_multidomain_restart!() to reject restart for subdomains using
# material models with internal state that the restart snapshot can't carry
# (currently just `j2 plasticity`), without having to open the mesh or build
# the full model. Mirrors the block -> material_name -> model resolution in
# SolidMechanics() (model.jl) and create_material() (constitutive.jl).
function _material_models_used(model_params::Parameters)
    material_params = get(model_params, "material", Parameters())
    blocks = get(material_params, "blocks", Parameters())
    return [
        get(get(material_params, material_name, Parameters()), "model", "")
        for material_name in values(blocks)
    ]
end

# Resolve and validate restart for a multi-domain (Schwarz) simulation. Run
# once, at the very start of MultiDomainSimulation() construction, before the
# top-level controller is built.
#
# Unlike single-domain restart, Schwarz restart is specified *once*, in the
# top-level (multi-domain) file's own `restart: index: N` block — not
# per-subdomain — because every subdomain must resume the coupled iteration
# at the same physical time, and duplicating the same index into every
# subdomain file would just invite drift (a forgotten file, a typo'd index).
# Subdomains keep supplying their own restart checkpoint the normal way,
# through `input mesh file` (still subdomain-local, since each subdomain has
# its own checkpoint mesh); only the `index` into that checkpoint is shared.
#
# This function peeks at each subdomain's checkpoint mesh just far enough to
# read the snapshot time at that shared index (as a consistency check —
# every subdomain should land on the same physical time, since Exodus output
# interval is forced equal across all subdomains) and requires those times
# to all agree. It then overwrites the top-level params["initial time"] with
# that restart time, and records the shared index on params so
# MultiDomainSimulation()'s per-domain loop can hand each subdomain a
# `restart: index: N` of its own — reusing process_restart!() unmodified for
# the actual per-subdomain work (reading displacement/velocity, the
# `j2 plasticity` check, building `restart_info`) once real subdomain
# construction happens.
#
# This has to happen before create_controller(params) builds the top-level
# Schwarz controller, because sync_control_time() unconditionally
# re-broadcasts the top-level controller's time fields (initial_time, time,
# prev_time, num_stops, stop) onto every subdomain at the start of every
# control step — so the top-level controller's initial time is the only one
# that actually has to be correct; whatever a subdomain's own controller
# computes is overwritten regardless.
#
# No-op when the top-level file has no `restart:` block. Aborts if a
# subdomain file has its own `restart:` block instead (restart belongs at
# the top level only, for Schwarz), if subdomains' checkpoints disagree on
# what time the shared index lands on, or if a restarting subdomain uses the
# `j2 plasticity` material model, for the same internal-state-variable reason
# monolithic restart rejects it (see SolidMechanics() in model.jl). Both FOM
# (`solid mechanics`) and ROM subdomains may restart, including mixed
# FOM-ROM pairings.
function process_multidomain_restart!(params::Parameters)
    domain_paths = params["domains"]
    for domain_path in domain_paths
        subparams = YAML.load_file(domain_path; dicttype=Parameters)
        if haskey(subparams, "restart")
            norma_abort(
                "`restart:` must be specified once, in the top-level multi-domain file, " *
                "not in individual subdomain files (found one in '$domain_path'). Every " *
                "subdomain restarts from the same physical time; move the `restart:` " *
                "block there and remove it from '$domain_path'.",
            )
        end
    end
    haskey(params, "restart") || return nothing
    restart_params = params["restart"]
    haskey(restart_params, "index") || norma_abort("`restart:` block must specify an `index`.")
    restart_index = _validate_restart_index(restart_params["index"])
    restart_time = nothing
    for domain_path in domain_paths
        subparams = YAML.load_file(domain_path; dicttype=Parameters)
        model_params = get(subparams, "model", Parameters())
        if "j2 plasticity" in _material_models_used(model_params)
            norma_abortf(
                "Schwarz restart does not support the `j2 plasticity` material model " *
                "(subdomain '%s'), for the same reason monolithic restart does not: the " *
                "restart snapshot only stores nodal displacement and velocity fields; J2 " *
                "plasticity's internal state variables (e.g. plastic strain, back stress) " *
                "are not written to or read from the restart file, so resuming would " *
                "silently discard the accumulated plastic history. Remove the `restart:` " *
                "block, or switch to a material model without internal state variables, " *
                "until restart support for internal variables is implemented.",
                domain_path,
            )
        end
        haskey(subparams, "input mesh file") ||
            norma_abort("Subdomain '$domain_path' has no `input mesh file`.")
        input_mesh_file = subparams["input mesh file"]
        input_mesh = Exodus.ExodusDatabase(input_mesh_file, "r")
        domain_time = try
            num_steps = Exodus.read_number_of_time_steps(input_mesh)
            if num_steps < 1
                norma_abort("Restart mesh '$input_mesh_file' (subdomain '$domain_path') contains no time steps.")
            end
            # Negative indices count back from the last snapshot, Python-style:
            # -1 is the final snapshot, -2 the one before it, and so on. Resolved
            # per subdomain against that subdomain's own checkpoint mesh, since
            # `restart_index` (possibly negative) is what gets propagated to each
            # subdomain's own `restart:` block later — see the comment on
            # `_multidomain_restart_index` below.
            domain_restart_index = restart_index < 0 ? num_steps + restart_index + 1 : restart_index
            if domain_restart_index < 1 || domain_restart_index > num_steps
                norma_abortf(
                    "Restart index %d is out of range for mesh '%s' (subdomain '%s'), " *
                    "which has %d time step(s) (valid range is 1:%d, or -1:-%d counting " *
                    "back from the last snapshot).",
                    restart_index,
                    input_mesh_file,
                    domain_path,
                    num_steps,
                    num_steps,
                    num_steps,
                )
            end
            Exodus.read_time(input_mesh, domain_restart_index)
        finally
            Exodus.close(input_mesh)
        end
        if restart_time === nothing
            restart_time = domain_time
        elseif !isapprox(domain_time, restart_time; atol=1.0e-9)
            norma_abortf(
                "Schwarz restart: subdomain '%s' checkpoint at index %d is at time %.6e, " *
                "which does not match the time %.6e read from an earlier subdomain's " *
                "checkpoint at the same index. All subdomains must restart from the same " *
                "physical time.",
                domain_path,
                restart_index,
                domain_time,
                restart_time,
            )
        end
    end
    # A restart at or past `final time` leaves nothing to integrate; see the
    # matching check in process_restart!() for why this must be rejected
    # rather than silently allowed to run one control step past final_time.
    final_time = get(params, "final time", nothing)
    if final_time !== nothing
        final_time = Float64(final_time)
        if restart_time > final_time || isapprox(restart_time, final_time; rtol=1e-9, atol=1e-12)
            norma_abortf(
                "Restart time %.6e (index %d) is at or past top-level `final time` %.6e; " *
                "there is nothing left to integrate. Choose an earlier restart index, or " *
                "increase `final time`.",
                restart_time,
                restart_index,
                final_time,
            )
        end
    end
    requested_initial_time = get(params, "initial time", nothing)
    params["initial time"] = restart_time
    params["_multidomain_restart_index"] = restart_index
    norma_log(0, :restart, "Restarting Schwarz simulation from snapshot data.")
    norma_logf(0, :restart, "Restart index:  %d", restart_index)
    norma_logf(0, :restart, "Restart time:   %.6e", restart_time)
    if requested_initial_time !== nothing && !isapprox(Float64(requested_initial_time), restart_time; atol=1e-12)
        norma_logf(
            0,
            :restart,
            "Overriding top-level `initial time` (%.6e) with restart time (%.6e).",
            Float64(requested_initial_time),
            restart_time,
        )
    end
    return nothing
end

# Warn (but do not abort) when a Schwarz restart is combined with
# non-overlapping Schwarz coupling (`Schwarz DN nonoverlap`, `Schwarz RR
# nonoverlap`, or `Schwarz impedance nonoverlap`). Nothing in process_restart!()
# or process_multidomain_restart!() is specific to overlapping vs.
# non-overlapping coupling — the positional nodal-field read is unaffected
# either way — but this combination has not been exercised by the test suite,
# so flag it for the user's attention rather than silently proceeding.
# Called from create_simulation() after create_bcs(), since boundary
# conditions are not populated until then.
function warn_restart_with_nonoverlap_schwarz(sim::MultiDomainSimulation)
    haskey(sim.params, "restart") || return nothing
    # Check the geometry-axis abstract type (overlap vs. non-overlap), not the
    # concrete `SolidMechanicsNonOverlapSchwarzBoundaryCondition` (the "Schwarz
    # DN nonoverlap" type): "Schwarz RR nonoverlap" and "Schwarz impedance
    # nonoverlap" both build SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition,
    # a sibling concrete type under the same
    # SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition abstract type that
    # the earlier concrete-type check never matched.
    has_nonoverlap = any(
        any(bc isa SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition for bc in sub.model.boundary_conditions)
        for sub in sim.subsims
    )
    if has_nonoverlap
        norma_log(
            0,
            :warning,
            "Restart is being used with non-overlapping Schwarz coupling. This should " *
            "work but has not been tested; please use with caution.",
        )
    end
    return nothing
end

warn_restart_with_nonoverlap_schwarz(::SingleDomainSimulation) = nothing

# The restart refinement pass (see the Mfc correction in time_integrator.jl)
# runs before detect_contact() re-establishes active_contact for the
# restarted step, so apply_bc skips every Schwarz contact BC during that
# pass and it silently keeps whatever stale interface acceleration was in
# the checkpoint. Restart + Schwarz contact has not been tested against
# this gap, so reject the combination outright rather than risk a quietly
# wrong result.
function reject_restart_with_contact(sim::MultiDomainSimulation)
    haskey(sim.params, "restart") || return nothing
    if sim.controller.schwarz_contact
        norma_abort(
            "Restart is not supported in combination with `Schwarz contact` boundary " *
            "conditions: the restart refinement pass runs before contact is " *
            "re-detected for the restarted step, so it would silently reuse the stale " *
            "interface acceleration recorded in the checkpoint instead of the correct " *
            "contact state. This combination is untested; remove the `restart:` block " *
            "or remove the `Schwarz contact` boundary conditions.",
        )
    end
    return nothing
end

reject_restart_with_contact(::SingleDomainSimulation) = nothing

function create_simulation(params::Parameters)
    sim_type = params["type"]
    if sim_type == "single"
        sim = SingleDomainSimulation(params)
        create_bcs(sim)
        initialize_storage(sim)
        return sim
    elseif sim_type == "multi"
        sim = MultiDomainSimulation(params)
        create_bcs(sim)
        validate_swap_criteria(sim)
        warn_restart_with_nonoverlap_schwarz(sim)
        reject_restart_with_contact(sim)
        initialize_storage(sim)
        return sim
    else
        norma_abort("Unknown type of simulation: $sim_type")
    end
end

function initialize_storage(sim::SingleDomainSimulation)
    model = sim.model
    integrator = sim.integrator
    solver = sim.solver
    model isa SolidMechanics || return nothing
    num_dof = length(model.free_dofs)
    integrator.displacement = unsafe_wrap(Vector{Float64}, pointer(model.displacement), num_dof)
    integrator.velocity = unsafe_wrap(Vector{Float64}, pointer(model.velocity), num_dof)
    integrator.acceleration = unsafe_wrap(Vector{Float64}, pointer(model.acceleration), num_dof)
    if solver isa ExplicitSolver
        solver.solution = unsafe_wrap(Vector{Float64}, pointer(model.acceleration), num_dof)
    else
        solver.solution = unsafe_wrap(Vector{Float64}, pointer(model.displacement), num_dof)
    end
    return nothing
end

function initialize_storage(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        initialize_storage(subsim)
    end
    return nothing
end

function SingleDomainSimulation(params::Parameters)
    t_setup = time()
    basename = params["name"]
    _reject_restart_with_swaps!(params, basename)
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    # `output mesh file` is unconditionally rm'd and then reopened fresh below
    # (as a copy of `input mesh file`) before the checkpoint is ever read. If
    # the two paths coincide -- e.g. a restart run accidentally configured to
    # write back over its own input checkpoint -- that rm destroys the only
    # copy of the data this run needs to read, with nothing left to restore
    # from on failure. Abort instead of running in place.
    if abspath(input_mesh_file) == abspath(output_mesh_file)
        norma_abortf(
            "`input mesh file` and `output mesh file` resolve to the same path ('%s'). " *
            "This would delete the input mesh before it can be read (it is rm'd and " *
            "recreated as a copy of itself) -- most dangerous on restart, where the " *
            "input mesh is often the only copy of the checkpoint. Use a different " *
            "`output mesh file`.",
            abspath(output_mesh_file),
        )
    end
    norma_log(0, :setup, "Input:  $input_mesh_file")
    norma_log(0, :setup, "Output: $output_mesh_file")
    rm(output_mesh_file; force=true)
    input_mesh = Exodus.ExodusDatabase(input_mesh_file, "r")
    local output_mesh
    try
        Exodus.copy(input_mesh, output_mesh_file)
        output_mesh = Exodus.ExodusDatabase(output_mesh_file, "rw")
    catch
        # Close the input handle so a subsequent simulation can re-open the
        # mesh files cleanly even if this construction failed mid-flight.
        try; Exodus.close(input_mesh); catch; end
        rethrow()
    end
    params["output_mesh"] = output_mesh
    params["input_mesh"] = input_mesh
    try
        n_nodes = Exodus.num_nodes(input_mesh.init)
        n_elems = Exodus.num_elements(input_mesh.init)
        norma_logf(0, :setup, "Mesh:   %d nodes, %d elements", n_nodes, n_elems)
        process_restart!(params, input_mesh, basename)
        controller = create_controller(params)
        norma_logf(0, :setup, "Time:   [%.2e, %.2e], Δt = %.2e, %d steps",
                   controller.initial_time, controller.final_time,
                   controller.time_step, controller.num_stops - 1)
        model = create_model(params)
        # The restart snapshot's displacement/velocity fields (if any) were
        # copied into model.displacement/model.velocity -- and, for a ROM
        # model, again into model.fom_model.displacement/velocity -- inside
        # SolidMechanics() above (model.jl); nothing reads them from
        # params["restart_info"] from here on. Only the fact that this run
        # *is* a restart still matters downstream: apply_ics(::RomModel, ...)
        # (opinf_ics_bcs.jl) checks `restart_info !== nothing` to decide
        # whether to project the restored FOM state onto the reduced basis.
        # Replace the two full 3xN snapshot matrices with a lightweight
        # marker instead of letting them sit alive, unused, in `params` (and
        # therefore in `sim.params`) for the rest of the run.
        if params["restart_info"] !== nothing
            params["restart_info"] = (index=params["restart_info"].index, time=params["restart_info"].time)
        end
        _log_materials(model, input_mesh)
        integrator = create_time_integrator(params, model)
        solver = create_solver(params, model)
        norma_logf(0, :setup, "Solver: %s, %s", _integrator_name(integrator), _solver_name(solver))
        norma_log(0, :setup, "Setup complete ($(format_time(time() - t_setup)))")
        failed = false
        swaps = parse_swap_plans(params)
        sim = SingleDomainSimulation(basename, params, controller, integrator, solver, model, failed, nothing, nothing, swaps)
        validate_swap_plans(sim)
        return sim
    catch
        # Close both Exodus handles on any failure after they were opened.
        try; Exodus.close(input_mesh); catch; end
        try; Exodus.close(output_mesh); catch; end
        rethrow()
    end
end

# True when this subsim belongs to a MultiDomainSimulation (i.e., is coupled).
is_coupled(sim::SingleDomainSimulation) = sim.parent !== nothing

# Resolve a domain name to its stable handle within a MultiDomainSimulation.
handle(sim::MultiDomainSimulation, name::String) = sim.handle_by_name[name]

_integrator_name(::QuasiStatic) = "Quasi-static"
_integrator_name(::Newmark) = "Newmark"
_integrator_name(::CentralDifference) = "Central difference"
_integrator_name(ig) = replace(string(typeof(ig)), r"^.*\." => "")

_solver_name(::HessianMinimizer) = "Newton (Hessian minimizer)"
_solver_name(::SteepestDescent) = "Steepest descent"
_solver_name(::ExplicitSolver) = "Explicit"
_solver_name(s) = replace(string(typeof(s)), r"^.*\." => "")

function log_dof_counts(model::SolidMechanics)
    n_dofs = length(model.free_dofs)
    n_free = count(model.free_dofs)
    norma_logf(0, :setup, "DOFs:   %d total, %d free, %d constrained",
               n_dofs, n_free, n_dofs - n_free)
end

log_dof_counts(::Model) = nothing

_material_name(m) = replace(string(typeof(m)), r"^.*\." => "")

# Format a single material struct field as "name = value".  Floats get
# scientific notation for readability; integers (e.g. SethHill's m, n)
# print as-is.
function _format_material_field(name::Symbol, value)
    if value isa AbstractFloat
        return @sprintf("%s = %.3e", String(name), value)
    else
        return string(name) * " = " * string(value)
    end
end

function _log_materials(model::SolidMechanics, input_mesh)
    block_names = Exodus.read_names(input_mesh, Block)
    for (block_name, material) in zip(block_names, model.materials)
        props = join((_format_material_field(f, getfield(material, f))
                      for f in fieldnames(typeof(material))), ", ")
        norma_log(0, :setup, "Block \"$block_name\": $(_material_name(material)) : $props")
    end
end

_log_materials(::Model, _) = nothing

function MultiDomainSimulation(params::Parameters)
    basename = params["name"]
    domain_paths = params["domains"]
    process_multidomain_restart!(params)
    _reject_restart_with_swaps!(params, basename)
    subsims = Vector{SingleDomainSimulation}()
    controller = create_controller(params)
    initial_time = controller.initial_time
    final_time = controller.final_time
    integrator_dt = params["time step"]
    exodus_interval = Float64(get(params, "Exodus output interval", integrator_dt))
    csv_interval = Float64(get(params, "CSV output interval", 0.0))
    handle_by_name = Dict{String,DomainHandle}()
    name_by_handle = String[]
    subsim_index = 1
    try
        for domain_path in domain_paths
            domain_name = stripped_name(domain_path)
            norma_log(4, :domain, domain_name)
            subparams = YAML.load_file(domain_path; dicttype=Parameters)
            subparams["name"] = domain_name
            if haskey(params, "_multidomain_restart_index")
                subparams["restart"] = Parameters("index" => params["_multidomain_restart_index"])
            end
            integrator_params = subparams["time integrator"]
            subsim_time_step = get(integrator_params, "time step", integrator_dt)
            subsim_time_step = min(subsim_time_step, integrator_dt)
            integrator_params["initial time"] = initial_time
            integrator_params["final time"] = final_time
            integrator_params["time step"] = subsim_time_step
            subparams["Exodus output interval"] = exodus_interval
            subparams["CSV output interval"] = csv_interval
            subsim = SingleDomainSimulation(subparams)
            params[domain_name] = subsim.params
            push!(subsims, subsim)
            handle_by_name[domain_name] = DomainHandle(subsim_index)
            push!(name_by_handle, domain_name)
            subsim_index += 1
        end
    catch
        # Close any Exodus handles opened by already-constructed subsims so
        # the failure does not leak file descriptors into later simulations.
        for already in subsims
            try; finalize_writing(already); catch; end
        end
        rethrow()
    end
    num_domains = length(subsims)
    swaps = parse_swap_plans(params)
    failed = false
    sim = MultiDomainSimulation(
        basename, params, controller, num_domains, subsims,
        handle_by_name, name_by_handle, swaps, failed,
    )
    for (i, subsim) in enumerate(sim.subsims)
        subsim.parent = sim
        subsim.handle = DomainHandle(i)
    end
    validate_swap_plans(sim)
    return sim
end

function SolidMultiDomainTimeController(params::Parameters)
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["absolute tolerance"]
    relative_tolerance = params["relative tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    # Multi-domain: control step = integrator dt (Schwarz coupling each step).
    # Output timing is handled by _is_output_time in write_stop.
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
    absolute_error = relative_error = 0.0
    time = prev_time = initial_time
    stop = 0
    converged = false
    iteration_number = 0
    stop_disp = [Vector{Float64}[] for _ in 1:num_domains]
    stop_velo = [Vector{Float64}[] for _ in 1:num_domains]
    stop_acce = [Vector{Float64}[] for _ in 1:num_domains]
    stop_∂Ω_f = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_disp = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_velo = [Vector{Float64}[] for _ in 1:num_domains]
    schwarz_acce = [Vector{Float64}[] for _ in 1:num_domains]
    time_hist = [Float64[] for _ in 1:num_domains]
    disp_hist = [Vector{Float64}[] for _ in 1:num_domains]
    velo_hist = [Vector{Float64}[] for _ in 1:num_domains]
    acce_hist = [Vector{Float64}[] for _ in 1:num_domains]
    ∂Ω_f_hist = [Vector{Float64}[] for _ in 1:num_domains]
    has_relaxation_key = haskey(params, "relaxation")
    has_relaxation_parameter = haskey(params, "relaxation parameter")
    has_aitken_N0_parameter = haskey(params, "aitken N0 parameter")
    relaxation_method = :fixed
    relaxation_parameter = 1.0
    aitken_N0 = 1
    if has_relaxation_key
        relaxation_value = params["relaxation"]
        relaxation_string = relaxation_value isa AbstractString ? lowercase(relaxation_value) : ""
        if relaxation_string == "aitken recursive"
            # Recursive Irons-Tuck Aitken, carrying θ^(n-1) forward.
            relaxation_method = :aitken_recursive
            if has_aitken_N0_parameter
                aitken_N0 = Int(params["aitken N0 parameter"])
            end
        elseif relaxation_string == "aitken secant"
            # Aitken secant, non-recursive (Sambataro-Tezaur eq. 9 / Deparis-
            # Discacciati-Quarteroni), the original-paper relaxation factor.
            relaxation_method = :aitken_secant
            if has_aitken_N0_parameter
                aitken_N0 = Int(params["aitken N0 parameter"])
            end
        else
            norma_abort("Schwarz controller: unsupported `relaxation: $(relaxation_value)` (only `aitken recursive` and `aitken secant` are recognized).")
        end
    end
    if has_relaxation_parameter
        relaxation_parameter = Float64(params["relaxation parameter"])
    end
    # Echo the effective relaxation settings once, so a misspelled or misplaced
    # YAML key (silently ignored, like all unrecognized keys) is visible at startup.
    if relaxation_method === :fixed
        norma_logf(0, :schwarz, "Relaxation: %s, θ = %.4e",
            relaxation_method_name(relaxation_method), relaxation_parameter)
    else
        norma_logf(0, :schwarz, "Relaxation: %s, θ₀ = %.4e, N0 = %d",
            relaxation_method_name(relaxation_method), relaxation_parameter, aitken_N0)
    end
    naive_stabilized = get(params, "naive stabilized", false)
    # Keyed per interface and filled on demand, as the interfaces are known to
    # the boundary conditions, not to the controller (see RelaxationKey).
    lambda_time = Dict{RelaxationKey,Vector{Float64}}()
    lambda_disp = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    lambda_velo = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    lambda_acce = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    aitken_prev_residual_disp = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    aitken_prev_residual_velo = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    aitken_prev_residual_acce = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    aitken_theta_disp = Dict{RelaxationKey,Vector{Float64}}()
    aitken_prev_lambda_disp = Dict{RelaxationKey,Vector{Vector{Float64}}}()
    is_schwarz = true
    schwarz_contact = false
    active_contact = false
    contact_hist = Vector{Bool}[]
    schwarz_iters = zeros(Int64, num_stops - 1)
    use_interface_predictor = get(params, "interface predictor", false)
    predictor_disp = [Vector{Float64}() for _ in 1:num_domains]
    predictor_velo = [Vector{Float64}() for _ in 1:num_domains]
    predictor_acce = [Vector{Float64}() for _ in 1:num_domains]
    predictor_∂Ω_f = [Vector{Float64}() for _ in 1:num_domains]
    prev_stop_disp = [Vector{Float64}() for _ in 1:num_domains]
    prev_stop_∂Ω_f = [Vector{Float64}() for _ in 1:num_domains]

    return SolidMultiDomainTimeController(
        minimum_iterations,
        maximum_iterations,
        absolute_tolerance,
        relative_tolerance,
        absolute_error,
        relative_error,
        initial_time,
        final_time,
        time_step,
        time,
        prev_time,
        num_stops,
        stop,
        converged,
        iteration_number,
        stop_disp,
        stop_velo,
        stop_acce,
        stop_∂Ω_f,
        schwarz_disp,
        schwarz_velo,
        schwarz_acce,
        time_hist,
        disp_hist,
        velo_hist,
        acce_hist,
        ∂Ω_f_hist,
        relaxation_parameter,
        relaxation_method,
        aitken_N0, 
        naive_stabilized,
        lambda_time,
        lambda_disp,
        lambda_velo,
        lambda_acce,
        aitken_prev_residual_disp,
        aitken_prev_residual_velo,
        aitken_prev_residual_acce,
        aitken_theta_disp,
        aitken_prev_lambda_disp,
        is_schwarz,
        schwarz_contact,
        active_contact,
        contact_hist,
        schwarz_iters,
        use_interface_predictor,
        predictor_disp,
        predictor_velo,
        predictor_acce,
        predictor_∂Ω_f,
        prev_stop_disp,
        prev_stop_∂Ω_f,
        false,   # relaxation_frozen, set per sweep
    )
end

function SolidSingleDomainTimeController(params::Parameters)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    # Control step = integrator dt always. Output timing is handled
    # by _is_output_time in write_stop, not by controller stops.
    num_stops = max(round(Int64, (final_time - initial_time) / time_step) + 1, 2)
    time = prev_time = initial_time
    stop = 0
    return SolidSingleDomainTimeController(initial_time, final_time, time_step, time, prev_time, num_stops, stop)
end

function create_controller(params::Parameters)
    sim_type = params["type"]
    if sim_type == "single"
        return SolidSingleDomainTimeController(params)
    elseif sim_type == "multi"
        return SolidMultiDomainTimeController(params)
    else
        norma_abort("Unknown type of simulation: $sim_type")
    end
end

function create_bcs(sim::SingleDomainSimulation)
    boundary_conditions = _create_bcs(sim)
    return sim.model.boundary_conditions = boundary_conditions
end

function create_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        create_bcs(subsim)
    end
    pair_schwarz_bcs(sim)
    return nothing
end

function evolve(sim::Simulation)
    sync_control_time(sim)
    initialize(sim)
    initialize_writing(sim)
    write_stop(sim)
    t_batch = time()
    while true
        advance_control_time(sim)
        maybe_apply_swaps!(sim)
        sync_control_time(sim)
        advance_control(sim)
        write_stop(sim; wall_time=time() - t_batch)
        t_batch = time()
        if stop_evolve(sim) == true
            break
        end
    end
    return nothing
end

function stop_evolve(sim::Simulation)
    time = sim.controller.time
    final_time = sim.controller.final_time
    return isapprox(time, final_time; rtol=1e-9, atol=1e-12) || time > final_time
end

function solve_contact(sim::MultiDomainSimulation)
    if sim.controller.active_contact == true
        schwarz(sim)
    else
        advance_independent(sim)
    end
end

function decrease_time_step(sim::SingleDomainSimulation)
    norma_log(0, :recover, "Attempting to recover by decreasing time step...")
    time_step = sim.integrator.time_step
    decrease_factor = sim.integrator.decrease_factor
    if decrease_factor == 1.0
        norma_abortf(
            "Cannot adapt time step %.2e because decrease factor is %.2e. Enable adaptive time stepping.",
            time_step,
            decrease_factor,
        )
    end
    new_time_step = decrease_factor * time_step
    minimum_time_step = sim.integrator.minimum_time_step
    if new_time_step < minimum_time_step
        norma_abortf("Cannot adapt time step to %.2e because minimum is %.2e.", new_time_step, minimum_time_step)
    end
    sim.integrator.time_step = new_time_step
    return nothing
end

function increase_time_step(sim::SingleDomainSimulation)
    increase_factor = sim.integrator.increase_factor
    if increase_factor > 1.0
        time_step = sim.integrator.time_step
        maximum_time_step = sim.integrator.maximum_time_step
        new_time_step = min(increase_factor * time_step, maximum_time_step)
        if new_time_step > time_step
            sim.integrator.time_step = new_time_step
        end
    end
    return nothing
end

function advance_one_step(sim::SingleDomainSimulation)
    is_explicit = sim.integrator isa ExplicitDynamicTimeIntegrator
    while true
        prev_time = sim.integrator.prev_time
        time = sim.integrator.time
        time_step = sim.integrator.time_step
        # For implicit, print each step's time interval (Newton iters follow).
        # For explicit, suppress — only output stops are printed.
        if !is_explicit
            norma_logf(4, :advance, "Time = [%.2e, %.2e] : Δt = %.2e", prev_time, time, time_step)
        end
        apply_bcs(sim)
        # Defense-in-depth: constitutive model evaluate() catches math errors
        # and sets model.failed.  This outer catch handles the same errors if
        # they escape from any other point in the solve path.
        try
            solve(sim)
        catch e
            e isa _MATH_ERRORS || rethrow()
            norma_logf(4, :solve, "Caught %s during solve — treating as step failure", typeof(e))
            sim.failed = sim.model.failed = sim.solver.failed = true
        end
        if sim.failed == false
            increase_time_step(sim)
            save_curr_state(sim)
            break
        end
        restore_prev_state(sim)
        decrease_time_step(sim)
        advance_time(sim)
        sim.failed = sim.model.failed = sim.solver.failed = false
    end
    return nothing
end

function advance_control(sim::SingleDomainSimulation)
    subcycle(sim)
    return nothing
end

function advance_control(sim::MultiDomainSimulation)
    if sim.controller.schwarz_contact == false
        schwarz(sim)
        return nothing
    end
    save_stop_state(sim)
    solve_contact(sim)
    was_in_contact = sim.controller.active_contact
    detect_contact(sim)
    if sim.controller.active_contact ≠ was_in_contact
        if was_in_contact == true
            norma_log(0, :contact, "Released — reattempting control step.")
        else
            norma_log(0, :contact, "Initiated — reattempting control step.")
        end
        restore_stop_state(sim)
        solve_contact(sim)
    end
    return nothing
end

function apply_ics(sim::SingleDomainSimulation)
    apply_ics(sim.params, sim.model, sim.integrator, sim.solver)
    return nothing
end

function apply_ics(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_ics(subsim)
    end
end

function apply_bcs(sim::SingleDomainSimulation)
    apply_bcs(sim.model)
    return nothing
end

function apply_bcs(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        apply_bcs(subsim)
    end
end

function initialize(sim::SingleDomainSimulation)
    apply_ics(sim)
    apply_bcs(sim)
    log_dof_counts(sim.model)
    initialize(sim.integrator, sim.solver, sim.model)
    save_curr_state(sim)
    return nothing
end

function initialize(sim::MultiDomainSimulation)
    initialize_bc_projectors(sim)
    apply_ics(sim)
    # This first snapshot (used only as apply_bcs()'s fallback while every
    # subdomain still has a placeholder zero acceleration) is pushed for FOM
    # subdomains only, matching existing non-restart behavior; ROM neighbors
    # fall back to zero here exactly as they always have (see apply_bc() in
    # schwarz.jl). The restart-only refinement pass below is what actually
    # gives every subdomain — FOM and ROM alike — a real Schwarz-coupled
    # snapshot to interpolate from.
    for (subsim_index, subsim) in enumerate(sim.subsims)
        if subsim.model isa SolidMechanics
            save_history_snapshot(sim.controller, subsim, subsim_index)
        end
    end
    apply_bcs(sim)
    for subsim in sim.subsims
        log_dof_counts(subsim.model)
        initialize(subsim.integrator, subsim.solver, subsim.model)
        save_curr_state(subsim)
    end
    # Restart-only refinement: the pass above computes each subdomain's own
    # initial acceleration while every Schwarz-coupled boundary DOF still
    # holds the placeholder (zero) value apply_bcs() left before any
    # subdomain had a real acceleration — see the comment in
    # initialize(::Newmark, ...) (time_integrator.jl) for why that's the
    # correct thing to do on a fresh, at-rest-ish start. On a restart,
    # though, the simulation is resuming mid-motion and the true
    # acceleration at a Schwarz interface generally isn't small, so that
    # placeholder can be badly wrong right where it matters most. This
    # applies equally to FOM and ROM subdomains: a ROM subdomain's own
    # initial acceleration is solved on its internal fom_model (see
    # initialize(::RomNewmark, ...) / initialize(::RomCentralDifference, ...)
    # in opinf_time_integrator.jl), which is restored from the restart
    # snapshot exactly like a standalone FOM model (is_restarted() below
    # checks subsim.model.fom_model.restarted for a ROM subsim).
    #
    # apply_bc() for a Schwarz-coupled BC (schwarz.jl) does not read the
    # coupled side's *live* integrator state — it interpolates from
    # controller.{time,disp,velo,acce,∂Ω_f}_hist[coupled_index], populated by
    # save_history_snapshot(). The very first entry in that history was
    # pushed above, before any subdomain had a real acceleration, so it's
    # exactly as stale as the placeholder apply_bcs() left; simply calling
    # apply_bcs() again would still interpolate from that same one stale
    # entry (interpolate() with a single history point ignores the query
    # time and just returns it) and be a no-op. So: reset each subdomain's
    # history and re-push a fresh snapshot reflecting its now-real,
    # pass-1-solved state first, *then* refresh apply_bcs and re-solve, this
    # time trusting the Schwarz-coupled DOFs.
    #
    # Run this reset/apply/re-solve sequence twice (not once) when a ROM
    # subdomain is involved. For a FOM neighbor, model.displacement and
    # integrator.displacement are the same aliased array, so as soon as
    # apply_bc() writes a real value anywhere, every reader sees it
    # immediately, in any subdomain order. A ROM neighbor's full-order
    # fields are a *separate* reconstruction from its reduced state
    # (reconstruct_fom_fields!() in opinf_solver.jl), and that
    # reconstruction skips any DOF where the ROM's own fom_model.free_dofs
    # is false — including the DOFs where the ROM's own Schwarz coupling
    # prescribes values from a *different* neighbor. Whether those DOFs
    # hold real or still-placeholder data by the time some other subdomain
    # reads them within a single apply_bcs(sim) sweep depends on whether
    # that ROM's own coupling BC happened to already run earlier in the
    # same sweep — i.e. on sim.subsims order, not physics. A second
    # sweep removes that order-dependence: every subdomain enters it with
    # a real (round-1-refined) acceleration already in place, so every
    # reconstruction this time reads real data no matter which subdomain
    # is processed first.
    # A fresh (non-restart) start whose initial conditions put nonzero
    # displacement or velocity on an impedance-Schwarz-coupled interface has
    # the same consistency problem as a restart: the pass-1 accelerations
    # were solved against placeholder (zero) partner data, so the interface
    # force balance at t = 0 is one-sided, and the unbalanced Robin term
    # α W u produces a one-step impulsive kick (measured: interface-
    # localized kinetic energy scaling as α²Δt² from a stress-free rotated
    # state). The same refinement pass restores the balance; it is a no-op
    # for the common at-rest start, which is skipped to preserve existing
    # behavior exactly. Fresh starts of other Schwarz couplings (DN,
    # overlap) are excluded — see has_initial_interface_motion below.
    if any(is_restarted(subsim.model) for subsim in sim.subsims) ||
        any(has_initial_interface_motion(subsim) for subsim in sim.subsims)
        for _ in 1:2
            for (subsim_index, subsim) in enumerate(sim.subsims)
                reset_history(sim.controller, subsim_index)
                save_history_snapshot(sim.controller, subsim, subsim_index)
            end
            apply_bcs(sim)
            for subsim in sim.subsims
                initialize(subsim.integrator, subsim.solver, subsim.model; trust_schwarz=true)
                save_curr_state(subsim)
            end
        end
    end
    detect_contact(sim)
    return nothing
end

# True when any impedance-Schwarz-coupled DOF of this subdomain carries
# nonzero initial displacement or velocity — the condition under which the
# t = 0 acceleration refinement pass in initialize(sim::MultiDomainSimulation)
# is needed for a fresh start. Restricted to the impedance nonoverlap
# exchange, the only coupling the pass has been validated on: the unbalanced
# Robin term α W u is what produces the one-step kick the pass removes, and
# the DN d-form exchange demonstrably does not tolerate a trust-Schwarz
# re-solve from a moving initial state (element inversion on the DN
# cantilever with a parabolic initial bend).
function has_initial_interface_motion(subsim::SingleDomainSimulation)
    model = get_fom_model(subsim)
    model isa SolidMechanics || return false
    for bc in model.boundary_conditions
        bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition || continue
        for i_global in bc.global_from_local_map
            for comp in 1:3
                if model.displacement[comp, i_global] != 0.0 ||
                    model.velocity[comp, i_global] != 0.0
                    return true
                end
            end
        end
    end
    return false
end

function solve(sim::SingleDomainSimulation)
    solve(sim.integrator, sim.solver, sim.model)
    return sim.failed = sim.model.failed
end

function initialize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        initialize_writing(subsim)
    end
end

function finalize_writing(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        finalize_writing(subsim)
    end
end

function sync_control_time(sim::SingleDomainSimulation)
    set_initial_subcycle_time(sim)
    return nothing
end

function sync_control_time(sim::MultiDomainSimulation)
    set_initial_subcycle_time(sim)
    controller = sim.controller
    initial_time = controller.initial_time
    final_time = controller.final_time
    time_step = controller.time_step
    time = controller.time
    prev_time = controller.prev_time
    num_stops = controller.num_stops
    stop = controller.stop
    for subsim in sim.subsims
        subsim.controller.initial_time = initial_time
        subsim.controller.final_time = final_time
        subsim.controller.time_step = time_step
        subsim.controller.time = time
        subsim.controller.prev_time = prev_time
        subsim.controller.num_stops = num_stops
        subsim.controller.stop = stop
    end
    return nothing
end

function set_initial_subcycle_time(sim::SingleDomainSimulation)
    time = sim.controller.prev_time
    sim.model.time = sim.integrator.time = time
    return nothing
end

function set_initial_subcycle_time(sim::MultiDomainSimulation)
    time = sim.controller.prev_time
    for subsim in sim.subsims
        subsim.model.time = subsim.integrator.time = time
    end
    return nothing
end

function get_adjusted_timestep(t::Float64, dt::Float64, t_stop::Float64)::Float64
    gap = t_stop - t
    gap <= 0.0 && return 0.0
    t_next = t + dt
    return t_stop - t_next <= 0.5 * dt ? gap : dt
end

function advance_time(sim::SingleDomainSimulation)
    controller = sim.controller
    integrator = sim.integrator
    time = integrator.time
    time_stop = controller.time
    time_step = get_adjusted_timestep(time, integrator.time_step, time_stop)
    next_time = time + time_step
    integrator.prev_time = time
    integrator.time = sim.model.time = next_time
    integrator.time_step = time_step
    return nothing
end

function advance_control_time(sim::Simulation)
    controller = sim.controller
    controller.prev_time = controller.time
    stop = controller.stop + 1
    initial_time = controller.initial_time
    time_step = controller.time_step
    next_time = stop * time_step + initial_time
    controller.time = next_time
    controller.stop = stop
    return nothing
end

function advance_independent(sim::MultiDomainSimulation)
    sim.controller.iteration_number = 0
    sim.controller.is_schwarz = false
    save_stop_state(sim)
    set_initial_subcycle_time(sim)
    subcycle(sim)
    return nothing
end

function schwarz(sim::MultiDomainSimulation)
    iteration_number = 0
    sim.controller.is_schwarz = true
    save_stop_state(sim)
    save_schwarz_state(sim)
    reset_histories(sim)
    reset_relaxation_state!(sim.controller)
    swap_swappable_bcs(sim)

    if sim.controller.use_interface_predictor
        compute_interface_predictor!(sim)
    end

    # Interface-jump gate for adjoint-paired impedance interfaces: their slow
    # jump mode contributes almost nothing to ΔU while still far from the
    # fixed point, and the dashpot dissipates whatever jump the iteration
    # leaves behind (measured: −16% of a wave packet crossing a conforming
    # interface at the 1.0e-8 default tolerance) — see paired_impedance_jump
    # (schwarz.jl). The gate holds convergence only while the jump is actually
    # contracting: a jump mode with no dashpot authority (e.g. a quiescent
    # interface) can stall above the tolerance, and holding then just rides
    # the sweep cap without improving the answer, so a stalled jump is
    # accepted with a warning instead.
    prev_jump_rel = -1.0
    while true
        norma_log(0, :schwarz, "Iteration [$iteration_number]")
        sim.controller.iteration_number = iteration_number
        # Recorded by the coupling BCs as they apply their relaxation factors.
        sim.controller.relaxation_frozen = false
        set_initial_subcycle_time(sim)
        subcycle(sim)
        ΔU, Δu = update_schwarz_convergence_criterion(sim)
        if sim.controller.converged && sim.controller.relaxation_frozen
            # A relaxation factor of essentially zero left an interface iterate
            # unchanged, so every subdomain re-solved against the data it
            # already had and returned the solution it already had. The small
            # update that follows measures the frozen blend, not the interface
            # residual, and must not be read as convergence.
            norma_log(
                0,
                :schwarz,
                "Relaxation factor near zero froze an interface iterate; update is not convergence.",
            )
            sim.controller.converged = false
        end
        if sim.controller.converged
            jump_rel = paired_impedance_jump(sim)
            if jump_rel > sim.controller.relative_tolerance
                if prev_jump_rel ≥ 0.0 && jump_rel > 0.95 * prev_jump_rel
                    # A stalled jump is a representability floor of the pair's
                    # trace spaces: no further sweeping reduces it, and the
                    # dashpot dissipates it at a rate the tolerance no longer
                    # bounds. The default accepts it loudly; "abort" is for
                    # runs where that dissipation is unacceptable and the
                    # remedy (mesh conformity or refinement) is on the user.
                    stall_action = get(sim.params, "stalled interface jump action", "warn")
                    if stall_action == "abort"
                        norma_abort(
                            "Impedance interface jump $(jump_rel) stalled above " *
                            "tolerance $(sim.controller.relative_tolerance) and " *
                            "`stalled interface jump action: abort` is set. The " *
                            "jump is a representability limit of the interface " *
                            "trace spaces (nonconforming meshes); the dashpot " *
                            "dissipates it every stop. Refine toward conformity " *
                            "or accept the loss with the default `warn`.",
                        )
                    elseif stall_action != "warn"
                        norma_abort(
                            "Unknown `stalled interface jump action: $(stall_action)`. " *
                            "Valid values are `warn` (default) and `abort`.",
                        )
                    end
                    norma_logf(
                        0,
                        :schwarz,
                        "Impedance interface jump %.2e stalled above tolerance %.2e; accepting.",
                        jump_rel,
                        sim.controller.relative_tolerance,
                    )
                else
                    norma_logf(
                        0,
                        :schwarz,
                        "Impedance interface jump %.2e > %.2e holds convergence.",
                        jump_rel,
                        sim.controller.relative_tolerance,
                    )
                    sim.controller.converged = false
                end
            end
            prev_jump_rel = jump_rel
        end
        if iteration_number == 0
            # Initial Schwarz pass: there is no prior iterate to compare against,
            # so the relative criterion is not yet a Schwarz convergence measure.
            # Report the absolute update for reference and skip the relative test.
            sim.controller.converged = false
            norma_logf(
                0,
                :schwarz,
                "Criterion [%d] |ΔU| = %.2e : |ΔU|/|U| = %8s : %s",
                iteration_number,
                ΔU,
                "—",
                colored_status("[INITIAL]"),
            )
            # Early exit: if the initial absolute update already meets the absolute
            # tolerance, no further Schwarz iterations are needed. The relative test
            # still cannot be applied on iteration 0 (no prior iterate). Paired
            # impedance interfaces must also clear the jump gate here, since this
            # path bypasses controller.converged entirely.
            if ΔU ≤ sim.controller.absolute_tolerance &&
                paired_impedance_jump(sim) ≤ sim.controller.relative_tolerance
                norma_log(0, :schwarz, "Performed 0 Schwarz Iterations")
                sim.controller.schwarz_iters[sim.controller.stop] = 0
                break
            end
        else
            raw_status = sim.controller.converged ? "[CONVERGED]" : "[CONVERGING]"
            status = colored_status(raw_status)
            norma_logf(
                0,
                :schwarz,
                "Criterion [%d] |ΔU| = %.2e : |ΔU|/|U| = %.2e : %s",
                iteration_number,
                ΔU,
                Δu,
                status,
            )
            if stop_schwarz(sim, iteration_number + 1) == true
                plural = iteration_number == 1 ? "" : "s"
                norma_log(0, :schwarz, "Performed $iteration_number Schwarz Iteration" * plural)
                sim.controller.schwarz_iters[sim.controller.stop] = iteration_number
                sim.controller.converged || report_unconverged_step(sim, iteration_number)
                report_overlap_l2_errors(sim)
                break
            end
        end
        iteration_number += 1
        save_schwarz_state(sim)
        restore_stop_state(sim)
    end
    report_impedance_interface_work(sim)
end

function report_overlap_l2_errors(sim::MultiDomainSimulation)
    overlap_rows = Vector{Vector{Any}}()
    for subsim in sim.subsims
        for bc in subsim.model.boundary_conditions
            if !stores_overlap_l2_error(bc)
                continue
            end
            overlap_l2_error = compute_overlap_l2_error!(bc)
            field = overlap_l2_error_field(bc)
            write_overlap_l2_error_screen(subsim.name, bc.name, field, overlap_l2_error)
            push!(overlap_rows, Any[subsim.name, bc.name, overlap_l2_error, field])
        end
    end
    write_overlap_l2_error_csv(sim, overlap_rows)
    return nothing
end

function write_overlap_l2_error_screen(domain_name::String, side_set_name::String, field::String, overlap_l2_error::Float64)
    field_label = field == "disp" ? "displacement" : field == "velo" ? "velocity" : "acceleration"
    norma_logf(
        0,
        :summary,
        "Overlap L2 %s error [%s:%s] = %.6e",
        field_label,
        domain_name,
        side_set_name,
        overlap_l2_error,
    )
    flush(stdout)
    return nothing
end


function save_curr_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    integrator.prev_disp = copy(integrator.displacement)
    integrator.prev_velo = copy(integrator.velocity)
    integrator.prev_acce = copy(integrator.acceleration)
    integrator.prev_∂Ω_f = copy(get_internal_force(sim.model))
    # Save per-QP material state so restore_prev_state can roll back internal
    # variables (e.g. plastic strain) when a step fails and is retried with a
    # smaller time step.  Without this, model.state_old is left at the
    # non-converged end-of-failed-solve values and the retry starts from a
    # physically inconsistent internal-variable state.
    # Guard: only SolidMechanics has state_old; ROM models do not.
    if sim.model isa SolidMechanics
        sim.model.prev_state_old = deepcopy(sim.model.state_old)
    end
    return nothing
end

function restore_prev_state(sim::SingleDomainSimulation)
    integrator = sim.integrator
    integrator.time = integrator.prev_time
    integrator.displacement .= integrator.prev_disp
    integrator.velocity .= integrator.prev_velo
    integrator.acceleration .= integrator.prev_acce
    set_internal_force!(sim.model, copy(integrator.prev_∂Ω_f))
    # Restore the per-QP material state saved by save_curr_state so the retry
    # begins from the correct internal-variable state at the start of the step.
    # Guard: only SolidMechanics has prev_state_old; ROM models do not.  This
    # mirrors the guard in save_curr_state -- ROM subdomains have no
    # per-quadrature-point state of their own to roll back (their fom_model is
    # a reconstruction target for output/BCs only, not part of the reduced
    # time integration), so there is nothing to restore for them here.
    if sim.model isa SolidMechanics && !isempty(sim.model.prev_state_old)
        sim.model.state_old = deepcopy(sim.model.prev_state_old)
    end
    return nothing
end

function save_stop_state(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    subsims = sim.subsims
    if controller.use_interface_predictor
        for i in 1:num_domains
            if !isempty(controller.stop_disp[i])
                controller.prev_stop_disp[i] = controller.stop_disp[i]
            end
            if !isempty(controller.stop_∂Ω_f[i])
                controller.prev_stop_∂Ω_f[i] = controller.stop_∂Ω_f[i]
            end
        end
    end
    for i in 1:num_domains
        subsim = subsims[i]
        controller.stop_disp[i] = copy(subsim.integrator.displacement)
        controller.stop_velo[i] = copy(subsim.integrator.velocity)
        controller.stop_acce[i] = copy(subsim.integrator.acceleration)
        controller.stop_∂Ω_f[i] = copy(get_internal_force(subsim.model))
    end
end

function restore_stop_state(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    subsims = sim.subsims
    for i in 1:num_domains
        subsim = subsims[i]
        subsim.integrator.time = controller.prev_time
        subsim.integrator.displacement .= controller.stop_disp[i]
        subsim.integrator.velocity .= controller.stop_velo[i]
        subsim.integrator.acceleration .= controller.stop_acce[i]
        set_internal_force!(subsim.model, copy(controller.stop_∂Ω_f[i]))
        if subsim.model isa RomModel
            reconstruct_fom_fields!(subsim.integrator, subsim.solver, subsim.model)
        end
    end
end

function save_schwarz_state(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = sim.num_domains
    for i in 1:num_domains
        subsim = subsims[i]
        # Note: Integrator values are not rotated due to inclined support as the
        # schwarz variables are only used for Schwarz convergence which are compared
        # to simulation integrator values.
        controller.schwarz_disp[i] = copy(subsim.integrator.displacement)
        controller.schwarz_velo[i] = copy(subsim.integrator.velocity)
        controller.schwarz_acce[i] = copy(subsim.integrator.acceleration)
    end
end

function swap_swappable_bcs(sim::MultiDomainSimulation)
    has_rom = any(subsim -> subsim.model isa RomModel, sim.subsims)
    if has_rom
        for subsim in sim.subsims
            for bc in subsim.model.boundary_conditions
                if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition && bc.swap_bcs == true
                    norma_abort("swap BC types not supported with RomModel nonoverlapping Schwarz")
                end
            end
        end
    end
    for subsim in sim.subsims
        swap_swappable_bcs(subsim)
    end
end

function swap_swappable_bcs(sim::SingleDomainSimulation)
    for bc in sim.model.boundary_conditions
        if is_swappable_dn_schwarz(bc)
            if (bc.swap_bcs == true)
                bc.is_dirichlet = !bc.is_dirichlet
            end
        end
    end
end

function subcycle(sim::MultiDomainSimulation)
    controller = sim.controller
    num_domains = sim.num_domains
    for subsim_index in 1:num_domains
        subsim = sim.subsims[subsim_index]
        norma_log(4, :domain, subsim.name)
        reset_history(controller, subsim_index)
        # advance_time() clips the last substep to land exactly on the stop
        # time and stores the clipped value in integrator.time_step (the
        # integrator updates need the actual step taken). Without a reset the
        # clipped remainder becomes the nominal step of every later pass, and
        # the Schwarz re-integration of each stop amplifies the floating-point
        # remainder geometrically (×3 per iteration for a 4:1 subcycle) until
        # the step collapses to a sliver. Restore the nominal step before each
        # pass; genuine adaptive stepping (minimum < maximum) keeps its
        # adapted value.
        integrator = subsim.integrator
        if integrator.minimum_time_step == integrator.maximum_time_step
            integrator.time_step = integrator.maximum_time_step
        end
        # Anchor the history with the stop-start state (the converged
        # previous stop, just restored by restore_stop_state). Snapshots are
        # otherwise pushed only after each substep, so without this anchor a
        # non-subcycling subdomain's history is a single end-of-stop entry
        # and a finer-stepping partner receives constant end-state data
        # instead of the interpolant in time between the stop endpoints.
        save_history_snapshot(controller, subsim, subsim_index)
        while true
            advance_time(subsim)
            advance_one_step(subsim)
            if controller.active_contact == true && controller.naive_stabilized == true
                apply_naive_stabilized_bcs(subsim)
            end
            save_history_snapshot(controller, subsim, subsim_index)
            if stop_subcyle(subsim) == true
                break
            end
        end
    end
    return nothing
end

function subcycle(sim::SingleDomainSimulation)
    is_explicit = sim.integrator isa ExplicitDynamicTimeIntegrator
    t_last_log = time()
    while true
        advance_time(sim)
        advance_one_step(sim)
        if is_explicit && time() - t_last_log >= 1.0
            norma_logf(4, :progress, "Time = %.2e : Δt = %.2e",
                       sim.integrator.time, sim.integrator.time_step)
            t_last_log = time()
        end
        if stop_subcyle(sim) == true
            break
        end
    end
    return nothing
end

function stop_subcyle(sim::SingleDomainSimulation)
    return isapprox(sim.integrator.time, sim.controller.time; rtol=1.0e-06, atol=1.0e-12)
end

function reset_histories(sim::MultiDomainSimulation)
    controller = sim.controller
    if controller.is_schwarz == false
        return nothing
    end
    num_domains = sim.num_domains
    for subsim_index in 1:num_domains
        reset_history(controller, subsim_index)
    end
    return nothing
end

# Clear the per-interface, per-time-slot relaxation state at the start of every
# controller stop. Slots are keyed by substep time within the stop, so state
# left over from a previous stop is never a valid previous iterate. Dropping the
# interface keys as well is equivalent: they are re-created on demand.
function reset_relaxation_state!(controller::MultiDomainTimeController)
    empty!(controller.lambda_time)
    empty!(controller.lambda_disp)
    empty!(controller.lambda_velo)
    empty!(controller.lambda_acce)
    empty!(controller.aitken_prev_residual_disp)
    empty!(controller.aitken_prev_residual_velo)
    empty!(controller.aitken_prev_residual_acce)
    empty!(controller.aitken_theta_disp)
    empty!(controller.aitken_prev_lambda_disp)
    return nothing
end

function reset_history(controller::MultiDomainTimeController, subsim_index::Int64)
    if controller.is_schwarz == false
        return nothing
    end
    empty!(controller.time_hist[subsim_index])
    empty!(controller.disp_hist[subsim_index])
    empty!(controller.velo_hist[subsim_index])
    empty!(controller.acce_hist[subsim_index])
    empty!(controller.∂Ω_f_hist[subsim_index])
    return nothing
end

function save_history_snapshot(controller::MultiDomainTimeController, sim::SingleDomainSimulation, subsim_index::Int64)
    if controller.is_schwarz == false
        return nothing
    end
    push!(controller.time_hist[subsim_index], sim.integrator.time)
    push!(controller.disp_hist[subsim_index], copy(sim.integrator.displacement))
    push!(controller.velo_hist[subsim_index], copy(sim.integrator.velocity))
    push!(controller.acce_hist[subsim_index], copy(sim.integrator.acceleration))
    push!(controller.∂Ω_f_hist[subsim_index], copy(get_internal_force(sim.model)))
    return nothing
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    controller = sim.controller
    subsims = sim.subsims
    num_domains = sim.num_domains
    norms_pos = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i in 1:num_domains
        Δt = controller.time_step
        u_prev = controller.schwarz_disp[i] + Δt * controller.schwarz_velo[i]
        u_curr = subsims[i].integrator.displacement + Δt * subsims[i].integrator.velocity
        if subsims[i].model isa SolidMechanics
            X = vec(subsims[i].model.reference)
            norms_pos[i] = norm(X + u_curr)
        else
            norms_pos[i] = norm(u_curr)
        end
        norms_diff[i] = norm(u_curr - u_prev)
    end
    norm_pos = norm(norms_pos)
    norm_diff = norm(norms_diff)
    controller.absolute_error = norm_diff
    controller.relative_error = norm_pos > 0.0 ? norm_diff / norm_pos : norm_diff
    conv_abs = controller.absolute_error ≤ controller.absolute_tolerance
    conv_rel = controller.relative_error ≤ controller.relative_tolerance
    controller.converged = conv_abs || conv_rel
    return controller.absolute_error, controller.relative_error
end

# A step that leaves the Schwarz loop without meeting either tolerance has
# exhausted `maximum iterations`, and until this was reported the run said
# nothing: the log ended in "Simulation Complete" and the simulation was not
# marked failed, so an interface still far from converged was indistinguishable
# from a converged one. Announce it, and let the input file promote it to an
# abort, mirroring `stalled interface jump action` for the impedance condition.
function report_unconverged_step(sim::MultiDomainSimulation, iteration_number::Int64)
    controller = sim.controller
    action = get(sim.params, "unconverged step action", "warn")
    if action == "abort"
        norma_abortf(
            "Schwarz did not converge at stop %d in %d iterations: |ΔU| = %.2e against " *
            "absolute tolerance %.2e, |ΔU|/|U| = %.2e against relative tolerance %.2e, and " *
            "`unconverged step action: abort` is set. Raise `maximum iterations`, loosen the " *
            "tolerances, or change the relaxation.",
            controller.stop, iteration_number, controller.absolute_error,
            controller.absolute_tolerance, controller.relative_error, controller.relative_tolerance,
        )
    elseif action != "warn"
        norma_abort(
            "Unknown `unconverged step action: $(action)`. Valid values are `warn` (default) and `abort`.",
        )
    end
    norma_logf(
        0,
        :warning,
        "Schwarz did not converge at stop %d in %d iterations: |ΔU| = %.2e > %.2e and " *
        "|ΔU|/|U| = %.2e > %.2e; continuing.",
        controller.stop, iteration_number, controller.absolute_error,
        controller.absolute_tolerance, controller.relative_error, controller.relative_tolerance,
    )
    return nothing
end

function stop_schwarz(sim::MultiDomainSimulation, iteration_number::Int64)
    if sim.controller.absolute_error == 0.0
        return true
    end
    exceeds_minimum_iterations = iteration_number > sim.controller.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations = iteration_number > sim.controller.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return sim.controller.converged
end

function check_overlap(model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition)
    distance_tol = 0.0
    parametric_tol = 1.0e-06
    overlap = false
    unique_node_indices = unique(bc.side_set_node_indices)
    coupled_model = coupled_subsim_of(bc).model
    coupled_bc = coupled_model.boundary_conditions[bc.coupled_bc_index]
    coupled_side_set_id = coupled_bc.side_set_id

    # Tolerance value to calculate closest point projection, to avoid projection failure
    # Set to 10 times the minimum characteristic edge length of the side set nodes
    # TODO (BRP): perhaps calculate this "characteristic size" a priori, and only recalculate for finite kinematics
    nodal_points = model.reference[:, unique_node_indices] + model.displacement[:, unique_node_indices]
    # Compute minimum distances between the nodes in the side set
    max_distance_tolerance =
        10 * minimum([
            minimum([norm(p1 - p2) for (j, p2) in enumerate(eachcol(nodal_points)) if i != j]) for
            (i, p1) in enumerate(eachcol(nodal_points))
        ])

    for node_index in unique_node_indices
        point = model.reference[:, node_index] + model.displacement[:, node_index]
        # Precompute the face node distance
        face_nodes, face_node_indices, min_distance = closest_face_to_point(point, coupled_model, coupled_side_set_id)
        if min_distance > max_distance_tolerance
            continue
        end
        _, ξ, _, coupled_face_node_indices, _, distance = project_point_to_side_set(
            point, face_nodes, face_node_indices
        )
        num_nodes_coupled_side = length(coupled_face_node_indices)
        parametric_dim = length(ξ)
        element_type = get_element_type(parametric_dim, num_nodes_coupled_side)
        overlap = distance ≤ distance_tol && is_inside_parametric(element_type, ξ, parametric_tol)
        if overlap == true
            break
        end
    end
    return overlap
end

function check_compression(
    mesh::ExodusDatabase, model::SolidMechanics, bc::SolidMechanicsContactSchwarzBoundaryCondition
)
    compression_tol = 0.0
    compression = false
    nodal_forces = get_dst_force(bc)
    normals = compute_normal(mesh, bc.side_set_id, model)
    global_from_local_map = bc.global_from_local_map
    for local_node in eachindex(global_from_local_map)
        local_range = (3 * (local_node - 1) + 1):(3 * local_node)
        nodal_force = nodal_forces[local_range]
        normal = normals[:, local_node]
        normal_force = dot(nodal_force, normal)
        compressive_force = normal_force ≤ compression_tol
        if compressive_force == true
            compression = true
            break
        end
    end
    return compression
end

function initialize_bc_projectors(sim::MultiDomainSimulation)
    for subsim in sim.subsims
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition
                compute_impedance_schwarz_projectors!(subsim.model, bc)
            elseif bc isa SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
                compute_impedance_overlap_schwarz_projectors!(subsim.model, bc)
            elseif bc isa SolidMechanicsOverlapSchwarzBoundaryCondition && bc.use_weak
                coupled_model = get_fom_model(coupled_subsim_of(bc))
                fom_model = get_fom_model(subsim)
                W = get_square_projection_matrix(fom_model, bc)
                L = get_overlap_rectangular_projection_matrix(
                    fom_model, bc, coupled_model, bc.coupled_block_name, bc.search_tolerance
                )
                bc.dirichlet_projector = (W \ I) * L
            elseif is_swappable_dn_schwarz(bc)
                compute_dirichlet_projector(subsim.model, bc)
                compute_neumann_projector(subsim.model, bc)
                if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                    fom_model = get_fom_model(subsim)
                    bc.square_projector = get_square_projection_matrix(fom_model, bc)
                end
            end
        end
    end
end

function detect_contact(sim::MultiDomainSimulation)
    if sim.controller.schwarz_contact == false
        return nothing
    end
    num_domains = sim.num_domains
    persistence = sim.controller.active_contact
    contact_domain = falses(num_domains)
    for domain in 1:num_domains
        subsim = sim.subsims[domain]
        mesh = subsim.params["input_mesh"]
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SolidMechanicsContactSchwarzBoundaryCondition
                if persistence == true
                    compression = check_compression(mesh, subsim.model, bc)
                    contact_domain[domain] = compression == true
                else
                    overlap = check_overlap(subsim.model, bc)
                    contact_domain[domain] = overlap == true
                end
            end
        end
    end

    sim.controller.active_contact = any(contact_domain)
    for domain in 1:num_domains
        subsim = sim.subsims[domain]
        bcs = subsim.model.boundary_conditions
        for bc in bcs
            if bc isa SolidMechanicsContactSchwarzBoundaryCondition
                bc.active_contact = sim.controller.active_contact
            end
        end
    end
    if sim.controller.active_contact == true
        norma_log(0, :contact, "Detected")
    end
    resize!(sim.controller.contact_hist, sim.controller.stop + 1)
    sim.controller.contact_hist[sim.controller.stop + 1] = sim.controller.active_contact
    write_schwarz_params_csv(sim)
    return nothing
end

"""
    get_phys_state_for_predictor(state, subsim) -> Vector{Float64}

For a FOM subdomain, returns `state` unchanged (already a physical DOF vector indexed
as [3*(node-1)+d] for node index `node` and component d ∈ {1,2,3}).

For a ROM subdomain, `state` holds ROM latent coordinates (length = num_modes).
This function reconstructs the full physical DOF vector by applying the model basis:

    u_phys[3*(j-1)+d] = basis[d, j, :] ⋅ state    for each node j, component d

Only free DOFs are set; constrained DOFs default to zero. This is safe because
`extract_interface_state` only ever reads interface nodes, which are always free DOFs.
"""
function get_phys_state_for_predictor(state::Vector{Float64}, subsim::Simulation)
    model = subsim.model
    model isa RomModel || return state
    num_nodes = size(model.basis, 2)
    phys = zeros(3 * num_nodes)
    @inbounds for j in 1:num_nodes
        for d in 1:3
            phys[3 * (j - 1) + d] = model.basis[d, j, :]' * state
        end
    end
    return phys
end

function extract_interface_state(global_state::Vector{Float64}, bc::SolidMechanicsSchwarzBoundaryCondition)
    n_local = length(bc.global_from_local_map)
    local_mat = Matrix{Float64}(undef, 3, n_local)
    for (li, gi) in enumerate(bc.global_from_local_map)
        @inbounds local_mat[:, li] = global_state[(3gi - 2):(3gi)]
    end
    return local_mat
end

function write_interface_state!(global_state::Vector{Float64}, local_mat::Matrix{Float64}, bc::SolidMechanicsSchwarzBoundaryCondition)
    for (li, gi) in enumerate(bc.global_from_local_map)
        @inbounds global_state[(3gi - 2):(3gi)] = local_mat[:, li]
    end
    return nothing
end

function compute_interface_predictor!(sim::MultiDomainSimulation)
    controller = sim.controller
    processed = Set{Tuple{Int,Int}}()

    for (dom_k, subsim_k) in enumerate(sim.subsims)
        for bc_k in subsim_k.model.boundary_conditions
            bc_k isa SolidMechanicsNonOverlapCouplingSchwarzBoundaryCondition || continue

            dom_j = bc_k.coupled_handle.id
            pair = minmax(dom_k, dom_j)
            pair ∈ processed && continue
            push!(processed, pair)

            subsim_j = sim.subsims[dom_j]
            bc_j = subsim_j.model.boundary_conditions[bc_k.coupled_bc_index]

            # Projectors: P_kj maps from domain j's interface space to domain k's
            #             P_jk maps from domain k's interface space to domain j's
            P_kj = bc_k.dirichlet_projector
            P_jk = bc_j.dirichlet_projector

            # Start predictors from the current stop state (t_n)
            pred_disp_k = copy(controller.stop_disp[dom_k])
            pred_velo_k = copy(controller.stop_velo[dom_k])
            pred_acce_k = copy(controller.stop_acce[dom_k])
            pred_disp_j = copy(controller.stop_disp[dom_j])
            pred_velo_j = copy(controller.stop_velo[dom_j])
            pred_acce_j = copy(controller.stop_acce[dom_j])

            both_static = is_static(subsim_k.integrator) && is_static(subsim_j.integrator)

            if both_static
                # Quasi-static: velocities and accelerations are always zero.
                # Displacement predictor: linear extrapolation in time when prior step exists.
                isempty(controller.prev_stop_disp[dom_k]) && continue
                isempty(controller.prev_stop_disp[dom_j]) && continue

                u_k = extract_interface_state(get_phys_state_for_predictor(controller.stop_disp[dom_k], subsim_k), bc_k)
                u_k_prev = extract_interface_state(get_phys_state_for_predictor(controller.prev_stop_disp[dom_k], subsim_k), bc_k)
                u_j = extract_interface_state(get_phys_state_for_predictor(controller.stop_disp[dom_j], subsim_j), bc_j)
                u_j_prev = extract_interface_state(get_phys_state_for_predictor(controller.prev_stop_disp[dom_j], subsim_j), bc_j)

                # Each domain extrapolates its own interface displacement independently.
                # This correctly handles non-conforming meshes: domain k predicts its own
                # interface motion, and domain j predicts its own interface motion.
                u_k_pred = @. 2 * u_k - u_k_prev
                u_j_pred = @. 2 * u_j - u_j_prev

                # For ROM subdomains, pred_disp is the ROM latent state (not a physical
                # DOF vector), so write_interface_state! must not be called on it.
                # The t_n latent state is already a safe predictor for ROM domains.
                subsim_k.model isa RomModel || write_interface_state!(pred_disp_k, u_k_pred, bc_k)
                subsim_j.model isa RomModel || write_interface_state!(pred_disp_j, u_j_pred, bc_j)
            else
                # Dynamic: velocity predictor via perfectly inelastic collision using
                # the interface mass matrices W (row-sum lumping).
                # W = ∫ N Nᵀ dΓ is purely geometric; it works for both explicit
                # (lumped_mass) and implicit (consistent mass) time integrators.
                W_k = bc_k.square_projector    # (n_k × n_k)
                W_j = bc_j.square_projector    # (n_j × n_j)

                # Lumped interface masses: row sums of W
                m_k = vec(sum(W_k; dims=2))    # (n_k,)
                m_j = vec(sum(W_j; dims=2))    # (n_j,)

                # Interface velocities at t_n from stop state.
                # For ROM subdomains stop_velo holds latent coordinates; reconstruct
                # the physical DOF vector first so extract_interface_state can index it.
                v_k = extract_interface_state(get_phys_state_for_predictor(controller.stop_velo[dom_k], subsim_k), bc_k)
                v_j = extract_interface_state(get_phys_state_for_predictor(controller.stop_velo[dom_j], subsim_j), bc_j)

                # Project domain j's mass and velocity into domain k's interface space
                m_j_in_k = P_kj * m_j    # (n_k,)
                v_j_in_k = Matrix{Float64}(undef, 3, size(P_kj, 1))
                for d in 1:3
                    v_j_in_k[d, :] = P_kj * v_j[d, :]
                end

                # Perfectly inelastic collision: conserves linear momentum
                m_total = m_k .+ m_j_in_k
                v_pred_k = Matrix{Float64}(undef, 3, length(m_k))
                for d in 1:3
                    @. v_pred_k[d, :] = (m_k * v_k[d, :] + m_j_in_k * v_j_in_k[d, :]) / m_total
                end

                # Project predicted velocity from domain k's space back to domain j's space
                v_pred_j = Matrix{Float64}(undef, 3, size(P_jk, 1))
                for d in 1:3
                    v_pred_j[d, :] = P_jk * v_pred_k[d, :]
                end

                # Displacement predictor: velocity-consistent linear extrapolation
                # u_pred = u_n + Δt * v_pred  (each domain uses its own u_n as baseline)
                Δt = controller.time_step
                u_k = extract_interface_state(get_phys_state_for_predictor(controller.stop_disp[dom_k], subsim_k), bc_k)
                u_j = extract_interface_state(get_phys_state_for_predictor(controller.stop_disp[dom_j], subsim_j), bc_j)
                u_pred_k = @. u_k + Δt * v_pred_k
                u_pred_j = @. u_j + Δt * v_pred_j

                # Acceleration predictor: zero at the interface (standard Newmark prediction)
                # pred_acce_k and pred_acce_j already contain zero-valued stop_acce (from t_n
                # corrected acceleration); no modification needed here.

                # For ROM subdomains, pred_disp/pred_velo hold ROM latent coordinates.
                # Calling write_interface_state! on them with a physical-space matrix
                # would use FOM node indices that are out of bounds for the latent vector.
                # Instead, leave ROM predictors as the t_n latent state (a safe fallback):
                # schwarz.jl copies predictor_disp into coupled_integrator.displacement
                # and then calls reconstruct_fom_fields!, so latent coordinates are correct.
                if !(subsim_k.model isa RomModel)
                    write_interface_state!(pred_disp_k, u_pred_k, bc_k)
                    write_interface_state!(pred_velo_k, v_pred_k, bc_k)
                end
                if !(subsim_j.model isa RomModel)
                    write_interface_state!(pred_disp_j, u_pred_j, bc_j)
                    write_interface_state!(pred_velo_j, v_pred_j, bc_j)
                end
            end

            # Force predictor: linear temporal extrapolation f_pred = 2*f_n - f_{n-1}.
            # Which domain's force is read by the coupling:
            #   - Non-overlap Neumann (is_dirichlet == false) reads force FROM the coupled domain.
            #   - Robin reads force from BOTH coupled domains.
            # So set predictor_∂Ω_f[x] when domain x's force is read by the other domain.
            is_robin = bc_k isa SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition

            pred_∂Ω_f_k = copy(controller.stop_∂Ω_f[dom_k])
            pred_∂Ω_f_j = copy(controller.stop_∂Ω_f[dom_j])

            if (is_robin || !bc_j.is_dirichlet) && !isempty(controller.prev_stop_∂Ω_f[dom_k])
                @. pred_∂Ω_f_k = 2 * controller.stop_∂Ω_f[dom_k] - controller.prev_stop_∂Ω_f[dom_k]
            end
            if (is_robin || !bc_k.is_dirichlet) && !isempty(controller.prev_stop_∂Ω_f[dom_j])
                @. pred_∂Ω_f_j = 2 * controller.stop_∂Ω_f[dom_j] - controller.prev_stop_∂Ω_f[dom_j]
            end

            controller.predictor_disp[dom_k] = pred_disp_k
            controller.predictor_velo[dom_k] = pred_velo_k
            controller.predictor_acce[dom_k] = pred_acce_k
            controller.predictor_∂Ω_f[dom_k] = pred_∂Ω_f_k
            controller.predictor_disp[dom_j] = pred_disp_j
            controller.predictor_velo[dom_j] = pred_velo_j
            controller.predictor_acce[dom_j] = pred_acce_j
            controller.predictor_∂Ω_f[dom_j] = pred_∂Ω_f_j
        end
    end
    return nothing
end

function write_schwarz_params_csv(sim::MultiDomainSimulation)
    stop = sim.controller.stop
    csv_interval = get(sim.params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        index_string = "-" * string(stop; pad=4)
        contact_filename = "contact" * index_string * ".csv"
        writedlm(contact_filename, sim.controller.active_contact, '\n')
        iters_filename = "iterations" * index_string * ".csv"
        writedlm(iters_filename, sim.controller.iteration_number, '\n')
    end
end

function write_overlap_l2_error_csv(sim::MultiDomainSimulation, overlap_rows::Vector{Vector{Any}})
    isempty(overlap_rows) && return nothing
    time         = sim.controller.time
    initial_time = sim.controller.initial_time
    csv_interval = Float64(get(sim.params, "CSV output interval", 0.0))
    if !_is_output_time(time, initial_time, csv_interval)
        return nothing
    end
    stop = sim.controller.stop
    index_string = "-" * string(stop; pad=4)
    # Group rows by field type ("disp", "velo", "acce") and write a separate CSV for each.
    field_label_map = Dict("disp" => "disp", "velo" => "velo", "acce" => "acce")
    rows_by_field = Dict{String,Vector{Vector{Any}}}()
    for row in overlap_rows
        field = row[4]::String
        push!(get!(rows_by_field, field, Vector{Vector{Any}}()), row)
    end
    for (field, rows) in rows_by_field
        filename = "overlap-l2-" * field_label_map[field] * "-rel-errors" * index_string * ".csv"
        open(filename, "w") do io
            write(io, "domain,side_set,overlap_l2_error\n")
            for row in rows
                @printf(io, "%s,%s,%.16e\n", row[1], row[2], row[3])
            end
        end
    end
    return nothing
end
