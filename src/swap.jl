# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Swap a sub-simulation for a replacement at run time.  See examples/overlap/
# static-same-step/cuboids-swap for a proof-of-concept configuration.

# ---------------------------------------------------------------------------
# YAML parsing
# ---------------------------------------------------------------------------

function parse_swap_plans(params::Parameters)
    haskey(params, "swaps") || return SwapPlan[]
    return [SwapPlan(p) for p in params["swaps"]]
end

function SwapPlan(params::Parameters)
    # `subsim` is required for multi-domain swaps (selects which subsim to
    # replace) and absent for single-domain swaps (there is no choice).
    return SwapPlan(
        get(params, "subsim", nothing),
        create_swap_criterion(params["criterion"]),
        params["replacement"],
        false,
    )
end

function create_swap_criterion(params::Parameters)
    criterion_type = params["type"]
    if criterion_type == "time"
        return TimeSwapCriterion(Float64(params["t_swap"]))
    else
        norma_abort("Unknown swap criterion type: $criterion_type")
    end
end

function validate_swap_plans(sim::MultiDomainSimulation)
    for plan in sim.swaps
        plan.subsim_name === nothing &&
            norma_abort("Multi-domain swap plan must specify `subsim:` to select the subsim to replace")
        haskey(sim.handle_by_name, plan.subsim_name) ||
            norma_abort("Swap plan references unknown subsim: $(plan.subsim_name)")
        isfile(plan.replacement_file) ||
            norma_abort("Swap replacement file not found: $(plan.replacement_file)")
    end
    return nothing
end

function validate_swap_plans(sim::SingleDomainSimulation)
    for plan in sim.swaps
        plan.subsim_name === nothing ||
            norma_abort("Single-domain swap plan must not specify `subsim:` (there is no choice of target)")
        isfile(plan.replacement_file) ||
            norma_abort("Swap replacement file not found: $(plan.replacement_file)")
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Criterion dispatch
# ---------------------------------------------------------------------------

function should_swap(c::TimeSwapCriterion, sim::Simulation)
    return sim.controller.time > c.t_swap
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

# No-op for non-multi sims so the evolve loop stays generic.
maybe_apply_swaps!(::Simulation) = nothing

function maybe_apply_swaps!(sim::MultiDomainSimulation)
    isempty(sim.swaps) && return nothing
    for plan in sim.swaps
        plan.applied && continue
        should_swap(plan.criterion, sim) || continue
        slot = sim.handle_by_name[plan.subsim_name].id
        apply_swap!(sim, slot, plan)
    end
    return nothing
end

function maybe_apply_swaps!(sim::SingleDomainSimulation)
    isempty(sim.swaps) && return nothing
    for plan in sim.swaps
        plan.applied && continue
        should_swap(plan.criterion, sim) || continue
        apply_single_domain_swap!(sim, plan)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# The swap itself
# ---------------------------------------------------------------------------

function apply_swap!(sim::MultiDomainSimulation, slot::Int64, plan::SwapPlan)
    old = sim.subsims[slot]
    norma_logf(0, :swap, "Swapping subsim %s → %s at t = %.6e",
               old.name, plan.replacement_file, sim.controller.time)

    new = build_replacement_subsim(sim, slot, plan)
    copy_model_state!(new.model, old.model)
    align_replacement_time!(new, sim, old)

    finalize_writing(old)
    initialize_writing(new)

    sim.subsims[slot] = new

    # Cached coupled_bc_index and is_dirichlet flags on partner BCs may now
    # point into the old BC list — rebuild them.
    pair_schwarz_bcs(sim)
    # Projectors held on the replacement's BCs need their first fill.
    initialize_bc_projectors(sim)

    plan.applied = true
    return nothing
end

function build_replacement_subsim(sim::MultiDomainSimulation, slot::Int64, plan::SwapPlan)
    subparams = YAML.load_file(plan.replacement_file; dicttype=Parameters)
    subparams["name"] = stripped_name(plan.replacement_file)

    integrator_params = subparams["time integrator"]
    integrator_params["initial time"] = sim.controller.initial_time
    integrator_params["final time"] = sim.controller.final_time
    requested_dt = get(integrator_params, "time step", sim.controller.time_step)
    integrator_params["time step"] = min(requested_dt, sim.controller.time_step)

    subparams["Exodus output interval"] = Float64(
        get(sim.params, "Exodus output interval", sim.controller.time_step))
    subparams["CSV output interval"] = Float64(
        get(sim.params, "CSV output interval", 0.0))

    new = SingleDomainSimulation(subparams)
    # Wire parent/handle BEFORE create_bcs so Schwarz BC factories can
    # resolve peers via sim.handle_by_name in _create_bcs.
    new.parent = sim
    new.handle = DomainHandle(slot)
    create_bcs(new)
    initialize_storage(new)
    return new
end

function align_replacement_time!(new::SingleDomainSimulation,
                                 sim::MultiDomainSimulation,
                                 old::SingleDomainSimulation)
    nc = new.controller
    pc = sim.controller
    nc.initial_time = pc.initial_time
    nc.final_time = pc.final_time
    nc.time_step = pc.time_step
    nc.num_stops = pc.num_stops
    nc.stop = pc.stop
    nc.prev_time = pc.prev_time
    nc.time = pc.time
    new.integrator.prev_time = old.integrator.prev_time
    new.integrator.time = old.integrator.time
    new.integrator.time_step = old.integrator.time_step
    new.model.time = old.model.time
    return nothing
end

# ---------------------------------------------------------------------------
# Single-domain swap
# ---------------------------------------------------------------------------

# Swap a top-level SingleDomainSimulation in place: replace the model,
# integrator, solver, and controller fields with those of a freshly-built
# replacement, then transfer state across.  The outer evolve() loop's
# reference to `sim` stays valid because we only mutate fields, not the
# object identity.
function apply_single_domain_swap!(sim::SingleDomainSimulation, plan::SwapPlan)
    norma_logf(0, :swap, "Swapping single-domain sim %s → %s at t = %.6e",
               sim.name, plan.replacement_file, sim.controller.time)

    new = build_replacement_single_domain_sim(sim, plan)
    copy_model_state!(new.model, sim.model)
    align_replacement_time_single!(new, sim)

    finalize_writing(sim)
    initialize_writing(new)

    sim.name = new.name
    sim.params = new.params
    sim.controller = new.controller
    sim.integrator = new.integrator
    sim.solver = new.solver
    sim.model = new.model
    sim.swaps = new.swaps

    plan.applied = true
    return nothing
end

function build_replacement_single_domain_sim(sim::SingleDomainSimulation, plan::SwapPlan)
    subparams = YAML.load_file(plan.replacement_file; dicttype=Parameters)
    subparams["name"] = stripped_name(plan.replacement_file)

    integrator_params = subparams["time integrator"]
    integrator_params["initial time"] = sim.controller.initial_time
    integrator_params["final time"] = sim.controller.final_time
    requested_dt = get(integrator_params, "time step", sim.controller.time_step)
    integrator_params["time step"] = min(requested_dt, sim.controller.time_step)

    subparams["Exodus output interval"] = Float64(
        get(sim.params, "Exodus output interval", sim.controller.time_step))
    subparams["CSV output interval"] = Float64(
        get(sim.params, "CSV output interval", 0.0))

    new = SingleDomainSimulation(subparams)
    create_bcs(new)
    initialize_storage(new)
    return new
end

function align_replacement_time_single!(new::SingleDomainSimulation, old::SingleDomainSimulation)
    nc = new.controller
    oc = old.controller
    nc.initial_time = oc.initial_time
    nc.final_time = oc.final_time
    nc.time_step = oc.time_step
    nc.num_stops = oc.num_stops
    nc.stop = oc.stop
    nc.prev_time = oc.prev_time
    nc.time = oc.time
    new.integrator.prev_time = old.integrator.prev_time
    new.integrator.time = old.integrator.time
    new.integrator.time_step = old.integrator.time_step
    new.model.time = old.model.time
    return nothing
end

# ---------------------------------------------------------------------------
# State transfer
# ---------------------------------------------------------------------------

copy_model_state!(::Model, ::Model) = nothing

# Public entry: dispatches to the cheap same-mesh direct copy when meshes
# and per-QP state layouts match, otherwise to the cross-mesh L2 path.
function copy_model_state!(dst::SolidMechanics, src::SolidMechanics)
    if _meshes_match(dst, src) && _state_layouts_match(dst, src)
        return _copy_model_state_direct!(dst, src)
    end
    _ensure_recovery_for_swap!(dst)
    if _has_ivs_to_transfer(dst, src)
        _ensure_recovery_for_swap!(src)
    end
    return _copy_model_state_l2!(dst, src)
end

# Same-mesh, same per-QP state layout: direct copy.  Cheap, exact, and does
# not require any recovery infrastructure to be enabled.
function _copy_model_state_direct!(dst::SolidMechanics, src::SolidMechanics)
    dst.displacement .= src.displacement
    dst.velocity .= src.velocity
    dst.acceleration .= src.acceleration
    dst.internal_force .= src.internal_force
    dst.boundary_force .= src.boundary_force
    dst.time = src.time
    dst.strain_energy = src.strain_energy
    dst.state_old = deepcopy(src.state_old)
    dst.state = deepcopy(src.state)
    dst.stress = deepcopy(src.stress)
    dst.stored_energy = deepcopy(src.stored_energy)
    return nothing
end

# Cross-mesh (or same mesh with different per-QP state layout): L2-projected
# transfer of nodal kinematic state and per-QP internal variables.  Stress,
# stored_energy, internal_force, and boundary_force are recomputed by the
# next evaluate, so we don't transfer them.  The `.=` writes preserve the
# integrator's unsafe_wrap views into the model's displacement / velocity /
# acceleration buffers.
function _copy_model_state_l2!(dst::SolidMechanics, src::SolidMechanics)
    dst.displacement .= transfer_field(src, src.displacement, dst)
    dst.velocity     .= transfer_field(src, src.velocity,     dst)
    dst.acceleration .= transfer_field(src, src.acceleration, dst)
    _transfer_qp_internal_variables!(dst, src)
    dst.time = src.time
    dst.strain_energy = src.strain_energy
    return nothing
end

function _meshes_match(dst::SolidMechanics, src::SolidMechanics)
    size(dst.reference) == size(src.reference) || return false
    return isapprox(dst.reference, src.reference; atol=1.0e-12)
end

function _state_layouts_match(dst::SolidMechanics, src::SolidMechanics)
    length(dst.state) == length(src.state) || return false
    @inbounds for b in eachindex(dst.state)
        length(dst.state[b]) == length(src.state[b]) || return false
        for e in eachindex(dst.state[b])
            length(dst.state[b][e]) == length(src.state[b][e]) || return false
            for q in eachindex(dst.state[b][e])
                length(dst.state[b][e][q]) == length(src.state[b][e][q]) || return false
            end
        end
    end
    return true
end

function _has_ivs_to_transfer(dst::SolidMechanics, src::SolidMechanics)
    src_names = collect_internal_variable_names(src.materials)
    dst_names = collect_internal_variable_names(dst.materials)
    return !isempty(intersect(Set(src_names), Set(dst_names)))
end

# Build a consistent recovery on demand if the user didn't enable one.  The
# freshly-built recovery_data lives on the model afterward; the next
# initialize_writing(new) sees it and emits nodal stress output as a
# byproduct of the swap (small surprise, but transparent).
function _ensure_recovery_for_swap!(model::SolidMechanics)
    if model.recovery_data isa NoRecovery
        norma_log(0, :swap, "Building consistent recovery on demand for state transfer")
        model.recovery_data = build_recovery_data(
            :consistent, model.mesh, model.reference, model.num_int_pts,
        )
        model.recovered_stress = zeros(6, size(model.reference, 2))
    end
    return nothing
end
