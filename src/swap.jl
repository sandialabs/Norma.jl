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
    elseif criterion_type == "stress recovery"
        tolerance = Float64(get(params, "tolerance", 1.0e-2))
        haskey(params, "direction") || norma_abort(
            "Stress recovery swap criterion requires 'direction:' — use 'refine' to swap " *
            "when the relative lumped/consistent stress difference exceeds the tolerance, " *
            "or 'coarsen' to swap when it falls below the tolerance",
        )
        direction = Symbol(params["direction"])
        direction in (:refine, :coarsen) || norma_abort(
            "Unknown stress recovery swap direction: '$direction' (expected 'refine' or 'coarsen')",
        )
        return StressRecoverySwapCriterion(tolerance, direction)
    elseif criterion_type == "elastic to plastic transition"
        tolerance = Float64(get(params, "tolerance", 0.05))
        tolerance > 0.0 || norma_abort(
            "Elastic to plastic transition swap criterion requires a positive 'tolerance' (relative fraction of " *
            "the yield stress); got $tolerance",
        )
        return ElasticToPlasticTransitionSwapCriterion(tolerance)
    elseif criterion_type == "overlap l2 error"
        tolerance = Float64(get(params, "tolerance", 1.0e-6))
        haskey(params, "direction") || norma_abort(
            "Overlap L2 error swap criterion requires 'direction:' — use 'refine' to swap " *
            "when the overlap L2 displacement error exceeds the tolerance, " *
            "or 'coarsen' to swap when it falls below the tolerance",
        )
        direction = Symbol(params["direction"])
        direction in (:refine, :coarsen) || norma_abort(
            "Unknown overlap L2 error swap direction: '$direction' (expected 'refine' or 'coarsen')",
        )
        return SchwarzOverlapL2SwapCriterion(tolerance, direction)
    else
        norma_abort("Unknown swap criterion type: $criterion_type")
    end
end

function validate_swap_plans(sim::MultiDomainSimulation)
    isempty(sim.swaps) && return nothing
    # Swapping is not supported when any subsim runs an OpInf or NNOpInf ROM.
    # ROM state spaces differ from FOM state spaces, so the state-transfer step
    # that occurs at swap time is not well-defined in that context.
    for subsim in sim.subsims
        if subsim.model isa RomModel
            norma_abort(
                "Swap criteria cannot be used when any subsim is an OpInf or NNOpInf ROM. " *
                "Subsim '$(subsim.name)' uses a ROM model.  Swapping is only supported for " *
                "full-order model (FOM) subsims.",
            )
        end
    end
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
    isempty(sim.swaps) && return nothing
    # Swapping is not supported for OpInf or NNOpInf ROM simulations.
    if sim.model isa RomModel
        norma_abort(
            "Swap criteria cannot be used with an OpInf or NNOpInf ROM simulation.  " *
            "Swapping is only supported for full-order model (FOM) simulations.",
        )
    end
    for plan in sim.swaps
        plan.subsim_name === nothing ||
            norma_abort("Single-domain swap plan must not specify `subsim:` (there is no choice of target)")
        isfile(plan.replacement_file) ||
            norma_abort("Swap replacement file not found: $(plan.replacement_file)")
        validate_swap_criterion(plan.criterion, sim)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Per-criterion validation
# ---------------------------------------------------------------------------

# Default: no additional checks required.
validate_swap_criterion(::SwapCriterion, ::Simulation) = nothing

# SchwarzOverlapL2SwapCriterion is only meaningful for a MultiDomainSimulation
# with overlapping Schwarz coupling.  Single-domain simulations have no
# Schwarz interface at all, so the criterion can never produce a meaningful
# value there.
function validate_swap_criterion(::SchwarzOverlapL2SwapCriterion, ::SingleDomainSimulation)
    norma_abort(
        "SchwarzOverlapL2SwapCriterion is only valid for overlapping Schwarz multi-domain " *
        "simulations.  It cannot be used with a monolithic (single-domain) formulation, " *
        "non-overlapping Schwarz, or Schwarz contact.",
    )
end

# For a multi-domain simulation verify two things:
#   1. At least one subsim has a SolidMechanicsOverlapSchwarzBoundaryCondition
#      (ruling out non-overlapping Schwarz and contact-only setups).
#   2. At least one of those BCs has 'compute overlap L2 error: true'
#      (the criterion cannot evaluate without this flag).
#
# NOTE: this is called from validate_swap_criteria (below), NOT from
# validate_swap_plans, because boundary_conditions is populated by create_bcs
# which runs after the MultiDomainSimulation constructor returns.
function validate_swap_criterion(::SchwarzOverlapL2SwapCriterion, sim::MultiDomainSimulation)
    has_overlap_bc = any(
        any(bc isa SolidMechanicsOverlapSchwarzBoundaryCondition
            for bc in sub.model.boundary_conditions)
        for sub in sim.subsims
    )
    if !has_overlap_bc
        norma_abort(
            "SchwarzOverlapL2SwapCriterion requires overlapping Schwarz boundary conditions " *
            "(SolidMechanicsOverlapSchwarzBoundaryCondition) but none were found in any " *
            "subsim.  This criterion cannot be used with non-overlapping Schwarz or " *
            "Schwarz contact formulations.",
        )
    end
    has_l2_flag = any(
        any(stores_overlap_l2_error(bc) for bc in sub.model.boundary_conditions)
        for sub in sim.subsims
    )
    if !has_l2_flag
        norma_abort(
            "SchwarzOverlapL2SwapCriterion requires 'compute overlap L2 error: true' on at " *
            "least one Schwarz overlap boundary condition, but no such BC was found.  Add " *
            "'compute overlap L2 error: true' to the relevant 'Schwarz overlap:' entry in " *
            "the subsim input file.",
        )
    end
end

# ElasticToPlasticTransitionSwapCriterion is only meaningful for SolidMechanics models
# that have at least one block using a J2Plasticity material.  Validation is
# deferred to validate_swap_criterion (called after model construction) so that
# the material list is available.
function validate_swap_criterion(::ElasticToPlasticTransitionSwapCriterion, sim::SingleDomainSimulation)
    model = sim.model
    if !(model isa SolidMechanics)
        norma_abort(
            "ElasticToPlasticTransitionSwapCriterion is only valid for SolidMechanics models. " *
            "The current model type is $(typeof(model)).",
        )
    end
    _validate_elastic_to_plastic_transition_materials(model)
    return nothing
end

function validate_swap_criterion(::ElasticToPlasticTransitionSwapCriterion, sim::MultiDomainSimulation)
    for sub in sim.subsims
        if sub.model isa SolidMechanics
            _validate_elastic_to_plastic_transition_materials(sub.model)
        end
    end
    return nothing
end

function _validate_elastic_to_plastic_transition_materials(model::SolidMechanics)
    any(m isa J2Plasticity for m in model.materials) || norma_abort(
        "ElasticToPlasticTransitionSwapCriterion requires at least one block with a J2Plasticity " *
        "material, but none were found.  This criterion cannot be used with purely elastic " *
        "or other non-J2 material models.",
    )
    return nothing
end

# Called from create_simulation after create_bcs, so that boundary_conditions
# is fully populated when BC-dependent criterion checks run.
# Single-domain criteria that need BC checks are handled in validate_swap_plans
# (called earlier, inside SingleDomainSimulation()), so the single-domain
# overload here is a no-op.
validate_swap_criteria(::SingleDomainSimulation) = nothing

function validate_swap_criteria(sim::MultiDomainSimulation)
    for plan in sim.swaps
        validate_swap_criterion(plan.criterion, sim)
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
# StressRecoverySwapCriterion dispatch
# ---------------------------------------------------------------------------

# For a single-domain simulation apply the criterion to its one model.
function should_swap(c::StressRecoverySwapCriterion, sim::SingleDomainSimulation)
    return _stress_recovery_criterion_met(c, sim.model)
end

# For a multi-domain simulation return true only when ALL subsim models
# individually satisfy the criterion.
function should_swap(c::StressRecoverySwapCriterion, sim::MultiDomainSimulation)
    return all(_stress_recovery_criterion_met(c, sub.model) for sub in sim.subsims)
end

# Non-SolidMechanics models have no stress field; never trigger the swap so
# the criterion degrades gracefully in mixed-model setups.
_stress_recovery_criterion_met(::StressRecoverySwapCriterion, ::Model) = false

# Core evaluation for SolidMechanics models: compare the lumped and consistent
# L2-projected nodal stress fields against the tolerance in the chosen
# direction.
function _stress_recovery_criterion_met(c::StressRecoverySwapCriterion, model::SolidMechanics)
    # Fetch this model's lumped and consistent recovery operators, building
    # them on first use.  If the model already carries a BothRecovery, reuse
    # its operators rather than rebuilding them.
    lumped, consistent = get!(c._recovery, model) do
        rec = model.recovery_data
        if rec isa BothRecovery
            return (rec.lumped, rec.consistent)
        end
        norma_log(1, :swap, "StressRecoverySwapCriterion: building lumped and consistent recovery data")
        new_lumped = build_recovery_data(:lumped, model.mesh, model.reference, model.num_int_pts)::LumpedRecovery
        new_consistent =
            build_recovery_data(:consistent, model.mesh, model.reference, model.num_int_pts)::ConsistentRecovery
        return (new_lumped, new_consistent)
    end
    difference = _recovered_stress_difference(lumped, consistent, model)
    fired = c.direction === :refine ? difference > c.tolerance : difference < c.tolerance
    norma_logf(1, :swap,
        "StressRecoverySwapCriterion: relative lumped/consistent stress difference %.3e, tolerance %.3e, direction %s → %s",
        difference, c.tolerance, c.direction, fired ? "swap" : "no swap")
    return fired
end

# Function barrier: receiving the concrete LumpedRecovery / ConsistentRecovery
# types as arguments makes the projection and norm fully type-stable, even
# though the criterion's IdDict cache stores them with abstract element types.
function _recovered_stress_difference(
    lumped::LumpedRecovery, consistent::ConsistentRecovery, model::SolidMechanics
)
    # Assemble the shared L2 RHS once, copy it, then project each copy with
    # its respective mass operator.  The RHS copy is cheap (6 × n_nodes).
    n_nodes = size(model.reference, 2)
    σ_lumped = zeros(6, n_nodes)
    _assemble_l2_rhs_stress!(σ_lumped, model)
    σ_consistent = copy(σ_lumped)
    _apply_inverse_mass!(σ_lumped, lumped)
    _apply_inverse_mass!(σ_consistent, consistent)
    return _rel_frob_diff(σ_lumped, σ_consistent)
end

# Relative Frobenius-norm difference between two same-size matrices:
#   ‖A − B‖_F / max(‖A‖_F, ‖B‖_F)
# Returns 0.0 when both matrices are essentially zero, so a vanishing stress
# field is reported as zero difference rather than 0/0.
function _rel_frob_diff(A::Matrix{Float64}, B::Matrix{Float64})
    diff_sq = 0.0
    a_sq    = 0.0
    b_sq    = 0.0
    @inbounds for i in eachindex(A)
        d        = A[i] - B[i]
        diff_sq += d * d
        a_sq    += A[i] * A[i]
        b_sq    += B[i] * B[i]
    end
    ref = max(sqrt(a_sq), sqrt(b_sq))
    ref < eps(Float64) && return 0.0
    return sqrt(diff_sq) / ref
end

# ---------------------------------------------------------------------------
# ElasticToPlasticTransitionSwapCriterion dispatch
# ---------------------------------------------------------------------------

# For a single-domain simulation apply the criterion to its one model.
# Guard: skip the check before the first solve has completed.  In the evolve
# loop, maybe_apply_swaps! is called BEFORE advance_control (which runs the
# solve), so at the first iteration model.stress is still all zeros.
# prev_time is only advanced to the previous time value after a successful
# solve, so prev_time > initial_time reliably means at least one solve is done.
function should_swap(c::ElasticToPlasticTransitionSwapCriterion, sim::SingleDomainSimulation)
    sim.controller.prev_time > sim.controller.initial_time || return false
    return _elastic_to_plastic_transition_criterion_met(c, sim.model)
end

# For a multi-domain simulation return true only when ALL subsim models that
# carry J2Plasticity materials individually satisfy the criterion.  Models
# without any J2Plasticity block are skipped (they never block the swap).
# Same pre-solve guard as for SingleDomainSimulation.
function should_swap(c::ElasticToPlasticTransitionSwapCriterion, sim::MultiDomainSimulation)
    sim.controller.prev_time > sim.controller.initial_time || return false
    return all(_elastic_to_plastic_transition_criterion_met(c, sub.model) for sub in sim.subsims)
end

# Non-SolidMechanics models have no stress or material fields; skip them.
_elastic_to_plastic_transition_criterion_met(::ElasticToPlasticTransitionSwapCriterion, ::Model) = false

# Core evaluation: iterate over every block that uses J2Plasticity and check
# every integration point.  The criterion fires (returns true, triggering the
# swap) when ANY integration point has a von Mises stress that exceeds the
# upper edge of the yield band:
#   σ_vm > σy * (1 + tolerance)
# This means the elastic-to-plastic transition has already occurred at that
# point and the mesh should be swapped to capture it more accurately.
# Points below the lower edge (still elastic) or inside the band (transitioning)
# do not trigger the swap.
function _elastic_to_plastic_transition_criterion_met(c::ElasticToPlasticTransitionSwapCriterion, model::SolidMechanics)
    any_j2_block = false
    found_above_band = false
    for (block_index, material) in enumerate(model.materials)
        material isa J2Plasticity || continue
        any_j2_block = true
        σy = material.σy
        σy > 0.0 || continue   # degenerate case: zero yield stress, skip
        upper_band = σy * (1.0 + c.tolerance)
        block_stress = model.stress[block_index]
        for element_stress in block_stress
            for qp_stress in element_stress
                σ_vm = von_mises_from_voigt(qp_stress)
                if σ_vm > upper_band
                    found_above_band = true
                    @goto done_scanning   # early exit: one above-band point is enough
                end
            end
        end
    end
    @label done_scanning
    # If the model has no J2 blocks at all the criterion degrades gracefully.
    if !any_j2_block
        norma_log(
            1, :swap,
            "ElasticToPlasticTransitionSwapCriterion: no J2Plasticity blocks found in model; criterion will not fire",
        )
        return false
    end
    # Swap fires when at least one QP has exceeded the upper yield band edge.
    fired = found_above_band
    norma_logf(
        1, :swap,
        "ElasticToPlasticTransitionSwapCriterion: yield-band check (tolerance %.2f%%) → %s",
        100.0 * c.tolerance, fired ? "swap (QP(s) above yield band)" : "no swap (all QPs within or below yield band)",
    )
    return fired
end

# ---------------------------------------------------------------------------
# SchwarzOverlapL2SwapCriterion dispatch
# ---------------------------------------------------------------------------

# For a single-domain simulation scan the model's boundary conditions.
function should_swap(c::SchwarzOverlapL2SwapCriterion, sim::SingleDomainSimulation)
    return _overlap_l2_criterion_met(c, sim.model.boundary_conditions)
end

# For a multi-domain simulation aggregate over all BCs in all subsims.
function should_swap(c::SchwarzOverlapL2SwapCriterion, sim::MultiDomainSimulation)
    all_bcs = Iterators.flatten(sub.model.boundary_conditions for sub in sim.subsims)
    return _overlap_l2_criterion_met(c, all_bcs)
end

# Core logic: iterate over boundary conditions, call compute_overlap_l2_error!
# on each BC that has the computation enabled, take the maximum, and compare
# against the tolerance in the chosen direction.
#
# Returns false (never fires) when no BCs have the error computation enabled,
# so that a misconfigured input degrades gracefully.  A warning is logged in
# that case so the user can diagnose the issue.
function _overlap_l2_criterion_met(c::SchwarzOverlapL2SwapCriterion, boundary_conditions)
    max_error = -Inf
    found_any = false
    for bc in boundary_conditions
        stores_overlap_l2_error(bc) || continue
        err = compute_overlap_l2_error!(bc)
        isnan(err) && continue
        max_error = max(max_error, err)
        found_any = true
    end
    if !found_any
        norma_log(
            1, :swap,
            "SchwarzOverlapL2SwapCriterion: no boundary conditions with 'compute overlap L2 error: true' " *
            "found; criterion will not fire — add 'compute overlap L2 error: true' to the " *
            "relevant Schwarz overlap boundary condition",
        )
        return false
    end
    fired = c.direction === :refine ? max_error > c.tolerance : max_error < c.tolerance
    norma_logf(
        1, :swap,
        "SchwarzOverlapL2SwapCriterion: max overlap L2 error %.6e, tolerance %.6e, direction %s → %s",
        max_error, c.tolerance, c.direction, fired ? "swap" : "no swap",
    )
    return fired
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
# freshly-built recovery_data lives on the model afterward; output buffers
# stay empty so the swap does not silently enable nodal output the user did
# not request.
function _ensure_recovery_for_swap!(model::SolidMechanics)
    if model.recovery_data isa NoRecovery
        norma_log(0, :swap, "Building consistent recovery on demand for state transfer")
        model.recovery_data = build_recovery_data(
            :consistent, model.mesh, model.reference, model.num_int_pts,
        )
    end
    return nothing
end
