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
    elseif criterion_type == "overlap l2 relative error"
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
    # Only the concrete OpInfModel subtypes that have copy_model_state! overloads
    # support swapping.  Reject anything else (e.g. NeuralNetworkOpInfRom, or any
    # future RomModel subtype added without a corresponding state-transfer method)
    # with a clear message at parse time rather than a MethodError at swap time.
    for subsim in sim.subsims
        _assert_rom_swap_supported(subsim.model, "Subsim '$(subsim.name)'")
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
    _assert_rom_swap_supported(sim.model, "This simulation")
    for plan in sim.swaps
        plan.subsim_name === nothing ||
            norma_abort("Single-domain swap plan must not specify `subsim:` (there is no choice of target)")
        isfile(plan.replacement_file) ||
            norma_abort("Swap replacement file not found: $(plan.replacement_file)")
        validate_swap_criterion(plan.criterion, sim)
    end
    return nothing
end

# Whitelist of RomModel subtypes for which copy_model_state! overloads exist.
# SolidMechanics (FOM) is always supported and does not need to be listed here.
# Any RomModel subtype NOT in this list causes an early abort with a clear
# message, rather than a cryptic MethodError at the moment the swap fires.
const _SWAP_SUPPORTED_ROM_TYPES = Union{
    LinearOpInfRom,
    QuadraticOpInfRom,
    CubicOpInfRom,
    RBFKernelROM,
}

# FOM models are always swap-compatible; nothing to check.
_assert_rom_swap_supported(::SolidMechanics, ::String) = nothing

# Supported ROM subtypes: pass through.
_assert_rom_swap_supported(::_SWAP_SUPPORTED_ROM_TYPES, ::String) = nothing

# Catch-all: any other RomModel (including NeuralNetworkOpInfRom and any
# future subtype without a copy_model_state! overload) aborts here.
function _assert_rom_swap_supported(model::RomModel, context::String)
    norma_abortf(
        "%s uses a %s, which does not support mid-run swapping.  " *
        "Swap state transfer is implemented for: LinearOpInfRom, QuadraticOpInfRom, " *
        "CubicOpInfRom, RBFKernelROM.  To add support for a new ROM type, implement " *
        "copy_model_state!(dst::<YourRomType>, src::SolidMechanics), " *
        "copy_model_state!(dst::SolidMechanics, src::<YourRomType>), and " *
        "copy_model_state!(dst::<YourRomType>, src::<YourRomType>) in swap.jl, " *
        "then add <YourRomType> to _SWAP_SUPPORTED_ROM_TYPES.",
        context, typeof(model),
    )
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
#   2. At least one of those BCs has 'compute overlap L2 relative error: disp|velo|acce'
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
            "SchwarzOverlapL2SwapCriterion requires 'compute overlap L2 relative error: disp|velo|acce' on at " *
            "least one Schwarz overlap boundary condition, but no such BC was found.  Add " *
            "'compute overlap L2 relative error: disp' (or 'velo' or 'acce') to the relevant 'Schwarz overlap:' entry in " *
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
            "SchwarzOverlapL2SwapCriterion: no boundary conditions with 'compute overlap L2 relative error: disp|velo|acce' " *
            "found; criterion will not fire — add 'compute overlap L2 relative error: disp' (or 'velo' or 'acce') to the " *
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
        slot = sim.handle_by_name[plan.subsim_name].id
        # Evaluate the criterion against the specific subsim being replaced,
        # not the whole MultiDomainSimulation.  The multi-domain should_swap
        # dispatch aggregates over ALL subsims (e.g. with `all(...)`), which
        # is wrong here: only the target subsim's state is relevant for
        # deciding whether to swap it.
        should_swap(plan.criterion, sim.subsims[slot]) || continue
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

    # After state transfer, sync the new integrator's kinematic arrays from
    # the model.  ROM integrators hold reduced-space arrays that must be
    # copied from the model's reduced_state / reduced_velocity which
    # copy_model_state! just populated.
    _sync_integrator_from_model!(new.integrator, new.model)

    # Seed the new integrator's rollback buffers (same reason as single-domain).
    new.integrator.prev_disp  = copy(new.integrator.displacement)
    new.integrator.prev_velo  = copy(new.integrator.velocity)
    new.integrator.prev_acce  = copy(new.integrator.acceleration)
    new.integrator.prev_∂Ω_f = copy(get_internal_force(new.model))

    finalize_writing(old)
    initialize_writing(new)
    # Write the transferred state as the first frame of the new output file.
    # initialize_writing resets the per-file exodus_frame counter to 0, so
    # write_stop_exodus will write at time_index 1 regardless of the global
    # controller stop.  We set time to prev_time (the last converged stop) so
    # the frame carries the correct physical time stamp.
    let prev_time = sim.controller.prev_time
        new.controller.time = prev_time
        write_stop(new)
        new.controller.time = sim.controller.time  # restore for sync_control_time
    end

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
    # Apply BCs at the swap time before returning so that free_dofs is
    # correctly populated (constrained DOFs marked false) when the caller
    # subsequently invokes copy_model_state!.  Without this, free_dofs is
    # all-true (set in the SolidMechanics / OpInfModel constructor) and
    # _project_fom_to_rom! incorrectly includes constrained-node displacements
    # in the basis projection, producing a wrong initial reduced state and
    # negative Jacobians on the first FOM evaluate after a ROM→FOM swap.
    new.model.time = sim.controller.prev_time
    apply_bcs(new)
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

    # After state transfer, sync the new integrator's kinematic arrays from
    # the model.  For a FOM (SolidMechanics) the integrator holds unsafe_wrap
    # views into the model's arrays so no explicit sync is needed.  For a ROM
    # the integrator's displacement/velocity/acceleration are in the reduced
    # coordinate space and must be copied from the model's reduced_state /
    # reduced_velocity which copy_model_state! just populated.
    _sync_integrator_from_model!(new.integrator, new.model)

    # Seed the new integrator's "previous step" save-state buffers.
    # These are Float64[] after construction and normally filled by
    # save_curr_state at the end of the first successful step.  If that
    # first step fails (e.g. due to a difficult elastic-plastic transition),
    # restore_prev_state broadcasts into them and crashes with:
    #   DimensionMismatch: array could not be broadcast to match destination
    # Seeding with the just-transferred state gives a physically meaningful
    # fallback and prevents the crash.
    new.integrator.prev_disp  = copy(new.integrator.displacement)
    new.integrator.prev_velo  = copy(new.integrator.velocity)
    new.integrator.prev_acce  = copy(new.integrator.acceleration)
    new.integrator.prev_∂Ω_f = copy(get_internal_force(new.model))

    finalize_writing(sim)
    initialize_writing(new)
    # Write the transferred state as the first frame of the new output file
    # before overwriting sim's fields.  initialize_writing resets the per-file
    # exodus_frame counter to 0, so the frame lands at time_index 1.  We set
    # time to prev_time (the last converged stop) so the frame carries the
    # correct physical time stamp.
    let prev_time = sim.controller.prev_time
        new.controller.time = prev_time
        write_stop(new)
        new.controller.time = sim.controller.time  # restore
    end

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
    # Apply BCs at the swap time before returning so that free_dofs is
    # correctly populated (constrained DOFs marked false) when the caller
    # subsequently invokes copy_model_state!.  Without this, free_dofs is
    # all-true (set in the SolidMechanics / OpInfModel constructor) and
    # _project_fom_to_rom! incorrectly includes constrained-node displacements
    # in the basis projection, producing a wrong initial reduced state and
    # negative Jacobians on the first FOM evaluate after a ROM→FOM swap.
    new.model.time = sim.controller.prev_time
    apply_bcs(new)
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

# Safety net: if copy_model_state! is somehow called with a RomModel that has
# no specific overload (i.e. one not in _SWAP_SUPPORTED_ROM_TYPES), throw
# immediately with a clear message rather than silently doing nothing via the
# ::Model, ::Model no-op above.  validate_swap_plans should catch these at
# parse time, but this guards against future code paths that bypass validation.
function copy_model_state!(dst::Model, src::RomModel)
    src isa _SWAP_SUPPORTED_ROM_TYPES && return _copy_model_state_dispatch!(dst, src)
    error(
        "copy_model_state!: no state-transfer implementation for source ROM type $(typeof(src)).  " *
        "Add copy_model_state! overloads and register the type in _SWAP_SUPPORTED_ROM_TYPES.",
    )
end

function copy_model_state!(dst::RomModel, src::Model)
    dst isa _SWAP_SUPPORTED_ROM_TYPES && return _copy_model_state_dispatch!(dst, src)
    error(
        "copy_model_state!: no state-transfer implementation for destination ROM type $(typeof(dst)).  " *
        "Add copy_model_state! overloads and register the type in _SWAP_SUPPORTED_ROM_TYPES.",
    )
end

# Internal dispatch used by the catch-alls above once the type is confirmed
# supported.  Forwards to the concrete overloads defined further below.
_copy_model_state_dispatch!(dst::OpInfModel, src::SolidMechanics) = copy_model_state!(dst, src)
_copy_model_state_dispatch!(dst::SolidMechanics, src::OpInfModel) = copy_model_state!(dst, src)
_copy_model_state_dispatch!(dst::OpInfModel, src::OpInfModel)     = copy_model_state!(dst, src)
# Fallback for any other combination (e.g. two different non-SolidMechanics FOMs).
_copy_model_state_dispatch!(::Model, ::Model) = nothing

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

# ---------------------------------------------------------------------------
# ROM ↔ FOM and ROM ↔ ROM state transfer
# ---------------------------------------------------------------------------
#
# Design rationale
# ----------------
# An OpInf ROM stores kinematic state in the *reduced* coordinate q̂ ∈ ℝʳ
# (model.reduced_state / reduced_velocity) and keeps a "shadow" FOM in
# model.fom_model whose nodal fields are reconstructed from q̂ after every
# time step via reconstruct_fom_fields!.  Both the reduced and the full-order
# fields are therefore available at swap time.
#
# Four transitions are handled:
#
#   FOM  → FOM   (already covered by the SolidMechanics overload above)
#   FOM  → ROM   project the FOM kinematic fields onto the new ROM basis
#   ROM  → FOM   lift the ROM's shadow FOM fields into the destination FOM
#   ROM  → ROM   project the source shadow FOM fields onto the new ROM basis
#
# For FOM→ROM and ROM→ROM the projection uses the L2/Galerkin formula
#       q̂ₖ = Σᵢ Σⱼ Φ[j,i,k] * u_fom[j,i]
# where Φ has shape (n_var, n_node, n_mode) and u_fom has shape (n_var, n_node).
# This is the same inner product used in the ROM time integrators.
#
# When the source and destination share the *same* FOM mesh, the projection is
# exact up to machine precision.  When meshes differ (e.g. refining the FOM
# behind the ROM) we first L2-project the source full-order fields to the
# destination mesh before computing reduced coordinates, reusing the existing
# transfer_field machinery.
#
# The `time` field on the shadow FOM is also synced so that boundary-condition
# evaluations after the swap see the correct physical time.

# --- FOM → ROM ---------------------------------------------------------------
# Project a full-order displacement/velocity/acceleration into the ROM's
# reduced state and keep the ROM's shadow FOM consistent.
function copy_model_state!(dst::OpInfModel, src::SolidMechanics)
    norma_log(1, :swap, "State transfer: FOM → ROM (projecting FOM fields onto ROM basis)")
    dst_fom = dst.fom_model

    # Transfer full-order kinematic fields to the ROM's shadow FOM, using
    # mesh-to-mesh L2 projection when meshes differ, direct copy otherwise.
    if _meshes_match(dst_fom, src) && _state_layouts_match(dst_fom, src)
        _copy_model_state_direct!(dst_fom, src)
    else
        _ensure_recovery_for_swap!(dst_fom)
        _ensure_recovery_for_swap!(src)
        _copy_model_state_l2!(dst_fom, src)
    end

    # Project the (now up-to-date) shadow FOM fields onto the ROM basis.
    _project_fom_to_rom!(dst)
    return nothing
end

# --- ROM → FOM ---------------------------------------------------------------
# Lift the ROM's reconstructed full-order state into the destination FOM.
function copy_model_state!(dst::SolidMechanics, src::OpInfModel)
    norma_log(1, :swap, "State transfer: ROM → FOM (lifting ROM shadow FOM fields into FOM)")
    src_fom = src.fom_model

    # Use the shadow FOM of the source ROM as the origin.
    if _meshes_match(dst, src_fom) && _state_layouts_match(dst, src_fom)
        _copy_model_state_direct!(dst, src_fom)
    else
        _ensure_recovery_for_swap!(dst)
        _ensure_recovery_for_swap!(src_fom)
        _copy_model_state_l2!(dst, src_fom)
    end
    return nothing
end

# --- ROM → ROM ---------------------------------------------------------------
# Project the source ROM's reconstructed FOM fields onto the destination ROM basis.
function copy_model_state!(dst::OpInfModel, src::OpInfModel)
    norma_log(1, :swap, "State transfer: ROM → ROM (re-projecting onto new ROM basis)")
    src_fom = src.fom_model
    dst_fom = dst.fom_model

    # First, bring the destination shadow FOM up to date from the source.
    if _meshes_match(dst_fom, src_fom) && _state_layouts_match(dst_fom, src_fom)
        _copy_model_state_direct!(dst_fom, src_fom)
    else
        _ensure_recovery_for_swap!(dst_fom)
        _ensure_recovery_for_swap!(src_fom)
        _copy_model_state_l2!(dst_fom, src_fom)
    end

    # Project the destination shadow FOM fields onto the new ROM basis.
    _project_fom_to_rom!(dst)
    return nothing
end

# ---------------------------------------------------------------------------
# Core projection helper: FOM nodal fields → ROM reduced coordinates
# ---------------------------------------------------------------------------
#
# Computes the Galerkin projection
#   q̂ₖ = Σᵢ Σⱼ Φ[i,j,k] * u[i,j]
# for displacement and velocity, using the shadow FOM stored on `rom`.
# Acceleration is projected the same way.  The reduced coordinates on `rom`
# are updated in-place.  The shadow FOM fields are assumed to be current
# before this function is called.
function _project_fom_to_rom!(rom::OpInfModel)
    basis = rom.basis                   # (n_var, n_node_basis, n_mode)
    fom   = rom.fom_model
    _, _, n_mode = size(basis)

    u = fom.displacement  # (n_var, n_nodes_fom)
    v = fom.velocity
    a = fom.acceleration
    # Use the actual FOM node count, not the basis middle dimension.  The
    # basis array may have been generated with a different mesh or padded with
    # extra entries, so size(basis, 2) can exceed the number of FOM nodes and
    # must not be used to bound the node loop (it would produce out-of-bounds
    # accesses into fom.free_dofs).  This matches the pattern used everywhere
    # else in the codebase (e.g. reconstruct_fom_fields!).
    n_nodes_fom = size(u, 2)

    q̂_u = zeros(n_mode)
    q̂_v = zeros(n_mode)
    q̂_a = zeros(n_mode)

    for j in 1:n_nodes_fom
        x_dof = 3 * (j - 1) + 1
        y_dof = 3 * (j - 1) + 2
        z_dof = 3 * (j - 1) + 3
        for k in 1:n_mode
            if fom.free_dofs[x_dof]
                q̂_u[k] += basis[1, j, k] * u[1, j]
                q̂_v[k] += basis[1, j, k] * v[1, j]
                q̂_a[k] += basis[1, j, k] * a[1, j]
            end
            if fom.free_dofs[y_dof]
                q̂_u[k] += basis[2, j, k] * u[2, j]
                q̂_v[k] += basis[2, j, k] * v[2, j]
                q̂_a[k] += basis[2, j, k] * a[2, j]
            end
            if fom.free_dofs[z_dof]
                q̂_u[k] += basis[3, j, k] * u[3, j]
                q̂_v[k] += basis[3, j, k] * v[3, j]
                q̂_a[k] += basis[3, j, k] * a[3, j]
            end
        end
    end

    rom.reduced_state    .= q̂_u
    rom.reduced_velocity .= q̂_v
    # Sync the ROM's time from its shadow FOM.
    rom.time              = fom.time

    norma_logf(1, :swap,
        "_project_fom_to_rom!: projected %d FOM DOFs onto %d ROM modes",
        3 * n_nodes_fom, n_mode)
    return nothing
end

# ---------------------------------------------------------------------------
# Integrator ↔ model synchronisation after a swap
# ---------------------------------------------------------------------------
#
# For FOM (SolidMechanics) models the Newmark / CentralDifference integrators
# hold unsafe_wrap views into the model's displacement / velocity / acceleration
# arrays, so writing to the model arrays automatically updates the integrator.
# No explicit sync is needed and the FOM overload is a no-op.
#
# For ROM (OpInfModel) models the integrators hold *separate* reduced-space
# vectors of length n_modes that are NOT linked by view.  After copy_model_state!
# updates model.reduced_state and model.reduced_velocity we must copy those
# into the integrator so that the very first predict/correct step after the
# swap uses the transferred state.

# FOM: no action needed.
_sync_integrator_from_model!(::TimeIntegrator, ::SolidMechanics) = nothing

# ROM: copy reduced coordinates into the integrator's kinematic arrays.
function _sync_integrator_from_model!(integrator::TimeIntegrator, model::OpInfModel)
    n = length(model.reduced_state)
    # Guard against size mismatch (should not happen if the replacement YAML
    # specifies the same reduced dimension, but we abort with a clear message
    # rather than silently broadcasting wrong data).
    if length(integrator.displacement) != n
        norma_abortf(
            "_sync_integrator_from_model!: integrator displacement length %d does not " *
            "match ROM reduced_state length %d after swap.  Check that the replacement " *
            "ROM input file uses the same reduced dimension as the source model.",
            length(integrator.displacement), n,
        )
    end
    integrator.displacement .= model.reduced_state
    integrator.velocity     .= model.reduced_velocity
    # Acceleration is zeroed rather than projected: it would require a second
    # force evaluation on the new model which may not be valid yet.  The first
    # Newmark predict step will overwrite it with the correct value.
    fill!(integrator.acceleration, 0.0)
    return nothing
end
