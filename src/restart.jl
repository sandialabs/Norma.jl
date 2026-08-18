# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# In-memory restart: switch a subdomain between pre-built model variants
# (e.g. FOM <-> ROM) at any control step, decided on the fly.
#
# `swap.jl` does two separable things when its scheduled criterion fires: it
# transfers state (`copy_model_state!` plus the integrator and clock bookkeeping
# that must follow), and it rebuilds the subsim from disk. This file keeps the
# first and drops the second. Variants are built once, up front; a switch is then
# a state transfer and a pointer repoint -- no YAML re-read, no restart file, and
# the decision need not be known in advance.

"One entry per switch that actually changed a slot's occupant."
const SwitchRecord = NamedTuple{(:step, :slot, :from, :to),Tuple{Int64,Int64,String,String}}

"""
    InMemoryRestart(sim::MultiDomainSimulation)
    InMemoryRestart(input_file::String)

A resident simulation plus pre-built alternate variants for its subdomains.
`variants[slot]` maps a name to a constructed `SingleDomainSimulation`; exactly
one occupies `sim.subsims[slot]` at a time and the rest sit idle. The initial
occupant is registered under its own name, so it can be switched back to without
extra setup.

Construction runs Norma's start-of-`evolve` sequence (`sync_control_time`,
`initialize`), leaving the simulation ready to step.
"""
struct InMemoryRestart
    sim::MultiDomainSimulation
    variants::Vector{Dict{String,SingleDomainSimulation}}
    "Name of the variant currently occupying each slot."
    active::Vector{String}
    "Name each variant's input file asked for, before `uniquify_swap_output!` renamed it."
    intended::Vector{Dict{String,String}}
    "Whether this object opened the log file, and so is the one that must close it."
    owns_log::Bool
    history::Vector{SwitchRecord}
end

function InMemoryRestart(sim::MultiDomainSimulation; owns_log::Bool=false)
    sync_control_time(sim)
    initialize(sim)
    n = length(sim.subsims)
    variants = [Dict{String,SingleDomainSimulation}(sim.subsims[i].name => sim.subsims[i]) for i in 1:n]
    active = [sim.subsims[i].name for i in 1:n]
    intended = [Dict{String,String}(sim.subsims[i].name => sim.subsims[i].name) for i in 1:n]
    return InMemoryRestart(sim, variants, active, intended, owns_log, SwitchRecord[])
end

function InMemoryRestart(input_file::String)
    # `open_log_file` is a no-op under a session log or an already-open one, so
    # compare before and after rather than assume: only the call that actually
    # opened the file may close it in `close_restart!`.
    had_log = NORMA_LOG_FILE[] !== nothing
    open_log_file(input_file)   # as `Norma.run` does
    owns_log = !had_log && NORMA_LOG_FILE[] !== nothing
    return InMemoryRestart(create_simulation(input_file)::MultiDomainSimulation; owns_log=owns_log)
end

"""
    slot_of(r::InMemoryRestart, target) -> Int64

Resolve a subdomain from its slot index or the name of any subsim that has
occupied it. `switch!` keeps `sim.handle_by_name` current, so a subdomain stays
addressable by whichever variant name the caller finds convenient.
"""
function slot_of(r::InMemoryRestart, slot::Integer)
    1 <= slot <= length(r.sim.subsims) ||
        norma_abortf("In-memory restart: slot %d is out of range (1:%d)", slot, length(r.sim.subsims))
    return Int64(slot)
end

function slot_of(r::InMemoryRestart, name::AbstractString)
    key = String(name)
    haskey(r.sim.handle_by_name, key) ||
        norma_abort("In-memory restart: unknown subdomain '$key'")
    return r.sim.handle_by_name[key].id
end

"Name of the variant currently occupying `target`'s slot."
active_variant(r::InMemoryRestart, target) = r.active[slot_of(r, target)]

"Names registered for `target`, sorted, including the initial occupant."
variant_names(r::InMemoryRestart, target) = sort(collect(keys(r.variants[slot_of(r, target)])))

"""
    register_variant!(r::InMemoryRestart, target, name, yaml_file)

Build the subdomain described by `yaml_file` and register it for `target` under
`name`, without activating it. This is the expensive step (mesh read, model and
operator construction, BC wiring); doing it once is what makes switches cheap.

The variant must be a legal replacement for the slot: same mesh and degrees of
freedom, same Schwarz coupling, same integrator family, and for a ROM a
`copy_model_state!` path to and from whatever else may occupy it.

Only the last is checked here, by `_assert_rom_swap_supported`, since that would
otherwise raise a `MethodError` hundreds of steps into a run. A mesh or
degree-of-freedom mismatch is not caught until `copy_model_state!`.

Registering mid-march is well defined: `build_replacement_subsim` reads the
controller's current clock, so a variant built at step k is born aligned to it.
"""
function register_variant!(r::InMemoryRestart, target, name::AbstractString, yaml_file::AbstractString)
    slot = slot_of(r, target)
    key = String(name)
    haskey(r.variants[slot], key) &&
        norma_abort("In-memory restart: variant '$key' is already registered for slot $slot")
    isfile(yaml_file) ||
        norma_abort("In-memory restart: variant input file not found: $yaml_file")

    # Throwaway plan: `build_replacement_subsim` only reads `replacement_file`.
    plan = SwapPlan(r.sim.subsims[slot].name, TimeSwapCriterion(0.0), String(yaml_file), false)
    new = build_replacement_subsim(r.sim, slot, plan)
    _assert_rom_swap_supported(new.model, "In-memory restart variant '$key'")

    r.variants[slot][key] = new
    r.intended[slot][key] = stripped_name(yaml_file)
    norma_logf(0, :restart, "Registered in-memory variant '%s' for slot %d from %s", key, slot, yaml_file)
    return new
end

"""
    register_variants!(r::InMemoryRestart, target, pairs)

Register several variants for one subdomain, e.g.
`register_variants!(r, "sd1", ["rom" => "sd1rom.yaml"])`.
"""
function register_variants!(r::InMemoryRestart, target, pairs)
    for (name, file) in pairs
        register_variant!(r, target, name, file)
    end
    return r
end

"""
Register `name` as resolving to `slot`, refusing to steal it from another slot.
The in-memory counterpart of the collision guard in `apply_swap!`.
"""
function _claim_name!(sim::MultiDomainSimulation, slot::Int64, name::AbstractString)
    if haskey(sim.handle_by_name, name) && sim.handle_by_name[name].id != slot
        norma_abortf(
            "In-memory restart: switching slot %d would register subsim name '%s', but " *
            "that name is already claimed by slot %d. Two different subsims cannot share " *
            "a name while both are live; rename one of the variant input files.",
            slot, name, sim.handle_by_name[name].id,
        )
    end
    sim.handle_by_name[name] = DomainHandle(slot)
    return nothing
end

"""
    switch!(r::InMemoryRestart, target, name) -> Bool

Hand `target`'s state to its `name` variant and make that variant the occupant,
effective from the next solve. Returns `false` and does nothing if `name` is
already active, so it is safe to call unconditionally from a per-step hook.

This is the compute-only subset of `apply_swap!`: state transfer, clock and
integrator alignment, pointer repoint. Nothing is built or written.

Call it between `advance_control_time` and `sync_control_time` -- the only point
in the step where the controller knows the target time but no subdomain has
solved at it. `march!` calls its hook exactly there.
"""
function switch!(r::InMemoryRestart, target, name::AbstractString)
    slot = slot_of(r, target)
    key = String(name)
    haskey(r.variants[slot], key) ||
        norma_abort("In-memory restart: no variant '$key' registered for slot $slot " *
                    "(have: " * join(sort(collect(keys(r.variants[slot]))), ", ") * ")")
    sim = r.sim
    old = sim.subsims[slot]
    new = r.variants[slot][key]
    old === new && return false

    # A ROM's shadow FOM is refreshed only as a side effect of its Schwarz
    # partner reading it, so it can be one iteration stale. That lag is below
    # tolerance during ordinary stepping but not here, where it is treated as
    # authoritative -- reconstruct from the just-converged reduced state first.
    old.model isa RomModel && reconstruct_fom_fields!(old.integrator, old.solver, old.model)

    copy_model_state!(new.model, old.model)   # FOM->ROM projects, ROM->FOM lifts
    align_replacement_time!(new, sim, old)
    # ROM integrators hold reduced arrays that are not views into the model, so
    # the transferred state has to be pushed into them.
    _sync_integrator_from_model!(new.integrator, new.model)

    # Seed the rollback buffers: they are empty until the variant's first
    # successful step, and a step failure right after the switch would otherwise
    # broadcast into a zero-length array.
    new.integrator.prev_disp = copy(new.integrator.displacement)
    new.integrator.prev_velo = copy(new.integrator.velocity)
    new.integrator.prev_acce = copy(new.integrator.acceleration)
    new.integrator.prev_∂Ω_f = copy(get_internal_force(new.model))

    sim.subsims[slot] = new
    # Previous names survive as aliases (Schwarz partners resolve by slot id),
    # so `target` may keep being addressed by whichever name the caller started with.
    _claim_name!(sim, slot, new.name)
    sim.name_by_handle[slot] = new.name
    # If `uniquify_swap_output!` renamed the variant, keep the name its input
    # file asked for resolvable too.
    intended = get(r.intended[slot], key, new.name)
    intended == new.name || _claim_name!(sim, slot, intended)

    # Same rewiring `apply_swap!` does after re-pointing a slot: partner BCs
    # cache `coupled_bc_index` and `is_dirichlet` against the old BC list, and
    # the incoming variant's projectors have never been filled. Skipping this
    # leaves the neighbour coupled through stale indices -- measured as a 1e-4
    # disagreement against the scheduled-swap route on an otherwise identical
    # problem (see test/inmem-restart-vs-swap.jl).
    pair_schwarz_bcs(sim)
    initialize_bc_projectors(sim)

    from = r.active[slot]
    r.active[slot] = key
    push!(r.history, (step=sim.controller.stop, slot=slot, from=from, to=key))
    norma_logf(0, :restart, "Step %d: slot %d switched %s -> %s at t = %.6e",
               sim.controller.stop, slot, from, key, sim.controller.time)
    return true
end

"""
    step!(r::InMemoryRestart; on_step=nothing) -> Int64

Take one control step and return its index. Mirrors the body of Norma's `evolve`
loop with the scheduled-swap hook replaced by `on_step(r, step)`, called after
the controller has advanced to the target time but before any subdomain has
solved at it -- the window in which `switch!` is valid.

Kept separate from `evolve` so this feature does not touch the core loop; if
`evolve` gains a call, this copy must too. `test/inmem-restart-fidelity.jl`
compares a restart march against `Norma.run`, so a divergence surfaces there.

No Exodus output is written: each variant of a slot owns a separate output file,
and interleaving frames across them is the disk round-trip this avoids. Sample
the in-memory fields from the hook instead (see [`full_order_displacement`](@ref)).
"""
function step!(r::InMemoryRestart; on_step=nothing)
    sim = r.sim
    advance_control_time(sim)
    step = sim.controller.stop
    on_step === nothing || on_step(r, step)
    sync_control_time(sim)
    advance_control(sim)
    return step
end

"""
    march!(r::InMemoryRestart; on_step=nothing) -> InMemoryRestart

Step to the controller's final time, calling `on_step(r, step)` before each
step's physics. The hook may consult the current state (decision depends on the
solution), consult nothing (a fixed schedule), or be omitted (a plain march).
"""
function march!(r::InMemoryRestart; on_step=nothing)
    while true
        step!(r; on_step=on_step)
        stop_evolve(r.sim) && break
    end
    return r
end

"""
    full_order_displacement(subsim::SingleDomainSimulation) -> Matrix{Float64}

The `(3, n_node)` nodal displacement in full-order space, whichever variant is
active. For a ROM the shadow FOM is reconstructed from the current reduced state
first, so results are comparable across variants.

Returns the model's own buffer -- copy it before stepping again.
"""
function full_order_displacement(subsim::SingleDomainSimulation)
    if subsim.model isa RomModel
        reconstruct_fom_fields!(subsim.integrator, subsim.solver, subsim.model)
        return subsim.model.fom_model.displacement
    end
    return subsim.model.displacement
end

full_order_displacement(r::InMemoryRestart, target) = full_order_displacement(r.sim.subsims[slot_of(r, target)])

"""
    displacement_vector(r::InMemoryRestart) -> Vector{Float64}

Full-order displacements of every subdomain, concatenated in slot order. Length
and node ordering do not depend on which variants are active, which makes it
usable as a comparison vector between differently-switched runs.
"""
displacement_vector(r::InMemoryRestart) =
    vcat((vec(full_order_displacement(r.sim.subsims[i])) for i in 1:length(r.sim.subsims))...)

"""
    close_restart!(r::InMemoryRestart)

Close the Exodus handles of every variant, active or not, and the log file if
`InMemoryRestart(input_file)` was the call that opened it. Idle variants hold open
meshes just like active ones, so skipping this leaks file descriptors and a later
simulation in the same process will fail in `ex_create_int`. A close that fails is
logged rather than discarded: teardown must not throw, but must not hide a real error.
"""
function close_restart!(r::InMemoryRestart)
    for slot in eachindex(r.variants), subsim in values(r.variants[slot])
        try
            finalize_writing(subsim)
        catch e
            norma_logf(0, :warning, "In-memory restart: closing '%s' failed: %s",
                       subsim.name, sprint(showerror, e))
        end
    end
    r.owns_log && close_log_file()
    return nothing
end
