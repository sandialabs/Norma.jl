# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

using Base.Threads: @threads

# Whether `model` was constructed by resuming from a restart checkpoint (see
# process_restart!() / process_multidomain_restart!() in simulation.jl). For
# a ROM model this reflects whether the ROM's internal FOM model
# (model.fom_model) was restarted: ROM restart is layered on top of FOM
# restart — the FOM displacement/velocity fields are restored from the
# snapshot inside SolidMechanics() construction, then projected onto the
# reduced basis by apply_ics(::Parameters, ::RomModel, ...) in
# opinf_ics_bcs.jl — so the two always agree. Used by
# initialize(sim::MultiDomainSimulation) (simulation.jl) to decide whether
# the restart-only Schwarz refinement pass is needed, and forwarded as
# `include_schwarz_rows` into initialize(::TimeIntegrator, ::Solver, ::Model) /
# initialize(::RomNewmark, ...) / initialize(::RomCentralDifference, ...).
is_restarted(model::SolidMechanics) = model.restarted
is_restarted(model::RomModel) = model.fom_model.restarted

# The element blocks the model key `composite tet 10` switches to the
# composite tetrahedron: `true` for every TETRA10 block, or a list of block
# names, or a map of block names to booleans. Returns whether every TETRA10
# block is switched and the blocks named.
function composite_tet10_blocks(model_params, element_block_names::AbstractVector{<:AbstractString})
    spec = get(model_params, "composite tet 10", false)
    spec isa Bool && return spec, Set{String}()
    names = if spec isa AbstractVector
        Set{String}(string.(spec))
    elseif spec isa AbstractDict
        # A quoted "true" or "yes" would otherwise read as off.
        all(on -> on isa Bool, values(spec)) ||
            norma_abort("The values of a `composite tet 10` map of block names must be true or false.")
        Set{String}(string(name) for (name, on) in spec if on)
    else
        norma_abort("`composite tet 10` must be true, false, a list of block names or a map of block names to booleans.")
    end
    for name in names
        name ∈ element_block_names || norma_abortf("`composite tet 10` names block %s, which is not in the mesh.", name)
    end
    return false, names
end

# The element type of a block of the model: that of its mesh, except for the
# blocks switched to the composite tetrahedron.
function block_element_type(model::SolidMechanics, block_id::Integer)
    block_index = findfirst(block -> block.id == block_id, model.blocks)
    block_index === nothing && norma_abortf("No element block with id %d.", block_id)
    return model.blocks[block_index].element_type
end

block_element_types(model::SolidMechanics) = [block.element_type for block in model.blocks]

function SolidMechanics(params::Parameters)
    input_mesh = params["input_mesh"]
    model_params = params["model"]
    coords = read_coordinates(input_mesh)
    num_nodes = Exodus.num_nodes(input_mesh.init)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    restart_info = get(params, "restart_info", nothing)
    if restart_info === nothing
        displacement = zeros(3, num_nodes)
        velocity = zeros(3, num_nodes)
    else
        displacement = copy(restart_info.displacement)
        velocity = copy(restart_info.velocity)
    end
    acceleration = zeros(3, num_nodes)
    for node in 1:num_nodes
        reference[:, node] = coords[:, node]
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blocks_params = length(material_blocks)
    blocks = Exodus.read_sets(input_mesh, Block)
    num_blocks = length(blocks)
    if num_blocks_params ≠ num_blocks
        norma_abortf(
            "Number of element blocks in mesh '%s' (%d) does not match number of material blocks (%d).",
            input_mesh.file_name,
            num_blocks,
            num_blocks_params,
        )
    end
    element_block_names = Exodus.read_names(input_mesh, Block)
    materials = Vector{Solid}(undef, 0)
    kinematics = Undefined
    for element_block_name in element_block_names
        material_name = material_blocks[element_block_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        if kinematics == Undefined
            kinematics = get_kinematics(material_model)
        else
            if kinematics ≠ get_kinematics(material_model)
                norma_abortf(
                    "Material of type %s has inconsistent kinematics %s compared to previous materials of type %s.",
                    string(typeof(material_model)),
                    string(get_kinematics(material_model)),
                    string(kinematics),
                )
            end
        end
        push!(materials, material_model)
    end
    if restart_info !== nothing
        stateful = stateful_material_models(model_params)
        isempty(stateful) || norma_abort(
            "Restart is not currently supported for material models with internal state variables " *
            "($(join(("`" * model * "`" for model in stateful), ", "))). " *
            restart_internal_state_reason(),
        )
    end
    time = 0.0
    failed = false
    internal_force = zeros(3 * num_nodes)
    boundary_force = zeros(3 * num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(3 * num_nodes)
    stress = Vector{Vector{Vector{SVector{6,Float64}}}}()
    state_old = Vector{Vector{Vector{Vector{Float64}}}}()
    state = Vector{Vector{Vector{Vector{Float64}}}}()
    prev_state_old = Vector{Vector{Vector{Vector{Float64}}}}()  # empty until first save_curr_state
    stop_state_old = Vector{Vector{Vector{Vector{Float64}}}}()  # empty until first save_stop_state
    stored_energy = Vector{Vector{Float64}}()
    num_int_pts_overrides = get(model_params, "num integration points", Dict{String,Any}())
    all_composite, composite_blocks = composite_tet10_blocks(model_params, element_block_names)
    num_int_pts = Vector{Int}(undef, num_blocks)
    block_data = Vector{ElementBlockData}(undef, num_blocks)
    lattice_orientation = Vector{Vector{SMatrix{3,3,Float64,9}}}(undef, num_blocks)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string, num_block_elements, _, _, _, _ = Exodus.read_block_parameters(input_mesh, block_id)
        element_type = element_type_from_string(element_type_string)
        block_name = element_block_names[block_index]
        if block_name ∈ composite_blocks
            element_type == TETRA10 || norma_abortf(
                "Block %s is a %s block: the composite tet 10 needs TETRA10 elements.",
                block_name,
                element_type_string,
            )
            element_type = COMPOSITE_TETRA10
        elseif all_composite && element_type == TETRA10
            element_type = COMPOSITE_TETRA10
        end
        num_points = if haskey(num_int_pts_overrides, block_name)
            Int(num_int_pts_overrides[block_name])
        else
            default_num_int_pts(element_type)
        end
        num_int_pts[block_index] = num_points
        # Stored as one column per element (the flat Exodus order is element by
        # element), so an element's nodes are a contiguous column.
        raw_connectivity = get_block_connectivity(input_mesh, block_id)
        num_nodes_per_element = size(raw_connectivity, 2)
        connectivity = reshape(Int64.(vec(raw_connectivity)), (num_nodes_per_element, num_block_elements))
        N, dN, weights = isoparametric(element_type, num_points)
        mass_N, mass_dN, mass_weights = mass_quadrature(element_type, num_points)
        block_data[block_index] = ElementBlockData(
            Int64(block_id),
            element_type,
            num_points,
            num_block_elements,
            num_nodes_per_element,
            connectivity,
            N,
            dN,
            weights,
            mass_N,
            mass_dN,
            mass_weights,
        )
        material = materials[block_index]
        lattice_orientation[block_index] = if has_orientation(material)
            material_props = material_params[material_blocks[block_name]]
            block_lattice_orientations(material_props, input_mesh, block_id, num_block_elements)
        else
            SMatrix{3,3,Float64,9}[]
        end
        num_states = number_states(material)
        if (num_states > 0)
            point_init_state = initial_state(material)
        else
            point_init_state = Vector{Float64}()
        end
        block_stress = Vector{Vector{SVector{6,Float64}}}()
        block_state = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ in 1:num_block_elements
            element_stress = Vector{SVector{6,Float64}}()
            element_state = Vector{Vector{Float64}}()
            for _ in 1:num_points
                push!(element_stress, zero(SVector{6,Float64}))
                push!(element_state, copy(point_init_state))
            end
            push!(block_stress, element_stress)
            push!(block_state, element_state)
            element_stored_energy = 0.0
            push!(block_stored_energy, element_stored_energy)
        end
        push!(stress, block_stress)
        # state_old and state must be INDEPENDENT copies: state_old is the
        # converged state at the start of the step (what the constitutive
        # return mapping integrates from) and state is the trial state of the
        # current iterate. They used to alias the same nested arrays, so every
        # residual assembly committed its trial state instantly — each Newton
        # iteration (and each Schwarz iteration) restarted the return mapping from
        # the previous ITERATE's plastic state instead of the previous step's.
        # An early un-equilibrated iterate that fake-yields then corrupts
        # state_old and Newton stalls chasing a residual that moves with it.
        # The explicit commit lives in commit_state (simulation.jl), called
        # only when a step is accepted. A material without state has empty
        # vectors that nothing writes, shared as copy_state shares them: their
        # deep copy took 0.5 s of the setup on 110 592 hexahedra.
        push!(state_old, block_state)
        push!(state, num_states > 0 ? deepcopy(block_state) : block_state)
        push!(stored_energy, block_stored_energy)
    end
    if all_composite && !any(data -> data.element_type == COMPOSITE_TETRA10, block_data)
        norma_log(
            0,
            :warning,
            "`composite tet 10: true` switches no block: the mesh has no TETRA10 block, so the standard elements are used.",
        )
    end
    strain_energy = 0.0
    stiffness = spzeros(0, 0)
    mass = spzeros(0, 0)
    lumped_mass = Float64[]
    body_force = Float64[]
    compute_stiffness = true
    compute_mass = true
    compute_lumped_mass = true
    mesh_smoothing = get(params, "mesh smoothing", false)
    smooth_reference = get(model_params, "smooth reference", "")
    size_field = create_size_field(smooth_reference, get(model_params, "size field", nothing))
    metric_field = create_metric_field(smooth_reference, get(model_params, "metric field", nothing), input_mesh)
    for legacy in ("stress recovery", "recover internal variables")
        if haskey(model_params, legacy)
            norma_abort(
                "Legacy key '$legacy' is no longer supported. " *
                "Replace with a `nodal recovery:` block under `model:` " *
                "(method: lumped|consistent|both; stress|von mises stress|" *
                "internal variables|deformation gradient: true|false).",
            )
        end
    end
    recovery_params = get(model_params, "nodal recovery", nothing)
    recovery_kind = :none
    rec_stress = false
    rec_vm = false
    rec_iv = false
    rec_F = false
    if recovery_params !== nothing
        method = lowercase(string(get(recovery_params, "method", "")))
        if method == "lumped"
            recovery_kind = :lumped
        elseif method == "consistent"
            recovery_kind = :consistent
        elseif method == "both"
            recovery_kind = :both
        else
            norma_abort("nodal recovery: 'method' must be 'lumped', 'consistent', or 'both' " * "(got '$method').")
        end
        rec_stress = Bool(get(recovery_params, "stress", true))
        rec_vm = Bool(get(recovery_params, "von mises stress", false))
        rec_iv = Bool(get(recovery_params, "internal variables", false))
        rec_F = Bool(get(recovery_params, "deformation gradient", false))
    end
    recovery_data = build_recovery_data(
        recovery_kind, input_mesh, reference, num_int_pts, [block.element_type for block in block_data]
    )
    is_both = recovery_kind === :both
    n_iv = rec_iv ? length(collect_internal_variable_names(materials)) : 0
    # Single-mode buffers: allocated only when the quantity is enabled AND
    # the recovery is not BothRecovery (the latter uses the lumped_/consistent_
    # pairs below instead).
    recovered_stress = (rec_stress && !is_both) ? zeros(6, num_nodes) : zeros(0, 0)
    recovered_von_mises = (rec_vm && !is_both) ? zeros(1, num_nodes) : zeros(0, 0)
    recovered_F = (rec_F && !is_both) ? zeros(9, num_nodes) : zeros(0, 0)
    recovered_internal_variables = (n_iv > 0 && !is_both) ? zeros(n_iv, num_nodes) : zeros(0, 0)
    # BothRecovery buffers: lumped + consistent pair per enabled quantity.
    lumped_recovered_stress = (rec_stress && is_both) ? zeros(6, num_nodes) : zeros(0, 0)
    consistent_recovered_stress = (rec_stress && is_both) ? zeros(6, num_nodes) : zeros(0, 0)
    lumped_recovered_von_mises = (rec_vm && is_both) ? zeros(1, num_nodes) : zeros(0, 0)
    consistent_recovered_von_mises = (rec_vm && is_both) ? zeros(1, num_nodes) : zeros(0, 0)
    lumped_recovered_F = (rec_F && is_both) ? zeros(9, num_nodes) : zeros(0, 0)
    consistent_recovered_F = (rec_F && is_both) ? zeros(9, num_nodes) : zeros(0, 0)
    lumped_recovered_internal_variables = (n_iv > 0 && is_both) ? zeros(n_iv, num_nodes) : zeros(0, 0)
    consistent_recovered_internal_variables = (n_iv > 0 && is_both) ? zeros(n_iv, num_nodes) : zeros(0, 0)
    return SolidMechanics(
        input_mesh,
        materials,
        reference,
        displacement,
        velocity,
        acceleration,
        internal_force,
        boundary_force,
        boundary_conditions,
        state_old,
        state,
        prev_state_old,
        stop_state_old,
        stress,
        stored_energy,
        strain_energy,
        stiffness,
        mass,
        lumped_mass,
        body_force,
        free_dofs,
        time,
        compute_stiffness,
        compute_mass,
        compute_lumped_mass,
        failed,
        mesh_smoothing,
        smooth_reference,
        size_field,
        metric_field,
        kinematics,
        recovery_data,
        recovered_stress,
        recovered_von_mises,
        recovered_F,
        recovered_internal_variables,
        lumped_recovered_stress,
        consistent_recovered_stress,
        lumped_recovered_von_mises,
        consistent_recovered_von_mises,
        lumped_recovered_F,
        consistent_recovered_F,
        lumped_recovered_internal_variables,
        consistent_recovered_internal_variables,
        num_int_pts,
        block_data,
        lattice_orientation,
        create_element_value_buffers(block_data, num_nodes),
        nothing,
        restart_info !== nothing,
    )
end

# Maps a `model: type:` string to the Julia Model subtype used to construct
# it. The single source of truth create_model() (below) dispatches through;
# process_restart!() (simulation.jl) looks up the same mapping to resolve
# supports_restart() (model_types.jl) for a `model: type:` string before any
# model is actually constructed, instead of keeping a separate
# hand-maintained list of restart-capable type strings in sync with this
# function. Returns `nothing` for an unrecognized string. "mesh smoothing"
# maps to SolidMechanics too (it is a mode of that same model, selected via
# params["mesh smoothing"] in create_model() below, not a distinct type) --
# supports_restart() only sees the resolved Julia type, so process_restart!()
# still special-cases the "mesh smoothing" string directly wherever the
# distinction matters (mesh smoothing is not a stateful dynamic simulation
# and was never restart-capable).
function model_type_for(model_name::AbstractString)
    if model_name in ("solid mechanics", "mesh smoothing")
        return SolidMechanics
    elseif model_name in ("linear opinf rom", "linear kernel rom")
        return LinearOpInfRom
    elseif model_name in ("quadratic opinf rom", "quadratic kernel rom")
        return QuadraticOpInfRom
    elseif model_name in ("cubic opinf rom", "cubic kernel rom")
        return CubicOpInfRom
    elseif model_name == "neural network opinf rom"
        return NeuralNetworkOpInfRom
    elseif model_name == "rbf kernel rom"
        return RBFKernelROM
    else
        return nothing
    end
end

function create_model(params::Parameters)
    model_params = params["model"]
    model_name = model_params["type"]
    model_type = model_type_for(model_name)
    model_type === nothing && norma_abort("Unknown type of model : $model_name")
    if model_name == "mesh smoothing"
        params["mesh smoothing"] = true
    end
    return model_type(params)
end

# Compile a user-defined size field s(t, x, y, z) into a callable that returns
# the target edge length for the smooth reference element at a given location
# and time.  Returns `nothing` for smoothing modes that do not use a size field.
# Reuses the same Symbolics pipeline as boundary/initial conditions; the module
# variables t, x, y, z are declared in boundary_conditions.jl.
function create_size_field(smooth_reference::String, expression)
    if smooth_reference ∉ ("size field", "size field unrestricted")
        return nothing
    end
    if expression === nothing
        norma_abort(
            "smooth reference = \"$smooth_reference\" requires a \"size field\" expression under the model parameters"
        )
    end
    size_num = eval(Meta.parse(string(expression)))
    return eval(build_function(size_num, [t, x, y, z]; expression=Val(false)))
end

# Ideal reference element and metric factor of one element of a smoothing
# model, from the element's coordinates in the mesh that samples the target
# (the original mesh during a smoothing phase, the current mesh at a topology
# phase, docs/notes/ems-adaptivity) and its node indices, which a nodal
# metric needs; `metric` overrides the model's field, for a proposal with a
# trial node.  Without a metric field the factor is the identity and the
# ideal element follows the scalar rule.
function smoothing_reference(
    model::SolidMechanics,
    element_type::ElementType,
    sample_pos::AbstractMatrix{Float64};
    node_indices=nothing,
    metric::Union{MetricField,Nothing}=model.metric_field,
)
    element_type == TETRA4 || norma_abort("Mesh smoothing supports four-node tetrahedra only")
    X = SMatrix{3,4,Float64,12}(sample_pos)
    if metric !== nothing
        return create_metric_reference(metric, X, model.time; node_indices=node_indices)
    end
    I3 = SMatrix{3,3,Float64,9}(I)
    return smoothing_size(model.smooth_reference, X, model.size_field, model.time) * UNIT_TETRAHEDRON, I3, I3
end

# Volume of the ideal reference element of every element of a block, for the
# energy density written with the smoothing output.
function ideal_volumes(model::SolidMechanics, block_index::Int)
    block = model.blocks[block_index]
    # Behind a barrier on the untyped shape functions, as element_energies is.
    return ideal_volumes(model, block.element_type, block.connectivity, block.N)
end

function ideal_volumes(model::SolidMechanics, element_type::ElementType, connectivity::Matrix{<:Integer}, N::SMatrix)
    volumes = Vector{Float64}(undef, size(connectivity, 2))
    for element in 1:size(connectivity, 2)
        node_indices = view(connectivity, :, element)
        X, _, _ = smoothing_reference(model, element_type, gather_nodal(model.reference, node_indices, N); node_indices)
        volumes[element] = abs(tetrahedron_volume(X))
    end
    return volumes
end

# Energies of a set of elements given by their connectivity into `positions`,
# evaluated exactly as the assembly evaluates them (same shape functions,
# quadrature, ideal element, and metric-space deformation gradient), for the
# acceptance test of the adaptivity loop: the elements need not exist in the
# mesh.  `sample_positions` holds the coordinates at which the target is
# sampled; at a topology phase it is the same array as `positions`, during a
# smoothing phase it is the original mesh.  An inverted or degenerate element
# gets an infinite energy.
const DEGENERATE_JACOBIAN = 1.0e-12

function element_energies(
    model::SolidMechanics,
    block_index::Int,
    connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64};
    sample_positions::AbstractMatrix{Float64}=positions,
    metric::Union{MetricField,Nothing}=model.metric_field,
)
    block = model.blocks[block_index]
    material = model.materials[block_index]
    Q = adaptivity_orientation(model, material, block_index)
    # The shape functions are stored untyped on the block, so the loop runs
    # behind a barrier that specializes on their concrete types; without it
    # every operation in it is dispatched at run time and allocates.
    return element_energies(
        model,
        block.element_type,
        material,
        Q,
        block.N,
        block.dN,
        block.weights,
        connectivity,
        positions,
        sample_positions,
        metric,
    )
end

# The lattice orientation with which the adaptivity loop evaluates the
# energies of a block, whose new elements exist in no mesh: the block's one
# orientation, or the identity for a material without one. The energy of a
# material with internal state variables needs a state that new elements do
# not have, and such a material has no stored energy method.
function adaptivity_orientation(model::SolidMechanics, material::Material, block_index::Int)
    number_states(material) > 0 && norma_abortf(
        "The adaptivity of mesh smoothing evaluates the stored energy of new elements, which have no internal " *
        "state, so the material of block %d, which has internal state variables, cannot be used with it.",
        block_index,
    )
    has_orientation(material) || return SMatrix{3,3,Float64,9}(I)
    orientations = model.lattice_orientation[block_index]
    isempty(orientations) && return SMatrix{3,3,Float64,9}(I)
    all(==(orientations[1]), orientations) || norma_abortf(
        "The adaptivity of mesh smoothing creates new elements, so block %d needs one lattice orientation, " *
        "not one per element.",
        block_index,
    )
    return orientations[1]
end

function element_energies(
    model::SolidMechanics,
    element_type::ElementType,
    material::Material,
    Q::SMatrix{3,3,Float64,9},
    N::SMatrix,
    dN::SArray,
    ip_weights::AbstractVector,
    connectivity::AbstractMatrix{<:Integer},
    positions::AbstractMatrix{Float64},
    sample_positions::AbstractMatrix{Float64},
    metric::Union{MetricField,Nothing},
)
    num_points = size(N, 2)
    num_elements = size(connectivity, 2)
    energies = zeros(num_elements)
    for element in 1:num_elements
        node_indices = view(connectivity, :, element)
        sample = tetrahedron_coordinates(sample_positions, node_indices)
        X, F_M, F_M_inv = smoothing_reference(model, element_type, sample; node_indices, metric)
        element_reference_position = gather_nodal(X, N)
        element_current_position = tetrahedron_coordinates(positions, node_indices)
        energy = 0.0
        for point in 1:num_points
            dNdξ = dN[:, :, point]
            # The orientation is checked on the element itself: the
            # equal-volume rule would build a reflected ideal element for an
            # inverted element and hide the inversion in a positive Jacobian.
            dxdξ = dNdξ * element_current_position'
            dXdξ = dNdξ * element_reference_position'
            dNdX = dXdξ \ dNdξ
            F = element_current_position * dNdX'
            # Both determinants are checked: for a nearly flat element they
            # can differ in sign by roundoff, and the energy takes fractional
            # powers of det(F).  The energy computes det(F) again, and the
            # fused multiply-add contraction of the kernels can round the two
            # evaluations differently, so a Jacobian within roundoff of zero
            # counts as degenerate rather than being passed on.
            J = det(F)
            if det(dxdξ) ≤ 0.0 || J ≤ DEGENERATE_JACOBIAN || isfinite(J) == false
                energy = Inf
                break
            end
            F_U = F_M * F * F_M_inv
            # An oriented material is evaluated in its lattice frame, as the
            # assembly evaluates it.
            W = has_orientation(material) ? strain_energy(material, Q' * F_U * Q) : strain_energy(material, F_U)
            energy += W * det(dXdξ) * ip_weights[point]
        end
        energies[element] = energy
    end
    return energies
end

# Pull the tangent of the energy in the metric space back to the reference
# configuration: with F_U = F_M F F_M⁻¹ and P = F_M' P_U F_M⁻ᵀ,
# A_iJkL = (F_M)_pi (F_M⁻¹)_Jq A^U_pqrs (F_M)_rk (F_M⁻¹)_Ls.
function pull_back_tangent(
    AA::SArray{Tuple{3,3,3,3},Float64,4,81}, F_M::SMatrix{3,3,Float64,9}, F_M_inv::SMatrix{3,3,Float64,9}
)
    B = SArray{Tuple{3,3,3,3},Float64,4,81}(
        sum(AA[p, q, r, s] * F_M[r, k] * F_M_inv[L, s] for r in 1:3, s in 1:3) for
        p in 1:3, q in 1:3, k in 1:3, L in 1:3
    )
    return SArray{Tuple{3,3,3,3},Float64,4,81}(
        sum(F_M[p, i] * F_M_inv[J, q] * B[p, q, k, L] for p in 1:3, q in 1:3) for i in 1:3, J in 1:3, k in 1:3, L in 1:3
    )
end

# Target size of the ideal element of a TETRA4 under a smoothing rule, from
# the element's coordinates in the mesh that samples the target.
function smoothing_size(
    smooth_reference::String,
    element_ref_pos::SMatrix{3,4,Float64,12},
    size_field::Union{Function,Nothing},
    time::Float64,
)
    u = element_ref_pos[:, 2] - element_ref_pos[:, 1]
    v = element_ref_pos[:, 3] - element_ref_pos[:, 1]
    w = element_ref_pos[:, 4] - element_ref_pos[:, 1]
    if smooth_reference == "equal volume"
        return equal_volume_tet_h(u, v, w)
    elseif smooth_reference == "average edge length"
        return avg_edge_length_tet_h(u, v, w)
    elseif smooth_reference == "max"
        return max(avg_edge_length_tet_h(u, v, w), equal_volume_tet_h(u, v, w))
    elseif smooth_reference == "size field"
        # Target edge length from the user-defined size field at the element
        # reference centroid, combined with the volume criterion (max) to
        # anchor the reference size and avoid sliver pathologies.
        return max(size_field_tet_h(size_field, element_ref_pos, time), equal_volume_tet_h(u, v, w))
    elseif smooth_reference == "size field unrestricted"
        # Target edge length from the user-defined size field at the element
        # reference centroid
        return size_field_tet_h(size_field, element_ref_pos, time)
    elseif smooth_reference == "metric field" || smooth_reference == "metric field unrestricted"
        norma_abort("A metric field reference is built by create_metric_reference, not create_smooth_reference")
    end
    return norma_abort("Unknown type of mesh smoothing reference : $smooth_reference")
end

function create_smooth_reference(
    smooth_reference::String,
    element_type::ElementType,
    element_ref_pos::AbstractMatrix{Float64},
    size_field::Union{Function,Nothing}=nothing,
    time::Float64=0.0,
)
    element_type == TETRA4 || norma_abort("Unknown element type")
    X = SMatrix{3,4,Float64,12}(element_ref_pos)
    return smoothing_size(smooth_reference, X, size_field, time) * UNIT_TETRAHEDRON
end

function equal_volume_tet_h(u::SVector{3,Float64}, v::SVector{3,Float64}, w::SVector{3,Float64})
    return cbrt(sqrt(2.0) * dot(u, cross(v, w)))
end

function avg_edge_length_tet_h(u::SVector{3,Float64}, v::SVector{3,Float64}, w::SVector{3,Float64})
    return (norm(u) + norm(v) + norm(w) + norm(u - v) + norm(u - w) + norm(v - w)) / 6.0
end

# The compiled size field called through a function barrier, so that the
# value is a Float64 to the caller although the field is stored untyped.
evaluate_size_field(size_field::Function, args::NTuple{4,Float64})::Float64 = Float64(size_field(args))

# Target edge length from the user-defined size field, evaluated at the element
# reference centroid.
function size_field_tet_h(size_field::Union{Function,Nothing}, element_ref_pos::SMatrix{3,4,Float64,12}, time::Float64)
    if size_field === nothing
        norma_abort(
            "smooth reference = \"size field\" or \"size field unrestricted\" selected but no size field was compiled"
        )
    end
    centroid = (element_ref_pos[:, 1] + element_ref_pos[:, 2] + element_ref_pos[:, 3] + element_ref_pos[:, 4]) / 4.0
    h = evaluate_size_field(size_field, (time, centroid[1], centroid[2], centroid[3]))
    if !isfinite(h) || h ≤ 0.0
        norma_abort(
            "Size field must be strictly positive and finite; got $h at centroid " *
            "($(centroid[1]), $(centroid[2]), $(centroid[3])) and time $time",
        )
    end
    return h
end

# The per-point material state is copied at several points of a step. For
# materials without internal variables the nested vectors are all empty and
# nothing ever writes into them, so the copy is skipped and the arrays are
# shared. The copies were a fifth of an explicit step for elastic materials.
has_material_state(model::SolidMechanics) = any(m -> number_states(m) > 0, model.materials)

function copy_state(model::SolidMechanics, state::Vector{Vector{Vector{Vector{Float64}}}})
    return has_material_state(model) ? deepcopy(state) : state
end

# Cap the integrator's nominal step at the stable step. advance_time applies
# the same cap to every substep it takes, so the kinematic update and the
# time stamp always use the same step and a stop longer than the stable step
# is reached by subcycling.
function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    stable = stable_time_step(integrator, model)
    if stable < integrator.time_step
        norma_logf(0, :warning, "Δt = %.3e exceeds stable Δt = %.3e; using stable step.", integrator.time_step, stable)
    end
    integrator.time_step = min(stable, integrator.time_step)
    return nothing
end

function voigt_cauchy_from_stress(_::Solid, P::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9}, J::Float64)
    σ = F * P' ./ J
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function voigt_cauchy_from_stress(
    _::Union{Linear_Elastic,AnisotropicLinearElastic}, σ::SMatrix{3,3,Float64,9}, _::SMatrix{3,3,Float64,9}, _::Float64
)
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

@generated function create_element_matrix(::Type{T}, ::Val{N}) where {T,N}
    dof_per_node = 3
    total_dofs = dof_per_node * N
    quote
        MMatrix{$total_dofs,$total_dofs,$T}(undef)
    end
end

@generated function create_element_vector(::Type{T}, ::Val{N}) where {T,N}
    dof_per_node = 3
    total_dofs = dof_per_node * N
    quote
        MVector{$total_dofs,$T}(undef)
    end
end

# Internal force from the first Piola-Kirchhoff stress: for node a and
# direction i, f[(a,i)] += P[i,k] dN_a/dX_k dV, written as loops over the
# static arrays so that no operator matrix is formed.
function add_internal_force!(Fi::MVector{M,T}, dNdX::SMatrix{3,N,T}, P::SMatrix{3,3,T,9}, dV::T) where {M,N,T}
    @inbounds for a in 1:N
        base = 3 * (a - 1)
        for i in 1:3
            s = P[i, 1] * dNdX[1, a] + P[i, 2] * dNdX[2, a] + P[i, 3] * dNdX[3, a]
            Fi[base + i] += s * dV
        end
    end
    return nothing
end

# Element stiffness K[(a,i),(b,j)] += dN_a/dX_k A[i,k,j,l] dN_b/dX_l dV,
# accumulated in place. The previous form built the 9 x 3N gradient operator
# and the product B' A B as static matrices, which for eight nodes is a
# 24 x 24 result beyond the size StaticArrays keeps on the stack: 3.7 µs and
# 4.7 kB of allocation per integration point against about 0.5 µs here.
function add_stiffness!(K::MMatrix{M,M,T}, dNdX::SMatrix{3,N,T}, AA::SArray{Tuple{3,3,3,3},T}, dV::T) where {M,N,T}
    # A[i,k] is the 3 x 3 matrix AA[i,k,:,:], extracted once per point.
    A = ntuple(m -> begin
        i = (m - 1) % 3 + 1
        k = (m - 1) ÷ 3 + 1
        SMatrix{3,3,T,9}(ntuple(n -> AA[i, k, (n - 1) % 3 + 1, (n - 1) ÷ 3 + 1], Val(9)))
    end, Val(9))
    @inbounds for a in 1:N
        ra = 3 * (a - 1)
        g1, g2, g3 = dNdX[1, a], dNdX[2, a], dNdX[3, a]
        for i in 1:3
            # W[j,l] = dN_a/dX_k A[i,k,j,l]; row block Kb[j,b] = W[j,l] dN_b/dX_l
            W = g1 * A[i] + g2 * A[i + 3] + g3 * A[i + 6]
            Kb = W * dNdX
            for b in 1:N
                rb = 3 * (b - 1)
                K[ra + i, rb + 1] += Kb[1, b] * dV
                K[ra + i, rb + 2] += Kb[2, b] * dV
                K[ra + i, rb + 3] += Kb[3, b] * dV
            end
        end
    end
    return nothing
end

# Consistent mass: M[(a,i),(b,i)] += ρ N_a N_b dV for each direction i.
function add_mass!(Me::MMatrix{M,M,T}, Np::SVector{N,T}, density::T, dV::T) where {M,N,T}
    @inbounds for b in 1:N
        rb = 3 * (b - 1)
        for a in 1:N
            ra = 3 * (a - 1)
            m = Np[a] * Np[b] * density * dV
            Me[ra + 1, rb + 1] += m
            Me[ra + 2, rb + 2] += m
            Me[ra + 3, rb + 3] += m
        end
    end
    return nothing
end

function add_lumped_mass!(M::MVector{R,T}, Nξ::SVector{N,T}, density::T, dV::T) where {R,N,T}
    @assert R == 3N
    s = sum(Nξ)
    w = (density * dV) .* Nξ .* s

    @inbounds for a in 1:N
        idx = 3 * (a - 1) + 1
        m = w[a]
        M[idx] += m
        M[idx + 1] += m
        M[idx + 2] += m
    end
    return nothing
end

# Degrees of freedom of an element, three per node in node order.
@inline function fill_dofs!(dofs::AbstractVector{Int64}, connectivity::Matrix{Int64}, element::Integer)
    @inbounds for (a, node) in enumerate(view(connectivity, :, element))
        base = 3 * (a - 1)
        dofs[base + 1] = 3 * node - 2
        dofs[base + 2] = 3 * node - 1
        dofs[base + 3] = 3 * node
    end
    return nothing
end

function create_element_value_buffers(blocks::Vector{ElementBlockData}, num_nodes::Integer)
    force = [zeros(3 * block.num_nodes_per_element * block.num_elements) for block in blocks]
    lumped_mass = [zeros(3 * block.num_nodes_per_element * block.num_elements) for block in blocks]
    dof_target = Vector{Vector{Int64}}(undef, length(blocks))
    for (block_index, block) in enumerate(blocks)
        nd = 3 * block.num_nodes_per_element
        dofs = Vector{Int64}(undef, nd)
        targets = Vector{Int64}(undef, nd * block.num_elements)
        k = 0
        for element in 1:block.num_elements
            fill_dofs!(dofs, block.connectivity, element)
            @inbounds for i in 1:nd
                k += 1
                targets[k] = dofs[i]
            end
        end
        dof_target[block_index] = targets
    end
    dof_ptr, dof_block, dof_source = inverse_map(dof_target, 3 * Int64(num_nodes))
    return ElementValueBuffers(force, lumped_mass, dof_target, dof_ptr, dof_block, dof_source)
end

# The sparsity pattern is that of the sum over elements of the element
# matrices; it is built once from the connectivity. For every element entry
# the position in the value array is found by a search in its column.
function build_matrix_pattern(model::SolidMechanics)
    blocks = model.blocks
    num_dofs = 3 * size(model.reference, 2)
    total = sum(block.num_elements * (3 * block.num_nodes_per_element)^2 for block in blocks; init=0)
    rows = Vector{Int64}(undef, total)
    cols = Vector{Int64}(undef, total)
    k = 0
    for block in blocks
        nd = 3 * block.num_nodes_per_element
        dofs = Vector{Int64}(undef, nd)
        for element in 1:block.num_elements
            fill_dofs!(dofs, block.connectivity, element)
            @inbounds for j in 1:nd, i in 1:nd
                k += 1
                rows[k] = dofs[i]
                cols[k] = dofs[j]
            end
        end
    end
    structure = sparse(rows, cols, ones(total), num_dofs, num_dofs)
    colptr = structure.colptr
    rowval = structure.rowval
    slots = Vector{Vector{Int64}}(undef, length(blocks))
    for (block_index, block) in enumerate(blocks)
        nd = 3 * block.num_nodes_per_element
        dofs = Vector{Int64}(undef, nd)
        block_slots = Vector{Int64}(undef, block.num_elements * nd * nd)
        k = 0
        for element in 1:block.num_elements
            fill_dofs!(dofs, block.connectivity, element)
            @inbounds for j in 1:nd
                col = dofs[j]
                lo = colptr[col]
                hi = colptr[col + 1] - 1
                column_rows = view(rowval, lo:hi)
                for i in 1:nd
                    k += 1
                    block_slots[k] = lo - 1 + searchsortedfirst(column_rows, dofs[i])
                end
            end
        end
        slots[block_index] = block_slots
    end
    stiffness = [zeros(length(block_slots)) for block_slots in slots]
    # The mass values are sized on first use: quasi-static runs never need them.
    mass = [Float64[] for _ in slots]
    entry_ptr, entry_block, entry = inverse_map(slots, length(rowval))
    return MatrixPattern(colptr, rowval, slots, stiffness, mass, entry_ptr, entry_block, entry)
end

# Copy the element arrays of one element into its value segments. Segments of
# different elements are disjoint, so this is safe from any thread.
function store_element_values!(
    arrays::ElementArrays,
    block_element_index::Int64,
    force_values::Vector{Float64},
    lumped_mass_values::Vector{Float64},
    stiffness_values::Vector{Float64},
    mass_values::Vector{Float64},
    flags::EvaluationFlags,
)
    nd = length(arrays.dofs)
    offset = (block_element_index - 1) * nd
    fe = arrays.internal_force
    @inbounds for i in 1:nd
        force_values[offset + i] = fe[i]
    end
    if flags.compute_lumped_mass == true
        me = arrays.lumped_mass
        @inbounds for i in 1:nd
            lumped_mass_values[offset + i] = me[i]
        end
    end
    n2 = nd * nd
    offset2 = (block_element_index - 1) * n2
    if flags.compute_stiffness == true
        Ke = arrays.stiffness
        @inbounds for k in 1:n2
            stiffness_values[offset2 + k] = Ke[k]
        end
    end
    if flags.compute_mass == true
        Me = arrays.mass
        @inbounds for k in 1:n2
            mass_values[offset2 + k] = Me[k]
        end
    end
    return nothing
end

# A vector of length n, zeroed: the given one when it has that length, else a
# new one.
function sized_vector!(v::Vector{Float64}, n::Int64)
    length(v) == n || return zeros(n)
    fill!(v, 0.0)
    return v
end

# The global vector from the element value segments. On one thread the values
# are scattered into the vector in element order. On several, every degree of
# freedom gathers its sources (dof_source) in element order, so the sums are
# the same as those of the scatter, and the degrees of freedom are divided
# among the threads without conflicts. The gather reads the values out of
# order, which on one thread takes twice the time of the scatter.
function reduce_element_vector(buffers::ElementValueBuffers, values::Vector{Vector{Float64}}, num_dofs::Int64)
    return reduce_element_vector!(Vector{Float64}(undef, num_dofs), buffers, values)
end

function reduce_element_vector!(
    global_vector::Vector{Float64}, buffers::ElementValueBuffers, values::Vector{Vector{Float64}}
)
    num_dofs = length(global_vector)
    ptr = buffers.dof_ptr
    block = buffers.dof_block
    source = buffers.dof_source
    if Threads.nthreads() == 1
        fill!(global_vector, 0.0)
        @inbounds for block_index in eachindex(values)
            block_values = values[block_index]
            target = buffers.dof_target[block_index]
            for k in eachindex(block_values)
                global_vector[target[k]] += block_values[k]
            end
        end
        return global_vector
    end
    @threads :static for dof in 1:num_dofs
        s = 0.0
        @inbounds for p in ptr[dof]:(ptr[dof + 1] - 1)
            s += values[block[p]][source[p]]
        end
        @inbounds global_vector[dof] = s
    end
    return global_vector
end

# The global matrix from the element value segments, as the vector above: on
# one thread the scatter through the slots, on several every nonzero gathers
# its entries in element order. The values are written into `A` when it has
# the structure of the pattern (the matrix of the previous evaluation); a
# fresh vector of 3.4 million values per evaluation was 27 MB on the 24³ cube.
function reduce_element_matrix(pattern::MatrixPattern, values::Vector{Vector{Float64}}, num_dofs::Int64)
    return reduce_element_matrix!(spzeros(0, 0), pattern, values, num_dofs)
end

function reduce_element_matrix!(
    A::SparseMatrixCSC{Float64,Int64}, pattern::MatrixPattern, values::Vector{Vector{Float64}}, num_dofs::Int64
)
    num_nonzeros = length(pattern.rowval)
    reuse = A.colptr === pattern.colptr && A.rowval === pattern.rowval && length(A.nzval) == num_nonzeros
    nzval = reuse ? A.nzval : Vector{Float64}(undef, num_nonzeros)
    ptr = pattern.entry_ptr
    block = pattern.entry_block
    entry = pattern.entry
    if Threads.nthreads() == 1
        fill!(nzval, 0.0)
        @inbounds for block_index in eachindex(values)
            block_values = values[block_index]
            block_slots = pattern.slots[block_index]
            for k in eachindex(block_values)
                nzval[block_slots[k]] += block_values[k]
            end
        end
        return SparseMatrixCSC(num_dofs, num_dofs, pattern.colptr, pattern.rowval, nzval)
    end
    @threads :static for j in 1:num_nonzeros
        s = 0.0
        @inbounds for p in ptr[j]:(ptr[j + 1] - 1)
            s += values[block[p]][entry[p]]
        end
        @inbounds nzval[j] = s
    end
    return SparseMatrixCSC(num_dofs, num_dofs, pattern.colptr, pattern.rowval, nzval)
end

# K + M / β / Δt² for the Newmark Hessian. Matrices built from the same pattern
# share their structure arrays, so the sum is a pass over the values with the
# same elementwise arithmetic as the generic expression. The values are written
# into `H` when it shares the structure of K and holds values of its own (the
# Hessian of the previous Newton iteration): the vector was allocated at every
# iteration before (28 MB on 46 875 degrees of freedom, with the collection of
# those bytes 14% of a coupled run on 12 threads). The reduced system and the
# direct solver recognize the structure by identity, so the matrix returned
# keeps the arrays of K.
function newmark_hessian!(
    H::SparseMatrixCSC{Float64,Int64},
    K::SparseMatrixCSC{Float64,Int64},
    M::SparseMatrixCSC{Float64,Int64},
    β::Float64,
    Δt::Float64,
)
    if K.colptr === M.colptr && K.rowval === M.rowval
        K_values = K.nzval
        M_values = M.nzval
        reuse =
            H.colptr === K.colptr &&
            H.rowval === K.rowval &&
            H.nzval !== K_values &&
            H.nzval !== M_values &&
            length(H.nzval) == length(K_values)
        H_values = reuse ? H.nzval : Vector{Float64}(undef, length(K_values))
        chunk_results(Nothing, () -> nothing, length(H_values); chunk=VECTOR_CHUNK) do range, _
            @inbounds for k in range
                H_values[k] = K_values[k] + M_values[k] / β / Δt / Δt
            end
            return nothing
        end
        return reuse ? H : SparseMatrixCSC(K.m, K.n, K.colptr, K.rowval, H_values)
    end
    return K + M / β / Δt / Δt
end

function newmark_hessian(K::SparseMatrixCSC{Float64,Int64}, M::SparseMatrixCSC{Float64,Int64}, β::Float64, Δt::Float64)
    return newmark_hessian!(spzeros(0, 0), K, M, β, Δt)
end

function compute_flags(model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    is_implicit_dynamic = integrator isa Newmark
    is_explicit_dynamic = integrator isa CentralDifference
    is_implicit_static = integrator isa QuasiStatic
    is_dynamic = is_implicit_dynamic || is_explicit_dynamic
    is_implicit = is_implicit_dynamic || is_implicit_static
    is_hessian_opt = solver isa HessianMinimizer
    is_matrix_free = solver isa SteepestDescent
    need_lumped_mass = is_explicit_dynamic || (is_implicit_dynamic && is_matrix_free)
    need_stiffness = is_implicit && is_hessian_opt
    need_mass = is_dynamic && is_hessian_opt
    compute_lumped_mass = need_lumped_mass && model.compute_lumped_mass
    compute_stiffness = need_stiffness && model.compute_stiffness
    compute_mass = need_mass && model.compute_mass
    mesh_smoothing = model.mesh_smoothing

    return EvaluationFlags(
        is_dynamic,
        is_implicit,
        is_hessian_opt,
        is_matrix_free,
        need_lumped_mass,
        need_stiffness,
        need_mass,
        compute_lumped_mass,
        compute_stiffness,
        compute_mass,
        mesh_smoothing,
    )
end

# The element arrays of one task of the assembly: the degrees of freedom, the
# force, and the matrices of the element being evaluated. Arrays that the flags
# do not ask for have size zero.
function create_element_arrays(valN::Val{N}, flags::EvaluationFlags) where {N}
    dofs = create_element_vector(Int64, valN)
    internal_force = create_element_vector(Float64, valN)
    lumped_mass = create_element_vector(Float64, flags.compute_lumped_mass ? valN : Val(0))
    stiffness = create_element_matrix(Float64, flags.compute_stiffness ? valN : Val(0))
    mass = create_element_matrix(Float64, flags.compute_mass ? valN : Val(0))
    return ElementArrays(dofs, internal_force, lumped_mass, stiffness, mass)
end

function reset_element_arrays!(
    arrays::ElementArrays,
    element_block_connectivity::Matrix{<:Integer},
    block_element_index::Integer,
    flags::EvaluationFlags,
)
    node_indices = view(element_block_connectivity, :, block_element_index)
    fill_dofs!(arrays.dofs, element_block_connectivity, block_element_index)
    fill!(arrays.internal_force, 0.0)
    if flags.compute_lumped_mass == true
        fill!(arrays.lumped_mass, 0.0)
    end
    if flags.compute_stiffness == true
        fill!(arrays.stiffness, 0.0)
    end
    if flags.compute_mass == true
        fill!(arrays.mass, 0.0)
    end
    return node_indices
end

function accumulate_element_arrays!(
    arrays::ElementArrays,
    Np::SVector{N,T},
    dNdX::SMatrix{3,N,T},
    P::SMatrix{3,3,T,9},
    AA::SArray{Tuple{3,3,3,3},T},
    density::T,
    dvol::T,
    flags::EvaluationFlags,
    mass_at_points::Bool,
) where {T,N}
    add_internal_force!(arrays.internal_force, dNdX, P, dvol)
    if flags.compute_lumped_mass == true && mass_at_points
        add_lumped_mass!(arrays.lumped_mass, Np, density, dvol)
    end
    if flags.compute_stiffness == true
        add_stiffness!(arrays.stiffness, dNdX, AA, dvol)
    end
    if flags.compute_mass == true && mass_at_points
        add_mass!(arrays.mass, Np, density, dvol)
    end
    return nothing
end

# The mass of an element integrated with a rule of its own (the composite
# tetrahedron), on the reference configuration.
function add_element_mass!(
    arrays::ElementArrays,
    element_reference_position::SMatrix{3,NN,T},
    mass_N::SMatrix{NN,NM,T},
    mass_dN::SArray{Tuple{3,NN,NM},T},
    mass_weights::AbstractVector{T},
    density::T,
    flags::EvaluationFlags,
) where {NN,NM,T}
    for point in 1:NM
        Np = mass_N[:, point]
        dXdξ = mass_dN[:, :, point] * element_reference_position'
        dvol = det(dXdξ) * mass_weights[point]
        if flags.compute_lumped_mass == true
            add_lumped_mass!(arrays.lumped_mass, Np, density, dvol)
        end
        if flags.compute_mass == true
            add_mass!(arrays.mass, Np, density, dvol)
        end
    end
    return nothing
end

function evaluate(model::SolidMechanics, integrator::TimeIntegrator, solver::Solver)
    flags = compute_flags(model, integrator, solver)
    if (flags.compute_stiffness == true || flags.compute_mass == true) && model.matrix_pattern === nothing
        model.matrix_pattern = build_matrix_pattern(model)
    end
    # Quantities computed once are marked done here; a failed evaluation
    # switches them back on.
    if flags.compute_lumped_mass == true
        model.compute_lumped_mass = false
    end
    if flags.compute_stiffness == true && model.kinematics == Infinitesimal
        model.compute_stiffness = false
    end
    if flags.compute_mass == true
        model.compute_mass = false
    end
    strain_energy = 0.0
    Δt = material_time_step(integrator)
    materials = model.materials
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    body_force_vector = sized_vector!(model.body_force, num_dofs)
    buffers = model.value_buffers
    pattern = model.matrix_pattern
    for (block_index, block) in enumerate(model.blocks)
        material = materials[block_index]
        density = material.ρ
        element_type = block.element_type
        N, dN, ip_weights = block.N, block.dN, block.weights
        mass_N, mass_dN, mass_weights = block.mass_N, block.mass_dN, block.mass_weights
        element_block_connectivity = block.connectivity
        force_values = buffers.force[block_index]
        lumped_mass_values = flags.compute_lumped_mass ? buffers.lumped_mass[block_index] : Float64[]
        stiffness_values = flags.compute_stiffness ? pattern.stiffness[block_index] : Float64[]
        mass_values = flags.compute_mass ? pattern.mass[block_index] : Float64[]
        if flags.compute_mass == true && length(mass_values) != length(pattern.slots[block_index])
            resize!(mass_values, length(pattern.slots[block_index]))
        end
        # Function barrier: the shape function tables are static arrays whose
        # type depends on the element type, so the loop must be compiled for
        # the concrete types to avoid dynamic dispatch on every operation.
        ok, block_energy = evaluate_block!(
            model,
            force_values,
            lumped_mass_values,
            stiffness_values,
            mass_values,
            block_index,
            material,
            density,
            element_type,
            N,
            dN,
            ip_weights,
            mass_N,
            mass_dN,
            mass_weights,
            element_block_connectivity,
            flags,
            Δt,
        )
        ok || return nothing
        strain_energy += block_energy
    end
    model.strain_energy = strain_energy
    # The global vectors are written in place: a vector held elsewhere is a
    # copy (set_internal_force!, the saved states), so no one else sees the
    # write, and the allocation and first touch of a fresh vector were the
    # cost of the reduction on 12 threads.
    model.internal_force = reduce_element_vector!(sized_vector!(model.internal_force, num_dofs), buffers, buffers.force)
    if flags.compute_lumped_mass == true
        model.lumped_mass = reduce_element_vector!(
            sized_vector!(model.lumped_mass, num_dofs), buffers, buffers.lumped_mass
        )
    end
    # The matrices too: the quasi-static Hessian aliases the stiffness, which
    # the readers of the previous values (the tangent predictor) copy.
    if flags.compute_stiffness == true
        model.stiffness = reduce_element_matrix!(model.stiffness, pattern, pattern.stiffness, num_dofs)
    end
    if flags.compute_mass == true
        model.mass = reduce_element_matrix!(model.mass, pattern, pattern.mass, num_dofs)
    end
    model.body_force = body_force_vector
    return nothing
end

# Gather the nodal values of an element into a static 3 x NN matrix; the node
# count comes from the type of the shape function table.
@inline function gather_nodal(
    A::Matrix{Float64}, node_indices::AbstractVector{<:Integer}, ::SMatrix{NN,NP,T}
) where {NN,NP,T}
    return SMatrix{3,NN,T}(view(A, :, node_indices))
end
@inline function gather_nodal(A::AbstractMatrix, ::SMatrix{NN,NP,T}) where {NN,NP,T}
    return SMatrix{3,NN,T}(A)
end

# Lattice orientation of an element, the identity for materials without one.
@inline function element_lattice_orientation(
    model::SolidMechanics, material::Material, block_index::Int64, block_element_index::Int64
)
    has_orientation(material) || return SMatrix{3,3,Float64,9}(I)
    return model.lattice_orientation[block_index][block_element_index]
end

# Time increment over which rate-dependent materials integrate: the current
# step, or zero before the first step (the initial equilibrium solve and the
# initial acceleration), while the previous time of the integrator is unset.
material_time_step(integrator::TimeIntegrator) = isfinite(integrator.prev_time) ? integrator.time_step : 0.0

# The constitutive call is guarded against domain errors from the material
# models, but a try block inside the threaded loop body defeats type inference
# for every variable live across it, so the guard lives in this small function
# with a concrete return type. On failure the returned values are zeros and the
# first entry is false.
function guarded_constitutive!(
    model::SolidMechanics,
    material::Material,
    F::SMatrix{3,3,Float64,9},
    block_index::Int64,
    block_element_index::Int64,
    point::Int64,
    need_tangent::Bool,
    Δt::Float64,
)
    try
        Q = element_lattice_orientation(model, material, block_index, block_element_index)
        state = model.state_old[block_index][block_element_index][point]
        W, P, AA, state_new = constitutive_update(material, F, Q, state, Δt; need_tangent=need_tangent)
        if number_states(material) > 0
            model.state[block_index][block_element_index][point] = state_new
        end
        return true, W, P, AA
    catch e
        e isa _MATH_ERRORS || rethrow()
        norma_logf(4, :solve, "evaluate: caught %s in constitutive model", typeof(e))
        return false, 0.0, zero(SMatrix{3,3,Float64,9}), zero(SArray{Tuple{3,3,3,3},Float64,4,81})
    end
end

# Per-quadrature-point weights of the assembly of a subdomain of an overlap
# matching coupling (overlap_energy.jl): weights[block][element][point], one
# half inside the partner and one outside, multiplying the volume element of
# every assembled quantity. Keyed by model, as the Arlequin weight cache is;
# a model absent from the cache is assembled with unit weights.
const ASSEMBLY_WEIGHT_CACHE = IdDict{SolidMechanics,Vector{Vector{Vector{Float64}}}}()

# The same weights at the facet integration points of the traction and
# pressure conditions of such a subdomain: weights[side][point], keyed by
# condition. A condition absent from the cache is integrated with unit weights.
const TRACTION_WEIGHT_CACHE = IdDict{BoundaryCondition,Vector{Vector{Float64}}}()

function assembly_block_weights(model::SolidMechanics, block_index::Int64)
    weights = get(ASSEMBLY_WEIGHT_CACHE, model, nothing)
    return weights === nothing ? nothing : weights[block_index]
end

function evaluate_element!(
    model::SolidMechanics,
    arrays::ElementArrays,
    force_values::Vector{Float64},
    lumped_mass_values::Vector{Float64},
    stiffness_values::Vector{Float64},
    mass_values::Vector{Float64},
    block_index::Int64,
    block_element_index::Int64,
    material::Material,
    density::Float64,
    element_type::ElementType,
    N::SMatrix,
    dN::SArray,
    ip_weights::AbstractVector,
    mass_N::SMatrix,
    mass_dN::SArray,
    mass_weights::AbstractVector,
    element_block_connectivity::Matrix{<:Integer},
    flags::EvaluationFlags,
    num_points::Int64,
    Δt::Float64,
    block_weights::Union{Nothing,Vector{Vector{Float64}}},
)
    node_indices = reset_element_arrays!(arrays, element_block_connectivity, block_element_index, flags)
    # Anisotropic smoothing: the ideal element comes from the metric field and
    # the energy is evaluated on the deformation gradient in the metric space,
    # F_U = F_M F F_M⁻¹, with stress and tangent pulled back to the reference.
    anisotropic = flags.mesh_smoothing == true && model.metric_field !== nothing
    F_M = SMatrix{3,3,Float64,9}(I)
    F_M_inv = F_M
    if flags.mesh_smoothing == true
        smooth, F_M, F_M_inv = smoothing_reference(
            model, element_type, gather_nodal(model.reference, node_indices, N); node_indices
        )
        element_reference_position = gather_nodal(smooth, N)
    else
        element_reference_position = gather_nodal(model.reference, node_indices, N)
    end
    # The current position always uses the mesh reference, not the smoothed one.
    element_current_position =
        gather_nodal(model.reference, node_indices, N) + gather_nodal(model.displacement, node_indices, N)
    mass_at_points = element_type != COMPOSITE_TETRA10
    if !mass_at_points && (flags.compute_mass == true || flags.compute_lumped_mass == true)
        add_element_mass!(arrays, element_reference_position, mass_N, mass_dN, mass_weights, density, flags)
    end
    element_energy = 0.0
    for point in 1:num_points
        Np = N[:, point]
        dNdξ = dN[:, :, point]
        dXdξ = dNdξ * element_reference_position'
        dNdX = dXdξ \ dNdξ
        F = element_current_position * dNdX'
        J = det(F)
        if J ≤ 0.0 || isfinite(J) == false
            model.failed = true
            model.compute_mass = model.compute_lumped_mass = true
            model.compute_stiffness = true
            norma_log(0, :error, "Non-positive Jacobian detected! This may indicate element distortion.")
            norma_logf(4, :warning, "det(F) = %.3e", J)
            log_matrix(4, :info, "Reference Configuration", Matrix(element_reference_position))
            log_matrix(4, :info, "Current Configuration", Matrix(element_current_position))
            return 0.0
        end
        if anisotropic
            F_U = F_M * F * F_M_inv
            ok, W, P_U, AA_U = guarded_constitutive!(
                model, material, F_U, block_index, block_element_index, point, flags.compute_stiffness, Δt
            )
            P = F_M' * P_U * F_M_inv'
            AA = flags.compute_stiffness ? pull_back_tangent(AA_U, F_M, F_M_inv) : AA_U
        else
            ok, W, P, AA = guarded_constitutive!(
                model, material, F, block_index, block_element_index, point, flags.compute_stiffness, Δt
            )
        end
        if ok == false
            model.failed = true
            return 0.0  # skip remaining integration points for this element
        end
        ip_weight = ip_weights[point]
        det_dXdξ = det(dXdξ)
        dvol = det_dXdξ * ip_weight
        block_weights === nothing || (dvol *= block_weights[block_element_index][point])
        element_energy += W * dvol
        accumulate_element_arrays!(arrays, Np, dNdX, P, AA, density, dvol, flags, mass_at_points)
        voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
        model.stress[block_index][block_element_index][point] = voigt_cauchy
    end
    model.stored_energy[block_index][block_element_index] = element_energy
    store_element_values!(
        arrays, block_element_index, force_values, lumped_mass_values, stiffness_values, mass_values, flags
    )
    return element_energy
end

# The elements of a chunk, by the task that owns `arrays`; the sum of their
# strain energies.
function evaluate_chunk!(
    model::SolidMechanics,
    arrays::ElementArrays,
    force_values::Vector{Float64},
    lumped_mass_values::Vector{Float64},
    stiffness_values::Vector{Float64},
    mass_values::Vector{Float64},
    block_index::Int64,
    material::Material,
    density::Float64,
    element_type::ElementType,
    N::SMatrix,
    dN::SArray,
    ip_weights::AbstractVector,
    mass_N::SMatrix,
    mass_dN::SArray,
    mass_weights::AbstractVector,
    element_block_connectivity::Matrix{<:Integer},
    flags::EvaluationFlags,
    num_points::Int64,
    Δt::Float64,
    block_weights::Union{Nothing,Vector{Vector{Float64}}},
    elements::UnitRange{Int64},
)
    energy = 0.0
    for block_element_index in elements
        model.failed && break
        energy += evaluate_element!(
            model,
            arrays,
            force_values,
            lumped_mass_values,
            stiffness_values,
            mass_values,
            block_index,
            block_element_index,
            material,
            density,
            element_type,
            N,
            dN,
            ip_weights,
            mass_N,
            mass_dN,
            mass_weights,
            element_block_connectivity,
            flags,
            num_points,
            Δt,
            block_weights,
        )
    end
    return energy
end

# One block: whether every element evaluated, and the strain energy of the
# block. The elements are evaluated in chunks of fixed size by one task per
# thread, each task with element arrays of its own (chunk_results): no two
# tasks write one cache line, and the sum of the chunk energies, taken in
# chunk order, is the same on every thread count.
function evaluate_block!(
    model::SolidMechanics,
    force_values::Vector{Float64},
    lumped_mass_values::Vector{Float64},
    stiffness_values::Vector{Float64},
    mass_values::Vector{Float64},
    block_index::Int64,
    material::Material,
    density::Float64,
    element_type::ElementType,
    N::SMatrix{NN,NP,T},
    dN::SArray{Tuple{3,NN,NP},T},
    ip_weights::AbstractVector{T},
    mass_N::SMatrix,
    mass_dN::SArray,
    mass_weights::AbstractVector,
    element_block_connectivity::Matrix{<:Integer},
    flags::EvaluationFlags,
    Δt::Float64,
) where {NN,NP,T}
    num_block_elements = size(element_block_connectivity, 2)
    num_points = size(N, 2)
    block_weights = assembly_block_weights(model, block_index)
    # No static parameter (NN, NP, T) is used inside the closures below: a
    # captured static parameter is a runtime value there, which would make
    # every static array constructor a dynamic call. The element arrays are
    # created from this value instead.
    valN = Val(NN)
    init() = create_element_arrays(valN, flags)
    energies = chunk_results(Float64, init, num_block_elements) do elements, arrays
        return evaluate_chunk!(
            model,
            arrays,
            force_values,
            lumped_mass_values,
            stiffness_values,
            mass_values,
            block_index,
            material,
            density,
            element_type,
            N,
            dN,
            ip_weights,
            mass_N,
            mass_dN,
            mass_weights,
            element_block_connectivity,
            flags,
            num_points,
            Δt,
            block_weights,
            elements,
        )
    end
    return !model.failed, sum(energies)
end

# The connectivity of a block as Exodus stores it, element by element, shaped
# num_elements × num_nodes: the nodes of element e are the entries
# (e - 1) * num_nodes + 1 to e * num_nodes in linear order, not row e.
function get_block_connectivity(mesh::ExodusDatabase, block_id::Integer)
    _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(mesh, Int32(block_id))
    conn = Exodus.read_block_connectivity(mesh, Int32(block_id), num_elements * num_nodes)
    return reshape(conn, (num_elements, num_nodes))
end
