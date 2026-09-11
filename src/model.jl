# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.


using Base.Threads: @threads, threadid, nthreads, maxthreadid

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
# `trust_schwarz` into initialize(::TimeIntegrator, ::Solver, ::Model) /
# initialize(::RomNewmark, ...) / initialize(::RomCentralDifference, ...).
is_restarted(model::SolidMechanics) = model.restarted
is_restarted(model::RomModel) = model.fom_model.restarted

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
    if restart_info !== nothing && any(material isa J2Plasticity for material in materials)
        norma_abort(
            "Restart is not currently supported for the `j2 plasticity` material model. " *
            "The restart snapshot only stores nodal displacement and velocity fields; " *
            "J2 plasticity's internal state variables (e.g. plastic strain, back stress) " *
            "are not written to or read from the restart file, so resuming would silently " *
            "discard the accumulated plastic history. Remove the `restart:` block, or switch " *
            "to a material model without internal state variables, until restart support for " *
            "internal variables is implemented.",
        )
    end
    time = 0.0
    failed = false
    internal_force = zeros(3 * num_nodes)
    boundary_force = zeros(3 * num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(3 * num_nodes)
    stress = Vector{Vector{Vector{Vector{Float64}}}}()
    state_old = Vector{Vector{Vector{Vector{Float64}}}}()
    state = Vector{Vector{Vector{Vector{Float64}}}}()
    prev_state_old = Vector{Vector{Vector{Vector{Float64}}}}()  # empty until first save_curr_state
    stop_state_old = Vector{Vector{Vector{Vector{Float64}}}}()  # empty until first save_stop_state
    stored_energy = Vector{Vector{Float64}}()
    num_int_pts_overrides = get(model_params, "num integration points", Dict{String,Any}())
    num_int_pts = Vector{Int}(undef, num_blocks)
    block_data = Vector{ElementBlockData}(undef, num_blocks)
    for (block_index, block) in enumerate(blocks)
        block_id = block.id
        element_type_string, num_block_elements, _, _, _, _ = Exodus.read_block_parameters(input_mesh, block_id)
        element_type = element_type_from_string(element_type_string)
        block_name = element_block_names[block_index]
        num_points = haskey(num_int_pts_overrides, block_name) ?
            Int(num_int_pts_overrides[block_name]) : default_num_int_pts(element_type)
        num_int_pts[block_index] = num_points
        # Stored as one column per element (the flat Exodus order is element by
        # element), so an element's nodes are a contiguous column.
        raw_connectivity = get_block_connectivity(input_mesh, block_id)
        num_nodes_per_element = size(raw_connectivity, 2)
        connectivity = reshape(Int64.(vec(raw_connectivity)), (num_nodes_per_element, num_block_elements))
        N, dN, weights = isoparametric(element_type, num_points)
        block_data[block_index] = ElementBlockData(
            Int64(block_id), element_type, num_points, num_block_elements, num_nodes_per_element, connectivity, N, dN, weights
        )
        material = materials[block_index]
        num_states = number_states(material)
        if (num_states > 0)
            point_init_state = initial_state(material)
        else
            point_init_state = Vector{Float64}()
        end
        block_stress = Vector{Vector{Vector{Float64}}}()
        block_state = Vector{Vector{Vector{Float64}}}()
        block_stored_energy = Vector{Float64}()
        for _ in 1:num_block_elements
            element_stress = Vector{Vector{Float64}}()
            element_state = Vector{Vector{Float64}}()
            for _ in 1:num_points
                push!(element_stress, zeros(6))
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
        # only when a step is accepted.
        push!(state_old, block_state)
        push!(state, deepcopy(block_state))
        push!(stored_energy, block_stored_energy)
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
            norma_abort(
                "nodal recovery: 'method' must be 'lumped', 'consistent', or 'both' " *
                "(got '$method').",
            )
        end
        rec_stress = Bool(get(recovery_params, "stress", true))
        rec_vm = Bool(get(recovery_params, "von mises stress", false))
        rec_iv = Bool(get(recovery_params, "internal variables", false))
        rec_F = Bool(get(recovery_params, "deformation gradient", false))
    end
    # The impedance-overlap Schwarz BC evaluates the partner's traction at its
    # Schwarz boundary by interpolating nodal-recovered stress (see
    # apply_bc_detail in schwarz.jl), so a model that participates in such a
    # coupling must carry a recovered stress field even if the input file does
    # not request nodal recovery for output.
    bc_params = get(params, "boundary conditions", nothing)
    if bc_params !== nothing && haskey(bc_params, "Schwarz impedance overlap")
        if recovery_kind === :none
            recovery_kind = :consistent
        end
        rec_stress = true
    end
    recovery_data = build_recovery_data(recovery_kind, input_mesh, reference, num_int_pts)
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
        create_element_value_buffers(block_data),
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
            "smooth reference = \"$smooth_reference\" requires a \"size field\" expression under the model parameters",
        )
    end
    size_num = eval(Meta.parse(string(expression)))
    return eval(build_function(size_num, [t, x, y, z]; expression=Val(false)))
end

function create_smooth_reference(
    smooth_reference::String,
    element_type::ElementType,
    element_ref_pos::Matrix{Float64},
    size_field::Union{Function,Nothing}=nothing,
    time::Float64=0.0,
)::Matrix{Float64}
    if element_type == TETRA4
        u = element_ref_pos[:, 2] - element_ref_pos[:, 1]
        v = element_ref_pos[:, 3] - element_ref_pos[:, 1]
        w = element_ref_pos[:, 4] - element_ref_pos[:, 1]

        if smooth_reference == "equal volume"
            h = equal_volume_tet_h(u, v, w)
        elseif smooth_reference == "average edge length"
            h = avg_edge_length_tet_h(u, v, w)
        elseif smooth_reference == "max"
            h = max(avg_edge_length_tet_h(u, v, w), equal_volume_tet_h(u, v, w))
        elseif smooth_reference == "size field"
            # Target edge length from the user-defined size field at the element
            # reference centroid, combined with the volume criterion (max) to
            # anchor the reference size and avoid sliver pathologies.
            h = max(size_field_tet_h(size_field, element_ref_pos, time), equal_volume_tet_h(u, v, w))
        elseif smooth_reference == "size field unrestricted"
            # Target edge length from the user-defined size field at the element
            # reference centroid
            h = size_field_tet_h(size_field, element_ref_pos, time)
        else
            norma_abort("Unknown type of mesh smoothing reference : $smooth_reference")
        end

        c = h * 0.5 / sqrt(2.0)
        A = [
            1 -1 -1 1
            1 -1 1 -1
            1 1 -1 -1
        ]
        return c * A
    else
        norma_abort("Unknown element type")
    end
end

function equal_volume_tet_h(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    h = cbrt(sqrt(2.0) * dot(u, cross(v, w)))
    return h
end

function avg_edge_length_tet_h(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    h = (norm(u) + norm(v) + norm(w) + norm(u - v) + norm(u - w) + norm(v - w)) / 6.0
    return h
end

# Target edge length from the user-defined size field, evaluated at the element
# reference centroid and the current time.  The field must be strictly positive
# and finite to yield a valid (non-degenerate) reference element.
function size_field_tet_h(size_field::Union{Function,Nothing}, element_ref_pos::Matrix{Float64}, time::Float64)
    if size_field === nothing
        norma_abort("smooth reference = \"size field\" or \"size field unrestricted\" selected but no size field was compiled")
    end
    centroid = vec(sum(element_ref_pos; dims=2) / size(element_ref_pos, 2))
    h = size_field((time, centroid[1], centroid[2], centroid[3]))
    if !isfinite(h) || h ≤ 0.0
        norma_abort(
            "Size field must be strictly positive and finite; got $h at centroid " *
            "($(centroid[1]), $(centroid[2]), $(centroid[3])) and time $time",
        )
    end
    return h
end

function characteristic_element_length_centroid(nodal_coordinates::Matrix{Float64})::Float64
    centroid = sum(nodal_coordinates; dims=2) / size(nodal_coordinates, 2)
    total = 0.0
    @inbounds for i in 1:size(nodal_coordinates, 2)
        δ = nodal_coordinates[:, i] - centroid
        total += norm(δ)
    end
    return 2 * total / size(nodal_coordinates, 2)  # Approximate diameter
end

# Same measure on a static nodal coordinate matrix, without allocation.
function characteristic_element_length_centroid(X::SMatrix{3,N,T})::T where {N,T}
    centroid = zero(SVector{3,T})
    @inbounds for a in 1:N
        centroid += X[:, a]
    end
    centroid = centroid / N
    total = zero(T)
    @inbounds for a in 1:N
        total += norm(X[:, a] - centroid)
    end
    return 2 * total / N
end

# Minimum characteristic length over the elements of a block at the current
# configuration. Threaded over elements with a per-thread minimum; the shape
# function table only supplies the static node count through the barrier.
function minimum_characteristic_length(model::SolidMechanics, block::ElementBlockData, ::SMatrix{NN,NP,T}) where {NN,NP,T}
    connectivity = block.connectivity
    minima = fill(T(Inf), maxthreadid())
    @threads for element in 1:block.num_elements
        node_indices = view(connectivity, :, element)
        X = SMatrix{3,NN,T}(view(model.reference, :, node_indices)) + SMatrix{3,NN,T}(view(model.displacement, :, node_indices))
        h = characteristic_element_length_centroid(X)
        t = threadid()
        @inbounds minima[t] = min(minima[t], h)
    end
    return minimum(minima)
end

# The per-point material state is copied at several points of a step. For
# materials without internal variables the nested vectors are all empty and
# nothing ever writes into them, so the copy is skipped and the arrays are
# shared. The copies were a fifth of an explicit step for elastic materials.
has_material_state(model::SolidMechanics) = any(m -> number_states(m) > 0, model.materials)

function copy_state(model::SolidMechanics, state::Vector{Vector{Vector{Vector{Float64}}}})
    return has_material_state(model) ? deepcopy(state) : state
end

function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    materials = model.materials
    stable_time_step = Inf
    for (block_index, block) in enumerate(model.blocks)
        material = materials[block_index]
        ρ = material.ρ
        M = get_p_wave_modulus(material)
        wave_speed = sqrt(M / ρ)
        minimum_block_characteristic_length = minimum_characteristic_length(model, block, block.N)
        block_stable_time_step = integrator.CFL * minimum_block_characteristic_length / wave_speed
        stable_time_step = min(stable_time_step, block_stable_time_step)
    end
    if stable_time_step < integrator.time_step
        norma_logf(
            0,
            :warning,
            "Δt = %.3e exceeds stable Δt = %.3e — using stable step.",
            integrator.time_step,
            stable_time_step,
        )
    end
    integrator.time_step = min(stable_time_step, integrator.time_step)
    return nothing
end

function voigt_cauchy_from_stress(_::Solid, P::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9}, J::Float64)
    σ = F * P' ./ J
    return SVector{6,Float64}(σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2])
end

function voigt_cauchy_from_stress(_::Linear_Elastic, σ::SMatrix{3,3,Float64,9}, _::SMatrix{3,3,Float64,9}, _::Float64)
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











using Base.Threads

function create_threadlocal_element_matrices(::Type{T}, ::Val{N}) where {T,N}
    return [create_element_matrix(T, Val(N)) for _ in 1:Threads.maxthreadid()]
end

function create_threadlocal_element_vectors(::Type{T}, ::Val{N}) where {T,N}
    return [create_element_vector(T, Val(N)) for _ in 1:Threads.maxthreadid()]
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

function create_element_value_buffers(blocks::Vector{ElementBlockData})
    force = [zeros(3 * block.num_nodes_per_element * block.num_elements) for block in blocks]
    lumped_mass = [zeros(3 * block.num_nodes_per_element * block.num_elements) for block in blocks]
    return ElementValueBuffers(force, lumped_mass)
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
    mass = [zeros(length(block_slots)) for block_slots in slots]
    return MatrixPattern(colptr, rowval, slots, stiffness, mass)
end

# Copy the element arrays of one element into its value segments. Segments of
# different elements are disjoint, so this is safe from any thread.
function store_element_values!(
    element_arrays_tl::SMElementThreadLocalArrays,
    block_element_index::Int64,
    force_values::Vector{Float64},
    lumped_mass_values::Vector{Float64},
    stiffness_values::Vector{Float64},
    mass_values::Vector{Float64},
    flags::EvaluationFlags,
)
    t = threadid()
    nd = length(element_arrays_tl.dofs[t])
    offset = (block_element_index - 1) * nd
    fe = element_arrays_tl.internal_force[t]
    @inbounds for i in 1:nd
        force_values[offset + i] = fe[i]
    end
    if flags.compute_lumped_mass == true
        me = element_arrays_tl.lumped_mass[t]
        @inbounds for i in 1:nd
            lumped_mass_values[offset + i] = me[i]
        end
    end
    n2 = nd * nd
    offset2 = (block_element_index - 1) * n2
    if flags.compute_stiffness == true
        Ke = element_arrays_tl.stiffness[t]
        @inbounds for k in 1:n2
            stiffness_values[offset2 + k] = Ke[k]
        end
    end
    if flags.compute_mass == true
        Me = element_arrays_tl.mass[t]
        @inbounds for k in 1:n2
            mass_values[offset2 + k] = Me[k]
        end
    end
    return nothing
end

function reduce_element_vector(blocks::Vector{ElementBlockData}, values::Vector{Vector{Float64}}, num_dofs::Int64)
    global_vector = zeros(num_dofs)
    for (block_index, block) in enumerate(blocks)
        block_values = values[block_index]
        connectivity = block.connectivity
        num_nodes_per_element = block.num_nodes_per_element
        k = 0
        @inbounds for element in 1:block.num_elements, a in 1:num_nodes_per_element
            base = 3 * (connectivity[a, element] - 1)
            global_vector[base + 1] += block_values[k + 1]
            global_vector[base + 2] += block_values[k + 2]
            global_vector[base + 3] += block_values[k + 3]
            k += 3
        end
    end
    return global_vector
end

function reduce_element_matrix(pattern::MatrixPattern, values::Vector{Vector{Float64}}, num_dofs::Int64)
    nzval = zeros(length(pattern.rowval))
    for block_index in eachindex(values)
        block_values = values[block_index]
        block_slots = pattern.slots[block_index]
        @inbounds for k in eachindex(block_values)
            nzval[block_slots[k]] += block_values[k]
        end
    end
    return SparseMatrixCSC(num_dofs, num_dofs, pattern.colptr, pattern.rowval, nzval)
end

# K + M / β / Δt² for the Newmark Hessian. Matrices built from the same pattern
# share their structure arrays, so the sum is a pass over the values with the
# same elementwise arithmetic as the generic expression.
function newmark_hessian(K::SparseMatrixCSC{Float64,Int64}, M::SparseMatrixCSC{Float64,Int64}, β::Float64, Δt::Float64)
    if K.colptr === M.colptr && K.rowval === M.rowval
        return SparseMatrixCSC(K.m, K.n, K.colptr, K.rowval, K.nzval .+ M.nzval ./ β ./ Δt ./ Δt)
    end
    return K + M / β / Δt / Δt
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


function create_element_threadlocal_arrays(num_element_nodes::Int64, flags::EvaluationFlags)
    valN = Val(num_element_nodes)
    energy = zeros(maxthreadid())
    dofs = create_threadlocal_element_vectors(Int64, valN)
    internal_force = create_threadlocal_element_vectors(Float64, valN)
    lumped_mass = create_threadlocal_element_vectors(Float64, flags.compute_lumped_mass ? valN : Val(0))
    stiffness = create_threadlocal_element_matrices(Float64, flags.compute_stiffness ? valN : Val(0))
    mass = create_threadlocal_element_matrices(Float64, flags.compute_mass ? valN : Val(0))
    return SMElementThreadLocalArrays(energy, dofs, internal_force, lumped_mass, stiffness, mass)
end

function reset_element_threadlocal_arrays!(
    element_arrays_tl::SMElementThreadLocalArrays,
    element_block_connectivity::Matrix{<:Integer},
    block_element_index::Integer,
    flags::EvaluationFlags,
)
    t = threadid()
    node_indices = view(element_block_connectivity, :, block_element_index)
    fill_dofs!(element_arrays_tl.dofs[t], element_block_connectivity, block_element_index)
    element_arrays_tl.energy[t] = 0.0
    fill!(element_arrays_tl.internal_force[t], 0.0)
    if flags.compute_lumped_mass == true
        fill!(element_arrays_tl.lumped_mass[t], 0.0)
    end
    if flags.compute_stiffness == true
        fill!(element_arrays_tl.stiffness[t], 0.0)
    end
    if flags.compute_mass == true
        fill!(element_arrays_tl.mass[t], 0.0)
    end
    return node_indices
end

function compute_element_threadlocal_arrays!(
    element_arrays_tl::SMElementThreadLocalArrays,
    Np::SVector{N,T},
    dNdX::SMatrix{3,N,T},
    W::T,
    P::SMatrix{3,3,T,9},
    AA::SArray{Tuple{3,3,3,3},T},
    density::T,
    dvol::T,
    flags::EvaluationFlags,
) where {T,N}
    t = threadid()
    element_arrays_tl.energy[t] += W * dvol
    add_internal_force!(element_arrays_tl.internal_force[t], dNdX, P, dvol)
    if flags.compute_lumped_mass == true
        add_lumped_mass!(element_arrays_tl.lumped_mass[t], Np, density, dvol)
    end
    if flags.compute_stiffness == true
        add_stiffness!(element_arrays_tl.stiffness[t], dNdX, AA, dvol)
    end
    if flags.compute_mass == true
        add_mass!(element_arrays_tl.mass[t], Np, density, dvol)
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
    energy_tl = zeros(maxthreadid())
    materials = model.materials
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    body_force_vector = zeros(num_dofs)
    buffers = model.value_buffers
    pattern = model.matrix_pattern
    for (block_index, block) in enumerate(model.blocks)
        material = materials[block_index]
        density = material.ρ
        element_type = block.element_type
        N, dN, ip_weights = block.N, block.dN, block.weights
        element_block_connectivity = block.connectivity
        element_arrays_tl = create_element_threadlocal_arrays(block.num_nodes_per_element, flags)
        force_values = buffers.force[block_index]
        lumped_mass_values = flags.compute_lumped_mass ? buffers.lumped_mass[block_index] : Float64[]
        stiffness_values = flags.compute_stiffness ? pattern.stiffness[block_index] : Float64[]
        mass_values = flags.compute_mass ? pattern.mass[block_index] : Float64[]
        # Function barrier: the shape function tables are static arrays whose
        # type depends on the element type, so the loop must be compiled for
        # the concrete types to avoid dynamic dispatch on every operation.
        ok = evaluate_block!(
            model, energy_tl, element_arrays_tl, force_values, lumped_mass_values, stiffness_values, mass_values,
            block_index, material, density, element_type, N, dN, ip_weights, element_block_connectivity, flags,
        )
        ok || return nothing
    end
    model.strain_energy = sum(energy_tl)
    model.internal_force = reduce_element_vector(model.blocks, buffers.force, num_dofs)
    if flags.compute_lumped_mass == true
        model.lumped_mass = reduce_element_vector(model.blocks, buffers.lumped_mass, num_dofs)
    end
    if flags.compute_stiffness == true
        model.stiffness = reduce_element_matrix(pattern, pattern.stiffness, num_dofs)
    end
    if flags.compute_mass == true
        model.mass = reduce_element_matrix(pattern, pattern.mass, num_dofs)
    end
    model.body_force = body_force_vector
    return nothing
end

# Gather the nodal values of an element into a static 3 x NN matrix; the node
# count comes from the type of the shape function table.
@inline function gather_nodal(A::Matrix{Float64}, node_indices::AbstractVector{<:Integer}, ::SMatrix{NN,NP,T}) where {NN,NP,T}
    return SMatrix{3,NN,T}(view(A, :, node_indices))
end
@inline function gather_nodal(A::AbstractMatrix, ::SMatrix{NN,NP,T}) where {NN,NP,T}
    return SMatrix{3,NN,T}(A)
end

# The constitutive call is guarded against domain errors from the material
# models, but a try block inside the threaded loop body defeats type inference
# for every variable live across it, so the guard lives in this small function
# with a concrete return type. On failure the returned values are zeros and the
# first entry is false.
function guarded_constitutive!(
    model::SolidMechanics, material::Material, F::SMatrix{3,3,Float64,9},
    block_index::Int64, block_element_index::Int64, point::Int64, need_tangent::Bool,
)
    try
        if material isa Elastic
            W, P, AA = constitutive(material, F; need_tangent=need_tangent)
            return true, W, P, AA
        else
            state = model.state_old[block_index][block_element_index][point]
            W, P, AA, state_new = constitutive(material, F, state; need_tangent=need_tangent)
            model.state[block_index][block_element_index][point] = state_new
            return true, W, P, AA
        end
    catch e
        e isa _MATH_ERRORS || rethrow()
        norma_logf(4, :solve, "evaluate: caught %s in constitutive model", typeof(e))
        return false, 0.0, zero(SMatrix{3,3,Float64,9}), zero(SArray{Tuple{3,3,3,3},Float64,4,81})
    end
end

function evaluate_element!(
    model::SolidMechanics,
    energy_tl::Vector{Float64},
    element_arrays_tl::SMElementThreadLocalArrays,
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
    element_block_connectivity::Matrix{<:Integer},
    flags::EvaluationFlags,
    num_points::Int64,
)
    node_indices = reset_element_threadlocal_arrays!(
        element_arrays_tl, element_block_connectivity, block_element_index, flags
    )
    if flags.mesh_smoothing == true
        smooth = create_smooth_reference(
            model.smooth_reference, element_type, model.reference[:, node_indices], model.size_field, model.time
        )
        element_reference_position = gather_nodal(smooth, N)
    else
        element_reference_position = gather_nodal(model.reference, node_indices, N)
    end
    # The current position always uses the mesh reference, not the smoothed one.
    element_current_position =
        gather_nodal(model.reference, node_indices, N) + gather_nodal(model.displacement, node_indices, N)
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
            return nothing
        end
        ok, W, P, AA = guarded_constitutive!(model, material, F, block_index, block_element_index, point, flags.compute_stiffness)
        if ok == false
            model.failed = true
            return nothing  # skip remaining integration points for this element
        end
        ip_weight = ip_weights[point]
        det_dXdξ = det(dXdξ)
        dvol = det_dXdξ * ip_weight
        compute_element_threadlocal_arrays!(element_arrays_tl, Np, dNdX, W, P, AA, density, dvol, flags)
        voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
        model.stress[block_index][block_element_index][point] = voigt_cauchy
    end
    t = threadid()
    element_energy = element_arrays_tl.energy[t]
    model.stored_energy[block_index][block_element_index] = element_energy
    energy_tl[t] += element_energy
    store_element_values!(
        element_arrays_tl, block_element_index, force_values, lumped_mass_values, stiffness_values, mass_values, flags
    )
    return nothing
end

function evaluate_block!(
    model::SolidMechanics,
    energy_tl::Vector{Float64},
    element_arrays_tl::SMElementThreadLocalArrays,
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
    element_block_connectivity::Matrix{<:Integer},
    flags::EvaluationFlags,
)::Bool where {NN,NP,T}
    num_block_elements = size(element_block_connectivity, 2)
    num_points = size(N, 2)
    # No static parameter (NN, NP, T) is used inside the threaded loop body: it
    # is a closure, and captured static parameters are runtime values there,
    # which would make every static array constructor a dynamic call.
    @threads for block_element_index in 1:num_block_elements
        model.failed && continue
        evaluate_element!(
            model, energy_tl, element_arrays_tl, force_values, lumped_mass_values, stiffness_values, mass_values,
            block_index, block_element_index, material, density, element_type, N, dN, ip_weights,
            element_block_connectivity, flags, num_points,
        )
    end
    return !model.failed
end

function get_block_connectivity(mesh::ExodusDatabase, block_id::Integer)
    _, num_elements, num_nodes, _, _, _ = Exodus.read_block_parameters(mesh, Int32(block_id))
    conn = Exodus.read_block_connectivity(mesh, Int32(block_id), num_elements * num_nodes)
    return reshape(conn, (num_elements, num_nodes))
end
