# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using DelimitedFiles
using Format

# Exodus.jl (0.14 and 0.15) reads every block, set, and variable name into a
# buffer of 32 bytes, where the Exodus library writes the name padded to its
# read length plus a terminator, 33 bytes, whatever the length of the name
# (cmhamel/Exodus.jl#232).  The extra zero byte lands on the Julia heap next
# to the buffer; on macOS it zeroed part of the type tag of a live object,
# which crashed later in dispatch (the intermittent segmentation faults of
# the adaptivity tests on CI).  Its read-mode constructor reads the name of
# every block, set, and variable, so every open of a file did this.  Norma
# opens files and reads names through the functions below, which give the
# library a buffer of its full read length.

function exodus_name_buffer(exoid::Cint)
    read_length = @ccall Exodus.libexodus.ex_inquire_int(
        exoid::Cint, Exodus.EX_INQ_MAX_READ_NAME_LENGTH::Exodus.ex_inquiry
    )::Int64
    return zeros(UInt8, max(read_length, Int64(Exodus.MAX_STR_LENGTH)) + 1)
end

function name_from_buffer(buffer::Vector{UInt8})
    terminator = something(findfirst(iszero, buffer), length(buffer) + 1)
    return String(buffer[1:(terminator - 1)])
end

function read_exodus_name(exo::ExodusDatabase, ::Type{S}, id::Integer) where {S<:Exodus.AbstractExodusSet}
    exoid = Exodus.get_file_id(exo)
    buffer = exodus_name_buffer(exoid)
    error_code = @ccall Exodus.libexodus.ex_get_name(
        exoid::Cint, Exodus.entity_type(S)::Exodus.ex_entity_type, Int64(id)::Int64, buffer::Ptr{UInt8}
    )::Cint
    error_code < 0 && norma_abort("Cannot read the name of $S $id of Exodus file $(exo.file_name)")
    return name_from_buffer(buffer)
end

function read_exodus_name(exo::ExodusDatabase, ::Type{V}, index::Integer) where {V<:Exodus.AbstractExodusVariable}
    exoid = Exodus.get_file_id(exo)
    buffer = exodus_name_buffer(exoid)
    error_code = @ccall Exodus.libexodus.ex_get_variable_name(
        exoid::Cint, Exodus.entity_type(V)::Exodus.ex_entity_type, Cint(index)::Cint, buffer::Ptr{UInt8}
    )::Cint
    error_code < 0 && norma_abort("Cannot read the name of $V $index of Exodus file $(exo.file_name)")
    return name_from_buffer(buffer)
end

read_exodus_names(exo::ExodusDatabase, ::Type{S}) where {S<:Exodus.AbstractExodusSet} =
    [read_exodus_name(exo, S, id) for id in Exodus.read_ids(exo, S)]

read_exodus_names(exo::ExodusDatabase, ::Type{V}) where {V<:Exodus.AbstractExodusVariable} =
    [read_exodus_name(exo, V, i) for i in 1:Exodus.read_number_of_variables(exo, V)]

# The counts of an Exodus file, read with room for the full title, which the
# library writes with up to 80 characters and a terminator.
function read_exodus_initialization(exoid::Cint, ::Type{B}) where {B}
    title = zeros(UInt8, Exodus.MAX_LINE_LENGTH + 1)
    dimensions, nodes, elements = Ref{B}(0), Ref{B}(0), Ref{B}(0)
    blocks, node_sets, side_sets = Ref{B}(0), Ref{B}(0), Ref{B}(0)
    error_code = @ccall Exodus.libexodus.ex_get_init(
        exoid::Cint, title::Ptr{UInt8}, dimensions::Ptr{B}, nodes::Ptr{B}, elements::Ptr{B}, blocks::Ptr{B},
        node_sets::Ptr{B}, side_sets::Ptr{B},
    )::Cint
    error_code < 0 && norma_abort("Cannot read the initialization of an Exodus file")
    return Exodus.Initialization{B}(dimensions[], nodes[], elements[], blocks[], node_sets[], side_sets[])
end

# Open an existing Exodus file for reading ("r") or reading and writing
# ("rw"), as the constructor ExodusDatabase(file, mode) of Exodus.jl does,
# with its name dictionaries filled through read_exodus_name.
function open_exodus_database(file_name::AbstractString, mode::AbstractString="r")
    mode in ("r", "rw") || norma_abort("An existing Exodus file opens in mode \"r\" or \"rw\", not \"$mode\"")
    isfile(file_name) || norma_abort("Exodus file $file_name does not exist")
    exoid = Exodus.open_exodus_file(String(file_name), String(mode))
    return open_exodus_database(
        exoid, String(file_name), String(mode), Exodus.map_int_mode(exoid), Exodus.id_int_mode(exoid),
        Exodus.bulk_int_mode(exoid), Exodus.float_mode(exoid),
    )
end

function open_exodus_database(
    exoid::Cint, file_name::String, mode::String, ::Type{M}, ::Type{I}, ::Type{B}, ::Type{F}
) where {M,I,B,F}
    init = read_exodus_initialization(exoid, B)
    D = Dict{String,I}
    exo = Exodus.ExodusDatabase{M,I,B,F}(exoid, mode, file_name, init, D(), D(), D(), D(), D(), D(), D(), D())
    counts = (Exodus.num_element_blocks(init), Exodus.num_node_sets(init), Exodus.num_side_sets(init))
    unnamed_prefixes = ("unnamed_block_", "unnamed_nset_", "unnamed_sset_")
    for (S, count, unnamed) in zip((Block, NodeSet, SideSet), counts, unnamed_prefixes)
        count == 0 && continue
        for id in Exodus.read_ids(exo, S)
            name = read_exodus_name(exo, S, id)
            Exodus.set_name_dict(exo, S)[isempty(name) ? unnamed * string(id) : name] = id
        end
    end
    for V in (ElementVariable, GlobalVariable, NodalVariable, NodeSetVariable, SideSetVariable)
        for (index, name) in enumerate(read_exodus_names(exo, V))
            Exodus.var_name_dict(exo, V)[name] = index
        end
    end
    return exo
end

# Create an Exodus database for writing, as the write-mode constructor of
# Exodus.jl does, but with a title.  That constructor passes an uninitialized
# buffer as the title, so the file gets up to 80 bytes of stale memory, and
# its reader gives ex_get_init a buffer of 80 bytes where the library writes
# the title and its terminator, up to 81.  A file written by the constructor
# can thus overwrite one byte of the Julia heap each time it is opened, which
# crashed the adaptivity tests on macOS (a corrupted type tag found during
# dispatch).  With a title shorter than 80 characters the reader is safe.
function create_exodus_database(file_name::AbstractString, init::Exodus.Initialization; title::AbstractString="Norma")
    length(codeunits(title)) < Exodus.MAX_LINE_LENGTH || norma_abort("Exodus title must be shorter than 80 bytes")
    isfile(file_name) && rm(file_name; force=true)
    exo = @ccall Exodus.libexodus.ex_create_int(
        String(file_name)::Cstring, Exodus.EX_WRITE::Cint, Exodus.cpu_word_size::Ref{Cint},
        sizeof(Float64)::Ref{Cint}, Exodus.EX_API_VERS_NODOT::Cint,
    )::Cint
    exo < 0 && norma_abort("Cannot create Exodus file $file_name")
    error_code = @ccall Exodus.libexodus.ex_put_init(
        exo::Cint, String(title)::Cstring, Exodus.num_dimensions(init)::Clonglong,
        Exodus.num_nodes(init)::Clonglong, Exodus.num_elements(init)::Clonglong,
        Exodus.num_element_blocks(init)::Clonglong, Exodus.num_node_sets(init)::Clonglong,
        Exodus.num_side_sets(init)::Clonglong,
    )::Cint
    error_code < 0 && norma_abort("Cannot write the initialization of Exodus file $file_name")
    D = Dict{String,Int32}
    return Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(
        exo, "w", String(file_name), init, D(), D(), D(), D(), D(), D(), D(), D()
    )
end

function _is_output_time(time::Float64, initial_time::Float64, interval::Float64; tol::Float64=1e-10)
    interval <= 0.0 && return false
    elapsed = time - initial_time
    n = round(elapsed / interval)
    return abs(elapsed - n * interval) < tol * interval
end

function collect_internal_variable_names(materials::Vector{Solid})
    all_names = String[]
    for material in materials
        for name in internal_variable_names(material)
            if !(name in all_names)
                push!(all_names, name)
            end
        end
    end
    return all_names
end

function initialize_writing(sim::SingleDomainSimulation)
    params = sim.params
    integrator = sim.integrator
    output_mesh = params["output_mesh"]
    # Per-file sequential frame counter.  write_stop_exodus increments this
    # and uses it as time_index, so the file always gets sequential frames
    # regardless of the global controller stop value (which may be offset when
    # a subsim is swapped in mid-run).
    params["exodus_frame"] = 0

    # setup nodal variables
    num_node_vars = 6
    node_var_names = ["refe_x", "refe_y", "refe_z", "disp_x", "disp_y", "disp_z"]
    if is_dynamic(integrator) == true
        num_node_vars += 6
        append!(node_var_names, ["velo_x", "velo_y", "velo_z", "acce_x", "acce_y", "acce_z"])
    end
    # For RomModel, nodal recovery is performed on the underlying fom_model after
    # reconstructing the full displacement field.  Use fom_model as the source of
    # recovery metadata so that the Exodus variable names are registered correctly
    # regardless of whether the top-level model is a ROM or a FOM.
    solid_model = sim.model isa RomModel ? sim.model.fom_model : sim.model
    if !(solid_model.recovery_data isa NoRecovery)
        F_comps = ("F_xx", "F_yx", "F_zx", "F_xy", "F_yy", "F_zy", "F_xz", "F_yz", "F_zz")
        if solid_model.recovery_data isa BothRecovery
            if size(solid_model.lumped_recovered_stress, 1) > 0
                num_node_vars += 12
                append!(node_var_names, ["sigma_xx_cons_n", "sigma_yy_cons_n", "sigma_zz_cons_n", "sigma_yz_cons_n", "sigma_xz_cons_n", "sigma_xy_cons_n"])
                append!(node_var_names, ["sigma_xx_lump_n", "sigma_yy_lump_n", "sigma_zz_lump_n", "sigma_yz_lump_n", "sigma_xz_lump_n", "sigma_xy_lump_n"])
            end
            if size(solid_model.lumped_recovered_von_mises, 1) > 0
                num_node_vars += 2
                append!(node_var_names, ["von_mises_stress_cons_n", "von_mises_stress_lump_n"])
            end
            if size(solid_model.lumped_recovered_F, 1) > 0
                num_node_vars += 18
                append!(node_var_names, [c * "_cons_n" for c in F_comps])
                append!(node_var_names, [c * "_lump_n" for c in F_comps])
            end
            if size(solid_model.lumped_recovered_internal_variables, 1) > 0
                iv_names = collect_internal_variable_names(solid_model.materials)
                num_node_vars += 2 * length(iv_names)
                append!(node_var_names, [name * "_cons_n" for name in iv_names])
                append!(node_var_names, [name * "_lump_n" for name in iv_names])
            end
        else
            if size(solid_model.recovered_stress, 1) > 0
                num_node_vars += 6
                append!(node_var_names, ["sigma_xx_n", "sigma_yy_n", "sigma_zz_n", "sigma_yz_n", "sigma_xz_n", "sigma_xy_n"])
            end
            if size(solid_model.recovered_von_mises, 1) > 0
                num_node_vars += 1
                push!(node_var_names, "von_mises_stress_n")
            end
            if size(solid_model.recovered_F, 1) > 0
                num_node_vars += 9
                append!(node_var_names, [c * "_n" for c in F_comps])
            end
            if size(solid_model.recovered_internal_variables, 1) > 0
                iv_names = collect_internal_variable_names(solid_model.materials)
                num_node_vars += length(iv_names)
                append!(node_var_names, [name * "_n" for name in iv_names])
            end
        end
    end
    # Energetic mesh smoothing driven by a size field: expose the target
    # edge-length field, sampled at each node, as a nodal variable so it can be
    # visualized on the (smoothed) mesh.
    if solid_model isa SolidMechanics && solid_model.mesh_smoothing == true && solid_model.size_field !== nothing
        num_node_vars += 1
        push!(node_var_names, "size")
    end
    # Anisotropic smoothing: the principal sizes and the rotation vector of the
    # metric field at each node, and the metric and target tensors built from
    # them (see nodal_metric_output and nodal_metric_tensors).
    if solid_model isa SolidMechanics && solid_model.mesh_smoothing == true && solid_model.metric_field !== nothing
        num_node_vars += 3
        append!(node_var_names, ["size_1", "size_2", "size_3"])
        if metric_output_has_rotation(solid_model.metric_field.source)
            num_node_vars += 3
            append!(node_var_names, ["rotation_1", "rotation_2", "rotation_3"])
        end
        num_node_vars += 12
        append!(node_var_names, METRIC_TENSOR_NAMES)
        append!(node_var_names, TARGET_TENSOR_NAMES)
    end
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_names(output_mesh, NodalVariable, node_var_names)

    # get maximum number of quadrature points (per-block override on the model)
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for (block_index, block) in enumerate(blocks)
        if sim.model isa RomModel
            block_id = block.id
            element_type_string = Exodus.read_block_parameters(output_mesh, block_id)[1]
            element_type = element_type_from_string(element_type_string)
            num_points = default_num_int_pts(element_type)
        else
            num_points = sim.model.num_int_pts[block_index]
        end
        max_num_int_points = max(max_num_int_points, num_points)
    end

    all_iv_names = sim.model isa RomModel ? String[] : collect_internal_variable_names(sim.model.materials)
    num_element_vars = (7 + length(all_iv_names)) * max_num_int_points + 1
    # Energetic mesh smoothing: the energy per unit ideal volume, the quantity
    # the adaptivity thresholds are stated on (docs/notes/ems-adaptivity).
    smoothing_density = sim.model isa SolidMechanics && sim.model.mesh_smoothing == true
    if smoothing_density
        num_element_vars += 1
    end
    Exodus.write_number_of_variables(output_mesh, ElementVariable, num_element_vars)

    el_var_names = String[]
    for point in 1:max_num_int_points
        ip_str = "_" * string(point)
        push!(el_var_names, "stress_xx" * ip_str)
        push!(el_var_names, "stress_yy" * ip_str)
        push!(el_var_names, "stress_zz" * ip_str)
        push!(el_var_names, "stress_yz" * ip_str)
        push!(el_var_names, "stress_xz" * ip_str)
        push!(el_var_names, "stress_xy" * ip_str)
        push!(el_var_names, "von_mises_stress" * ip_str)
        for iv_name in all_iv_names
            push!(el_var_names, iv_name * ip_str)
        end
    end
    push!(el_var_names, "stored_energy")
    if smoothing_density
        push!(el_var_names, "energy_density")
    end
    Exodus.write_names(output_mesh, ElementVariable, el_var_names)
    return nothing
end

function finalize_writing(sim::SingleDomainSimulation)
    input_mesh = sim.params["input_mesh"]
    Exodus.close(input_mesh)
    output_mesh = sim.params["output_mesh"]
    Exodus.close(output_mesh)
    return nothing
end

function writedlm_nodal_array(filename::String, nodal_array::Matrix{Float64})
    open(filename, "w") do io
        for col in 1:size(nodal_array, 2)
            # Write each column as a comma-separated line
            println(io, join(nodal_array[:, col], ","))
        end
    end
    return nothing
end

function get_umax(model::SolidMechanics)
  u_max = maximum(abs, model.displacement)
  return u_max
end

function get_umax(model::RomModel)
  u_max = maximum(abs, model.fom_model.displacement)
  return u_max
end

function write_stop(sim::SingleDomainSimulation; wall_time::Float64=0.0)
    params = sim.params
    stop = sim.controller.stop
    num_steps = sim.controller.num_stops - 1
    time = sim.controller.time
    name = sim.name
    model = sim.model
    is_explicit = sim.integrator isa ExplicitDynamicTimeIntegrator
    dt_default = params["time integrator"]["time step"]
    exodus_interval = Float64(get(params, "Exodus output interval", dt_default))
    csv_interval    = Float64(get(params, "CSV output interval", 0.0))
    initial_time = sim.controller.initial_time
    is_exodus_step = _is_output_time(time, initial_time, exodus_interval)
    is_csv_step    = _is_output_time(time, initial_time, csv_interval)
    is_output_step = is_exodus_step || is_csv_step

    # For explicit dynamics, only print at output steps (suppresses per-step noise).
    # For implicit/quasi-static, print every step (Newton iterations provide context).
    if !is_coupled(sim)
        if !is_explicit || is_output_step
            percent = 100 * stop / num_steps
            digits = max(0, Int64(ceil(log10(num_steps))) - 2)
            u_max = get_umax(model)
            if is_output_step && wall_time > 0.01
                norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.2e : |U|_max = %.2e : wall = %s",
                           stop, num_steps, percent, time, u_max, format_time(wall_time))
            elseif wall_time > 0.01
                norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.2e : wall = %s",
                           stop, num_steps, percent, time, format_time(wall_time))
            else
                norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.2e",
                           stop, num_steps, percent, time)
            end
        end
    end
    if is_exodus_step
        norma_log(0, :output, "Exodus II Database for $name [EXO]")
        write_stop_exodus(sim, sim.model)
    end
    if is_csv_step
        norma_log(0, :output, "Comma Separated Values for $name [CSV]")
        write_stop_csv(sim, sim.model)
        if haskey(params, "CSV write sidesets") == true
            write_sideset_stop_csv(sim, sim.model)
        end
    end
    return nothing
end

function write_stop(sim::MultiDomainSimulation; wall_time::Float64=0.0)
    stop = sim.controller.stop
    num_steps = sim.controller.num_stops - 1
    time = sim.controller.time
    percent = 100 * stop / num_steps
    digits = max(0, Int64(ceil(log10(num_steps))) - 2)
    if wall_time > 0.01
        norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.2e : wall = %s",
                   stop, num_steps, percent, time, format_time(wall_time))
    else
        norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.2e",
                   stop, num_steps, percent, time)
    end
    for subsim in sim.subsims
        write_stop(subsim)
    end
    if get(sim.params, "blended energy output", false) == true
        write_blended_energy_csv(sim)
    end
end

function write_stop_csv(sim::SingleDomainSimulation, model::SolidMechanics)
    stop = sim.controller.stop
    integrator = sim.integrator
    index_string = "-" * string(stop; pad=4)
    prefix = sim.name * "-"
    if stop == 0
        refe_filename = prefix * "refe" * ".csv"
        writedlm_nodal_array(refe_filename, model.reference)
    end
    free_dofs_filename = prefix * "free-dofs" * index_string * ".csv"
    writedlm(free_dofs_filename, model.free_dofs)
    curr_filename = prefix * "curr" * index_string * ".csv"
    writedlm_nodal_array(curr_filename, model.reference .+ model.displacement)
    disp_filename = prefix * "disp" * index_string * ".csv"
    writedlm_nodal_array(disp_filename, model.displacement)
    time_filename = prefix * "time" * index_string * ".csv"
    writedlm(time_filename, integrator.time, '\n')
    potential_filename = prefix * "potential" * index_string * ".csv"
    writedlm(potential_filename, integrator.stored_energy, '\n')
    if is_dynamic(integrator) == true
        velo_filename = prefix * "velo" * index_string * ".csv"
        writedlm_nodal_array(velo_filename, model.velocity)
        acce_filename = prefix * "acce" * index_string * ".csv"
        writedlm_nodal_array(acce_filename, model.acceleration)
        kinetic_filename = prefix * "kinetic" * index_string * ".csv"
        writedlm(kinetic_filename, integrator.kinetic_energy, '\n')
        total_filename = prefix * "total-energy" * index_string * ".csv"
        writedlm(total_filename, integrator.stored_energy + integrator.kinetic_energy, '\n')
    end
    return nothing
end

function write_sideset_stop_csv(sim::SingleDomainSimulation, model::SolidMechanics)
    stop = sim.controller.stop
    integrator = sim.integrator
    index_string = "-" * string(stop; pad=4)
    prefix = sim.name * "-"
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsDirichletBoundaryCondition
            node_set_name = bc.name
            offset = bc.offset
            if offset == 1
                offset_name = "x"
            end
            if offset == 2
                offset_name = "y"
            end
            if offset == 3
                offset_name = "z"
            end
            curr_filename = prefix * node_set_name * "-" * offset_name * "-curr" * index_string * ".csv"
            disp_filename = prefix * node_set_name * "-" * offset_name * "-disp" * index_string * ".csv"
            velo_filename = prefix * node_set_name * "-" * offset_name * "-velo" * index_string * ".csv"
            acce_filename = prefix * node_set_name * "-" * offset_name * "-acce" * index_string * ".csv"
            writedlm(curr_filename, model.reference[bc.offset, bc.node_set_node_indices] + model.displacement[bc.offset, bc.node_set_node_indices])
            writedlm(velo_filename, model.velocity[bc.offset, bc.node_set_node_indices])
            writedlm(acce_filename, model.acceleration[bc.offset, bc.node_set_node_indices])
            writedlm(disp_filename, model.displacement[bc.offset, bc.node_set_node_indices])
        elseif bc isa SolidMechanicsCouplingSchwarzBoundaryCondition
            # Every coupling Schwarz BC (overlap/non-overlap, DBC/impedance)
            # exposes name and side_set_node_indices; the force block below is
            # specific to the non-overlap DN variant.
            side_set_name = bc.name
            curr_filename = prefix * side_set_name * "-curr" * index_string * ".csv"
            disp_filename = prefix * side_set_name * "-disp" * index_string * ".csv"
            velo_filename = prefix * side_set_name * "-velo" * index_string * ".csv"
            acce_filename = prefix * side_set_name * "-acce" * index_string * ".csv"
            unique_indices = unique(bc.side_set_node_indices)
            writedlm_nodal_array(curr_filename, model.reference[:, unique_indices] .+ model.displacement[:, unique_indices])
            writedlm_nodal_array(velo_filename, model.velocity[:, unique_indices])
            writedlm_nodal_array(acce_filename, model.acceleration[:, unique_indices])
            writedlm_nodal_array(disp_filename, model.displacement[:, unique_indices])
            if bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
                force_filename = prefix * side_set_name * "-force" * index_string * ".csv"
                # See schwarz.get_dst_force() for this
                # Will get projected with Neumann projector onto destination simulation boundary
                force_global = model.internal_force
                force = extract_local_vector(bc, force_global, 3)
                num_nodes = length(bc.global_from_local_map)
                force_out = reshape(force, (3, num_nodes))
                writedlm_nodal_array(force_filename, force_out)
            end
        end
    end
    return nothing
end

function write_stop_exodus(sim::SingleDomainSimulation, model::SolidMechanics)
    params = sim.params
    integrator = sim.integrator
    # A run that continues a sequence of files writes its times after those
    # of the file before it (see run_adaptive).
    time = sim.controller.time + get(params, "exodus_time_offset", 0.0)
    output_mesh = params["output_mesh"]
    # Increment the per-file frame counter and use it as the Exodus time_index.
    # This ensures sequential writes (1, 2, 3, …) even when the global
    # controller stop is offset (e.g. after a mid-run subsim swap).
    params["exodus_frame"] = params["exodus_frame"] + 1
    time_index = params["exodus_frame"]
    Exodus.write_time(output_mesh, time_index, time)

    displacement = model.displacement
    refe_x = model.reference[1, :] + model.displacement[1, :]
    refe_y = model.reference[2, :] + model.displacement[2, :]
    refe_z = model.reference[3, :] + model.displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_x", refe_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_y", refe_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_z", refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_x", disp_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_y", disp_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_z", disp_z)
    # Energetic mesh smoothing size field, sampled at each node's current
    # position X = reference + displacement (see initialize_writing).
    if model isa SolidMechanics && model.mesh_smoothing == true && model.size_field !== nothing
        size_vals = [
            model.size_field((time, refe_x[n], refe_y[n], refe_z[n])) for n in 1:length(refe_x)
        ]
        Exodus.write_values(output_mesh, NodalVariable, time_index, "size", size_vals)
    end
    if model isa SolidMechanics && model.mesh_smoothing == true && model.metric_field !== nothing
        positions = Matrix{Float64}(vcat(refe_x', refe_y', refe_z'))
        sizes, rotation = nodal_metric_output(model, positions, time)
        for i in 1:3
            Exodus.write_values(output_mesh, NodalVariable, time_index, "size_$i", sizes[i, :])
        end
        if rotation !== nothing
            for i in 1:3
                Exodus.write_values(output_mesh, NodalVariable, time_index, "rotation_$i", rotation[i, :])
            end
        end
        metric, target = nodal_metric_tensors(model, positions, time)
        for (i, name) in enumerate(METRIC_TENSOR_NAMES)
            Exodus.write_values(output_mesh, NodalVariable, time_index, name, metric[i, :])
        end
        for (i, name) in enumerate(TARGET_TENSOR_NAMES)
            Exodus.write_values(output_mesh, NodalVariable, time_index, name, target[i, :])
        end
    end
    if is_dynamic(integrator) == true
        velocity = model.velocity
        velo_x = velocity[1, :]
        velo_y = velocity[2, :]
        velo_z = velocity[3, :]
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_x", velo_x)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_y", velo_y)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_z", velo_z)
        acceleration = model.acceleration
        acce_x = acceleration[1, :]
        acce_y = acceleration[2, :]
        acce_z = acceleration[3, :]
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_x", acce_x)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_y", acce_y)
        Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_z", acce_z)
    end
    if !(model.recovery_data isa NoRecovery)
        if model.recovery_data isa BothRecovery
            if size(model.lumped_recovered_stress, 1) > 0
                recover_stress!(model)
                nodal_sigma_c = model.consistent_recovered_stress
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xx_cons_n", nodal_sigma_c[1, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yy_cons_n", nodal_sigma_c[2, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_zz_cons_n", nodal_sigma_c[3, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yz_cons_n", nodal_sigma_c[4, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xz_cons_n", nodal_sigma_c[5, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xy_cons_n", nodal_sigma_c[6, :])
                nodal_sigma_l = model.lumped_recovered_stress
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xx_lump_n", nodal_sigma_l[1, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yy_lump_n", nodal_sigma_l[2, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_zz_lump_n", nodal_sigma_l[3, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yz_lump_n", nodal_sigma_l[4, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xz_lump_n", nodal_sigma_l[5, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xy_lump_n", nodal_sigma_l[6, :])
            end
            if size(model.lumped_recovered_von_mises, 1) > 0
                recover_von_mises_stress!(model)
                Exodus.write_values(output_mesh, NodalVariable, time_index, "von_mises_stress_cons_n", model.consistent_recovered_von_mises[1, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "von_mises_stress_lump_n", model.lumped_recovered_von_mises[1, :])
            end
            if size(model.lumped_recovered_F, 1) > 0
                recover_deformation_gradient!(model)
                nodal_F_c = model.consistent_recovered_F
                nodal_F_l = model.lumped_recovered_F
                for (c, comp) in enumerate(("F_xx", "F_yx", "F_zx", "F_xy", "F_yy", "F_zy", "F_xz", "F_yz", "F_zz"))
                    Exodus.write_values(output_mesh, NodalVariable, time_index, comp * "_cons_n", nodal_F_c[c, :])
                    Exodus.write_values(output_mesh, NodalVariable, time_index, comp * "_lump_n", nodal_F_l[c, :])
                end
            end
            if size(model.lumped_recovered_internal_variables, 1) > 0
                iv_names = collect_internal_variable_names(model.materials)
                recover_internal_variables!(model, iv_names)
                nodal_iv_c = model.consistent_recovered_internal_variables
                nodal_iv_l = model.lumped_recovered_internal_variables
                for (k, name) in enumerate(iv_names)
                    Exodus.write_values(output_mesh, NodalVariable, time_index, name * "_cons_n", nodal_iv_c[k, :])
                    Exodus.write_values(output_mesh, NodalVariable, time_index, name * "_lump_n", nodal_iv_l[k, :])
                end
            end
        else
            if size(model.recovered_stress, 1) > 0
                recover_stress!(model)
                nodal_sigma = model.recovered_stress
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xx_n", nodal_sigma[1, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yy_n", nodal_sigma[2, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_zz_n", nodal_sigma[3, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_yz_n", nodal_sigma[4, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xz_n", nodal_sigma[5, :])
                Exodus.write_values(output_mesh, NodalVariable, time_index, "sigma_xy_n", nodal_sigma[6, :])
            end
            if size(model.recovered_von_mises, 1) > 0
                recover_von_mises_stress!(model)
                Exodus.write_values(output_mesh, NodalVariable, time_index, "von_mises_stress_n", model.recovered_von_mises[1, :])
            end
            if size(model.recovered_F, 1) > 0
                recover_deformation_gradient!(model)
                nodal_F = model.recovered_F
                for (c, name) in enumerate(
                    ("F_xx_n", "F_yx_n", "F_zx_n", "F_xy_n", "F_yy_n", "F_zy_n", "F_xz_n", "F_yz_n", "F_zz_n"),
                )
                    Exodus.write_values(output_mesh, NodalVariable, time_index, name, nodal_F[c, :])
                end
            end
            if size(model.recovered_internal_variables, 1) > 0
                iv_names = collect_internal_variable_names(model.materials)
                recover_internal_variables!(model, iv_names)
                nodal_iv = model.recovered_internal_variables
                for (k, name) in enumerate(iv_names)
                    Exodus.write_values(output_mesh, NodalVariable, time_index, name * "_n", nodal_iv[k, :])
                end
            end
        end
    end
    stress = model.stress
    stored_energy = model.stored_energy
    state = model.state
    all_iv_names = collect_internal_variable_names(model.materials)
    blocks = Exodus.read_sets(output_mesh, Block)
    for (block_index, (block, block_stress, block_stored_energy, block_state)) in
        enumerate(zip(blocks, stress, stored_energy, state))
        block_id = block.id
        element_type_string, num_block_elements, _, _, _, _ = Exodus.read_block_parameters(output_mesh, block_id)
        element_type = element_type_from_string(element_type_string)
        num_points = model.num_int_pts[block_index]
        stress_xx = zeros(num_block_elements, num_points)
        stress_yy = zeros(num_block_elements, num_points)
        stress_zz = zeros(num_block_elements, num_points)
        stress_yz = zeros(num_block_elements, num_points)
        stress_xz = zeros(num_block_elements, num_points)
        stress_xy = zeros(num_block_elements, num_points)
        von_mises = zeros(num_block_elements, num_points)
        for block_element_index in 1:num_block_elements
            element_stress = block_stress[block_element_index]
            for point in 1:num_points
                point_stress = element_stress[point]
                stress_xx[block_element_index, point] = point_stress[1]
                stress_yy[block_element_index, point] = point_stress[2]
                stress_zz[block_element_index, point] = point_stress[3]
                stress_yz[block_element_index, point] = point_stress[4]
                stress_xz[block_element_index, point] = point_stress[5]
                stress_xy[block_element_index, point] = point_stress[6]
                von_mises[block_element_index, point] = von_mises_from_voigt(point_stress)
            end
        end
        mat_iv_names = isempty(all_iv_names) ? String[] : internal_variable_names(model.materials[block_index])
        for point in 1:num_points
            ip_str = "_" * string(point)
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_xx" * ip_str, stress_xx[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_yy" * ip_str, stress_yy[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_zz" * ip_str, stress_zz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_yz" * ip_str, stress_yz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_xz" * ip_str, stress_xz[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "stress_xy" * ip_str, stress_xy[:, point]
            )
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "von_mises_stress" * ip_str, von_mises[:, point]
            )
            for iv_name in all_iv_names
                iv_idx = findfirst(==(iv_name), mat_iv_names)
                if iv_idx !== nothing
                    values = [block_state[elem][point][iv_idx] for elem in 1:num_block_elements]
                else
                    values = zeros(num_block_elements)
                end
                Exodus.write_values(
                    output_mesh, ElementVariable, time_index, Int64(block_id), iv_name * ip_str, values
                )
            end
        end
        Exodus.write_values(
            output_mesh, ElementVariable, time_index, Int64(block_id), "stored_energy", block_stored_energy
        )
        if model isa SolidMechanics && model.mesh_smoothing == true
            Exodus.write_values(
                output_mesh, ElementVariable, time_index, Int64(block_id), "energy_density",
                block_stored_energy ./ ideal_volumes(model, block_index),
            )
        end
    end
    return nothing
end

include("opinf/opinf_io.jl")
