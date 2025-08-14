# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using DelimitedFiles
using Format

function initialize_writing(sim::SingleDomainSimulation)
    params = sim.params
    integrator = sim.integrator
    output_mesh = params["output_mesh"]

    # global variables
    num_global_vars = Exodus.read_number_of_variables(output_mesh, GlobalVariable)
    num_global_vars += 1  # adjust as necessary
    Exodus.write_number_of_variables(output_mesh, GlobalVariable, num_global_vars)

    # setup nodal variables
    num_node_vars = 6
    node_var_names = ["refe_x", "refe_y", "refe_z", "disp_x", "disp_y", "disp_z"]
    if is_dynamic(integrator) == true
        num_node_vars += 6
        append!(node_var_names, ["velo_x", "velo_y", "velo_z", "acce_x", "acce_y", "acce_z"])
    end
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_names(output_mesh, NodalVariable, node_var_names)

    # get maximum number of quadrature points
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for block in blocks
        block_id = block.id
        element_type_string = Exodus.read_block_parameters(output_mesh, block_id)[1]
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end

    num_element_vars = 6 * max_num_int_points + 1
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
    end
    push!(el_var_names, "stored_energy")
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

function write_stop(sim::SingleDomainSimulation)
    params = sim.params
    stop = sim.controller.stop
    num_steps = sim.controller.num_stops - 1
    time = sim.controller.time
    name = sim.name
    if haskey(params, "parent_simulation") == false
        percent = 100 * stop / num_steps
        digits = max(0, Int64(ceil(log10(num_steps))) - 2)
        norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.4e", stop, num_steps, percent, time)
    end
    exodus_interval = get(params, "Exodus output interval", 1)::Int64
    if exodus_interval > 0 && stop % exodus_interval == 0
        norma_log(0, :output, "Exodus II Database for $name [EXO]")
        write_stop_exodus(sim, sim.model)
    end
    csv_interval = get(params, "CSV output interval", 0)::Int64
    if csv_interval > 0 && stop % csv_interval == 0
        norma_log(0, :output, "Comma Separated Values for $name [CSV]")
        write_stop_csv(sim, sim.model)
        if haskey(params, "CSV write sidesets") == true
            write_sideset_stop_csv(sim, sim.model)
        end
    end
    return nothing
end

function write_stop(sim::MultiDomainSimulation)
    stop = sim.controller.stop
    num_steps = sim.controller.num_stops - 1
    time = sim.controller.time
    percent = 100 * stop / num_steps
    digits = max(0, Int64(ceil(log10(num_steps))) - 2)
    norma_logf(0, :stop, "[%d/%d, %.$(digits)f%%] : Time = %.4e", stop, num_steps, percent, time)
    for subsim in sim.subsims
        write_stop(subsim)
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
    free_dofs_filename = prefix * "free_dofs" * index_string * ".csv"
    writedlm(free_dofs_filename, model.free_dofs)
    curr_filename = prefix * "curr" * index_string * ".csv"
    writedlm_nodal_array(curr_filename, model.current)
    disp_filename = prefix * "disp" * index_string * ".csv"
    writedlm_nodal_array(disp_filename, model.current - model.reference)
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
            writedlm(curr_filename, model.current[bc.offset, bc.node_set_node_indices])
            writedlm(velo_filename, model.velocity[bc.offset, bc.node_set_node_indices])
            writedlm(acce_filename, model.acceleration[bc.offset, bc.node_set_node_indices])
            writedlm(
                disp_filename,
                model.current[bc.offset, bc.node_set_node_indices] -
                model.reference[bc.offset, bc.node_set_node_indices],
            )
        elseif bc isa SolidMechanicsOverlapSchwarzBoundaryCondition ||
            bc isa SolidMechanicsNonOverlapSchwarzBoundaryCondition
            side_set_name = bc.name
            curr_filename = prefix * side_set_name * "-curr" * index_string * ".csv"
            disp_filename = prefix * side_set_name * "-disp" * index_string * ".csv"
            velo_filename = prefix * side_set_name * "-velo" * index_string * ".csv"
            acce_filename = prefix * side_set_name * "-acce" * index_string * ".csv"
            unique_indices = unique(bc.side_set_node_indices)
            writedlm_nodal_array(curr_filename, model.current[:, unique_indices])
            writedlm_nodal_array(velo_filename, model.velocity[:, unique_indices])
            writedlm_nodal_array(acce_filename, model.acceleration[:, unique_indices])
            writedlm_nodal_array(disp_filename, model.current[:, unique_indices] - model.reference[:, unique_indices])
        end
    end
    return nothing
end


function write_stop_exodus(sim::SingleDomainSimulation, model::SolidMechanics)
    params = sim.params
    integrator = sim.integrator
    time = sim.controller.time
    stop = sim.controller.stop
    time_index = stop + 1
    output_mesh = params["output_mesh"]
    Exodus.write_time(output_mesh, time_index, time)

    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_x", refe_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_y", refe_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_z", refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_x", disp_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_y", disp_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_z", disp_z)
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
    stress = model.stress
    stored_energy = model.stored_energy
    blocks = Exodus.read_sets(output_mesh, Block)
    for (block, block_stress, block_stored_energy) in zip(blocks, stress, stored_energy)
        block_id = block.id
        element_type_string, num_block_elements, _, _, _, _ = Exodus.read_block_parameters(output_mesh, block_id)
        element_type = element_type_from_string(element_type_string)
        num_points = default_num_int_pts(element_type)
        stress_xx = zeros(num_block_elements, num_points)
        stress_yy = zeros(num_block_elements, num_points)
        stress_zz = zeros(num_block_elements, num_points)
        stress_yz = zeros(num_block_elements, num_points)
        stress_xz = zeros(num_block_elements, num_points)
        stress_xy = zeros(num_block_elements, num_points)
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
            end
        end
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
        end
        Exodus.write_values(
            output_mesh, ElementVariable, time_index, Int64(block_id), "stored_energy", block_stored_energy
        )
    end
    return nothing
end

include("opinf/opinf_io.jl")

