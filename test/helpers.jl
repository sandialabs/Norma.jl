# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
function average_components(V::Vector{Float64})
    x = mean(V[1:3:end])
    y = mean(V[2:3:end])
    z = mean(V[3:3:end])
    return [x y z]
end

function maximum_components(V::Vector{Float64})
    x = maximum(V[1:3:end])
    y = maximum(V[2:3:end])
    z = maximum(V[3:3:end])
    return [x y z]
end

function minimum_components(V::Vector{Float64})
    x = minimum(V[1:3:end])
    y = minimum(V[2:3:end])
    z = minimum(V[3:3:end])
    return [x y z]
end

function average_components(stress::Vector{Vector{Vector{Vector{Float64}}}})
    xx = 0.0
    yy = 0.0
    zz = 0.0
    yz = 0.0
    xz = 0.0
    xy = 0.0
    num_stress = 0
    for block_index in 1:length(stress)
        for block_element_index in 1:length(stress[block_index])
            for point in 1:length(stress[block_index][block_element_index])
                xx += stress[block_index][block_element_index][point][1]
                yy += stress[block_index][block_element_index][point][2]
                zz += stress[block_index][block_element_index][point][3]
                yz += stress[block_index][block_element_index][point][4]
                xz += stress[block_index][block_element_index][point][5]
                xy += stress[block_index][block_element_index][point][6]
                num_stress += 1
            end
        end
    end
    return [xx yy zz yz xz xy] ./ num_stress
end

using Exodus
using Symbolics
@variables t, x, y, z

function create_force(expression::String, mesh::ExodusDatabase, side_set_id::Int64, time::Float64)
    force_num = eval(Meta.parse(expression))
    force_fun = eval(build_function(force_num, [t, x, y, z]; expression=Val(false)))
    coords = read_coordinates(mesh)
    num_nodes = size(coords)[2]
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    local_from_global_map = Norma.get_side_set_local_from_global_map(mesh, side_set_id)
    num_nodes = length(local_from_global_map)
    force = zeros(num_nodes)
    ss_node_index = 1
    for side in num_nodes_sides
        side_nodes = side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = coords[:, side_nodes]
        nodal_force_component = Norma.get_side_set_nodal_forces(side_coordinates, force_fun, time)
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        force[local_indices] += nodal_force_component
        ss_node_index += side
    end
    return force
end

function create_displacement(expression::String, mesh::ExodusDatabase, side_set_id::Int64, time::Float64)
    disp_num = eval(Meta.parse(expression))
    coords = read_coordinates(mesh)
    num_nodes = size(coords)[2]
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)
    local_from_global_map = Norma.get_side_set_local_from_global_map(mesh, side_set_id)
    num_nodes = length(local_from_global_map)
    displacements = zeros(num_nodes)
    ss_node_index = 1
    for side in num_nodes_sides
        side_nodes = side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        for side_node in side_nodes
            node_coordinates = coords[:, side_node]
            values = Dict(t => time, x => node_coordinates[1], y => node_coordinates[2], z => node_coordinates[3])
            disp_sym = substitute(disp_num, values)
            disp_val = Norma.extract_value(disp_sym)
            local_index = get.(Ref(local_from_global_map), side_node, 0)
            displacements[local_index] = disp_val
        end
        ss_node_index += side
    end
    return displacements
end
