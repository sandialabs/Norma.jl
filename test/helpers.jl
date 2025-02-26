# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
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
    for blk_index ∈ 1:length(stress)
        for blk_elem_index ∈ 1:length(stress[blk_index])
            for point ∈ 1:length(stress[blk_index][blk_elem_index])
                xx += stress[blk_index][blk_elem_index][point][1]
                yy += stress[blk_index][blk_elem_index][point][2]
                zz += stress[blk_index][blk_elem_index][point][3]
                yz += stress[blk_index][blk_elem_index][point][4]
                xz += stress[blk_index][blk_elem_index][point][5]
                xy += stress[blk_index][blk_elem_index][point][6]
                num_stress += 1
            end
        end
    end
    return [xx yy zz yz xz xy] ./ num_stress
end

using Exodus
using Symbolics
@variables t

function get_boundary_traction_force(mesh::ExodusDatabase, side_set_id::Int64)
    expression = "1.0 * t"
    t = 1.0
    traction_num = eval(Meta.parse(expression))
    coords = read_coordinates(mesh)
    num_nodes = size(coords)[2]
    local_from_global_map, num_nodes_sides, side_set_node_indices =
        Norma.get_side_set_local_from_global_map(mesh, side_set_id)
    num_nodes = length(local_from_global_map)
    boundary_tractions_force = zeros(num_nodes)
    ss_node_index = 1
    for side ∈ num_nodes_sides
        side_nodes = side_set_node_indices[ss_node_index:ss_node_index+side-1]
        side_coordinates = coords[:, side_nodes]
        nodal_force_component =
            Norma.get_side_set_nodal_forces(side_coordinates, traction_num, t)
        local_indices = get.(Ref(local_from_global_map), side_nodes, 0)
        boundary_tractions_force[local_indices] += nodal_force_component
        ss_node_index += side
    end
    return boundary_tractions_force
end
