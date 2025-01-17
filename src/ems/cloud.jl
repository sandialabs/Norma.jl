# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using YAML

@enum Geometry Cube Sphere

function run(input_file::String)
    params = create_params(input_file)
    create_point_cloud(params)
end

function create_params(input_file::String)
    println("Reading simulation file: ", input_file)
    params = YAML.load_file(input_file; dicttype=Dict{String,Any})
    return params
end

function create_geometry(params::Dict{String,Any})
    geometry_string = params["geometry"]
    if geometry_string == "cube"
        return Cube
    elseif geometry_string == "sphere"
        return Sphere
    else
        error("Unrecognized geometry, options are cube or sphere")
    end
end

function create_point_cloud(params::Dict{String,Any})
    geometry = create_geometry(params)
    number_points = params["number points"]
    radius = params["radius"]
    points = zeros(3, number_points)
    point_index = 1
    while point_index ≤ number_points
        trial_point = 2.0 * rand(3) .- 1.0
        if geometry == Cube && norm(trial_point, Inf) > radius
            continue
        elseif geometry == Sphere && norm(trial_point, 2) > radius
            continue
        else
            error("Unrecognized geometry, options are cube or sphere")
        end
        points[:, point_index] = trial_point
        point_index += 1
    end
    return points
end