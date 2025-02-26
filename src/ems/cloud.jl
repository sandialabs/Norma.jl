# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using YAML

abstract type Geometry end

struct Cube <: Geometry
    side_length::Float64
end

struct Sphere <: Geometry
    radius::Float64
end

abstract type Potential end

struct Lennard_Jones <: Potential
    ε::Float64
    σ::Float64
end

struct Repulsive <: Potential
    k::Float64
    c::Float64
end

struct SymmetricSethHill <: Potential
    k::Float64
    h::Float64
    m::Int64
    c::Float64
end

struct Polynomial <: Potential
    k::Float64
    h::Float64
    m::Int64
    c::Float64
end

abstract type Solver end

struct SteepestDescent <: Solver
    step_length::Float64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    line_search_iterations::Int64
    backtrack_factor::Float64
    decrease_factor::Float64
    output_file_root::String
    output_interval::Int64
end

function run(input_file::String)
    params = create_params(input_file)
    geometry, points = create_point_cloud(params)
    potential = create_potential(params)
    solver = create_solver(params)
    solve(solver, geometry, potential, points)
end

function create_params(input_file::String)
    println("Reading simulation file: ", input_file)
    params = YAML.load_file(input_file; dicttype = Parameters)
    return params
end

function create_geometry(params::Parameters)
    geometry_params = params["geometry"]
    type = geometry_params["type"]
    if type == "cube"
        side_length = geometry_params["side length"]
        return Cube(side_length)
    elseif type == "sphere"
        radius = geometry_params["radius"]
        return Sphere(radius)
    else
        error("Unrecognized geometry type: ", type, ", options are cube or sphere")
    end
end

function create_potential(params::Parameters)
    potential_params = params["potential"]
    type = potential_params["type"]
    if type == "Lennard-Jones"
        ε = potential_params["ε"]
        σ = potential_params["σ"]
        return Lennard_Jones(ε, σ)
    elseif type == "repulsive"
        k = potential_params["k"]
        c = potential_params["cutoff"]
        return Repulsive(k, c)
    elseif type == "symmetric Seth-Hill"
        k = potential_params["k"]
        h = potential_params["h"]
        m = potential_params["m"]
        c = potential_params["cutoff"]
        return SymmetricSethHill(k, h, m, c)
    elseif type == "polynomial"
        k = potential_params["k"]
        h = potential_params["h"]
        m = potential_params["m"]
        c = potential_params["cutoff"]
        return Polynomial(k, h, m, c)
    else
        error(
            "Unrecognized potential type: ",
            type,
            ", options are Lennard-Jones, repulsive, or symmetric Seth-Hill",
        )
    end
end

function create_solver(params::Parameters)
    solver_params = params["solver"]
    type = solver_params["type"]
    if type == "steepest descent"
        step_length = solver_params["step length"]
        maximum_iterations = solver_params["maximum iterations"]
        absolute_tolerance = solver_params["absolute tolerance"]
        relative_tolerance = solver_params["relative tolerance"]
        line_search_iterations = solver_params["maximum line search iterations"]
        backtrack_factor = solver_params["backtrack factor"]
        decrease_factor = solver_params["decrease factor"]
        output_file_root = solver_params["output file root"]
        output_interval = solver_params["output interval"]
        return SteepestDescent(
            step_length,
            maximum_iterations,
            absolute_tolerance,
            relative_tolerance,
            line_search_iterations,
            backtrack_factor,
            decrease_factor,
            output_file_root,
            output_interval,
        )
    else
        error("Unrecognized solver type: ", type, ", options are steepest descent")
    end
end

function is_inside(cube::Cube, point::Vector{Float64})
    return norm(point, Inf) ≤ 0.5 * cube.side_length
end

function is_inside(sphere::Sphere, point::Vector{Float64})
    return norm(point, 2) ≤ sphere.radius
end

function rescale(cube::Cube, point::Vector{Float64})
    new_point = clamp.(point, -0.5 * cube.side_length, 0.5 * cube.side_length)
    return new_point
end

function rescale(sphere::Sphere, point::Vector{Float64})
    n = norm(point, 2)
    new_point = (sphere.radius / n) * point
    return new_point
end

function return_map!(geometry::Geometry, points::Matrix{Float64})
    num_points = size(points, 2)
    for point_index ∈ 1:num_points
        point = points[:, point_index]
        if is_inside(geometry, point) == false
            new_point = rescale(geometry, point)
            points[:, point_index] = new_point
        end
    end
end

function create_point_cloud(params::Parameters)
    geometry = create_geometry(params)
    number_points = params["number points"]
    points = zeros(3, number_points)
    point_index = 1
    while point_index ≤ number_points
        trial_point = 2.0 * rand(3) .- 1.0
        if is_inside(geometry, trial_point) == false
            continue
        end
        points[:, point_index] = trial_point
        point_index += 1
    end
    if number_points == 2
        points = [-0.05 0.05; 0.0 0.0; 0.0 0.0]
    end
    if number_points == 4
        points = [-0.05 0.05 0.05 -0.05; -0.05 -0.05 0.05 0.05; 0.0 0.0 0.0 0.0]
    end
    return geometry, points
end

function compute_energy_force!(
    potential::Lennard_Jones,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    ε = potential.ε
    σ = potential.σ
    number_points = size(points, 2)
    energy = 0.0
    fill!(force, 0.0)
    for i ∈ 1:number_points-1
        for j ∈ i+1:number_points
            r_ij = points[:, j] - points[:, i]
            r = max(norm(r_ij), 1.0e-06) # Avoid singularity
            energy += 4 * ε * ((σ / r)^12 - (σ / r)^6)
            magnitude = 24 * ε * (2 * (σ / r)^12 - (σ / r)^6) / r^2
            f_ij = magnitude * r_ij
            force[:, i] .+= f_ij
            force[:, j] .-= f_ij
        end
    end
    return energy
end

function compute_energy_force!(
    potential::Repulsive,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    k = potential.k
    c = potential.c
    number_points = size(points, 2)
    energy = 0.0
    fill!(force, 0.0)
    for i ∈ 1:number_points-1
        for j ∈ i+1:number_points
            r_ij = points[:, j] - points[:, i]
            r = max(norm(r_ij), 1.0e-06) # Avoid singularity
            if r > c
                continue
            end
            energy += k / r
            f_ij = k / r^3 * r_ij
            force[:, i] .-= f_ij
            force[:, j] .+= f_ij
        end
    end
    return energy
end

function compute_energy_force!(
    potential::SymmetricSethHill,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    k = potential.k
    h = potential.h
    m = potential.m
    c = potential.c
    number_points = size(points, 2)
    energy = 0.0
    fill!(force, 0.0)
    for i ∈ 1:number_points-1
        for j ∈ i+1:number_points
            r_ij = points[:, j] - points[:, i]
            r = max(norm(r_ij), 1.0e-06) # Avoid singularity
            if r > c
                continue
            end
            γ = r / h
            κ = γ^m
            α = κ - 1.0
            β = 1.0 / κ - 1.0
            energy += (0.5 * k / m / m) * (α * α + β * β)
            f_ij = (-k / m) * (α * κ / γ - β / κ / γ) * (r_ij / r)
            force[:, i] .-= f_ij
            force[:, j] .+= f_ij
        end
    end
    return energy
end

function compute_energy_force!(
    potential::Polynomial,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    k = potential.k
    h = potential.h
    m = potential.m
    c = potential.c
    number_points = size(points, 2)
    energy = 0.0
    fill!(force, 0.0)
    for i ∈ 1:number_points-1
        for j ∈ i+1:number_points
            r_ij = points[:, j] - points[:, i]
            r = max(norm(r_ij), 1.0e-06) # Avoid singularity
            if r > c
                continue
            end
            y = r / h
            x = log(y)
            xm = x^m
            energy += 0.5 * k * xm
            f_ij = (-0.5 * k * m) * (xm / x / y) * (r_ij / r)
            force[:, i] .-= f_ij
            force[:, j] .+= f_ij
        end
    end
    return energy
end

function solve(
    solver::SteepestDescent,
    geometry::Geometry,
    potential::Potential,
    points::Matrix{Float64},
)
    maximum_iterations = solver.maximum_iterations
    absolute_tolerance = solver.absolute_tolerance
    relative_tolerance = solver.relative_tolerance
    output_file_root = solver.output_file_root
    output_interval = solver.output_interval
    force = zeros(size(points))
    energy = compute_energy_force!(potential, points, force)
    initial_norm = norm(force)
    output_file = output_file_root * "-" * lpad(string(0), 6, '0') * ".vtk"
    write_points_vtk(output_file, points)
    for iter ∈ 1:maximum_iterations
        absolute_error = norm(force)
        relative_error = initial_norm > 0.0 ? absolute_error / initial_norm : absolute_error
        println(
            "Iteration : ",
            iter,
            ", energy : ",
            energy,
            ", absolute error : ",
            absolute_error,
            ", relative error : ",
            relative_error,
        )
        if absolute_error < absolute_tolerance || relative_error < relative_tolerance
            println("Converged")
            println(
                "Absolute tolerance : ",
                absolute_tolerance,
                ", relative tolerance : ",
                relative_tolerance,
            )
            break
        end
        points, energy = backtrack_line_search!(solver, geometry, potential, points, force)
        if iter % output_interval == 0
            output_file = output_file_root * "-" * lpad(string(iter), 6, '0') * ".vtk"
            write_points_vtk(output_file, points)
        end
    end
end

function backtrack_line_search!(
    solver::Solver,
    geometry::Geometry,
    potential::Potential,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    backtrack_factor = solver.backtrack_factor
    decrease_factor = solver.decrease_factor
    max_iters = solver.line_search_iterations
    energy = compute_energy_force!(potential, points, force)
    direction = force
    step_length = solver.step_length
    initial_solution = 1.0 * points
    initial_energy = energy
    for _ ∈ 1:max_iters
        step = step_length * direction
        points = initial_solution + step
        return_map!(geometry, points)
        energy = compute_energy_force!(potential, points, force)
        if energy <= initial_energy + decrease_factor * step_length * dot(-force, direction)
            return points, energy
        end
        step_length *= backtrack_factor
    end
    @warn "Line search did nor converge in ", max_iters, " iterations, taking full step"
    # Exceeded line search iterations. Take full step and hope for the best.
    step = solver.step_length * direction
    points = initial_solution + step
    return_map!(geometry, points)
    energy, _ = compute_energy_force!(potential, points, force)
    return points, energy
end

function backtrack_line_search_orig!(
    solver::Solver,
    geometry::Geometry,
    potential::Potential,
    points::Matrix{Float64},
    force::Matrix{Float64},
)
    backtrack_factor = solver.backtrack_factor
    decrease_factor = solver.decrease_factor
    max_iters = solver.line_search_iterations
    energy, search_distance = compute_energy_force!(potential, points, force)
    direction = search_distance * force / norm(force)
    merit = energy
    merit_prime = -2.0 * merit
    step_length = solver.step_length
    initial_solution = 1.0 * points
    for _ ∈ 1:max_iters
        merit_old = merit
        step = step_length * direction
        points = initial_solution + step
        return_map!(geometry, points)
        energy, search_distance = compute_energy_force!(potential, points, force)
        merit = energy
        if merit ≤ (1.0 - 2.0 * decrease_factor * step_length) * merit_old
            break
        end
        step_length =
            max(
                backtrack_factor * step_length,
                -0.5 * step_length * step_length * merit_prime,
            ) / (merit - merit_old - step_length * merit_prime)
    end
    return points, energy
end

function write_points_vtk(filename::String, points::Matrix{Float64})
    # Ensure the input is 3xN
    if size(points, 1) != 3
        error("The input matrix must have size 3xN, where N is the number of points.")
    end

    N = size(points, 2)  # Number of points

    # Open the file for writing
    open(filename, "w") do io
        # Write VTK header
        println(io, "# vtk DataFile Version 3.0")
        println(io, "3D points")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        println(io, "POINTS $N float")

        # Write point data
        for i = 1:N
            println(io, join(points[:, i], " "))
        end
    end
end

for input_file ∈ ARGS
    run(input_file)
end
