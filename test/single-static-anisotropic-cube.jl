# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.
using LinearAlgebra
using StaticArrays
using YAML

# A copy of `source` with the given element variables on its single block at
# time step 1.
function write_element_variables_mesh(source::String, target::String, variables::Vector{Pair{String,Vector{Float64}}})
    input = Exodus.ExodusDatabase(source, "r")
    Exodus.copy(input, target)
    close(input)
    mesh = Exodus.ExodusDatabase(target, "rw")
    Exodus.write_number_of_variables(mesh, ElementVariable, length(variables))
    Exodus.write_names(mesh, ElementVariable, first.(variables))
    Exodus.write_time(mesh, 1, 0.0)
    for (name, values) in variables
        Exodus.write_values(mesh, ElementVariable, 1, 1, name, values)
    end
    close(mesh)
    return nothing
end

voigt_to_matrix(v) = @SMatrix [v[1] v[6] v[5]; v[6] v[2] v[4]; v[5] v[4] v[3]]

# Largest deviation of the stress at any integration point from σ_ref.
function stress_deviation(sim, σ_ref)
    return maximum(maximum(abs, voigt_to_matrix(p) - σ_ref) for element in sim.model.stress[1] for p in element)
end

@testset "Single Static Anisotropic Cube" begin
    cp("../examples/single/static-solid/cube/anisotropic/cube-shear.yaml", "cube-shear.yaml"; force=true)
    cp("../examples/single/static-solid/cube/anisotropic/cube.g", "cube.g"; force=true)
    base = YAML.load_file("cube-shear.yaml"; dicttype=Norma.Parameters)
    base["name"] = "cube-shear.yaml"
    function run_shear(orientation; mesh="cube.g", material_params=nothing)
        params = deepcopy(base)
        params["input mesh file"] = mesh
        material_params === nothing || (params["model"]["material"]["copper"] = deepcopy(material_params))
        material = params["model"]["material"]["copper"]
        if orientation === nothing
            delete!(material, "orientation")
        else
            material["orientation"] = orientation
        end
        return Norma.run(params)
    end
    F = @SMatrix [1.0 0.1 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    φ = (0.12, 0.03, 0.28)
    Q = Norma.orientation_from_bunge_angles(φ...)

    @testset "Aligned Lattice Matches LCM" begin
        sim = run_shear(nothing)
        σ_lcm = lcm_tensor(only(lcm_final_rows("ElasticShear3d")), "Cauchy_Stress")
        @test stress_deviation(sim, σ_lcm) < 1.0e-10 * maximum(abs, σ_lcm)
    end

    @testset "Rotated Lattice Matches LCM" begin
        # LCM receives the orientation as basis vectors rounded to eight
        # digits, hence the looser comparison with its gold.
        sim = run_shear(Dict("Bunge Euler angles" => collect(φ), "angle unit" => "radians"))
        σ_lcm = lcm_tensor(only(lcm_final_rows("ElasticShearRotated3d")), "Cauchy_Stress")
        @test stress_deviation(sim, σ_lcm) < 1.0e-7 * maximum(abs, σ_lcm)
        σ_ref = cubic_svk_cauchy(168.4, 121.4, 75.4, Q, F)
        @test stress_deviation(sim, σ_ref) < 1.0e-10 * maximum(abs, σ_ref)
    end

    @testset "Orientation From Mesh Variables" begin
        num_elements = 8
        # A different orientation for each element, and the same orientation
        # as Bunge angles and as the rotation matrix, row-major.
        varied = [(0.1 * e, 0.05 * e, -0.07 * e) for e in 1:num_elements]
        variables = Pair{String,Vector{Float64}}[
            "phi1" => [a[1] for a in varied],
            "Phi" => [a[2] for a in varied],
            "phi2" => [a[3] for a in varied],
            "uniform_phi1" => fill(φ[1], num_elements),
            "uniform_Phi" => fill(φ[2], num_elements),
            "uniform_phi2" => fill(φ[3], num_elements),
        ]
        for i in 1:3, j in 1:3
            push!(variables, "Q$i$j" => fill(Q[i, j], num_elements))
        end
        write_element_variables_mesh("cube.g", "cube-oriented.g", variables)
        sim_varied = run_shear(
            Dict("mesh variables" => ["phi1", "Phi", "phi2"], "angle unit" => "radians"); mesh="cube-oriented.g"
        )
        @test all(
            sim_varied.model.lattice_orientation[1][e] ≈ Norma.orientation_from_bunge_angles(varied[e]...) for
            e in 1:num_elements
        )
        sim_angles = run_shear(
            Dict("mesh variables" => ["uniform_phi1", "uniform_Phi", "uniform_phi2"], "angle unit" => "radians");
            mesh="cube-oriented.g",
        )
        sim_matrix = run_shear(Dict("mesh variables" => ["Q$i$j" for i in 1:3 for j in 1:3]); mesh="cube-oriented.g")
        σ_ref = cubic_svk_cauchy(168.4, 121.4, 75.4, Q, F)
        @test stress_deviation(sim_angles, σ_ref) < 1.0e-10 * maximum(abs, σ_ref)
        @test stress_deviation(sim_matrix, σ_ref) < 1.0e-10 * maximum(abs, σ_ref)
        rm("cube-oriented.g"; force=true)
    end

    @testset "Fiber-Reinforced Damage Matches The Material Point" begin
        # Fibers along [110] of the material frame, which the orientation
        # rotates, are stretched by the shear and damage as it grows. The
        # deformation is homogeneous, so every integration point follows the
        # material point driver along the same steps.
        ahd = Dict{String,Any}(
            "model" => "anisotropic hyperelastic damage",
            "elastic modulus" => 200.0,
            "Poisson's ratio" => 0.3,
            "density" => 1.0,
            "matrix volume fraction" => 0.6,
            "matrix maximum damage" => 0.2,
            "matrix damage saturation" => 0.05,
            "fibers" => [
                Dict(
                    "direction" => [1.0, 1.0, 0.0],
                    "k" => 50.0,
                    "q" => 5.0,
                    "volume fraction" => 0.4,
                    "maximum damage" => 0.9,
                    "damage saturation" => 0.5,
                ),
            ],
        )
        orientation = Dict("Bunge Euler angles" => collect(φ), "angle unit" => "radians")
        sim = run_shear(orientation; material_params=ahd)
        material = Norma.create_material(Norma.Parameters(ahd))
        times = collect(0.0:0.01:0.1)
        history = Norma.material_point_history(
            material, s -> SMatrix{3,3,Float64,9}([1.0 s 0; 0 1 0; 0 0 1]), times; orientation=Q
        )
        @test all(history.state[end][3:4] .> 0.0)
        @test stress_deviation(sim, history.σ[end]) < 1.0e-10 * maximum(abs, history.σ[end])
        @test all(p ≈ history.state[end] for element in sim.model.state_old[1] for p in element)
    end

    rm("cube-shear.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube-shear.e"; force=true)
end
