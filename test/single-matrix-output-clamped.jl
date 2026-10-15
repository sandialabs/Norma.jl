# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Matrix output: the stiffness, mass, external force and free-DOF map written
# by `matrix output` must round-trip exactly to the operators Norma used, and
# together with Norma's own solution they must satisfy the discrete equations
# of motion, which is what an external solver relies on.

using DelimitedFiles
using SparseArrays
using YAML

# Minimal Matrix Market reader (coordinate, real, general), independent of the
# writer under test.
function read_matrix_market_test(filename::String)
    lines = filter(line -> !startswith(line, "%"), readlines(filename))
    num_rows, num_cols, num_entries = parse.(Int64, split(lines[1]))
    I = Vector{Int64}(undef, num_entries)
    J = Vector{Int64}(undef, num_entries)
    V = Vector{Float64}(undef, num_entries)
    for k in 1:num_entries
        fields = split(lines[k + 1])
        I[k] = parse(Int64, fields[1])
        J[k] = parse(Int64, fields[2])
        V[k] = parse(Float64, fields[3])
    end
    return sparse(I, J, V, num_rows, num_cols)
end

read_free_dofs_test(filename::String) = BitVector(strip.(readlines(filename)) .== "true")

# Runs examples/ahead/single/clamped/dynamic-linear-matrix-output, the input used to
# generate data for an external solver, shortened to 10 steps.
@testset "Single Matrix Output Clamped" begin
    cp("../examples/ahead/single/clamped/dynamic-linear-matrix-output/clamped.yaml", "clamped.yaml"; force=true)
    cp("../examples/ahead/single/clamped/clamped.g", "../clamped.g"; force=true)
    params = YAML.load_file("clamped.yaml"; dicttype=Norma.Parameters)
    params["name"] = "clamped"
    params["time integrator"]["final time"] = 1.0e-05  # 10 steps of 1.0e-06; CSV output at stop 0
    params["matrix output"]["steps"] = [0, 10]
    simulation = Norma.run(params)
    model = simulation.model
    integrator = simulation.integrator
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes

    expected_files = [
        "clamped-stiffness-0000.mm",
        "clamped-stiffness-0010.mm",
        "clamped-mass.mm",
        "clamped-external-force-0000.csv",
        "clamped-external-force-0010.csv",
        "clamped-free-dofs-0000.csv",
        "clamped-free-dofs-0010.csv",
        "clamped-disp-0000.csv",
        "clamped-velo-0000.csv",
        "clamped-acce-0000.csv",
    ]
    for file in expected_files
        @test isfile(file)
    end
    # Only the requested stops are written, and the mass only once.
    @test !isfile("clamped-stiffness-0005.mm")
    @test !isfile("clamped-free-dofs-0005.csv")
    @test !isfile("clamped-lumped-mass.mm")

    read_nodal_csv(filename) = vec(readdlm(filename, ',', Float64)')
    K = read_matrix_market_test("clamped-stiffness-0000.mm")
    K10 = read_matrix_market_test("clamped-stiffness-0010.mm")
    M = read_matrix_market_test("clamped-mass.mm")
    f0 = read_nodal_csv("clamped-external-force-0000.csv")
    f10 = read_nodal_csv("clamped-external-force-0010.csv")
    free0 = read_free_dofs_test("clamped-free-dofs-0000.csv")
    free10 = read_free_dofs_test("clamped-free-dofs-0010.csv")
    u0 = read_nodal_csv("clamped-disp-0000.csv")
    v0 = read_nodal_csv("clamped-velo-0000.csv")
    a0 = read_nodal_csv("clamped-acce-0000.csv")

    for file in readdir()
        if startswith(file, "clamped-") && (endswith(file, ".csv") || endswith(file, ".mm"))
            rm(file; force=true)
        end
    end
    rm("clamped.yaml"; force=true)
    rm("../clamped.g"; force=true)
    rm("clamped.e"; force=true)

    # Exact round trip of what Norma assembled.
    @test size(K) == (num_dofs, num_dofs)
    @test size(M) == (num_dofs, num_dofs)
    @test K == model.stiffness
    @test M == model.mass
    @test K10 == K  # linear elastic: the stiffness does not change
    @test length(f0) == num_dofs
    @test f10 == model.body_force + model.boundary_force
    @test norm(f0) == 0.0
    @test free0 == model.free_dofs
    @test free10 == model.free_dofs

    # Operator properties.
    @test norm(K - K') <= 1.0e-12 * norm(K)
    @test norm(M - M') <= 1.0e-12 * norm(M)
    # Total mass: rho * volume per component, 1000 * (0.001 * 0.001 * 1.0).
    @test sum(M) ≈ 3 * 1000.0 * 1.0e-6 rtol = 1.0e-12
    # Rigid translation in z produces no internal force.
    rigid_z = zeros(num_dofs)
    rigid_z[3:3:end] .= 1.0
    @test norm(K * rigid_z) <= 1.0e-10 * norm(K) * norm(rigid_z)

    # Constrained DOFs: x on nsx±, y on nsy±, z on nsz±. Interior lines of the
    # 1 x 1 x N column keep only z free, end nodes have nothing free.
    @test count(.!free0) > 0
    @test count(free0) > 0
    @test all(.!free0[1:3:end])  # every node lies on nsx- or nsx+
    @test all(.!free0[2:3:end])  # every node lies on nsy- or nsy+

    # Initial conditions as an external solver reads them: the Gaussian pulse
    # in z, at rest, and the initial acceleration solves M a0 = f - K u0 on
    # the free DOFs (the prescribed accelerations are zero).
    z = model.reference[3, :]
    @test u0[3:3:end] ≈ 0.001 .* exp.(-z .^ 2 ./ (2 * 0.02^2)) atol = 1.0e-15
    @test norm(u0[1:3:end]) == 0.0
    @test norm(u0[2:3:end]) == 0.0
    @test norm(v0) == 0.0
    rhs0 = (f0 - K * u0)[free0]
    @test norm(M[free0, free0] * a0[free0] - rhs0) <= 1.0e-10 * norm(rhs0)
    @test norm(a0[.!free0]) == 0.0

    # The final state satisfies the Newmark equations of motion with the
    # dumped operators: K u + M a = f on the free DOFs.
    u = vec(integrator.displacement)
    a = vec(integrator.acceleration)
    residual = K * u + M * a - f10
    @test norm(residual[free10]) <= 1.0e-6 * norm((M * a)[free10])

    # Input parsing.
    Norma.NORMA_TEST_MODE[] = true
    parse_block(block) = Norma.parse_matrix_output(Norma.Parameters("matrix output" => block))
    default_options = parse_block(Norma.Parameters())
    @test default_options.matrices == ["stiffness", "mass", "external force"]
    @test Norma.is_matrix_output_stop(default_options, 0)
    @test !Norma.is_matrix_output_stop(default_options, 1)
    all_options = parse_block(Norma.Parameters("steps" => "all"))
    @test Norma.is_matrix_output_stop(all_options, 7)
    every_options = parse_block(Norma.Parameters("steps" => Norma.Parameters("every" => 5)))
    @test Norma.is_matrix_output_stop(every_options, 0)
    @test Norma.is_matrix_output_stop(every_options, 10)
    @test !Norma.is_matrix_output_stop(every_options, 7)
    list_options = parse_block(Norma.Parameters("matrices" => ["mass"], "steps" => [3, 8]))
    @test list_options.matrices == ["mass"]
    @test Norma.is_matrix_output_stop(list_options, 8)
    @test !Norma.is_matrix_output_stop(list_options, 4)
    @test_throws Norma.NormaAbortException parse_block(Norma.Parameters("matrices" => ["damping"]))
    @test_throws Norma.NormaAbortException parse_block(Norma.Parameters("steps" => "sometimes"))
    @test_throws Norma.NormaAbortException parse_block(Norma.Parameters("steps" => Norma.Parameters("every" => 0)))
    @test_throws Norma.NormaAbortException parse_block(Norma.Parameters("steps" => [-1]))
    @test_throws Norma.NormaAbortException parse_block("stiffness")

    # Schwarz (multi-domain) inputs are rejected.
    multi_params = Norma.Parameters("type" => "multi", "matrix output" => Norma.Parameters())
    @test_throws Norma.NormaAbortException Norma.reject_matrix_output_in_schwarz(multi_params, "multi.yaml")
    @test Norma.reject_matrix_output_in_schwarz(Norma.Parameters("type" => "single"), "single.yaml") === nothing
    Norma.NORMA_TEST_MODE[] = false
end
