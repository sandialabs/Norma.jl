# Metric targets of energetic mesh smoothing given as tensors or carried by
# the nodes (docs/notes/ems-anisotropic): the factorization of a tensor into
# sizes and a rotation with a frame that follows the reference, the frames
# followed through a mesh, the equivalence of the four sources, and the nodal
# data through an adaptivity loop.
using LinearAlgebra
using Random
using Test
using Exodus
using StaticArrays

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

identity_3 = SMatrix{3,3,Float64,9}(I)
rotation_angle(R) = acos(clamp((tr(R) - 1.0) / 2.0, -1.0, 1.0))

@testset "principal_of_tensor" begin
    # Distinct sizes: the factorization returns the sizes in the order of the
    # frame closest to the reference, and the tensor is reproduced.
    h = SVector(0.2, 0.5, 1.0)
    R = Norma.rt_of_rv(SVector(0.3, -0.2, 0.4))
    M = Norma.metric_tensor(h, R)
    hr, Rr = Norma.principal_of_tensor(M)
    @test hr ≈ h rtol = 1.0e-12
    @test Rr ≈ R atol = 1.0e-12
    @test Norma.metric_tensor(hr, Rr) ≈ M rtol = 1.0e-12
    # A large rotation is followed when the reference is close to it.
    R2 = Norma.rt_of_rv(SVector(0.0, 0.0, 2.0))
    M2 = Norma.metric_tensor(h, R2)
    hr2, Rr2 = Norma.principal_of_tensor(M2, Norma.rt_of_rv(SVector(0.0, 0.0, 1.9)))
    @test hr2 ≈ h rtol = 1.0e-12
    @test Rr2 ≈ R2 atol = 1.0e-12
    # With the identity as reference the same tensor takes an equivalent
    # frame closer to the axes: a different rotation, the same tensor.
    hr3, Rr3 = Norma.principal_of_tensor(M2)
    @test rotation_angle(Rr3) < rotation_angle(R2)
    @test Norma.metric_tensor(hr3, Rr3) ≈ M2 rtol = 1.0e-12
    # Repeated sizes: the axis is the odd principal direction, the in-plane
    # phase comes from the reference, and the tensor is reproduced.
    ht = SVector(0.25, 1.0, 1.0)
    Rt = Norma.rt_of_rv(SVector(0.0, 0.0, 0.4))
    Mt = Norma.metric_tensor(ht, Rt)
    hrt, Rrt = Norma.principal_of_tensor(Mt)
    @test hrt ≈ ht rtol = 1.0e-10
    @test Rrt ≈ Rt atol = 1.0e-8
    @test Norma.metric_tensor(hrt, Rrt) ≈ Mt rtol = 1.0e-10
    # Isotropic: the reference frame is kept.
    hi, Ri = Norma.principal_of_tensor(Norma.metric_tensor(SVector(0.3, 0.3, 0.3), R), R)
    @test hi ≈ SVector(0.3, 0.3, 0.3) rtol = 1.0e-12
    @test Ri == R
    # Random tensors and references: the tensor is always reproduced by a
    # proper rotation.
    for _ in 1:20
        hx = exp.(randn(3))
        Rx = Norma.rt_of_rv(SVector{3,Float64}(randn(3)))
        Mx = Norma.metric_tensor(SVector{3,Float64}(hx), Rx)
        ref = Norma.rt_of_rv(SVector{3,Float64}(randn(3)))
        hy, Ry = Norma.principal_of_tensor(Mx, ref)
        @test det(Ry) ≈ 1.0 rtol = 1.0e-12
        @test Ry' * Ry ≈ identity_3 atol = 1.0e-12
        @test Norma.metric_tensor(hy, Ry) ≈ Mx rtol = 1.0e-10
    end
    # Not positive definite.
    @test Norma.principal_of_tensor(SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0)) === nothing
    @test Norma.log_symmetric(SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)) === nothing
    L = Norma.log_symmetric(M)
    @test Norma.exp_symmetric(L) ≈ M rtol = 1.0e-12
    c = Norma.components_from_symmetric(M)
    @test Norma.symmetric_from_components(c) ≈ M rtol = 1.0e-15
end

@testset "principal_of_nodal_tensors" begin
    # A chain of nodes whose frame turns by 0.3 rad per node: every frame is
    # followed from its neighbor, so the rotation vectors are continuous
    # although the total turn exceeds a half turn, and no jump is reported.
    n = 12
    h = SVector(0.2, 0.5, 1.0)
    tensors = [Norma.metric_tensor(h, Norma.rt_of_rv(SVector(0.0, 0.0, 0.3 * (k - 1)))) for k in 1:n]
    neighbors = [Int[] for _ in 1:n]
    for k in 1:(n - 1)
        push!(neighbors[k], k + 1)
        push!(neighbors[k + 1], k)
    end
    sizes, rotations, jumps = Norma.principal_of_nodal_tensors(tensors, neighbors)
    @test jumps == 0
    @test all(sizes[:, k] ≈ h for k in 1:n)
    @test all(norm(rotations[:, k + 1] - rotations[:, k] - [0.0, 0.0, 0.3]) < 1.0e-8 for k in 1:(n - 1))
    # The frames are equivalent up to the 24 signed permutations of the
    # axes, so a turn between two nodes counts as a jump only when no
    # equivalent frame is within the threshold.  A generic turn beyond it is
    # found by search and reported once (the node is reached from one side).
    signed_permutations = [
        SMatrix{3,3,Float64,9}(hcat(s1 * identity_3[:, p[1]], s2 * identity_3[:, p[2]], s3 * identity_3[:, p[3]])) for
        p in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)) for s1 in (1, -1) for s2 in (1, -1) for
        s3 in (1, -1)
    ]
    proper = [P for P in signed_permutations if det(P) > 0.0]
    @test length(proper) == 24
    minimal_angle(Q) = minimum(rotation_angle(Q * P) for P in proper)
    turn = SVector(0.0, 0.0, 0.0)
    for _ in 1:200
        v = SVector{3,Float64}(randn(3))
        if minimal_angle(Norma.rt_of_rv(v)) > π / 4 + 0.05
            turn = v
            break
        end
    end
    @test minimal_angle(Norma.rt_of_rv(turn)) > π / 4
    tensors[6] = Norma.metric_tensor(h, Norma.rt_of_rv(SVector(0.0, 0.0, 1.5)) * Norma.rt_of_rv(turn))
    sizes, rotations, jumps = Norma.principal_of_nodal_tensors(tensors, neighbors)
    # The turned node is reached from one neighbor and reaches the other.
    @test jumps == 2
    # Every node still reproduces its tensor.
    @test all(
        Norma.metric_tensor(SVector{3,Float64}(sizes[:, k]), Norma.rt_of_rv(SVector{3,Float64}(rotations[:, k]))) ≈
        tensors[k] for k in 1:n
    )
end

# Smoothing model on the cube with a metric field block.  Every simulation
# gets its own output file, since creating one opens the file.
const nodal_output_counter = Ref(0)
function metric_model(mesh_file, metric_field; output=nothing, interval=0)
    if output === nothing
        nodal_output_counter[] += 1
        output = "metric-nodal-$(nodal_output_counter[]).e"
    end
    return Dict{String,Any}(
        "type" => "single",
        "name" => "metric_nodal",
        "input mesh file" => mesh_file,
        "output mesh file" => output,
        "Exodus output interval" => interval,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing",
            "smooth reference" => "metric field unrestricted",
            "metric field" => metric_field,
            "material" => Dict{String,Any}(
                "elastic" => Dict{String,Any}(
                    "model" => "seth-hill",
                    "m" => 2,
                    "n" => 2,
                    "bulk modulus" => 1.0,
                    "shear modulus" => 1.0,
                    "density" => 1.0,
                ),
                "blocks" => Dict{String,Any}("cube" => "elastic"),
            ),
        ),
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "boundary conditions" => Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => ns, "component" => c, "function" => "0.0") for
                (ns, c) in (("nsx-", "x"), ("nsx+", "x"), ("nsy-", "y"), ("nsy+", "y"), ("nsz-", "z"), ("nsz+", "z"))
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent",
            "step" => "lbfgs",
            "memory" => 10,
            "minimum iterations" => 1,
            "maximum iterations" => 4,
            "relative tolerance" => 1.0e-12,
            "absolute tolerance" => 1.0e-10,
            "step length" => 1.0e-3,
            "use line search" => true,
            "line search backtrack factor" => 0.5,
            "line search decrease factor" => 1.0e-04,
            "line search maximum iterations" => 16,
        ),
    )
end

function assembled_energy(params, displacement)
    params = copy(params)
    nodal_output_counter[] += 1
    params["output mesh file"] = "metric-nodal-$(nodal_output_counter[]).e"
    sim = Norma.create_simulation(params)
    model = sim.model
    model.displacement .= displacement
    Norma.evaluate(model, sim.integrator, sim.solver)
    @test model.failed == false
    return model.strain_energy, copy(model.internal_force), sim
end

# The metric of the cube examples: sizes (0.05, 0.1414, 0.1414) rotated by a
# quarter of a half turn about z, and its tensor.
const nodal_cube_sizes = ["0.05", "0.1414", "0.1414"]
const nodal_cube_rotation = ["0.0", "0.0", "0.7853981633974483"]
const nodal_cube_tensor = let
    M = Norma.metric_tensor(SVector(0.05, 0.1414, 0.1414), Norma.rt_of_rv(SVector(0.0, 0.0, π / 4)))
    [string(c) for c in Norma.components_from_symmetric(M)]
end

@testset "tensor_functions" begin
    # The tensor given by its components and by its sizes and rotation vector
    # give the same energy and force: the factor chosen for the tensor does
    # not matter.
    principal = metric_model(
        "../examples/ems/cube/cube.g",
        Dict{String,Any}("sizes" => nodal_cube_sizes, "rotation vector" => nodal_cube_rotation),
    )
    tensor = metric_model("../examples/ems/cube/cube.g", Dict{String,Any}("tensor" => nodal_cube_tensor))
    sim0 = Norma.create_simulation(principal)
    u = 0.003 * randn(size(sim0.model.reference))
    Ep, fp, _ = assembled_energy(principal, u)
    Et, ft, sim_t = assembled_energy(tensor, u)
    @test Et ≈ Ep rtol = 1.0e-10
    @test ft ≈ fp rtol = 1.0e-8
    @test sim_t.model.metric_field.source isa Norma.TensorMetricFunctions
    # The output recovers sizes and rotation vectors from the tensor.
    positions = copy(sim_t.model.reference)
    sizes, rotation = Norma.nodal_metric_output(sim_t.model, positions, 0.0)
    @test size(sizes) == size(positions) && rotation !== nothing
    @test all(sort(sizes[:, n]) ≈ [0.05, 0.1414, 0.1414] for n in 1:size(sizes, 2))
    @test all(
        Norma.metric_tensor(SVector{3,Float64}(sizes[:, n]), Norma.rt_of_rv(SVector{3,Float64}(rotation[:, n]))) ≈
        Norma.metric_tensor(SVector(0.05, 0.1414, 0.1414), Norma.rt_of_rv(SVector(0.0, 0.0, π / 4))) for
        n in 1:size(sizes, 2)
    )
    # Input errors: two forms, a rotation vector without sizes, a bad count.
    @test_throws Exception Norma.create_metric_field(
        "metric field", Dict("sizes" => nodal_cube_sizes, "tensor" => nodal_cube_tensor)
    )
    @test_throws Exception Norma.create_metric_field(
        "metric field", Dict("tensor" => nodal_cube_tensor, "rotation vector" => nodal_cube_rotation)
    )
    @test_throws Exception Norma.create_metric_field("metric field", Dict("tensor" => nodal_cube_tensor[1:5]))
    @test_throws Exception Norma.create_metric_field("metric field", Dict("nodal sizes" => ["a", "b", "c"]))
    indefinite = Norma.create_metric_field("metric field", Dict("tensor" => ["1", "-1", "1", "0", "0", "0"]))
    @test_throws Exception Norma.create_metric_reference(indefinite, sim0.model.reference[:, 1:4], 0.0)
end

# A copy of the cube mesh that carries nodal variables, written through the
# topology writer that the adaptivity loop uses.
function mesh_with_variables(model, file_name, variables)
    topology = Norma.build_topology(model)
    Norma.write_topology(topology, file_name; nodal_variables=variables)
    return file_name
end

@testset "nodal_sources" begin
    principal = metric_model(
        "../examples/ems/cube/cube.g",
        Dict{String,Any}("sizes" => nodal_cube_sizes, "rotation vector" => nodal_cube_rotation),
    )
    sim0 = Norma.create_simulation(principal)
    model0 = sim0.model
    num_nodes = size(model0.reference, 2)
    u = 0.003 * randn(size(model0.reference))
    Ep, fp, _ = assembled_energy(principal, u)
    # Constant sizes and rotation at the nodes: the mean over an element is
    # the constant, and the energy equals the analytic one exactly.
    variables = Dict{String,Vector{Float64}}(
        "h1" => fill(0.05, num_nodes),
        "h2" => fill(0.1414, num_nodes),
        "h3" => fill(0.1414, num_nodes),
        "r1" => zeros(num_nodes),
        "r2" => zeros(num_nodes),
        "r3" => fill(π / 4, num_nodes),
    )
    M = Norma.metric_tensor(SVector(0.05, 0.1414, 0.1414), Norma.rt_of_rv(SVector(0.0, 0.0, π / 4)))
    for (i, suffix) in enumerate(Norma.METRIC_COMPONENT_SUFFIXES)
        variables["m_$suffix"] = fill(Norma.components_from_symmetric(M)[i], num_nodes)
    end
    mesh_file = mesh_with_variables(model0, "metric-nodal-cube.g", variables)
    exo = ExodusDatabase(mesh_file, "r")
    @test Exodus.read_number_of_time_steps(exo) == 1
    @test Set(Exodus.read_names(exo, NodalVariable)) == Set(keys(variables))
    Exodus.close(exo)
    nodal_sizes = Dict{String,Any}(
        "nodal sizes" => ["h1", "h2", "h3"], "nodal rotation vector" => ["r1", "r2", "r3"], "time index" => 1
    )
    En, fn, sim_n = assembled_energy(metric_model(mesh_file, nodal_sizes), u)
    @test En ≈ Ep rtol = 1.0e-10
    @test fn ≈ fp rtol = 1.0e-8
    @test sim_n.model.metric_field.source isa Norma.NodalPrincipalMetric
    tensor_names = ["m_$s" for s in Norma.METRIC_COMPONENT_SUFFIXES]
    Et, ft, sim_t = assembled_energy(metric_model(mesh_file, Dict{String,Any}("nodal tensor" => tensor_names)), u)
    @test Et ≈ Ep rtol = 1.0e-10
    @test ft ≈ fp rtol = 1.0e-8
    @test sim_t.model.metric_field.source isa Norma.NodalPrincipalMetric
    @test sim_t.model.metric_field.source.tensor_names == tensor_names
    El, fl, sim_l = assembled_energy(
        metric_model(mesh_file, Dict{String,Any}("nodal tensor" => tensor_names, "interpolation" => "log-Euclidean")), u
    )
    @test El ≈ Ep rtol = 1.0e-10
    @test fl ≈ fp rtol = 1.0e-8
    @test sim_l.model.metric_field.source isa Norma.NodalTensorMetric
    # The sizes without a rotation are the sizes along the axes.
    axes_only = Dict{String,Any}("nodal sizes" => ["h1", "h2", "h3"])
    Ea, _, sim_a = assembled_energy(metric_model(mesh_file, axes_only), u)
    axes_params = metric_model("../examples/ems/cube/cube.g", Dict{String,Any}("sizes" => nodal_cube_sizes))
    Eax, _, _ = assembled_energy(axes_params, u)
    @test Ea ≈ Eax rtol = 1.0e-10
    @test isempty(sim_a.model.metric_field.source.rotation_names)
    # The nodal variables written for an adapted mesh reproduce the input.
    written = Norma.metric_nodal_variables(sim_n.model.metric_field)
    @test Set(keys(written)) == Set(["h1", "h2", "h3", "r1", "r2", "r3"])
    @test written["h1"] ≈ variables["h1"] && written["r3"] ≈ variables["r3"]
    written_t = Norma.metric_nodal_variables(sim_t.model.metric_field)
    @test Set(keys(written_t)) == Set(tensor_names)
    @test all(written_t[name] ≈ variables[name] for name in tensor_names)
    written_l = Norma.metric_nodal_variables(sim_l.model.metric_field)
    @test all(written_l[name] ≈ variables[name] for name in tensor_names)
    # A graded field carried by the nodes agrees with the analytic field to
    # the second order in the element size.
    graded = Dict{String,Any}(
        "sizes" => ["0.08 + 0.04 * x", "0.1414", "0.1414"], "rotation vector" => ["0.0", "0.0", "0.5 * y"]
    )
    Eg, _, sim_g = assembled_energy(metric_model("../examples/ems/cube/cube.g", graded), u)
    positions = copy(sim_g.model.reference)
    sizes, rotation = Norma.nodal_metric_output(sim_g.model, positions, 0.0)
    graded_variables = Dict{String,Vector{Float64}}(
        "h1" => sizes[1, :], "h2" => sizes[2, :], "h3" => sizes[3, :],
        "r1" => rotation[1, :], "r2" => rotation[2, :], "r3" => rotation[3, :],
    )
    graded_file = mesh_with_variables(model0, "metric-nodal-graded.g", graded_variables)
    Egn, _, _ = assembled_energy(metric_model(graded_file, nodal_sizes), u)
    @test Egn ≈ Eg rtol = 2.0e-2
    @test !isapprox(Egn, Eg; rtol=1.0e-8)
    # Errors: a missing variable, a bad time index, a nonpositive size.
    missing_variable = Dict{String,Any}("nodal sizes" => ["h1", "h2", "no"])
    @test_throws Exception Norma.create_simulation(metric_model(mesh_file, missing_variable))
    @test_throws Exception Norma.create_simulation(metric_model(mesh_file, merge(nodal_sizes, Dict("time index" => 2))))
    bad = copy(variables)
    bad["h1"] = zeros(num_nodes)
    bad_file = mesh_with_variables(model0, "metric-nodal-bad.g", bad)
    @test_throws Exception Norma.create_simulation(metric_model(bad_file, axes_only))
    rm(bad_file; force=true)
    rm(graded_file; force=true)
end

@testset "nodal_metric_adaptivity" begin
    # A graded nodal metric finer than the mesh through one topology phase:
    # the split nodes take the mean of the ends, the compacted data matches
    # the adapted mesh, which carries the variables under the input names, and
    # the loop runs again on it.
    principal = metric_model("../examples/ems/cube/cube.g", Dict{String,Any}("sizes" => nodal_cube_sizes))
    model0 = Norma.create_simulation(principal).model
    num_nodes = size(model0.reference, 2)
    x = model0.reference[1, :]
    variables = Dict{String,Vector{Float64}}(
        "h1" => 0.06 .+ 0.02 .* x, "h2" => 0.07 .+ 0.02 .* x, "h3" => 0.065 .+ 0.02 .* x
    )
    mesh_file = mesh_with_variables(model0, "metric-nodal-adapt.g", variables)
    nodal = Dict{String,Any}("nodal sizes" => ["h1", "h2", "h3"])
    params = metric_model(mesh_file, nodal; output="metric-nodal-adapt.e")
    params["adaptivity"] = Dict{String,Any}(
        "desired energy density" => 0.05,
        "maximum passes" => 1,
        "outer iterations" => 2,
        "adjacency layers" => 1,
        "swaps" => false,
        "collapses" => false,
    )
    sim = Norma.run(params)
    model = sim.model
    source = model.metric_field.source
    @test source isa Norma.NodalPrincipalMetric
    # The final model was created from the mesh of the second iteration.
    @test isfile("metric-nodal-adapt-adapted-1.g")
    adapted = "metric-nodal-adapt-adapted-2.g"
    @test isfile(adapted)
    exo = ExodusDatabase(adapted, "r")
    names = Exodus.read_names(exo, NodalVariable)
    @test Set(names) == Set(["h1", "h2", "h3"])
    n_adapted = Exodus.num_nodes(exo.init)
    h1 = Exodus.read_values(exo, NodalVariable, 1, "h1")
    Exodus.close(exo)
    @test n_adapted > num_nodes
    @test size(source.log_sizes, 2) == n_adapted
    @test exp.(source.log_sizes[1, :]) ≈ h1
    # Every carried size lies in the range of the input, as a mean of ends
    # must, and the original nodes keep their values.
    @test minimum(h1) ≥ minimum(variables["h1"]) - 1.0e-12
    @test maximum(h1) ≤ maximum(variables["h1"]) + 1.0e-12
    @test all(iszero, source.rotation)
    for file in (mesh_file, adapted, "metric-nodal-adapt-adapted-1.g", "metric-nodal-adapt.e")
        rm(file; force=true)
    end
    rm("metric-nodal-adapt-adapted-1.e"; force=true)
    rm("metric-nodal-adapt-adapted-2.e"; force=true)
    rm("metric-nodal-cube.g"; force=true)
    for k in 1:nodal_output_counter[]
        rm("metric-nodal-$k.e"; force=true)
    end
end
