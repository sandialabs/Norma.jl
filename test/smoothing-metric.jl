# Anisotropic energetic mesh smoothing with a metric tensor field (issue: EMS
# extension to metric tensors).  The ideal element is the unit regular
# tetrahedron mapped by F_M⁻¹, the energy is evaluated on the deformation
# gradient in the metric space, F_U = F_M F F_M⁻¹, and stress and tangent are
# pulled back to the reference configuration.
using LinearAlgebra
using Random
using Test
using StaticArrays

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

tet_edges(X) = [norm(X[:, i] - X[:, j]) for i in 1:4 for j in (i + 1):4]
tet_volume(X) = dot(X[:, 2] - X[:, 1], cross(X[:, 3] - X[:, 1], X[:, 4] - X[:, 1])) / 6.0
# Deformation gradient of a linear tetrahedron from reference X to current x.
tet_gradient(X, x) = SMatrix{3,3,Float64,9}((x[:, 2:4] .- x[:, 1]) / (X[:, 2:4] .- X[:, 1]))

a = 1.0
t3 = sqrt(3) / 2 * a
h3 = sqrt(2.0 / 3.0) * a
reg_tet_coords = [
    -t3/3 -t3/3 2*t3/3 0.0
    a/2 -a/2 0.0 0.0
    0.0 0.0 0.0 h3
]

@testset "metric_reference" begin
    # Principal sizes along the global axes: the ideal element is the unit
    # regular tetrahedron scaled by h_i, and F_M maps it back to unit edges.
    mf = Norma.create_metric_field("metric field unrestricted", Dict("sizes" => ["0.2", "0.5", "1.0"]))
    X, F_M, F_M_inv = Norma.create_metric_reference(mf, reg_tet_coords, 0.0)
    @test F_M * F_M_inv ≈ I atol = 1.0e-14
    @test F_M' * F_M ≈ Diagonal([1 / 0.2^2, 1 / 0.5^2, 1.0]) atol = 1.0e-12
    for e in tet_edges(F_M * X)
        @test e ≈ 1.0 atol = 1.0e-12
    end
    @test tet_volume(X) ≈ 0.2 * 0.5 * 1.0 / (6.0 * sqrt(2.0)) atol = 1.0e-14

    # A rotation vector carries the global axes onto the principal directions:
    # M = R diag(1/h_i^2) R' and the ideal element is R diag(h_i) Y.
    mf = Norma.create_metric_field(
        "metric field unrestricted", Dict("sizes" => ["0.2", "0.5", "1.0"], "rotation vector" => ["0.3", "-0.2", "0.5"])
    )
    Xr, F_Mr, F_Mr_inv = Norma.create_metric_reference(mf, reg_tet_coords, 0.0)
    R = Norma.rt_of_rv(SVector(0.3, -0.2, 0.5))
    @test F_Mr' * F_Mr ≈ R * Diagonal([1 / 0.2^2, 1 / 0.5^2, 1.0]) * R' atol = 1.0e-12
    @test Xr ≈ R * X atol = 1.0e-12
    for e in tet_edges(F_Mr * Xr)
        @test e ≈ 1.0 atol = 1.0e-12
    end

    # Isotropic sizes reproduce the isotropic size-field reference.
    mf = Norma.create_metric_field("metric field unrestricted", Dict("sizes" => ["0.7", "0.7", "0.7"]))
    Xi, F_Mi, _ = Norma.create_metric_reference(mf, reg_tet_coords, 0.0)
    sf = Norma.create_size_field("size field unrestricted", "0.7")
    Xs = Norma.create_smooth_reference("size field unrestricted", Norma.TETRA4, reg_tet_coords, sf, 0.0)
    @test sort(tet_edges(Xi)) ≈ sort(tet_edges(Xs)) atol = 1.0e-12
    @test F_Mi ≈ I / 0.7 atol = 1.0e-12

    # The restricted rule never asks for a volume below that of the original
    # element, and scales the sizes uniformly so the shape is unchanged.
    tiny = Norma.create_metric_field("metric field", Dict("sizes" => ["1.0e-3", "2.0e-3", "4.0e-3"]))
    Xt, F_Mt, _ = Norma.create_metric_reference(tiny, reg_tet_coords, 0.0)
    @test tet_volume(Xt) ≈ tet_volume(reg_tet_coords) rtol = 1.0e-12
    d = diag(F_Mt' * F_Mt)
    @test d[1] / d[2] ≈ 4.0 rtol = 1.0e-12
    @test d[2] / d[3] ≈ 4.0 rtol = 1.0e-12
    tiny_u = Norma.create_metric_field("metric field unrestricted", Dict("sizes" => ["1.0e-3", "2.0e-3", "4.0e-3"]))
    Xu, _, _ = Norma.create_metric_reference(tiny_u, reg_tet_coords, 0.0)
    @test tet_volume(Xu) ≈ 8.0e-9 / (6.0 * sqrt(2.0)) rtol = 1.0e-12

    # A spatially and temporally varying field is evaluated at the reference
    # centroid and the supplied time.
    mf = Norma.create_metric_field("metric field unrestricted", Dict("sizes" => ["2.0 * x", "1.0 + t", "1.0"]))
    coords = reg_tet_coords .+ [1.5; 0.0; 0.0]
    cx = sum(coords[1, :]) / 4
    Xv, F_Mv, _ = Norma.create_metric_reference(mf, coords, 3.0)
    @test diag(F_Mv' * F_Mv) ≈ [1 / (2cx)^2, 1 / 16.0, 1.0] rtol = 1.0e-12

    # Other modes compile no metric field; malformed input aborts.
    @test Norma.create_metric_field("max", nothing) === nothing
    @test Norma.create_metric_field("size field", Dict("sizes" => ["1", "1", "1"])) === nothing
    @test_throws Exception Norma.create_metric_field("metric field", nothing)
    @test_throws Exception Norma.create_metric_field("metric field", Dict("sizes" => ["1", "1"]))
    @test_throws Exception Norma.create_metric_field(
        "metric field", Dict("sizes" => ["1", "1", "1"], "rotation vector" => ["0"])
    )
    bad = Norma.create_metric_field("metric field", Dict("sizes" => ["x - 100.0", "1", "1"]))
    @test_throws Exception Norma.create_metric_reference(bad, reg_tet_coords, 0.0)
    @test_throws Exception Norma.create_smooth_reference("metric field", Norma.TETRA4, reg_tet_coords, nothing, 0.0)
end

@testset "metric_energy_orientation" begin
    material = Norma.SethHill(Dict{String,Any}("bulk modulus" => 1.0, "shear modulus" => 1.0, "m" => 2, "n" => 2))
    θ = 0.4
    mf = Norma.create_metric_field(
        "metric field unrestricted", Dict("sizes" => ["0.25", "1.0", "1.0"], "rotation vector" => ["0", "0", "$θ"])
    )
    X, F_M, F_M_inv = Norma.create_metric_reference(mf, reg_tet_coords, 0.0)
    energy(x) = Norma.strain_energy(material, F_M * tet_gradient(X, x) * F_M_inv)
    # The ideal element, in its prescribed orientation, has zero energy; the
    # same shape rotated by a quarter turn about the axis of the metric does
    # not: an isotropic energy on F alone could not tell the two apart.
    @test energy(X) ≈ 0.0 atol = 1.0e-14
    Q = Norma.rt_of_rv(SVector(0.0, 0.0, π / 2))
    @test energy(Q * X) > 1.0
    @test Norma.strain_energy(material, tet_gradient(X, Q * X)) ≈ 0.0 atol = 1.0e-14
    # A rotation about the axis of transverse isotropy of the metric (the
    # principal direction of the small size) leaves the energy at zero.
    Qa =
        Norma.rt_of_rv(SVector(0.0, 0.0, θ)) * Norma.rt_of_rv(SVector(0.9, 0.0, 0.0)) *
        Norma.rt_of_rv(SVector(0.0, 0.0, -θ))
    @test energy(Qa * X) ≈ 0.0 atol = 1.0e-12
    # The energy does not depend on the factor chosen for M: Q F_M gives the
    # same value for any rotation Q.
    x = reg_tet_coords .+ 0.1 * randn(3, 4)
    F = tet_gradient(X, x)
    Qf = Norma.rt_of_rv(SVector(0.7, -0.4, 1.1))
    @test Norma.strain_energy(material, (Qf * F_M) * F * inv(Qf * F_M)) ≈ energy(x) rtol = 1.0e-12
end

# Assembled energy, force, and stiffness on a small mesh under a metric field.
metric_params = Dict{String,Any}(
    "type" => "single",
    "name" => "metric_gradient",
    "input mesh file" => "../examples/ems/cube/cube.g",
    "output mesh file" => "metric_gradient.e",
    "Exodus output interval" => 0,
    "CSV output interval" => 0,
    "model" => Dict{String,Any}(
        "type" => "mesh smoothing",
        "smooth reference" => "metric field unrestricted",
        "metric field" => Dict{String,Any}(
            "sizes" => ["0.09 + 0.02 * x", "0.11", "0.1 + 0.02 * z"],
            "rotation vector" => ["0.3 * y", "0.1", "0.4 + 0.2 * x"],
        ),
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
    "solver" => Dict{String,Any}(
        "type" => "Hessian minimizer",
        "step" => "full Newton",
        "minimum iterations" => 1,
        "maximum iterations" => 1,
        "relative tolerance" => 1.0e-12,
        "absolute tolerance" => 1.0e-8,
    ),
)

@testset "metric_gradient_and_tangent" begin
    sim = Norma.create_simulation(metric_params)
    model = sim.model
    num_dofs = length(model.internal_force)
    # Perturb the nodes by a few percent of the mean edge length.
    conn = model.blocks[1].connectivity
    mean_edge =
        sum(norm(model.reference[:, conn[1, e]] - model.reference[:, conn[2, e]]) for e in 1:size(conn, 2)) /
        size(conn, 2)
    model.displacement .= 0.03 * mean_edge * randn(3, size(model.reference, 2))
    u0 = copy(model.displacement)
    function energy_force(u)
        model.displacement .= u
        model.compute_stiffness = true
        Norma.evaluate(model, sim.integrator, sim.solver)
        @test model.failed == false
        return model.strain_energy, copy(model.internal_force), copy(model.stiffness)
    end
    E0, f0, K0 = energy_force(u0)
    @test E0 > 0.0
    # The step balances the truncation error of the central difference against
    # the cancellation in E(u + ε) - E(u - ε), which is of the order of the
    # roundoff of the total energy; the tolerance is set by that cancellation,
    # measured against the largest force so that small entries do not demand
    # more than the difference can deliver.
    ε = 1.0e-4 * mean_edge
    f_scale = maximum(abs, f0)
    for _ in 1:12
        i = rand(1:num_dofs)
        e = zeros(size(u0))
        e[i] = ε
        Ep, _, _ = energy_force(u0 + e)
        Em, _, _ = energy_force(u0 - e)
        @test (Ep - Em) / (2ε) ≈ f0[i] atol = 1.0e-6 * max(abs(f0[i]), 1.0e-2 * f_scale)
    end
    δ = randn(size(u0))
    δ ./= norm(δ)
    _, fp, _ = energy_force(u0 + ε * δ)
    _, fm, _ = energy_force(u0 - ε * δ)
    @test (fp - fm) / (2ε) ≈ K0 * vec(δ) atol = 1.0e-5 * norm(K0 * vec(δ))
    Norma.finalize_writing(sim)
end

# A full smoothing run with a constant anisotropic metric at 45 degrees: the
# thought experiment of the notes.  The energy decreases and the edges of the
# mesh measured in the metric move toward unit length.
cube_params = Dict{String,Any}(
    "type" => "single",
    "name" => "metric_cube",
    "input mesh file" => "../examples/ems/cube/cube.g",
    "output mesh file" => "metric_cube.e",
    "Exodus output interval" => 0,
    "CSV output interval" => 0,
    "model" => Dict{String,Any}(
        "type" => "mesh smoothing",
        "smooth reference" => "metric field unrestricted",
        "metric field" => Dict{String,Any}(
            "sizes" => ["1.0", "1.0", "1.0"], "rotation vector" => ["0", "0", "0.7853981633974483"]
        ),
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
        "maximum iterations" => 40,
        "relative tolerance" => 1.0e-12,
        "absolute tolerance" => 1.0e-10,
        "step length" => 1.0e-3,
        "use line search" => true,
        "line search backtrack factor" => 0.5,
        "line search decrease factor" => 1.0e-04,
        "line search maximum iterations" => 16,
    ),
)

@testset "metric_smoothing_run" begin
    # Mean edge length of the cube mesh sets the scale H of the metric.
    probe_params = deepcopy(cube_params)
    probe_params["output mesh file"] = "metric_probe.e"
    probe = Norma.create_simulation(probe_params)
    conn = probe.model.blocks[1].connectivity
    X0 = copy(probe.model.reference)
    Norma.finalize_writing(probe)
    edges = Set{Tuple{Int,Int}}()
    for e in 1:size(conn, 2), i in 1:4, j in (i + 1):4
        push!(edges, minmax(conn[i, e], conn[j, e]))
    end
    H = sum(norm(X0[:, i] - X0[:, j]) for (i, j) in edges) / length(edges)
    params = deepcopy(cube_params)
    params["model"]["metric field"]["sizes"] = ["$(0.5 * H)", "$(1.5 * H)", "$(1.5 * H)"]
    R = Norma.rt_of_rv(SVector(0.0, 0.0, π / 4))
    F_M = Diagonal([1 / (0.5H), 1 / (1.5H), 1 / (1.5H)]) * R'
    metric_misfit(X) = sum(abs(norm(F_M * (X[:, i] - X[:, j])) - 1.0) for (i, j) in edges) / length(edges)
    misfit0 = metric_misfit(X0)
    sim = Norma.create_simulation(params)
    Norma.evaluate(sim.model, sim.integrator, sim.solver)
    E0 = sim.model.strain_energy
    Norma.run(sim)
    X1 = sim.model.reference + sim.model.displacement
    E1 = sim.model.strain_energy
    @test E1 < E0
    @test metric_misfit(X1) < misfit0
    # The volume of each element stays positive.
    for e in 1:size(conn, 2)
        @test tet_volume(X1[:, conn[:, e]]) > 0.0
    end
end
