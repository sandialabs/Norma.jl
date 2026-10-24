# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The composite tet 10: its shape functions and projected gradients against
# the sub-tetrahedra built from the geometry alone and against the tables of
# Intrepid2's unit test, its mass against LCM's analytic one, the patch test,
# convergence on a thick-walled cylinder and its behavior as the material
# becomes incompressible. The comparisons with LCM itself are in
# composite-tet10-lcm.jl and lcm-rubiks-cube.jl.
using Exodus
using LinearAlgebra
using Random
using StaticArrays

# The 12 sub-tetrahedra of the composite tet, built from the geometry alone:
# the 4 corner tets of the regular refinement and the cones from the
# centroid over the 8 faces of the inner octahedron. Each is the list of its
# 4 vertices, as (reference coordinates, nodal values of the 10 shape
# functions there).
function composite_tet10_subtets()
    vertices = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)]
    edges = [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)]
    e(i) = SVector{10,Float64}(j == i ? 1.0 : 0.0 for j in 1:10)
    corner(k) = (vertices[k], e(k))
    midside(m) = ((vertices[edges[m][1]] + vertices[edges[m][2]]) / 2, e(4 + m))
    center = (SVector(0.25, 0.25, 0.25), SVector{10,Float64}(j > 4 ? 1 / 6 : 0.0 for j in 1:10))
    touching(k) = [m for m in 1:6 if k in edges[m]]
    subtets = [[corner(k); midside.(touching(k))] for k in 1:4]
    for k in 1:4
        push!(subtets, [center; midside.(touching(k))])
        push!(subtets, [center; midside.(setdiff(1:6, touching(k)))])
    end
    return subtets
end

# A structured tetrahedral mesh of the unit cube mapped by φ: n[1] × n[2] ×
# n[3] cells cut into the 6 Kuhn tetrahedra along their main diagonal. The
# edge nodes are mapped by φ too (curved elements for a curved φ), or put at
# the middle of the mapped corners with `straight`. TETRA4 without
# `quadratic`. Node sets xm, xp, ym, yp, zm, zp on the faces of the cube.
function write_structured_tet_mesh(path, n, φ; quadratic=true, straight=false)
    nx, ny, nz = n
    vertex(i, j, k) = 1 + i + (nx + 1) * (j + (ny + 1) * k)
    P = [SVector(i / nx, j / ny, k / nz) for k in 0:nz for j in 0:ny for i in 0:nx]
    X = [SVector{3,Float64}(φ(p)) for p in P]
    edge_nodes = Dict{Tuple{Int,Int},Int}()
    connectivity = Vector{Vector{Int}}()
    for k in 0:(nz - 1), j in 0:(ny - 1), i in 0:(nx - 1)
        for permutation in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
            c = [i, j, k]
            tet = [vertex(c...)]
            for d in permutation
                c[d] += 1
                push!(tet, vertex(c...))
            end
            if det(hcat(X[tet[2]] - X[tet[1]], X[tet[3]] - X[tet[1]], X[tet[4]] - X[tet[1]])) < 0
                tet = tet[[1, 3, 2, 4]]
            end
            if quadratic
                for (a, b) in ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))
                    key = minmax(tet[a], tet[b])
                    if !haskey(edge_nodes, key)
                        middle = (P[tet[a]] + P[tet[b]]) / 2
                        push!(P, middle)
                        push!(X, straight ? (X[tet[a]] + X[tet[b]]) / 2 : SVector{3,Float64}(φ(middle)))
                        edge_nodes[key] = length(P)
                    end
                    push!(tet, edge_nodes[key])
                end
            end
            push!(connectivity, tet)
        end
    end
    faces = (
        ("xm", p -> p[1] < 1.0e-12),
        ("xp", p -> p[1] > 1 - 1.0e-12),
        ("ym", p -> p[2] < 1.0e-12),
        ("yp", p -> p[2] > 1 - 1.0e-12),
        ("zm", p -> p[3] < 1.0e-12),
        ("zp", p -> p[3] > 1 - 1.0e-12),
    )
    init = Exodus.Initialization{Int32}(
        Int32(3), Int32(length(X)), Int32(length(connectivity)), Int32(1), Int32(length(faces)), Int32(0)
    )
    rm(path; force=true)
    exo = Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}(path, "w", init)
    Exodus.write_coordinates(exo, Matrix(reduce(hcat, X)))
    Exodus.write_block(exo, 1, quadratic ? "TETRA10" : "TETRA4", Matrix{Int32}(reduce(hcat, connectivity)))
    Exodus.write_name(exo, Block, 1, "solid")
    for (s, (name, on_face)) in enumerate(faces)
        node_set = Exodus.NodeSet(Int32(s), Int32.(findall(on_face, P)))
        Exodus.write_set(exo, node_set)
        Exodus.write_name(exo, node_set, name)
    end
    close(exo)
    return nothing
end

# A static linear elastic run on a structured mesh: the composite tet 10, the
# standard TETRA10 or, on the corner nodes, the TETRA4.
function run_structured_elastic(name, mesh, ν, dirichlet; composite=true)
    model = Dict{String,Any}(
        "type" => "solid mechanics",
        "material" => Dict{String,Any}(
            "blocks" => Dict("solid" => "elastic"),
            "elastic" => Dict{String,Any}(
                "model" => "linear elastic", "elastic modulus" => 1.0, "Poisson's ratio" => ν, "density" => 1.0
            ),
        ),
    )
    composite && (model["composite tet 10"] = ["solid"])
    params = Norma.Parameters(
        "name" => name * ".yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => name * ".e",
        "Exodus output interval" => 1.0,
        "CSV output interval" => 0,
        "model" => model,
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "boundary conditions" => Dict(
            "Dirichlet" =>
                [Dict{String,Any}("node set" => n, "component" => c, "function" => f) for (n, c, f) in dirichlet],
        ),
        "solver" => Dict{String,Any}(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => "direct",
            "minimum iterations" => 1,
            "maximum iterations" => 3,
            "absolute tolerance" => 1.0e-14,
            "relative tolerance" => 1.0e-12,
        ),
    )
    sim = Norma.run(params)
    for file in (mesh, name * ".e", name * ".log")
        rm(file; force=true)
    end
    return sim.model
end

# Patch test: the displacement u = A x prescribed on the boundary of a
# warped cube. Returns the largest deviations of the displacement and of the
# stress from the exact ones.
function composite_tet10_patch(n; straight)
    warp(p) = SVector(
        p[1] + 0.15 * sin(2 * p[2]) * p[1] * (1 - p[1]) + 0.1 * p[3]^2,
        p[2] + 0.12 * sin(3 * p[3]) * p[1],
        p[3] + 0.1 * p[1] * p[2],
    )
    write_structured_tet_mesh("composite-patch.g", (n, n, n), warp; straight=straight)
    A = [0.01 0.002 -0.003; 0.004 -0.005 0.001; 0.002 0.003 0.006]
    f = ["$(A[i, 1]) * x + $(A[i, 2]) * y + $(A[i, 3]) * z" for i in 1:3]
    dirichlet = [(s, c, f[i]) for s in ("xm", "xp", "ym", "yp", "zm", "zp") for (i, c) in enumerate(("x", "y", "z"))]
    model = run_structured_elastic("composite-patch", "composite-patch.g", 0.3, dirichlet)
    ν = 0.3
    λ, μ = ν / ((1 + ν) * (1 - 2ν)), 1 / (2 * (1 + ν))
    ε = (A + A') / 2
    σ = λ * tr(ε) * I + 2μ * ε
    voigt = [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
    u_error = maximum(abs.(model.displacement - A * model.reference))
    σ_error = maximum(maximum(abs.(s - voigt)) for element in model.stress[1] for s in element)
    return u_error, σ_error
end

# Quarter of a thick-walled cylinder, radii 1 and 2, in plane strain: radial
# displacement 0.01 on the inner surface, free outer surface. Returns the
# largest nodal error of the displacement relative to 0.01 (Lamé's solution
# u_r = A r + B / r).
function composite_tet10_lame(n, ν; kind=:composite, straight=false)
    a, b, u_a = 1.0, 2.0, 0.01
    height = (b - a) / n
    φ(p) = (r=a + (b - a) * p[1]; θ=π / 2 * p[2]; SVector(r * cos(θ), r * sin(θ), height * p[3]))
    write_structured_tet_mesh("composite-lame.g", (n, 2n, 1), φ; quadratic=kind != :tetra4, straight=straight)
    dirichlet = [
        ("xm", "x", "$u_a * x / sqrt(x^2 + y^2)"),
        ("xm", "y", "$u_a * y / sqrt(x^2 + y^2)"),
        ("ym", "y", "0.0"),
        ("yp", "x", "0.0"),
        ("zm", "z", "0.0"),
        ("zp", "z", "0.0"),
    ]
    model = run_structured_elastic("composite-lame", "composite-lame.g", ν, dirichlet; composite=kind == :composite)
    λ, μ = ν / ((1 + ν) * (1 - 2ν)), 1 / (2 * (1 + ν))
    A = u_a / (a + (λ + μ) * b^2 / (μ * a))
    B = (λ + μ) * A * b^2 / μ
    error = 0.0
    for node in axes(model.reference, 2)
        x = model.reference[:, node]
        r = hypot(x[1], x[2])
        error = max(error, norm(model.displacement[:, node] - (A * r + B / r) * [x[1] / r, x[2] / r, 0.0]))
    end
    return error / u_a
end

@testset "Composite Tet 10" begin
    subtets = composite_tet10_subtets()
    volumes = [abs(det(hcat((v[1] - s[1][1] for v in s[2:4])...))) / 6 for s in subtets]
    @test sum(volumes) ≈ 1 / 6

    @testset "Shape functions" begin
        # The shape functions are the linear interpolants of the vertex values
        # in every sub-tetrahedron, with constant gradients there.
        gradients = map(subtets) do s
            X = hcat((v[1] - s[1][1] for v in s[2:4])...)
            V = hcat((v[2] - s[1][2] for v in s[2:4])...)
            return V / X  # 10 × 3, row i is ∇N_i on the sub-tetrahedron
        end
        for (s, G) in zip(subtets, gradients), weights in ((0.25, 0.25, 0.25, 0.25), (0.1, 0.2, 0.3, 0.4))
            ξ = sum(w * v[1] for (w, v) in zip(weights, s))
            N, _, _ = Norma.composite_tetra10_shape_functions(ξ)
            @test N ≈ sum(w * v[2] for (w, v) in zip(weights, s)) atol = 1.0e-14
        end
        nodes = Norma.COMPOSITE_TETRA10_POINTS
        @test hcat((Norma.composite_tetra10_values(SVector{3}(nodes[:, a])) for a in 1:10)...) ≈ I atol = 1.0e-15

        # The projected gradients are the L2 projection, over the parent tet,
        # of those piecewise constant gradients onto the linear polynomials: in
        # the barycentric basis λ, M c = b with M = ∫ λλᵀ and b = Σ ∫ λ ∇N.
        λ(ξ) = SVector(1 - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3])
        M = [(i == j ? 2.0 : 1.0) / 120 for i in 1:4, j in 1:4]
        b = sum(vol * λ(sum(v[1] for v in s) / 4) * vec(G')' for (s, G, vol) in zip(subtets, gradients, volumes))
        coefficients = M \ b  # 4 × 30, column 3(i - 1) + j is ∂N_i/∂ξ_j
        for G in (4, 5, 14)
            N, dN, w, ξ = Norma.isoparametric(Norma.COMPOSITE_TETRA10, G)
            @test sum(w) ≈ 1 / 6
            for p in 1:G
                @test dN[:, :, p] ≈ reshape(λ(ξ[:, p])' * coefficients, 3, 10) atol = 1.0e-13
                @test sum(N[:, p]) ≈ 1.0
            end
        end

        # An affine deformation is reproduced exactly on a distorted
        # straight-sided tet: the projection of a constant gradient is itself.
        X = [0.1 1.3 0.2 0.3; -0.2 0.1 1.1 0.2; 0.05 0.1 -0.1 0.9]
        X = hcat(X, [(X[:, i] + X[:, j]) / 2 for (i, j) in ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))]...)
        A = @SMatrix [1.01 0.02 -0.03; 0.04 0.97 0.01; -0.02 0.03 1.05]
        x = A * X .+ [0.1, -0.2, 0.3]
        N, dN, w, _ = Norma.isoparametric(Norma.COMPOSITE_TETRA10, 4)
        for p in 1:4
            J = dN[:, :, p] * X'
            @test x * (J \ dN[:, :, p])' ≈ A atol = 1.0e-13
        end
        @test sum(det(dN[:, :, p] * X') * w[p] for p in 1:4) ≈ abs(det(X[:, 2:4] .- X[:, 1])) / 6
    end

    @testset "Intrepid2 unit test" begin
        # The checks of Intrepid2's unit test of Basis_HGRAD_TET_COMP12_FEM
        # (Trilinos, packages/intrepid2/unit-test/Discretization/Basis/
        # HGRAD_TET_COMP12_FEM/test_01.hpp), with its tables, here to the
        # digits of Float64, and its tolerance of 100 ε.
        tol = 100 * eps()
        shape_values(ξ) = Norma.composite_tetra10_values(SVector{3,Float64}(ξ...))
        shape_gradients(ξ) = Norma.composite_tetra10_gradients(SVector{3,Float64}(ξ...))
        nodes = (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.5, 0.0, 0.0),
            (0.5, 0.5, 0.0),
            (0.0, 0.5, 0.0),
            (0.0, 0.0, 0.5),
            (0.5, 0.0, 0.5),
            (0.0, 0.5, 0.5),
        )
        kronecker = [Float64(a == b) for a in 1:10, b in 1:10]
        @test maximum(maximum(abs.(shape_values(x) - kronecker[:, a])) for (a, x) in enumerate(nodes)) ≤ tol
        # The values and the projected gradients at the points of the 5-point
        # rule and of the 4-point rule: a row of values per point, and a
        # 10 × 3 block of gradients per point, with ∂N_a/∂ξ in row a.
        a4, b4 = 0.1381966011250105, 0.5854101966249684
        points = (
            (0.25, 0.25, 0.25),
            (0.5, 1 / 6, 1 / 6),
            (1 / 6, 0.5, 1 / 6),
            (1 / 6, 1 / 6, 0.5),
            (1 / 6, 1 / 6, 1 / 6),
            (a4, a4, a4),
            (b4, a4, a4),
            (a4, b4, a4),
            (a4, a4, b4),
        )
        v1, v2 = 0.17082039324993692, 0.276393202250021
        point_values = [
            0 0 0 0 1/6 1/6 1/6 1/6 1/6 1/6
            0 0 0 0 1/3 1/3 0 0 1/3 0
            0 0 0 0 0 1/3 1/3 0 0 1/3
            0 0 0 0 0 0 0 1/3 1/3 1/3
            0 0 0 0 1/3 0 1/3 1/3 0 0
            v1 0 0 0 v2 0 v2 v2 0 0
            0 v1 0 0 v2 v2 0 0 v2 0
            0 0 v1 0 0 v2 v2 0 0 v2
            0 0 0 v1 0 0 0 v2 v2 v2
        ]
        g1, g2, g3, g4 = 1.0885254915624212, 0.029508497187473712, 1.3043729868748772, 0.5636610018750176
        g5, g6, g7, g8 = 0.37732200375003505, 0.18633899812498247, 1.868033988749895, 0.19098300562505258
        point_gradients = (
            [
                -1/4 -1/4 -1/4
                1/4 0 0
                0 1/4 0
                0 0 1/4
                0 -3/4 -3/4
                3/4 3/4 0
                -3/4 0 -3/4
                -3/4 -3/4 0
                3/4 0 3/4
                0 3/4 3/4
            ],
            [
                -1/24 -1/24 -1/24
                7/8 0 0
                0 1/24 0
                0 0 1/24
                -35/36 -19/12 -19/12
                11/18 19/12 0
                -17/36 0 -1/3
                -17/36 -1/3 0
                11/18 0 19/12
                -5/36 1/3 1/3
            ],
            [
                -1/24 -1/24 -1/24
                1/24 0 0
                0 7/8 0
                0 0 1/24
                0 -17/36 -1/3
                19/12 11/18 0
                -19/12 -35/36 -19/12
                -1/3 -17/36 0
                1/3 -5/36 1/3
                0 11/18 19/12
            ],
            [
                -1/24 -1/24 -1/24
                1/24 0 0
                0 1/24 0
                0 0 7/8
                0 -1/3 -17/36
                1/3 1/3 -5/36
                -1/3 0 -17/36
                -19/12 -19/12 -35/36
                19/12 0 11/18
                0 19/12 11/18
            ],
            [
                -7/8 -7/8 -7/8
                1/24 0 0
                0 1/24 0
                0 0 1/24
                35/36 -11/18 -11/18
                17/36 17/36 5/36
                -11/18 35/36 -11/18
                -11/18 -11/18 35/36
                17/36 5/36 17/36
                5/36 17/36 17/36
            ],
            [
                -g1 -g1 -g1
                -g2 0 0
                0 -g2 0
                0 0 -g2
                g3 -g4 -g4
                g5 g5 g6
                -g4 g3 -g4
                -g4 -g4 g3
                g5 g6 g5
                g6 g5 g5
            ],
            [
                g2 g2 g2
                g1 0 0
                0 -g2 0
                0 0 -g2
                -g3 -g7 -g7
                g4 g7 0
                -g5 0 -g8
                -g5 -g8 0
                g4 0 g7
                -g6 g8 g8
            ],
            [
                g2 g2 g2
                -g2 0 0
                0 g1 0
                0 0 -g2
                0 -g5 -g8
                g7 g4 0
                -g7 -g3 -g7
                -g8 -g5 0
                g8 -g6 g8
                0 g4 g7
            ],
            [
                g2 g2 g2
                -g2 0 0
                0 -g2 0
                0 0 g1
                0 -g8 -g5
                g8 g8 -g6
                -g8 0 -g5
                -g7 -g7 -g3
                g7 0 g4
                0 g7 g4
            ],
        )
        @test maximum(maximum(abs.(shape_values(x) - point_values[p, :])) for (p, x) in enumerate(points)) ≤ tol
        @test maximum(maximum(abs.(shape_gradients(x)' - point_gradients[p])) for (p, x) in enumerate(points)) ≤ tol
        # The partition of unity at 16384 random points of the parent.
        rng = Random.Xoshiro(1)
        random_points = NTuple{3,Float64}[]
        while length(random_points) < 16384
            ξ = (rand(rng), rand(rng), rand(rng))
            sum(ξ) ≤ 1 && push!(random_points, ξ)
        end
        @test maximum(abs(sum(shape_values(ξ)) - 1) for ξ in random_points) ≤ tol
        # Points of the face r + s + t = 1 from normalized weights, whose sum
        # can exceed 1 by an ulp, and points of the boundary moved 4 ulps
        # outward: the values still sum to 1, and they are continuous across
        # the boundary.
        face_points = map(1:4096) do _
            w = (rand(rng), rand(rng), rand(rng))
            return w ./ sum(w)
        end
        boundary = (
            (0.3, 0.0, 0.3),
            (0.0, 0.2, 0.2),
            (0.2, 0.3, 0.0),
            (0.25, 0.25, 0.5),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.5, 0.5, 0.0),
            (0.0, 0.5, 0.5),
        )
        outward = (
            (0.0, -1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, 0.0, -1.0),
            (1.0, 1.0, 1.0),
            (-1.0, -1.0, -1.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, 1.0),
        )
        outside = [x .+ 4 * eps() .* n for (x, n) in zip(boundary, outward)]
        @test maximum(abs(sum(shape_values(ξ)) - 1) for ξ in [face_points; collect(boundary); outside]) ≤ tol
        @test maximum(maximum(abs.(shape_values(o) - shape_values(x))) for (x, o) in zip(boundary, outside)) ≤ tol
    end

    @testset "Mass" begin
        # On a straight-sided tetrahedron, the exact mass of the sub-tetrahedra
        # is LCM's analytic composite tet mass (AnalyticMassResidual: Use
        # Analytic Mass), V/1440 times this matrix, and its row sums LCM's
        # lumped one, V (3, 14)/96.
        corner_rows = [18 0 0 0 9 0 9 9 0 0; 0 18 0 0 9 9 0 0 9 0; 0 0 18 0 0 9 9 0 0 9; 0 0 0 18 0 0 0 9 9 9]
        edge_rows = [
            9 9 0 0 80 26 26 26 26 8
            0 9 9 0 26 80 26 8 26 26
            9 0 9 0 26 26 80 26 8 26
            9 0 0 9 26 8 26 80 26 26
            0 9 0 9 26 26 8 26 80 26
            0 0 9 9 8 26 26 26 26 80
        ]
        lcm_mass = [corner_rows; edge_rows] / 1440
        N, dN, w, _ = Norma.mass_quadrature(Norma.COMPOSITE_TETRA10, 4)
        mass = 6 * sum(w[p] * N[:, p] * N[:, p]' for p in eachindex(w))
        @test mass ≈ lcm_mass atol = 1.0e-15
        @test vec(sum(mass; dims=2)) ≈ [fill(3 / 96, 4); fill(14 / 96, 6)] atol = 1.0e-15
        # Its gradients are the true ones: the volume of a distorted tet.
        X = [0.0 1.2 0.1 0.3; 0.0 0.1 0.9 0.2; 0.0 -0.1 0.2 1.1]
        X = hcat(X, [(X[:, i] + X[:, j]) / 2 for (i, j) in ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))]...)
        @test sum(det(dN[:, :, p] * X') * w[p] for p in eachindex(w)) ≈ det(X[:, 2:4]) / 6
    end

    @testset "Patch test" begin
        # Exact on distorted straight-sided tets. On curved ones the projected
        # gradients are not those of the geometry, so the patch test holds
        # only in the limit: the stress error halves with the element size.
        u_error, σ_error = composite_tet10_patch(3; straight=true)
        @test u_error < 1.0e-15
        @test σ_error < 1.0e-14
        _, coarse = composite_tet10_patch(2; straight=false)
        _, fine = composite_tet10_patch(4; straight=false)
        @test 1.0e-6 < fine < 2.0e-5  # the stress is about 0.01
        @test coarse / fine > 1.7
    end

    @testset "Thick-walled cylinder" begin
        # Second-order convergence of the nodal displacement, as for the
        # piecewise linear interpolation on the sub-tetrahedra.
        errors = [composite_tet10_lame(n, 0.3) for n in (2, 4, 8)]
        @test errors[3] < 2.5e-3
        @test all(log2(errors[i] / errors[i + 1]) > 1.7 for i in 1:2)
        # Nearly incompressible: at a fixed mesh the error of the composite tet
        # stays bounded as ν → 1/2, while that of the TETRA4 grows without
        # bound. It still grows 31-fold from ν = 0.3, the partial locking that
        # Ostien et al. (2016) report for this formulation.
        composite = [composite_tet10_lame(4, ν; straight=true) for ν in (0.4999, 0.499999)]
        tetra4 = [composite_tet10_lame(4, ν; kind=:tetra4, straight=true) for ν in (0.4999, 0.499999)]
        @test composite[2] < 0.3
        @test composite[2] / composite[1] < 2
        @test tetra4[2] > 3
        @test tetra4[2] / tetra4[1] > 4
    end
end

# A simulation of a structured mesh, built but not run, so that its meshes
# stay open for recovery and transfer: the composite tet 10 or, on a TETRA4
# mesh, the TETRA4, with a nodal recovery of every field when `recovery` is
# given. The material is linear elastic unless `material` is given.
function composite_tet10_simulation(
    name, mesh; composite=true, num_points=nothing, recovery=nothing, key=nothing, material=nothing
)
    elastic = Dict{String,Any}(
        "model" => "linear elastic", "elastic modulus" => 1.0, "Poisson's ratio" => 0.3, "density" => 1.0
    )
    model = Dict{String,Any}(
        "type" => "solid mechanics",
        "material" => Dict{String,Any}(
            "blocks" => Dict("solid" => "elastic"), "elastic" => material === nothing ? elastic : material
        ),
    )
    composite && (model["composite tet 10"] = ["solid"])
    key === nothing || (model["composite tet 10"] = key)
    num_points === nothing || (model["num integration points"] = Dict("solid" => num_points))
    if recovery !== nothing
        model["nodal recovery"] = Dict{String,Any}(
            "method" => recovery,
            "stress" => true,
            "von mises stress" => true,
            "deformation gradient" => true,
            "internal variables" => true,
        )
    end
    params = Norma.Parameters(
        "name" => name * ".yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => name * ".e",
        "Exodus output interval" => 1.0,
        "CSV output interval" => 0,
        "model" => model,
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "boundary conditions" => Dict(
            "Dirichlet" => [
                Dict{String,Any}("node set" => "xm", "component" => c, "function" => "0.0") for c in ("x", "y", "z")
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => "direct",
            "minimum iterations" => 1,
            "maximum iterations" => 3,
            "absolute tolerance" => 1.0e-14,
            "relative tolerance" => 1.0e-12,
        ),
    )
    return Norma.create_simulation(params)
end

function close_composite_tet10_simulation(sim, name)
    Exodus.close(sim.params["input_mesh"])
    Exodus.close(sim.params["output_mesh"])
    rm(name * ".e"; force=true)
    return nothing
end

# The recovery of composite blocks loads the integration point values through
# their linear fit and projects them with the exact mass on the
# sub-tetrahedra. The consistent mass from the block's own 4 or 5 points was
# singular (a PosDefException, or garbage on some meshes), and with the
# 5-point rule the corner functions vanish at every point, so the lumped
# recovery and the transfer left the corner nodes at zero.
@testset "Composite Tet 10 Recovery And Transfer" begin
    # Distorted but straight-sided tetrahedra, on which F is exact.
    shear(p) = SVector(p[1] + 0.2 * p[2], p[2] + 0.1 * p[3], p[3] + 0.15 * p[1])
    write_structured_tet_mesh("composite-recovery.g", (2, 2, 2), shear; straight=true)
    write_structured_tet_mesh("composite-source.g", (2, 2, 2), shear; quadratic=false)
    A = [0.01 0.002 -0.003; 0.004 -0.005 0.001; 0.002 0.003 0.006]
    σ = [1.0, 2.0, 3.0, 0.4, 0.5, 0.6]
    linear(x) = 1.0 + x[1] + 2.0 * x[2] - x[3]
    quadratic(x) = x[1]^2 + x[2] * x[3] + 0.5 * x[3]^2
    j2 = Dict{String,Any}(
        "model" => "j2 plasticity",
        "elastic modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "density" => 1.0,
        "yield stress" => 1.0e-3,
        "hardening modulus" => 0.1,
    )
    source = composite_tet10_simulation("composite-source", "composite-source.g"; composite=false)
    source_field = reshape([linear(source.model.reference[:, n]) for n in axes(source.model.reference, 2)], 1, :)
    for num_points in (4, 5, 14)
        # A constant stress, its von Mises stress and the F of an affine
        # displacement are recovered exactly at every node by every method.
        for method in ("lumped", "consistent", "both")
            sim = composite_tet10_simulation("composite-recovery", "composite-recovery.g"; num_points, recovery=method)
            model = sim.model
            for element in model.stress[1]
                fill!(element, SVector{6,Float64}(σ))
            end
            model.displacement .= A * model.reference
            Norma.recover_stress!(model)
            Norma.recover_von_mises_stress!(model)
            Norma.recover_deformation_gradient!(model)
            for prefix in (method == "both" ? ("lumped_", "consistent_") : ("",))
                field(name) = getfield(model, Symbol(prefix * "recovered_" * name))
                @test maximum(abs.(field("stress") .- σ)) < 1.0e-12
                @test maximum(abs.(field("von_mises") .- Norma.von_mises_from_voigt(σ))) < 1.0e-12
                @test maximum(abs.(field("F") .- vec(I + A))) < 1.0e-12
            end
            close_composite_tet10_simulation(sim, "composite-recovery")
            # A constant internal variable of every point is recovered
            # exactly too, through the same composite load.
            sim = composite_tet10_simulation(
                "composite-recovery", "composite-recovery.g"; num_points, recovery=method, material=j2
            )
            model = sim.model
            names = Norma.collect_internal_variable_names(model.materials)
            values = collect(1.0:length(names))
            for element in model.state[1], point_state in element
                point_state .= values
            end
            Norma.recover_internal_variables!(model, names)
            for prefix in (method == "both" ? ("lumped_", "consistent_") : ("",))
                recovered = getfield(model, Symbol(prefix * "recovered_internal_variables"))
                @test size(recovered, 1) == length(names)
                @test maximum(abs.(recovered .- values)) < 1.0e-12
            end
            close_composite_tet10_simulation(sim, "composite-recovery")
        end
        # A field of the composite space transferred onto the same mesh comes
        # back unchanged: the load samples the nodal source at the points of
        # the exact mass rule. With the block's own points it kept only the
        # linear fit of each element, off by 10-40 % of the field.
        transfer_source = composite_tet10_simulation(
            "composite-transfer", "composite-recovery.g"; num_points, recovery="consistent"
        )
        sim = composite_tet10_simulation(
            "composite-recovery", "composite-recovery.g"; num_points, recovery="consistent"
        )
        field = reshape([quadratic(sim.model.reference[:, n]) for n in axes(sim.model.reference, 2)], 1, :)
        out = Norma.transfer_field(transfer_source.model, field, sim.model)
        @test maximum(abs.(out .- field)) < 1.0e-12
        close_composite_tet10_simulation(sim, "composite-recovery")
        close_composite_tet10_simulation(transfer_source, "composite-transfer")
        # A linear field transferred from the TETRA4 mesh of the same cube:
        # exact with the consistent recovery, as the composite space holds the
        # linear functions, and a constant exact with the lumped one.
        sim = composite_tet10_simulation(
            "composite-recovery", "composite-recovery.g"; num_points, recovery="consistent"
        )
        out = Norma.transfer_field(source.model, source_field, sim.model)
        exact = [linear(sim.model.reference[:, n]) for n in axes(sim.model.reference, 2)]
        @test maximum(abs.(out[1, :] .- exact)) < 1.0e-12
        close_composite_tet10_simulation(sim, "composite-recovery")
        sim = composite_tet10_simulation("composite-recovery", "composite-recovery.g"; num_points, recovery="lumped")
        out = Norma.transfer_field(source.model, ones(1, size(source_field, 2)), sim.model)
        @test maximum(abs.(out .- 1.0)) < 1.0e-12
        # The blended kinetic energy of an overlap (all weights 1 here) with
        # a velocity on the corner nodes only: the lumped mass of the corner
        # nodes is 3/96 of the element volume (Ostien et al.), and the
        # consistent energy that of the exact element mass.
        model = sim.model
        block = model.blocks[1]
        corners = unique(vec(block.connectivity[1:4, :]))
        model.velocity .= 0.0
        model.velocity[1, corners] .= 1.0
        model.acceleration .= 0.0
        weights = [ones(block.num_points) for _ in 1:(block.num_elements)]
        volume = 1.0 * det([1.0 0.2 0.0; 0.0 1.0 0.1; 0.15 0.0 1.0])
        for (half_dt, expected) in ((0.5, 0.5 * volume / 8), (0.0, nothing))
            stored, kinetic = zeros(Threads.maxthreadid()), zeros(Threads.maxthreadid())
            Norma.blended_block_energy!(
                stored, kinetic, model, block, 1, model.materials[1], 1.0, weights, true, half_dt
            )
            if expected === nothing
                rule = Norma.composite_tetra10_mass_rule()
                expected = 0.0
                for e in 1:(block.num_elements)
                    nodes = block.connectivity[:, e]
                    M = Norma.composite_tetra10_element_mass(model.reference[:, nodes], rule)
                    v = model.velocity[1, nodes]
                    expected += 0.5 * dot(v, M * v)
                end
            end
            @test sum(kinetic) ≈ expected rtol = 1.0e-12
        end
        close_composite_tet10_simulation(sim, "composite-recovery")
    end
    close_composite_tet10_simulation(source, "composite-source")
    rm("composite-recovery.g"; force=true)
    rm("composite-source.g"; force=true)
end

# Composite blocks with the 4-point rule next to composite blocks with the
# 5-point rule: the Exodus output defines on each block the variables of its
# own points only. With every variable declared on every block, the 4-point
# block's stress_xx_5 and von_mises_stress_5 existed and read as zeros.
@testset "Composite Tet 10 Output Of Blocks With Different Rules" begin
    mesh = "composite-mixed.g"
    cp(joinpath("lcm-reference", "meshes", "RubiksCubeTet10.g"), mesh; force=true)
    names = mesh_block_names(mesh)
    params = Norma.Parameters(
        "name" => "composite-mixed.yaml",
        "type" => "single",
        "input mesh file" => mesh,
        "output mesh file" => "composite-mixed.e",
        "Exodus output interval" => 1.0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "solid mechanics",
            "material" => Dict{String,Any}(
                "blocks" => Dict(name => "elastic" for name in names),
                "elastic" => Dict{String,Any}(
                    "model" => "linear elastic",
                    "elastic modulus" => 1.0,
                    "Poisson's ratio" => 0.0,
                    "density" => 1.0,
                ),
            ),
            "composite tet 10" => true,
            "num integration points" => Dict(name => 5 for name in names[2:end]),
        ),
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "boundary conditions" => Dict(
            "Dirichlet" => [
                Dict{String,Any}("node set" => n, "component" => c, "function" => f) for (n, c, f) in (
                    ("nodelist_1", "x", "0.0"),
                    ("nodelist_2", "x", "1.0e-3 * t"),
                    ("nodelist_3", "y", "0.0"),
                    ("nodelist_5", "z", "0.0"),
                )
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "Hessian minimizer",
            "step" => "full Newton",
            "linear solver" => "direct",
            "minimum iterations" => 1,
            "maximum iterations" => 3,
            "absolute tolerance" => 1.0e-14,
            "relative tolerance" => 1.0e-12,
        ),
    )
    Norma.run(params)
    exo = Exodus.ExodusDatabase("composite-mixed.e", "r")
    try
        ids = Exodus.read_ids(exo, Block)
        output_names = Exodus.read_names(exo, Block)
        four_id = ids[findfirst(==(names[1]), output_names)]
        five_id = ids[findfirst(==(names[2]), output_names)]
        # Uniaxial stress of 1e-3 along x, exact with both rules.
        @test maximum(abs.(Exodus.read_values(exo, ElementVariable, 2, four_id, "stress_xx_4") .- 1.0e-3)) < 1.0e-12
        @test maximum(abs.(Exodus.read_values(exo, ElementVariable, 2, five_id, "stress_xx_5") .- 1.0e-3)) < 1.0e-12
        for variable in ("stress_xx_5", "von_mises_stress_5")
            @test_throws Exodus.ExodusError Exodus.read_values(exo, ElementVariable, 2, four_id, variable)
        end
        @test length(Exodus.read_values(exo, ElementVariable, 2, four_id, "stored_energy")) == 6
    finally
        close(exo)
        rm(mesh; force=true)
        rm("composite-mixed.e"; force=true)
    end
end

# The search for a point in a mesh inverts the map of each candidate element.
# For the composite tet it used Newton's method with the projected gradients,
# which are not the derivative of the map on curved elements, so on a strongly
# curved element it diverged, and the inside test, which checks only the sum
# of the barycentric coordinates, accepted the garbage point. The map is now
# inverted exactly, on the sub-tetrahedra.
@testset "Composite Tet 10 Point Inversion" begin
    # A Kuhn tetrahedron of the unit cube on a quarter turn, r = 1 + p1 and
    # θ = π p2 / 2, with its edge nodes on the arcs.
    φ(p) = SVector((1 + p[1]) * cos(π / 2 * p[2]), (1 + p[1]) * sin(π / 2 * p[2]), p[3])
    corners = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0), SVector(1.0, 1.0, 1.0)]
    edges = ((1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4))
    X = Matrix(hcat([φ(c) for c in corners]..., [φ((corners[a] + corners[b]) / 2) for (a, b) in edges]...))
    if det(hcat(X[:, 2] - X[:, 1], X[:, 3] - X[:, 1], X[:, 4] - X[:, 1])) < 0
        X = X[:, [1, 3, 2, 4, 7, 6, 5, 8, 10, 9]]
    end
    worst_ξ, worst_x, outside = 0.0, 0.0, 0
    for k in 1:500
        # Points spread over the parent tetrahedron, without a random number generator.
        ξ = SVector(mod(0.6180339887 * k, 1.0), mod(0.4142135624 * k, 1.0), mod(0.7320508076 * k, 1.0))
        sum(ξ) < 1.0 || continue
        x = X * Norma.composite_tetra10_values(ξ)
        found = Norma.map_to_parametric(Norma.COMPOSITE_TETRA10, X, Vector(x))
        worst_ξ = max(worst_ξ, norm(found - ξ))
        worst_x = max(worst_x, norm(X * Norma.composite_tetra10_values(SVector{3,Float64}(found)) - x))
        Norma.is_inside_parametric(Norma.COMPOSITE_TETRA10, found) || (outside += 1)
    end
    @test worst_ξ < 1.0e-12
    @test worst_x < 1.0e-12
    @test outside == 0
    # A point holds in the element only when it lies in one of the
    # sub-tetrahedra the geometry defines. Outside, the composite functions
    # extend from another sub-tetrahedron than the one the inverse map used,
    # so the point would map elsewhere; the sum-only test of the other
    # simplices accepted most such points.
    centroid = sum(X[:, a] for a in 5:10) / 6
    vertex(v) = v == 11 ? centroid : X[:, v]
    function containment(x)
        worst = Inf
        for vertices in Norma.COMPOSITE_TETRA10_SUBTETS
            X1 = vertex(vertices[1])
            λ = hcat(vertex(vertices[2]) - X1, vertex(vertices[3]) - X1, vertex(vertices[4]) - X1) \ (x - X1)
            worst = min(worst, max(-minimum(λ), sum(λ) - 1.0, 0.0))
        end
        return worst
    end
    low, high = vec(minimum(X; dims=2)), vec(maximum(X; dims=2))
    wrong, mapped_wrong, outside_points = 0, 0.0, 0
    for k in 1:3000
        u = SVector(mod(0.6180339887 * k, 1.0), mod(0.4142135624 * k, 1.0), mod(0.7320508076 * k, 1.0))
        x = low + (high - low) .* u
        violation = containment(x)
        1.0e-9 < violation < 1.0e-3 && continue
        ξ, inside = Norma.is_inside(Norma.COMPOSITE_TETRA10, X, Vector(x))
        outside_points += violation > 0.0
        inside == (violation == 0.0) || (wrong += 1)
        inside &&
            (mapped_wrong = max(mapped_wrong, norm(X * Norma.composite_tetra10_values(SVector{3,Float64}(ξ)) - x)))
    end
    @test outside_points > 1000
    @test wrong == 0
    @test mapped_wrong < 1.0e-12
end

# A linear field transferred from a mesh of curved composite tets (the unit
# cube with a bump inside, so its boundary stays flat) to a TETRA4 mesh of
# the same cube is exact, as the composite map interpolates the geometry.
# With the sum-only inside test, points were taken from neighboring elements
# that do not hold them, and the transferred field was off by up to 0.14.
@testset "Composite Tet 10 Transfer From Curved Elements" begin
    bump(p) = p + 0.15 * sin(π * p[1]) * sin(π * p[2]) * sin(π * p[3]) * SVector(1.0, -0.5, 0.8)
    write_structured_tet_mesh("composite-curved.g", (2, 2, 2), bump)
    write_structured_tet_mesh("composite-cube.g", (2, 2, 2), p -> p; quadratic=false)
    linear(x) = 1.0 + x[1] + 2.0 * x[2] - x[3]
    source = composite_tet10_simulation("composite-curved", "composite-curved.g"; recovery="consistent")
    destination = composite_tet10_simulation(
        "composite-cube", "composite-cube.g"; composite=false, recovery="consistent"
    )
    field = reshape([linear(source.model.reference[:, n]) for n in axes(source.model.reference, 2)], 1, :)
    out = Norma.transfer_field(source.model, field, destination.model)
    exact = [linear(destination.model.reference[:, n]) for n in axes(destination.model.reference, 2)]
    @test maximum(abs.(out[1, :] .- exact)) < 1.0e-12
    close_composite_tet10_simulation(source, "composite-curved")
    close_composite_tet10_simulation(destination, "composite-cube")
    rm("composite-curved.g"; force=true)
    rm("composite-cube.g"; force=true)
end

# The `composite tet 10` key: a map value that is not a Bool (a quoted "true"
# or "yes") was silently off, and `true` on a mesh without TETRA10 blocks
# silently ran the standard elements.
# The overlap couplings with composite blocks. The weak coupling integrates
# the transfer onto a facet over the sections of the source elements, those of
# the parent tetrahedra for a composite block: a linear field, which the
# composite functions reproduce, comes across exactly, and the row sums of L
# are those of W. The matching coupling weights the volume element at the
# integration points, which the exact mass rule of a composite block does not
# use, so it rejects composite blocks. The consistent-mass term of the rate
# invariant is that of the exact mass.
@testset "Composite Tet 10 Overlap Couplings" begin
    shear(p) = SVector(p[1] + 0.2 * p[2], p[2] + 0.1 * p[3], p[3] + 0.15 * p[1])
    write_structured_tet_mesh("composite-overlap.g", (2, 2, 2), shear; straight=true)
    sim = composite_tet10_simulation("composite-overlap", "composite-overlap.g")
    model = sim.model
    block = model.blocks[1]
    @test block.element_type == Norma.COMPOSITE_TETRA10
    # A slanted triangle inside the sheared cube, one node per column.
    facet = [0.3 0.8 0.5; 0.3 0.35 0.75; 0.45 0.5 0.6]
    L = zeros(3, size(model.reference, 2))
    W = zeros(3, 3)
    boxes = Norma.element_bounding_boxes(model, block)
    @test Norma.add_overlap_facet_segmented_quadrature!(L, W, [1, 2, 3], facet, model, block, boxes)
    linear(x) = 1.0 + x[1] + 2.0 * x[2] - x[3]
    source = [linear(model.reference[:, n]) for n in axes(model.reference, 2)]
    destination = [linear(facet[:, a]) for a in 1:3]
    @test maximum(abs.(L * source - W * destination)) < 1.0e-14
    @test vec(sum(L; dims=2)) ≈ vec(sum(W; dims=2)) rtol = 1.0e-13
    @test_throws Exception Norma.set_assembly_weights!(sim)
    model.acceleration .= randn(Random.Xoshiro(2), size(model.acceleration))
    _, _, kinetic = Norma.quadratic_energies(model)
    rule = Norma.composite_tetra10_mass_rule()
    expected = 0.0
    for e in 1:(block.num_elements)
        nodes = block.connectivity[:, e]
        M = Norma.composite_tetra10_element_mass(model.reference[:, nodes], rule)
        a = model.acceleration[:, nodes]
        expected += 0.5 * sum(dot(a[i, :], M * a[i, :]) for i in 1:3)
    end
    @test kinetic ≈ expected rtol = 1.0e-13
    close_composite_tet10_simulation(sim, "composite-overlap")
    rm("composite-overlap.g"; force=true)
end

@testset "Composite Tet 10 Key" begin
    names = ["solid"]
    Norma.NORMA_TEST_MODE[] = true
    try
        @test Norma.composite_tet10_blocks(Dict("composite tet 10" => Dict("solid" => true)), names) ==
            (false, Set(["solid"]))
        @test Norma.composite_tet10_blocks(Dict("composite tet 10" => Dict("solid" => false)), names) ==
            (false, Set{String}())
        for value in ("true", "yes", 1)
            @test_throws Norma.NormaAbortException Norma.composite_tet10_blocks(
                Dict("composite tet 10" => Dict("solid" => value)), names
            )
        end
    finally
        Norma.NORMA_TEST_MODE[] = false
    end
    write_structured_tet_mesh("composite-key.g", (1, 1, 1), p -> p; quadratic=false)
    sim = Ref{Any}(nothing)
    log_file, log = mktemp()
    redirect_stdout(log) do
        sim[] = composite_tet10_simulation("composite-key", "composite-key.g"; composite=false, key=true)
    end
    close(log)
    @test occursin("switches no block", read(log_file, String))
    rm(log_file; force=true)
    @test sim[].model.blocks[1].element_type == Norma.TETRA4
    close_composite_tet10_simulation(sim[], "composite-key")
    rm("composite-key.g"; force=true)
end

# The example of examples/element-types/composite-tet10/cube: the TETRA10 cube
# of tet10-static-solid-cube.jl with the composite tet, in the same uniaxial
# stress.
@testset "Composite Tet 10 Example Cube" begin
    for file in ("cube.yaml", "cube.g")
        cp(joinpath("../examples/element-types/composite-tet10/cube", file), file; force=true)
    end
    sim = Norma.run("cube.yaml")
    for file in ("cube.yaml", "cube.g", "cube.e")
        rm(file; force=true)
    end
    @test sim.model.blocks[1].element_type == Norma.COMPOSITE_TETRA10
    u = sim.integrator.displacement
    @test minimum(u[1:3:end]) ≈ -0.25 rtol = 1.0e-06
    @test minimum(u[2:3:end]) ≈ -0.25 rtol = 1.0e-06
    @test maximum(u[3:3:end]) ≈ 1.0 rtol = 1.0e-06
    stress = average_components(sim.model.stress)
    @test stress[3] ≈ 1.0e+09 rtol = 1.0e-06
    @test maximum(abs, stress[[1, 2, 4, 5, 6]]) < 1.0e-06
end
