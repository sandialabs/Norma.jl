using LinearAlgebra
using Random
using Exodus

# Under the test suite, Norma is already loaded by runtests.jl and its package
# extensions (e.g. NormaPyTorchExt, the neural-network OpInf backend) are active.
# Re-including the source would define a second `Main.Norma` without those
# extension methods, breaking later tests that rely on them. Only load from
# source when this file is run on its own.
if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

a = 1
t = sqrt(3) / 2 * a
h = sqrt(2.0 / 3.0) * a
reg_tet_coords = [
    -t/3 -t/3 2*t/3 0
    a/2 -a/2 0.0 0
    0.0 0.0 0.0 h
]
tet_conn = reshape(
    [
        1
        2
        3
        4
    ],
    4,
    1,
)
tet_base = [1, 2, 3]
tet_back_edge = [1, 2]
tet_first_node = [1]

@testset "smooth_reference" begin
    ref_tet_eqvol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, reg_tet_coords);
    @test size(ref_tet_eqvol) == (3, 4)
    a1 = norm(ref_tet_eqvol[:, 1] - ref_tet_eqvol[:, 2])
    a2 = norm(ref_tet_eqvol[:, 1] - ref_tet_eqvol[:, 3])
    a3 = norm(ref_tet_eqvol[:, 1] - ref_tet_eqvol[:, 4])
    a4 = norm(ref_tet_eqvol[:, 2] - ref_tet_eqvol[:, 3])
    a5 = norm(ref_tet_eqvol[:, 2] - ref_tet_eqvol[:, 4])
    a6 = norm(ref_tet_eqvol[:, 3] - ref_tet_eqvol[:, 4])
    @test a1 ≈ a atol = 1.0e-12
    @test a1 ≈ a2 atol = 1.0e-12
    @test a1 ≈ a3 atol = 1.0e-12
    @test a1 ≈ a4 atol = 1.0e-12
    @test a1 ≈ a5 atol = 1.0e-12
    @test a1 ≈ a6 atol = 1.0e-12

    ref_tet_eqedl = Norma.create_smooth_reference("average edge length", Norma.TETRA4, reg_tet_coords);
    @test size(ref_tet_eqedl) == (3, 4)
    a1 = norm(ref_tet_eqedl[:, 1] - ref_tet_eqedl[:, 2])
    a2 = norm(ref_tet_eqedl[:, 1] - ref_tet_eqedl[:, 3])
    a3 = norm(ref_tet_eqedl[:, 1] - ref_tet_eqedl[:, 4])
    a4 = norm(ref_tet_eqedl[:, 2] - ref_tet_eqedl[:, 3])
    a5 = norm(ref_tet_eqedl[:, 2] - ref_tet_eqedl[:, 4])
    a6 = norm(ref_tet_eqedl[:, 3] - ref_tet_eqedl[:, 4])
    @test a1 ≈ a atol = 1.0e-12
    @test a1 ≈ a2 atol = 1.0e-12
    @test a1 ≈ a3 atol = 1.0e-12
    @test a1 ≈ a4 atol = 1.0e-12
    @test a1 ≈ a5 atol = 1.0e-12
    @test a1 ≈ a6 atol = 1.0e-12

    rand_tests = 10
    for _ in 1:rand_tests
        random_tet = reg_tet_coords + Random.randn(3, 4) * 0.1 * a
        random_tet_vol = det(random_tet[:, 2:4] .- random_tet[:, 1]) / 6.0
        random_tet_edl =
            norm(random_tet[:, 1] - random_tet[:, 2]) +
            norm(random_tet[:, 1] - random_tet[:, 3]) +
            norm(random_tet[:, 1] - random_tet[:, 4]) +
            norm(random_tet[:, 2] - random_tet[:, 3]) +
            norm(random_tet[:, 3] - random_tet[:, 4]) +
            norm(random_tet[:, 4] - random_tet[:, 2])

        ref_tet_eqvol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, random_tet);
        a_eqvol = norm(ref_tet_eqvol[:, 1] - ref_tet_eqvol[:, 2])
        @test random_tet_vol ≈ a_eqvol^3 / 6.0 / sqrt(2.0) atol = 1.0e-12

        ref_tet_eqedl = Norma.create_smooth_reference("average edge length", Norma.TETRA4, random_tet);
        a_eqedl = norm(ref_tet_eqedl[:, 1] - ref_tet_eqedl[:, 2])
        @test random_tet_edl ≈ 6 * a_eqedl atol = 1.0e-12

        @test a_eqedl >= a_eqvol
    end
end

@testset "size_field_reference" begin
    # The size field prescribes the target edge length of the regular reference
    # element, combined with the volume criterion via max(size field, volume).
    tet_edges(ref) = [
        norm(ref[:, 1] - ref[:, 2]), norm(ref[:, 1] - ref[:, 3]), norm(ref[:, 1] - ref[:, 4]),
        norm(ref[:, 2] - ref[:, 3]), norm(ref[:, 2] - ref[:, 4]), norm(ref[:, 3] - ref[:, 4]),
    ]

    # A regular tet of edge ~1: a large constant size field dominates the volume
    # criterion, so the reference is a regular tet of that prescribed edge length.
    big = Norma.create_size_field("size field", "5.0")
    ref_big = Norma.create_smooth_reference("size field", Norma.TETRA4, reg_tet_coords, big, 0.0)
    for e in tet_edges(ref_big)
        @test e ≈ 5.0 atol = 1.0e-12
    end

    # A tiny constant size field is dominated by the volume criterion (max), so
    # the reference matches the equal-volume reference instead.
    tiny = Norma.create_size_field("size field", "1.0e-6")
    ref_tiny = Norma.create_smooth_reference("size field", Norma.TETRA4, reg_tet_coords, tiny, 0.0)
    ref_vol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, reg_tet_coords)
    @test norm(ref_tiny[:, 1] - ref_tiny[:, 2]) ≈ norm(ref_vol[:, 1] - ref_vol[:, 2]) atol = 1.0e-12

    # A spatially varying field is evaluated at the element reference centroid.
    sfx = Norma.create_size_field("size field", "10.0*x")
    for shift in (1.0, 2.5, 4.0)
        coords = reg_tet_coords .+ [shift; 0.0; 0.0]
        centroid_x = sum(coords[1, :]) / 4
        ref = Norma.create_smooth_reference("size field", Norma.TETRA4, coords, sfx, 0.0)
        for e in tet_edges(ref)
            @test e ≈ 10.0 * centroid_x atol = 1.0e-10
        end
    end

    # A time-varying field uses the supplied time argument.
    sft = Norma.create_size_field("size field", "2.0 + t")
    for time in (0.0, 1.0, 3.0)
        ref = Norma.create_smooth_reference("size field", Norma.TETRA4, reg_tet_coords, sft, time)
        @test norm(ref[:, 1] - ref[:, 2]) ≈ (2.0 + time) atol = 1.0e-12
    end

    # Non-size-field modes compile no size field.
    @test Norma.create_size_field("max", nothing) === nothing
    @test Norma.create_size_field("equal volume", "1.0") === nothing

    # Selecting the size-field mode without an expression is an error.
    @test_throws Exception Norma.create_size_field("size field", nothing)

    # A non-positive field yields an invalid reference and must abort.
    bad = Norma.create_size_field("size field", "x - 100.0")
    @test_throws Exception Norma.create_smooth_reference("size field", Norma.TETRA4, reg_tet_coords, bad, 0.0)
end

    big = Norma.create_size_field("size field restricted", "5.0")
    ref_big = Norma.create_smooth_reference("size field", Norma.TETRA4, reg_tet_coords, big, 0.0)
    for e in tet_edges(ref_big)
        @test e ≈ 5.0 atol = 1.0e-12
    end

    # A tiny constant size field is dominated by the volume criterion (max), so
    # the reference matches the equal-volume reference instead.
    tiny = Norma.create_size_field("size field restricted", "1.0e-6")
    ref_tiny = Norma.create_smooth_reference("size field restricted", Norma.TETRA4, reg_tet_coords, tiny, 0.0)
    ref_vol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, reg_tet_coords)
    @test norm(ref_tiny[:, 1] - ref_tiny[:, 2]) ≈ norm(ref_vol[:, 1] - ref_vol[:, 2]) atol = 1.0e-12

    # A spatially varying field is evaluated at the element reference centroid.
    sfx = Norma.create_size_field("size field restricted", "10.0*x")
    for shift in (1.0, 2.5, 4.0)
        coords = reg_tet_coords .+ [shift; 0.0; 0.0]
        centroid_x = sum(coords[1, :]) / 4
        ref = Norma.create_smooth_reference("size field restricted", Norma.TETRA4, coords, sfx, 0.0)
        for e in tet_edges(ref)
            @test e ≈ 10.0 * centroid_x atol = 1.0e-10
        end
    end

    # A time-varying field uses the supplied time argument.
    sft = Norma.create_size_field("size field restricted", "2.0 + t")
    for time in (0.0, 1.0, 3.0)
        ref = Norma.create_smooth_reference("size field restricted", Norma.TETRA4, reg_tet_coords, sft, time)
        @test norm(ref[:, 1] - ref[:, 2]) ≈ (2.0 + time) atol = 1.0e-12
    end


base_params = Dict{String,Any}(
    "type" => "single",
    "model" => Dict{String,Any}(
        "material" => Dict{String,Any}(
            "elastic" => Dict{String,Any}("model" => "seth-hill", "m" => 2, "n" => 2, "density" => 1.0),
            "blocks" => Dict{String,Any}("block" => "elastic"),
        ),
        "type" => "mesh smoothing",
    ),
    "Exodus output interval" => 0,
    "CSV output interval" => 0,
    "time integrator" =>
        Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0e-1),
    "solver" => Dict{String,Any}(
        "step" => "steepest descent",
        "type" => "steepest descent",
        "minimum iterations" => 1,
        "maximum iterations" => 64,
        "absolute tolerance" => 1.0e-8,
        "relative tolerance" => 1.0e-12,
        "step length" => 1.0e-3,
        "use line search" => true,
        "line search backtrack factor" => 0.5,
        "line search decrease factor" => 1.0e-04,
        "line search maximum iterations" => 16,
    ),
)

input_mesh_file = "tet_smoothing.g"
output_mesh_file = "tet_smoothing.e"

@testset "tet_smoothing" begin
    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by perturbating the xy coordinates of the top vertex,
    # corresponding to a pure shear deformation.
    # The resulting tetrahedron's volume is unchanged, so the smoothing
    # procedure should find the original coordinates of the top vertex when the
    # base is fixed, using a deviatoric energy term only.
    a = 1
    t = sqrt(3) / 2 * a
    h = sqrt(2.0 / 3.0) * a
    top_xy_disp = Random.rand(2) * a * 0.1
    tet_coords = reg_tet_coords + [
        0.0 0.0 0.0 top_xy_disp[1]
        0.0 0.0 0.0 top_xy_disp[2]
        0.0 0.0 0.0 0.0
    ]

    node_sets = Dict{String,Vector}("base" => tet_base)

    num_dim, num_nodes = size(tet_coords)
    num_elems = size(tet_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tet_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end

        tet_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tet_init)
        write_coordinates(tet_exo, tet_coords)
        write_block(tet_exo, 1, "TETRA4", Matrix{Int32}(tet_conn))
        write_name(tet_exo, Block, 1, "block")
        for (i, (node_set_name, tet_base)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(tet_base))
            write_set(tet_exo, node_set)
            write_name(tet_exo, node_set, node_set_name)
        end

        close(tet_exo)

        tet_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "shear_tet_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}("node set" => "base", "component" => "x", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "y", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "z", "function" => "0.0"),
                    ],
                ),
            ),
        )
        tet_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tet_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tet_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tet_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tet_coords - tet_coords) atol = a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end
    end

    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by perturbating the z coordinates of the top vertex,
    # corresponding to a uniaxial z deformation. An additional uniaxial x displacement
    # is applied to the top and front vertices to preserve the volume of the tetrahedron.
    top_z_disp = Random.rand() * a * 0.1
    front_x_disp = t * h / (h + top_z_disp) - t
    top_x_disp = front_x_disp / 3
    tet_coords = reg_tet_coords + [
        0.0 0.0 front_x_disp top_x_disp
        0.0 0.0 0.0 0
        0.0 0.0 0.0 top_z_disp
    ]

    node_sets = Dict{String,Vector}("base" => tet_base, "back edge" => tet_back_edge)

    num_dim, num_nodes = size(tet_coords)
    num_elems = size(tet_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tet_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end

        tet_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tet_init)
        write_coordinates(tet_exo, tet_coords)
        write_block(tet_exo, 1, "TETRA4", Matrix{Int32}(tet_conn))
        write_name(tet_exo, Block, 1, "block")
        for (i, (node_set_name, tet_base)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(tet_base))
            write_set(tet_exo, node_set)
            write_name(tet_exo, node_set, node_set_name)
        end

        close(tet_exo)

        tet_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "bulk_tet_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}("node set" => "back edge", "component" => "x", "function" => "0.0"),
                        Dict{String,Any}("node set" => "back edge", "component" => "y", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "z", "function" => "0.0"),
                    ],
                ),
            ),
        )
        tet_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tet_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tet_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tet_test_params)

        @test sim.integrator.displacement ≈ [0; 0; 0; 0; 0; 0; -front_x_disp; 0; 0; -top_x_disp; 0; -top_z_disp] atol =
            a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end
    end

    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by applying equal uniaxial xy deformation gradients
    # to all vertices and then adjusting the z coordinate of the top vertex to
    # preserve the total edge length of the tetrahedron.
    # Using the average edge length smoothing method, the smoothing procedure
    # should retrieve the original tetrahedron coordinates.

    xy_strain_init = Random.rand() * 0.1
    tet_coords =
        diagm([1+xy_strain_init, 1+xy_strain_init, 1]) * (reg_tet_coords .- reg_tet_coords[:, 1]) .+
        reg_tet_coords[:, 1]
    total_edge_l =
        norm(tet_coords[:, 1] - tet_coords[:, 2]) +
        norm(tet_coords[:, 1] - tet_coords[:, 3]) +
        norm(tet_coords[:, 1] - tet_coords[:, 4]) +
        norm(tet_coords[:, 2] - tet_coords[:, 3]) +
        norm(tet_coords[:, 2] - tet_coords[:, 4]) +
        norm(tet_coords[:, 3] - tet_coords[:, 4])
    xyz_scale = 6*a / total_edge_l
    tet_coords = xyz_scale * (tet_coords .- tet_coords[:, 1]) .+ tet_coords[:, 1]

    node_sets = Dict{String,Vector}("base" => tet_base, "back edge" => tet_back_edge, "first node" => tet_first_node)

    num_dim, num_nodes = size(tet_coords)
    num_elems = size(tet_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tet_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end

        tet_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tet_init)
        write_coordinates(tet_exo, tet_coords)
        write_block(tet_exo, 1, "TETRA4", Matrix{Int32}(tet_conn))
        write_name(tet_exo, Block, 1, "block")
        for (i, (node_set_name, tet_base)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(tet_base))
            write_set(tet_exo, node_set)
            write_name(tet_exo, node_set, node_set_name)
        end

        close(tet_exo)

        tet_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "bulk_tet_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}("node set" => "back edge", "component" => "x", "function" => "0.0"),
                        Dict{String,Any}("node set" => "first node", "component" => "y", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "z", "function" => "0.0"),
                    ],
                ),
            ),
        )
        tet_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tet_test_params["model"]["material"]["elastic"]["bulk modulus"] = 1.0
        tet_test_params["model"]["smooth reference"] = "average edge length"

        sim = Norma.run(tet_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tet_coords - tet_coords) atol = a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file; force=true)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file; force=true)
        end
    end
end

@testset "tet_smoothing_lbfgs" begin
    # Same pure-shear analytical case as the first sub-test of "tet_smoothing",
    # but solved with the L-BFGS step instead of steepest descent.  Both must
    # recover the regular tetrahedron (the exact energy minimizer), so this is a
    # correctness check on the quasi-Newton step against a known solution.
    a = 1
    top_xy_disp = Random.rand(2) * a * 0.1
    tet_coords = reg_tet_coords + [
        0.0 0.0 0.0 top_xy_disp[1]
        0.0 0.0 0.0 top_xy_disp[2]
        0.0 0.0 0.0 0.0
    ]
    node_sets = Dict{String,Vector}("base" => tet_base)
    num_dim, num_nodes = size(tet_coords)
    num_elems = size(tet_conn, 2)

    try
        tet_init = Initialization{Int32}(
            Int32(num_dim), Int32(num_nodes), Int32(num_elems), Int32(1), Int32(length(node_sets)), Int32(0)
        )
        rm(input_mesh_file; force=true)
        rm(output_mesh_file; force=true)
        tet_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tet_init)
        write_coordinates(tet_exo, tet_coords)
        write_block(tet_exo, 1, "TETRA4", Matrix{Int32}(tet_conn))
        write_name(tet_exo, Block, 1, "block")
        for (i, (node_set_name, ns)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(ns))
            write_set(tet_exo, node_set)
            write_name(tet_exo, node_set, node_set_name)
        end
        close(tet_exo)

        lbfgs_params = deepcopy(base_params)
        lbfgs_params["solver"]["step"] = "lbfgs"
        lbfgs_params["solver"]["memory"] = 10
        lbfgs_params["solver"]["step length"] = 1.0e-3
        tet_test_params = merge(
            lbfgs_params,
            Dict{String,Any}(
                "name" => "shear_tet_smoothing_lbfgs",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}("node set" => "base", "component" => "x", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "y", "function" => "0.0"),
                        Dict{String,Any}("node set" => "base", "component" => "z", "function" => "0.0"),
                    ],
                ),
            ),
        )
        tet_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tet_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tet_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tet_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tet_coords - tet_coords) atol = a*1.0e-6
    finally
        rm(input_mesh_file; force=true)
        rm(output_mesh_file; force=true)
    end
end

@testset "awful_cube_multielement" begin
    # Multi-element regression on the real, badly distorted awful-cube mesh
    # (6955 tets, slivers along the boundary).  Unlike the single-tet cases
    # above there is no closed-form minimizer, so we check the two robust
    # invariants of a working smoother: the total pseudo-strain-energy (the
    # global distortion measure being minimized) drops sharply, and the mesh
    # stays valid (every det(F) > 0, i.e. no element inverts — surfaced by
    # model.failed).  We run a single pseudo-time step with a small iteration
    # budget for both the steepest-descent and L-BFGS steps, and also confirm
    # L-BFGS reaches a lower energy than SD for the same budget.
    mesh = joinpath(@__DIR__, "..", "examples", "ems", "awful-cube", "awful-cube.g")
    out = joinpath(@__DIR__, "awful_cube_regression.e")

    function smooth_energy(step)
        params = Dict{String,Any}(
            "name" => "awful_cube_regression",
            "type" => "single",
            "input mesh file" => mesh,
            "output mesh file" => out,
            "Exodus output interval" => 0,
            "CSV output interval" => 0,
            "model" => Dict{String,Any}(
                "type" => "mesh smoothing",
                "smooth reference" => "max",
                "material" => Dict{String,Any}(
                    "blocks" => Dict{String,Any}("awful" => "elastic"),
                    "elastic" => Dict{String,Any}(
                        "model" => "seth-hill", "m" => 2, "n" => 2,
                        "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                    ),
                ),
            ),
            "time integrator" => Dict{String,Any}(
                "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0
            ),
            "boundary conditions" => Dict{String,Any}(
                "Dirichlet" => [
                    Dict{String,Any}("node set" => "nsx-", "component" => "x", "function" => "0.0"),
                    Dict{String,Any}("node set" => "nsy-", "component" => "y", "function" => "0.0"),
                    Dict{String,Any}("node set" => "nsz-", "component" => "z", "function" => "0.0"),
                    Dict{String,Any}("node set" => "nsx+", "component" => "x", "function" => "0.0"),
                    Dict{String,Any}("node set" => "nsy+", "component" => "y", "function" => "0.0"),
                    Dict{String,Any}("node set" => "nsz+", "component" => "z", "function" => "0.0"),
                ],
            ),
            "solver" => Dict{String,Any}(
                "type" => "steepest descent", "step" => step,
                "minimum iterations" => 1, "maximum iterations" => 30,
                "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
                "step length" => 1.0e-3, "use line search" => true,
                "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
                "line search maximum iterations" => 16,
            ),
        )
        step == "lbfgs" && (params["solver"]["memory"] = 10)
        rm(out; force=true)  # Exodus refuses to overwrite an existing output mesh
        sim = Norma.create_simulation(params)
        Norma.initialize(sim)
        Norma.evaluate(sim.integrator, sim.solver, sim.model)
        e0 = sim.model.strain_energy        # energy of the original mesh
        Norma.evolve(sim)
        ef = sim.model.strain_energy        # energy after smoothing
        Norma.finalize_writing(sim)         # close the output Exodus handle
        return e0, ef, sim.model.failed
    end

    try
        e0_sd, ef_sd, failed_sd = smooth_energy("steepest descent")
        @test !failed_sd                    # no inverted elements: min det(F) > 0
        @test ef_sd < 0.05 * e0_sd          # SD reduces distortion energy by >20x

        e0_lb, ef_lb, failed_lb = smooth_energy("lbfgs")
        @test !failed_lb                    # no inverted elements: min det(F) > 0
        @test ef_lb < 1.0e-3 * e0_lb        # L-BFGS reduces it by >1000x
        @test ef_lb < ef_sd                 # and outperforms SD at equal budget
    finally
        rm(out; force=true)
    end
end

# Smooth the perturbed tet cylinder (radius 1, height 2) with its lateral
# surface held to x^2 + y^2 = 1 and its end caps to z = ∓1 by `Surface` boundary
# conditions in the given enforcement mode.  Returns the lateral surface residual
# max|x^2 + y^2 - 1| before and after smoothing, the pseudo-energy before/after,
# and the failure flag.  Every boundary node starts off its surface (the mesh is
# jittered), so enforcement must both pull it back and let it slide within.
function _smooth_surface_cylinder(enforcement::String)
    mesh = joinpath(@__DIR__, "..", "examples", "ems", "cylinder", "cylinder-awful.g")
    out = joinpath(@__DIR__, "surface_cylinder_$(enforcement).e")
    surface(ss, f) = begin
        d = Dict{String,Any}("side set" => ss, "function" => f, "enforcement" => enforcement)
        enforcement == "penalty" && (d["penalty"] = 1.0e6)
        d
    end
    params = Dict{String,Any}(
        "name" => "surface_cylinder_$(enforcement)", "type" => "single",
        "input mesh file" => mesh, "output mesh file" => out,
        "Exodus output interval" => 0, "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing", "smooth reference" => "max",
            "material" => Dict{String,Any}(
                "blocks" => Dict{String,Any}("cylinder" => "elastic"),
                "elastic" => Dict{String,Any}(
                    "model" => "seth-hill", "m" => 2, "n" => 2,
                    "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                ),
            ),
        ),
        "time integrator" => Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 5.0, "time step" => 1.0
        ),
        "boundary conditions" => Dict{String,Any}(
            "Surface" => [
                surface("lateral", "x^2 + y^2 - 1.0"),
                surface("top", "z + 1.0"),
                surface("bottom", "z - 1.0"),
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent", "step" => "steepest descent",
            "minimum iterations" => 1, "maximum iterations" => 64,
            "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
            "step length" => 1.0e-3, "use line search" => true,
            "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
            "line search maximum iterations" => 16,
        ),
    )
    try
        rm(out; force=true)
        sim = Norma.create_simulation(params)
        # Read the constrained node list straight off the lateral Surface BC —
        # this also checks the BC was created and resolved its side set's nodes.
        lateral_bc = nothing
        for bc in sim.model.boundary_conditions
            if bc isa Norma.SolidMechanicsSurfaceBoundaryCondition && bc.name == "lateral"
                lateral_bc = bc
            end
        end
        lat = lateral_bc === nothing ? Int64[] : lateral_bc.node_indices
        surf(x, y) = x^2 + y^2 - 1.0
        ref = sim.model.reference
        g0 = maximum(abs(surf(ref[1, n], ref[2, n])) for n in lat)  # off-surface start
        Norma.initialize(sim)
        Norma.evaluate(sim.integrator, sim.solver, sim.model)
        e0 = sim.model.strain_energy
        Norma.evolve(sim)
        ef = sim.model.strain_energy
        cur = sim.model.reference .+ sim.model.displacement
        gf = maximum(abs(surf(cur[1, n], cur[2, n])) for n in lat)
        Norma.finalize_writing(sim)
        return (bc_found=lateral_bc !== nothing, g0=g0, gf=gf, e0=e0, ef=ef, failed=sim.model.failed)
    finally
        rm(out; force=true)
    end
end

@testset "surface_constraint_cylinder" begin
    # Energetic mesh smoothing on a general (non-box) geometry: analytic
    # level-set surface constraints on a curved boundary no Cartesian component
    # can express.  Both enforcement modes are exercised on the same case, so
    # this is also the penalty-vs-exact comparison: exact holds the surface to
    # the return-to-surface tolerance (near machine precision), independent of
    # any penalty parameter, while penalty leaves an O(1/κ) residual.
    exact = _smooth_surface_cylinder("exact")
    penalty = _smooth_surface_cylinder("penalty")

    for r in (exact, penalty)
        @test r.bc_found          # the Surface BC was created and found its side set
        @test !r.failed           # no inverted elements: min det(F) > 0
        @test r.g0 > 1.0e-2       # the boundary really started off the surface
        @test r.ef < r.e0         # and the mesh is smoothed
    end
    # Exact enforcement is surface-exact up to the return-to-surface tolerance,
    @test exact.gf < 1.0e-8
    # far tighter than the penalty residual, which scales as 1/κ,
    @test penalty.gf < 5.0e-2
    @test exact.gf < 1.0e-4 * penalty.gf
end

@testset "surface_constraint_reduces_to_cube" begin
    # Reduce-to-cube regression (design note test 1): on axis-aligned faces the
    # exact surface roller must reproduce the Dirichlet fixed-component
    # smoothing.  A linear level set g = coord - c has a Cartesian-basis
    # gradient, so the roller pins that component and frees the other two —
    # exactly what a component Dirichlet BC does; an edge/corner node accumulates
    # two/three level sets and loses the matching DOFs, as its two/three
    # component Dirichlets would.  We smooth the same distorted cube two ways
    # (six component Dirichlets vs six linear-level-set exact Surfaces) and
    # require the meshes to match to solver precision.  (Over very many steps the
    # two code paths drift as near-boundary slivers amplify floating-point
    # differences; a few steps keeps the comparison in the regime where the two
    # are mathematically identical.)
    mesh = joinpath(@__DIR__, "..", "examples", "ems", "cube", "cube.g")
    function smooth_cube(bcs, tag)
        out = joinpath(@__DIR__, "reduce_cube_$(tag).e")
        params = Dict{String,Any}(
            "name" => "reduce_cube_$(tag)", "type" => "single",
            "input mesh file" => mesh, "output mesh file" => out,
            "Exodus output interval" => 0, "CSV output interval" => 0,
            "model" => Dict{String,Any}(
                "type" => "mesh smoothing", "smooth reference" => "max",
                "material" => Dict{String,Any}(
                    "blocks" => Dict{String,Any}("cube" => "elastic"),
                    "elastic" => Dict{String,Any}(
                        "model" => "seth-hill", "m" => 2, "n" => 2,
                        "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                    ),
                ),
            ),
            "time integrator" => Dict{String,Any}(
                "type" => "quasi static", "initial time" => 0.0, "final time" => 0.03, "time step" => 1.0e-2
            ),
            "boundary conditions" => bcs,
            "solver" => Dict{String,Any}(
                "type" => "steepest descent", "step" => "steepest descent",
                "step length" => 5.0e-4, "minimum iterations" => 1, "maximum iterations" => 16,
                "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
                "use line search" => true, "line search backtrack factor" => 0.5,
                "line search decrease factor" => 1.0e-4, "line search maximum iterations" => 16,
            ),
        )
        try
            rm(out; force=true)
            sim = Norma.create_simulation(params)
            Norma.initialize(sim)
            Norma.evolve(sim)
            u = copy(sim.model.displacement)
            failed = sim.model.failed
            Norma.finalize_writing(sim)
            return u, failed
        finally
            rm(out; force=true)
        end
    end

    dirichlet = Dict{String,Any}(
        "Dirichlet" => [
            Dict{String,Any}("node set" => "nsx-", "component" => "x", "function" => "0.0"),
            Dict{String,Any}("node set" => "nsx+", "component" => "x", "function" => "0.0"),
            Dict{String,Any}("node set" => "nsy-", "component" => "y", "function" => "0.0"),
            Dict{String,Any}("node set" => "nsy+", "component" => "y", "function" => "0.0"),
            Dict{String,Any}("node set" => "nsz-", "component" => "z", "function" => "0.0"),
            Dict{String,Any}("node set" => "nsz+", "component" => "z", "function" => "0.0"),
        ],
    )
    surface = Dict{String,Any}(
        "Surface" => [
            Dict{String,Any}("side set" => "ssx-", "function" => "x + 0.5"),
            Dict{String,Any}("side set" => "ssx+", "function" => "x - 0.5"),
            Dict{String,Any}("side set" => "ssy-", "function" => "y + 0.5"),
            Dict{String,Any}("side set" => "ssy+", "function" => "y - 0.5"),
            Dict{String,Any}("side set" => "ssz-", "function" => "z + 0.5"),
            Dict{String,Any}("side set" => "ssz+", "function" => "z - 0.5"),
        ],
    )
    uD, failedD = smooth_cube(dirichlet, "dirichlet")
    uS, failedS = smooth_cube(surface, "surface")
    @test !failedD
    @test !failedS
    @test norm(uD) > 1.0e-3                    # the smoothing actually moved nodes
    @test norm(uD - uS) ≤ 1.0e-9 * norm(uD)    # exact roller ≡ Dirichlet on planar faces
end

# Structured cube of n^3 hex cells, each split into 6 tets (Kuhn split on the
# v0-v6 diagonal).  Returns coordinates, connectivity (num_elems x 4), and the
# boundary / interior node-id lists.
function _make_tet_grid(n)
    idx(i, j, k) = i + (n + 1) * j + (n + 1)^2 * k + 1
    coords = zeros(3, (n + 1)^3)
    boundary = Int[]
    interior = Int[]
    h = 1.0 / n
    for k in 0:n, j in 0:n, i in 0:n
        p = idx(i, j, k)
        coords[:, p] = [i * h, j * h, k * h]
        (i in (0, n) || j in (0, n) || k in (0, n)) ? push!(boundary, p) : push!(interior, p)
    end
    tets = NTuple{4,Int}[]
    for k in 0:n-1, j in 0:n-1, i in 0:n-1
        v = (idx(i, j, k), idx(i+1, j, k), idx(i+1, j+1, k), idx(i, j+1, k),
             idx(i, j, k+1), idx(i+1, j, k+1), idx(i+1, j+1, k+1), idx(i, j+1, k+1))
        for t in ((1,2,3,7), (1,3,4,7), (1,4,8,7), (1,8,5,7), (1,5,6,7), (1,6,2,7))
            push!(tets, (v[t[1]], v[t[2]], v[t[3]], v[t[4]]))
        end
    end
    conn = Matrix{Int}(undef, length(tets), 4)
    for (e, t) in enumerate(tets)
        conn[e, :] .= collect(t)
    end
    return coords, conn, boundary, interior, h
end

# Flip node order so every tet has positive signed volume (the smoothing ideal
# reference is positively oriented, so this keeps det(F) > 0 at the start).
function _orient_tets!(coords, conn)
    for e in 1:size(conn, 1)
        p = coords[:, conn[e, :]]
        if dot(p[:,2]-p[:,1], cross(p[:,3]-p[:,1], p[:,4]-p[:,1])) < 0
            tmp = conn[e, 3]; conn[e, 3] = conn[e, 4]; conn[e, 4] = tmp
        end
    end
end

@testset "torture_mesh_robustness" begin
    # Robustness guard: on a sliver-laden mesh, the L-BFGS step must be no less
    # robust than steepest descent — i.e. it must not invert an element where SD
    # would survive.  We build a small structured tet grid, fix the whole
    # boundary, and shove the interior nodes by 30% of the cell size (a fixed
    # pseudo-random perturbation that produces strongly distorted, near-sliver
    # elements while keeping a valid start).  Both solvers must finish without a
    # non-positive Jacobian (model.failed) and reduce the pseudo-energy.
    mesh = joinpath(@__DIR__, "torture_robustness.g")
    coords, conn, boundary, interior, h = _make_tet_grid(3)
    _orient_tets!(coords, conn)
    Random.seed!(123)
    for p in interior
        coords[:, p] .+= (2 .* Random.rand(3) .- 1) .* (0.30 * h)
    end

    function smooth(step)
        out = joinpath(@__DIR__, "torture_robustness_$(step).e")
        rm(out; force=true)
        params = Dict{String,Any}(
            "name" => "torture_robustness", "type" => "single",
            "input mesh file" => mesh, "output mesh file" => out,
            "Exodus output interval" => 0, "CSV output interval" => 0,
            "model" => Dict{String,Any}(
                "type" => "mesh smoothing", "smooth reference" => "max",
                "material" => Dict{String,Any}(
                    "blocks" => Dict{String,Any}("block" => "elastic"),
                    "elastic" => Dict{String,Any}(
                        "model" => "seth-hill", "m" => 2, "n" => 2,
                        "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                    ),
                ),
            ),
            "time integrator" => Dict{String,Any}(
                "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0
            ),
            "boundary conditions" => Dict{String,Any}(
                "Dirichlet" => [
                    Dict{String,Any}("node set" => "boundary", "component" => "x", "function" => "0.0"),
                    Dict{String,Any}("node set" => "boundary", "component" => "y", "function" => "0.0"),
                    Dict{String,Any}("node set" => "boundary", "component" => "z", "function" => "0.0"),
                ],
            ),
            "solver" => Dict{String,Any}(
                "type" => "steepest descent", "step" => step,
                "minimum iterations" => 1, "maximum iterations" => 60,
                "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
                "step length" => 1.0e-3, "use line search" => true,
                "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
                "line search maximum iterations" => 16,
            ),
        )
        step == "lbfgs" && (params["solver"]["memory"] = 10)
        sim = Norma.create_simulation(params)
        Norma.initialize(sim)
        Norma.evaluate(sim.integrator, sim.solver, sim.model)
        e0 = sim.model.strain_energy
        Norma.evolve(sim)
        ef = sim.model.strain_energy
        Norma.finalize_writing(sim)
        rm(out; force=true)
        return e0, ef, sim.model.failed
    end

    try
        init = Initialization{Int32}(
            Int32(3), Int32(size(coords, 2)), Int32(size(conn, 1)), Int32(1), Int32(1), Int32(0)
        )
        rm(mesh; force=true)
        exo = ExodusDatabase{Int32,Int32,Int32,Float64}(mesh, "w", init)
        write_coordinates(exo, coords)
        write_block(exo, 1, "TETRA4", Matrix{Int32}(permutedims(conn)))
        write_name(exo, Block, 1, "block")
        ns = NodeSet(1, Vector{Int32}(boundary))
        write_set(exo, ns)
        write_name(exo, ns, "boundary")
        close(exo)

        e0_sd, ef_sd, failed_sd = smooth("steepest descent")
        e0_lb, ef_lb, failed_lb = smooth("lbfgs")

        @test !failed_sd            # SD does not invert the mesh
        @test !failed_lb            # L-BFGS does not invert it either (the guard)
        @test ef_sd < e0_sd         # both actually smooth
        @test ef_lb < e0_lb
        @test ef_lb ≤ ef_sd         # and L-BFGS is no worse than SD here
    finally
        rm(mesh; force=true)
    end
end

@testset "size_field_end_to_end" begin
    # Exercise the full smoothing pipeline with a "size field" reference on a
    # multi-element mesh.  This runs create_smooth_reference inside the threaded
    # element loop (catching any world-age/thread-safety issue with the compiled
    # field) and confirms the smoother reduces energy without inverting.
    mesh = joinpath(@__DIR__, "size_field_e2e.g")
    coords, conn, boundary, interior, h = _make_tet_grid(3)
    _orient_tets!(coords, conn)
    Random.seed!(7)
    for p in interior
        coords[:, p] .+= (2 .* Random.rand(3) .- 1) .* (0.25 * h)
    end

    out = joinpath(@__DIR__, "size_field_e2e.e")
    params = Dict{String,Any}(
        "name" => "size_field_e2e", "type" => "single",
        "input mesh file" => mesh, "output mesh file" => out,
        "Exodus output interval" => 0, "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing",
            "smooth reference" => "size field",
            "size field" => "0.3 + 0.1 * x",     # spatially varying target edge length
            "material" => Dict{String,Any}(
                "blocks" => Dict{String,Any}("block" => "elastic"),
                "elastic" => Dict{String,Any}(
                    "model" => "seth-hill", "m" => 2, "n" => 2,
                    "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                ),
            ),
        ),
        "time integrator" => Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0
        ),
        "boundary conditions" => Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => "boundary", "component" => "x", "function" => "0.0"),
                Dict{String,Any}("node set" => "boundary", "component" => "y", "function" => "0.0"),
                Dict{String,Any}("node set" => "boundary", "component" => "z", "function" => "0.0"),
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent", "step" => "lbfgs", "memory" => 10,
            "minimum iterations" => 1, "maximum iterations" => 40,
            "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
            "step length" => 1.0e-3, "use line search" => true,
            "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
            "line search maximum iterations" => 16,
        ),
    )

    try
        init = Initialization{Int32}(
            Int32(3), Int32(size(coords, 2)), Int32(size(conn, 1)), Int32(1), Int32(1), Int32(0)
        )
        rm(mesh; force=true)
        exo = ExodusDatabase{Int32,Int32,Int32,Float64}(mesh, "w", init)
        write_coordinates(exo, coords)
        write_block(exo, 1, "TETRA4", Matrix{Int32}(permutedims(conn)))
        write_name(exo, Block, 1, "block")
        ns = NodeSet(1, Vector{Int32}(boundary))
        write_set(exo, ns)
        write_name(exo, ns, "boundary")
        close(exo)

        rm(out; force=true)
        sim = Norma.create_simulation(params)
        @test sim.model.size_field !== nothing      # field compiled and stored
        Norma.initialize(sim)
        Norma.evaluate(sim.integrator, sim.solver, sim.model)
        e0 = sim.model.strain_energy
        Norma.evolve(sim)
        ef = sim.model.strain_energy
        Norma.finalize_writing(sim)

        @test !sim.model.failed                     # no inverted elements
        @test ef < e0                               # the mesh is smoothed
    finally
        rm(mesh; force=true)
        rm(out; force=true)
    end
end

@testset "size_field_exodus_output" begin
    # A size-field smoothing run must write the target edge-length field, sampled
    # at each node, as a nodal variable "size" in the Exodus output — equal to the
    # field evaluated at each node's current position (reference + displacement).
    mesh = joinpath(@__DIR__, "size_field_out.g")
    out = joinpath(@__DIR__, "size_field_out.e")
    coords, conn, boundary, interior, h = _make_tet_grid(3)
    _orient_tets!(coords, conn)
    Random.seed!(11)
    for p in interior
        coords[:, p] .+= (2 .* Random.rand(3) .- 1) .* (0.2 * h)
    end
    params = Dict{String,Any}(
        "name" => "size_field_out", "type" => "single",
        "input mesh file" => mesh, "output mesh file" => out,
        "Exodus output interval" => 1, "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing", "smooth reference" => "size field",
            "size field" => "0.3 + 0.1 * x",
            "material" => Dict{String,Any}(
                "blocks" => Dict{String,Any}("block" => "elastic"),
                "elastic" => Dict{String,Any}(
                    "model" => "seth-hill", "m" => 2, "n" => 2,
                    "bulk modulus" => 1.0e3, "shear modulus" => 1.0e3, "density" => 1.0e3,
                ),
            ),
        ),
        "time integrator" => Dict{String,Any}(
            "type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0
        ),
        "boundary conditions" => Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => "boundary", "component" => "x", "function" => "0.0"),
                Dict{String,Any}("node set" => "boundary", "component" => "y", "function" => "0.0"),
                Dict{String,Any}("node set" => "boundary", "component" => "z", "function" => "0.0"),
            ],
        ),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent", "step" => "lbfgs", "memory" => 10,
            "minimum iterations" => 1, "maximum iterations" => 30,
            "absolute tolerance" => 1.0e-8, "relative tolerance" => 1.0e-12,
            "step length" => 1.0e-3, "use line search" => true,
            "line search backtrack factor" => 0.5, "line search decrease factor" => 1.0e-4,
            "line search maximum iterations" => 16,
        ),
    )
    try
        init = Initialization{Int32}(
            Int32(3), Int32(size(coords, 2)), Int32(size(conn, 1)), Int32(1), Int32(1), Int32(0)
        )
        rm(mesh; force=true); rm(out; force=true)
        exo = ExodusDatabase{Int32,Int32,Int32,Float64}(mesh, "w", init)
        write_coordinates(exo, coords)
        write_block(exo, 1, "TETRA4", Matrix{Int32}(permutedims(conn)))
        write_name(exo, Block, 1, "block")
        ns = NodeSet(1, Vector{Int32}(boundary))
        write_set(exo, ns); write_name(exo, ns, "boundary")
        close(exo)

        Norma.run(params)

        db = ExodusDatabase(out, "r")
        try
            names = Exodus.read_names(db, NodalVariable)
            @test "size" in names                       # the field is written
            last_step = Exodus.read_number_of_time_steps(db)
            # "size" must equal 0.3 + 0.1*x at each node's current x (= refe_x)
            for step in (1, last_step)
                sz = Vector{Float64}(Exodus.read_values(db, NodalVariable, step, "size"))
                rx = Vector{Float64}(Exodus.read_values(db, NodalVariable, step, "refe_x"))
                @test maximum(abs.(sz .- (0.3 .+ 0.1 .* rx))) < 1.0e-10
            end
        finally
            Exodus.close(db)
        end
    finally
        rm(mesh; force=true); rm(out; force=true)
    end
end
