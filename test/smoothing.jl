using LinearAlgebra
using Random
using Exodus

include("../src/Norma.jl")
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
    for _ = 1:rand_tests
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

        ref_tet_eqedl =
            Norma.create_smooth_reference("average edge length", Norma.TETRA4, random_tet);
        a_eqedl = norm(ref_tet_eqedl[:, 1] - ref_tet_eqedl[:, 2])
        @test random_tet_edl ≈ 6 * a_eqedl atol = 1.0e-12

        @test a_eqedl >= a_eqvol
    end
end

base_params = Dict{String,Any}(
    "type" => "single",
    "model" => Dict{String,Any}(
        "material" => Dict{String,Any}(
            "elastic" => Dict{String,Any}(
                "model" => "seth-hill",
                "m" => 2,
                "n" => 2,
                "density" => 1.0,
            ),
            "blocks" => Dict{String,Any}("block" => "elastic"),
        ),
        "type" => "mesh smoothing",
    ),
    "Exodus output interval" => 0,
    "CSV output interval" => 0,
    "time integrator" => Dict{String,Any}(
        "type" => "quasi static",
        "initial time" => 0.0,
        "final time" => 1.0,
        "time step" => 1.0e-1,
    ),
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
    tet_coords =
        reg_tet_coords + [
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
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
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
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
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
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end

    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by perturbating the z coordinates of the top vertex,
    # corresponding to a uniaxial z deformation. An additional uniaxial x displacement
    # is applied to the top and front vertices to preserve the volume of the tetrahedron.
    top_z_disp = Random.rand() * a * 0.1
    front_x_disp = t * h / (h + top_z_disp) - t
    top_x_disp = front_x_disp / 3
    tet_coords =
        reg_tet_coords + [
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
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
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
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        tet_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tet_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tet_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tet_test_params)

        @test sim.integrator.displacement ≈
              [0; 0; 0; 0; 0; 0; -front_x_disp; 0; 0; -top_x_disp; 0; -top_z_disp] atol =
            a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
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
    total_edge_l = norm(tet_coords[:, 1] - tet_coords[:, 2]) +
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
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
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
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "first node",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
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
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end
end
