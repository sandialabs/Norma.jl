# The shape operations of the adaptivity loop accepted on the cavity minimum
# of the scaled Jacobian (`shape criterion: scaled Jacobian`), the face swap,
# the boundary edge swap, and the consistency of the operations of one pass
# (docs/notes/ems-adaptivity).
using LinearAlgebra
using Random
using Test
using Exodus
using StaticArrays

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function quality_model(mesh_file, block_name, output; size_field="0.214", surfaces=false)
    model = Dict{String,Any}(
        "type" => "mesh smoothing",
        "smooth reference" => "size field unrestricted",
        "size field" => size_field,
        "material" => Dict{String,Any}(
            "elastic" => Dict{String,Any}(
                "model" => "seth-hill",
                "m" => 2,
                "n" => 2,
                "bulk modulus" => 1.0,
                "shear modulus" => 1.0,
                "density" => 1.0,
            ),
            "blocks" => Dict{String,Any}(block_name => "elastic"),
        ),
    )
    boundary = if surfaces
        Dict{String,Any}(
            "Surface" => [
                Dict{String,Any}("side set" => "ss$(axis)$(sign)", "function" => "$axis $(sign == "-" ? "+" : "-") 1.0")
                for axis in ("x", "y", "z") for sign in ("-", "+")
            ],
        )
    else
        Dict{String,Any}(
            "Dirichlet" => [
                Dict{String,Any}("node set" => ns, "component" => c, "function" => "0.0") for
                (ns, c) in (("nsx-", "x"), ("nsx+", "x"), ("nsy-", "y"), ("nsy+", "y"), ("nsz-", "z"), ("nsz+", "z"))
            ],
        )
    end
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "quality",
        "input mesh file" => mesh_file,
        "output mesh file" => output,
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => model,
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "boundary conditions" => boundary,
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
            "energy stagnation window" => 5,
            "energy stagnation tolerance" => 1.0e-03,
        ),
    )
    return Norma.create_simulation(params)
end

# Options with the energy test made unattainable (margin 10), so that only
# the scaled Jacobian criterion can accept a shape operation.
function by_quality(; kwargs...)
    return Norma.AdaptivityOptions(0.05, Inf, 0.0, 10.0, 2, 10, 2, true, true, true; shape_by_quality=true, kwargs...)
end
by_energy(; kwargs...) = Norma.AdaptivityOptions(0.05, Inf, 0.0, 1.0e-8, 2, 10, 2, true, true, true; kwargs...)

function quality_closed_boundary(topology)
    counts = Dict{Tuple{Int,Int},Int}()
    for face in Norma.boundary_faces(topology)
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            edge = Norma.sorted_edge(face[i], face[j])
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    return all(v == 2 for v in values(counts))
end

@testset "quality_options" begin
    o = Norma.AdaptivityOptions(Dict{String,Any}())
    @test !o.shape_by_quality && !o.face_swaps && !o.boundary_swaps
    o = Norma.AdaptivityOptions(
        Dict{String,Any}(
            "shape criterion" => "scaled Jacobian",
            "face swaps" => true,
            "boundary swaps" => true,
            "boundary swap angle" => 5.0,
        ),
    )
    @test o.shape_by_quality && o.face_swaps && o.boundary_swaps
    @test o.boundary_swap_angle ≈ deg2rad(5.0)
    @test Norma.AdaptivityOptions(Dict{String,Any}("boundary swaps" => true)).boundary_swap_angle ≈ deg2rad(20.0)
    @test Norma.AdaptivityOptions(Dict{String,Any}()).desired_quality == 0.9
    @test !Norma.AdaptivityOptions(Dict{String,Any}()).size_first
    @test Norma.AdaptivityOptions(Dict{String,Any}("size operations first" => true)).size_first
    @test Norma.AdaptivityOptions(Dict{String,Any}("desired scaled Jacobian" => 0.6)).desired_quality == 0.6
    @test_throws Exception Norma.AdaptivityOptions(Dict{String,Any}("shape criterion" => "quality"))
    # The ranking of a cavity: the energy of its elements, or under the
    # scaled Jacobian criterion the negative of its worst element, so that
    # the worst cavities sort first in decreasing order either way.
    positions = [0.0 1.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0; 0.0 0.0 0.0 1.0 0.05]
    conn = [1 2; 2 5; 3 3; 4 4]
    @test all(Norma.tetrahedron_volume(positions[:, conn[:, e]]) > 0.0 for e in 1:2)
    topology = Norma.build_topology(positions, conn)
    densities = [1.0, 3.0]
    energy_options = Norma.AdaptivityOptions(Dict{String,Any}())
    quality_options = Norma.AdaptivityOptions(Dict{String,Any}("shape criterion" => "scaled Jacobian"))
    @test Norma.cavity_rank(topology, densities, [1, 2], energy_options) == 4.0
    sj = Norma.scaled_jacobians(positions, conn)
    @test Norma.cavity_rank(topology, densities, [1, 2], quality_options) == -minimum(sj)
    @test Norma.cavity_rank(topology, densities, [1], quality_options) == -sj[1]
    # The minimum must rise by the relative tolerance.
    positions = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    conn = reshape([1, 2, 3, 4], 4, 1)
    topology = Norma.build_topology(positions, conn)
    better = copy(positions)
    better[:, 4] = [1 / 3, 1 / 3, sqrt(2 / 3)]
    @test Norma.raises_minimum_quality(topology, [1], conn, better)
    @test !Norma.raises_minimum_quality(topology, [1], conn, positions)
end

@testset "face_swap" begin
    sim = quality_model("../examples/ems/cube/cube.g", "cube", "quality-face.e"; size_field="0.1")
    model = sim.model
    h = 0.1
    # Two elements sharing the face (1, 2, 3), with apexes 4 and 5 on either
    # side: the swap replaces them by the three elements around the edge
    # (4, 5).  With the apexes close to the face the two elements are flat
    # and the three are better; the swap is accepted on the minimum scaled
    # Jacobian, and the face (1, 2, 3) is gone.
    r = 0.6h
    positions = hcat(
        [r, 0.0, 0.0],
        [-r / 2, r * sqrt(3) / 2, 0.0],
        [-r / 2, -r * sqrt(3) / 2, 0.0],
        [0.0, 0.0, 0.25h],
        [0.0, 0.0, -0.25h],
    )
    conn = [1 1; 2 3; 3 2; 4 5]
    topology = Norma.build_topology(positions, conn)
    @test all(Norma.tetrahedron_volume(positions[:, conn[:, e]]) > 0.0 for e in 1:2)
    face = Norma.sorted_face(1, 2, 3)
    @test length(topology.faces[face]) == 2
    before = minimum(Norma.scaled_jacobians(positions, conn))
    proposal = Norma.try_face_swap(model, topology, face, by_quality())
    @test proposal !== nothing
    @test size(proposal.new_connectivity) == (4, 3)
    @test all(count(==(n), proposal.new_connectivity) == 3 for n in (4, 5))
    @test minimum(Norma.scaled_jacobians(positions, proposal.new_connectivity)) > before
    Norma.apply!(topology, proposal)
    Norma.compact!(topology)
    @test Norma.num_alive_elements(topology) == 3
    @test !haskey(topology.faces, face)
    @test length(Norma.edge_elements(topology, 4, 5)) == 3
    @test all(1 <= length(v) <= 2 for v in values(topology.faces))
    # The inverse: the edge (4, 5) with three elements around it swaps back
    # to the two elements only when that raises the minimum, which it does
    # not here.
    @test Norma.try_edge_swap(model, topology, 4, 5, by_quality()) === nothing
    # Tall apexes: the two elements are the better ones, so the face swap is
    # refused under both criteria.
    tall = copy(positions)
    tall[3, 4], tall[3, 5] = 1.2h, -1.2h
    topology = Norma.build_topology(tall, conn)
    @test Norma.try_face_swap(model, topology, face, by_quality()) === nothing
    @test Norma.try_face_swap(model, topology, face, by_energy()) === nothing
    # Apexes on the same side of the face: the swap is invalid.
    same = copy(positions)
    same[3, 5] = 0.5h
    conn_same = copy(conn)
    conn_same[3, 2], conn_same[4, 2] = conn_same[4, 2], conn_same[3, 2]
    @test all(Norma.tetrahedron_volume(same[:, conn_same[:, e]]) > 0.0 for e in 1:2)
    topology = Norma.build_topology(same, conn_same)
    @test Norma.try_face_swap(model, topology, face, by_quality()) === nothing
    # A face whose apexes are already joined by an edge is refused: the
    # third element around (4, 5) exists.
    positions6 = hcat(positions, [-0.3h, 0.0, 0.0])
    conn6 = hcat(conn, [4, 5, 2, 6], [4, 5, 6, 3])
    for e in 3:4
        if Norma.tetrahedron_volume(positions6[:, conn6[:, e]]) < 0.0
            conn6[3, e], conn6[4, e] = conn6[4, e], conn6[3, e]
        end
    end
    topology = Norma.build_topology(positions6, conn6)
    @test haskey(topology.edges, Norma.sorted_edge(4, 5))
    @test Norma.try_face_swap(model, topology, face, by_quality()) === nothing
    Norma.finalize_writing(sim)
    rm("quality-face.e"; force=true)
end

@testset "boundary_edge_swap" begin
    sim = quality_model("../examples/ems/cube/cube.g", "cube", "quality-boundary.e"; size_field="0.1")
    model = sim.model
    h = 0.1
    # Two elements behind the flat boundary patch made of the faces
    # (a, b, p1) and (a, b, p3), which share the boundary edge (a, b) and
    # the interior node p2 below the patch.  The edge (a, b) is long and the
    # quadrilateral (a, p1, b, p3) is convex, so the swap to the edge
    # (p1, p3) is accepted and the two faces are replaced by (a, p1, p3) and
    # (b, p1, p3) in the side set.
    a, b, p1, p2, p3 = 1, 2, 3, 4, 5
    positions = hcat([-h, 0.0, 0.0], [h, 0.0, 0.0], [0.0, 0.5h, 0.0], [0.0, 0.0, -0.8h], [0.0, -0.5h, 0.0])
    conn = [a a; b b; p1 p2; p2 p3]
    for e in 1:2
        if Norma.tetrahedron_volume(positions[:, conn[:, e]]) < 0.0
            conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
        end
    end
    topology = Norma.build_topology(positions, conn)
    top = Set([Norma.sorted_face(a, b, p1), Norma.sorted_face(a, b, p3)])
    topology.side_sets[1] = copy(top)
    topology.side_set_names[1] = "top"
    topology.node_side_sets[1] = BitVector([true, true, true, false, true])
    @test Norma.is_boundary_edge(topology, a, b)
    chain = Norma.boundary_edge_chain(topology, a, b)
    @test chain !== nothing
    @test chain[2] in ([p1, p2, p3], [p3, p2, p1])
    options = by_quality(; boundary_swaps=true, boundary_swap_angle=deg2rad(20.0))
    proposal = Norma.try_boundary_edge_swap(model, topology, a, b, options)
    @test proposal !== nothing
    @test proposal.surface !== nothing
    @test proposal.surface.side_set == 1
    @test Set(proposal.surface.removed) == top
    @test Set(proposal.surface.added) == Set([Norma.sorted_face(a, p1, p3), Norma.sorted_face(b, p1, p3)])
    @test size(proposal.new_connectivity) == (4, 2)
    Norma.apply!(topology, proposal)
    Norma.compact!(topology)
    @test Norma.num_alive_elements(topology) == 2
    @test !haskey(topology.edges, Norma.sorted_edge(a, b))
    @test length(Norma.edge_elements(topology, p1, p3)) == 2
    @test topology.side_sets[1] == Set([Norma.sorted_face(a, p1, p3), Norma.sorted_face(b, p1, p3)])
    @test Set(Norma.boundary_faces(topology)) ⊇ topology.side_sets[1]
    @test quality_closed_boundary(topology)
    # The same patch folded by 30 degrees is refused at 20 degrees and
    # accepted at 45.
    folded = copy(positions)
    folded[3, p3] = 0.5h * tan(deg2rad(30.0))
    topology = Norma.build_topology(folded, conn)
    topology.side_sets[1] = copy(top)
    topology.side_set_names[1] = "top"
    topology.node_side_sets[1] = BitVector([true, true, true, false, true])
    @test Norma.try_boundary_edge_swap(model, topology, a, b, options) === nothing
    wide = by_quality(; boundary_swaps=true, boundary_swap_angle=deg2rad(45.0))
    @test Norma.try_boundary_edge_swap(model, topology, a, b, wide) !== nothing
    # Faces in different side sets are never swapped across.
    topology = Norma.build_topology(positions, conn)
    topology.side_sets[1] = Set([Norma.sorted_face(a, b, p1)])
    topology.side_sets[2] = Set([Norma.sorted_face(a, b, p3)])
    topology.side_set_names[1], topology.side_set_names[2] = "one", "two"
    topology.node_side_sets[1] = BitVector([true, true, true, false, false])
    topology.node_side_sets[2] = BitVector([true, true, false, false, true])
    @test Norma.try_boundary_edge_swap(model, topology, a, b, options) === nothing
    # Without side sets the angle alone decides, and the proposal carries no
    # side set change.
    topology = Norma.build_topology(positions, conn)
    proposal = Norma.try_boundary_edge_swap(model, topology, a, b, options)
    @test proposal !== nothing && proposal.surface === nothing
    Norma.finalize_writing(sim)
    rm("quality-boundary.e"; force=true)
end

@testset "created_entities" begin
    # An operation may not create an edge or a face that exists in the mesh
    # or was created earlier in the pass.
    positions = [0.0 1.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0; 0.0 0.0 0.0 1.0 1.0]
    conn = reshape([1, 2, 3, 4], 4, 1)
    topology = Norma.build_topology(positions, conn)
    created = Norma.CreatedEntities()
    @test Norma.faces_are_new(topology, [1], reshape([2, 3, 4, 5], 4, 1), created)
    @test !Norma.faces_are_new(topology, Int[], reshape([1, 2, 3, 5], 4, 1), created)
    Norma.record!(created, reshape([2, 3, 4, 5], 4, 1))
    @test Norma.sorted_edge(2, 5) in created.edges
    @test Norma.sorted_face(3, 4, 5) in created.faces
    @test !Norma.faces_are_new(topology, [1], reshape([1, 3, 4, 5], 4, 1), created)
end

@testset "quality_phase_on_cube" begin
    # The distorted cube with side sets and Surface conditions on its six
    # faces: one topology phase under the scaled Jacobian criterion with
    # boundary swaps raises the minimum scaled Jacobian more than the energy
    # criterion does, the boundary stays closed, covered by the side sets,
    # and on the planes, and the mesh is consistent after every operator.
    results = Dict{String,Float64}()
    for (name, quality, boundary) in (("energy", false, false), ("quality", true, true), ("size-first", true, true))
        sim = quality_model("../examples/ems/awful-cube/awful-cube.g", "awful", "quality-$name.e"; surfaces=true)
        Norma.run(sim)
        model = sim.model
        topology = Norma.build_topology(model)
        chi0 = Norma.euler_characteristic(topology)
        sj0 = minimum(Norma.scaled_jacobians(topology.positions, topology.connectivity))
        options = Norma.AdaptivityOptions(
            0.05, Inf, 0.15, 1.0e-8, 2, 4, 2, true, true, true;
            size_by_length=true,
            shape_by_quality=quality,
            boundary_swaps=boundary,
            size_first=(name == "size-first"),
        )
        accepted, _ = Norma.topology_phase!(model, topology, options)
        @test accepted > 0
        sj = Norma.scaled_jacobians(topology.positions, topology.connectivity)
        results[name] = minimum(sj)
        @test minimum(sj) > sj0
        @test Norma.euler_characteristic(topology) == chi0
        @test quality_closed_boundary(topology)
        @test all(1 <= length(v) <= 2 for v in values(topology.faces))
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
        @test Set(Norma.boundary_faces(topology)) == union(values(topology.side_sets)...)
        X = topology.positions
        for (id, faces) in topology.side_sets
            name_ss = topology.side_set_names[id]
            axis = Dict('x' => 1, 'y' => 2, 'z' => 3)[name_ss[3]]
            value = name_ss[4] == '-' ? -1.0 : 1.0
            for face in faces, n in face
                @test abs(X[axis, n] - value) < 1.0e-9
            end
        end
        rm("quality-$name.e"; force=true)
    end
    println(
        "minimum scaled Jacobian after one phase: energy criterion $(results["energy"]), ",
        "scaled Jacobian criterion with boundary swaps $(results["quality"]), ",
        "size operations first $(results["size-first"])",
    )
    @test results["quality"] > results["energy"]
    @test results["size-first"] > results["energy"]
end
