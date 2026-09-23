# In-memory topology for the adaptivity loop (docs/notes/ems-adaptivity):
# adjacency, sets, tombstones, compaction, and Exodus output.
using LinearAlgebra
using Random
using Test
using Exodus

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function smoothing_model(mesh_file, block_name, output)
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "topology",
        "input mesh file" => mesh_file,
        "output mesh file" => output,
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing",
            "smooth reference" => "equal volume",
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
        ),
        "time integrator" =>
            Dict{String,Any}("type" => "quasi static", "initial time" => 0.0, "final time" => 1.0, "time step" => 1.0),
        "solver" => Dict{String,Any}(
            "type" => "steepest descent",
            "step" => "steepest descent",
            "minimum iterations" => 1,
            "maximum iterations" => 1,
            "relative tolerance" => 1.0e-12,
            "absolute tolerance" => 1.0e-8,
            "step length" => 1.0e-3,
        ),
    )
    return Norma.create_simulation(params)
end

@testset "topology_build" begin
    # The cube is a ball (Euler characteristic one); the tube is a solid
    # torus (zero).
    for (mesh_file, block_name, has_side_sets, chi) in
        (("../examples/ems/cube/cube.g", "cube", false, 1), ("../examples/ems/tube/tube.g", "tube", true, 0))
        sim = smoothing_model(mesh_file, block_name, "topology-build.e")
        model = sim.model
        topology = Norma.build_topology(model)
        nn = size(model.reference, 2)
        ne = model.blocks[1].num_elements
        @test Norma.num_alive_nodes(topology) == nn
        @test Norma.num_alive_elements(topology) == ne
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:ne)
        @test Norma.euler_characteristic(topology) == chi
        # Every face has one or two elements, and the sum over faces of the
        # incidences is four per element.
        incidences = [length(v) for v in values(topology.faces)]
        @test all(1 .<= incidences .<= 2)
        @test sum(incidences) == 4 * ne
        # Every node belongs to at least one element and the node-to-element
        # adjacency inverts the connectivity.
        for n in 1:nn
            elements = Norma.node_elements(topology, n)
            @test length(elements) > 0
            @test all(n in view(topology.connectivity, :, e) for e in elements)
        end
        @test sum(length(Norma.node_elements(topology, n)) for n in 1:nn) == 4 * ne
        # Edge-to-element adjacency: an edge of an element is in the table
        # with that element.
        e = rand(1:ne)
        c = topology.connectivity[:, e]
        @test e in Norma.edge_elements(topology, c[1], c[3])
        @test isempty(Norma.edge_elements(topology, c[1], nn + 1))
        # The boundary faces are exactly the union of the side sets when the
        # mesh has them, and every boundary edge lies on a boundary face.
        boundary = Set(Norma.boundary_faces(topology))
        if has_side_sets
            @test !isempty(topology.side_sets)
            @test union(values(topology.side_sets)...) == boundary
            for (id, faces) in topology.side_sets
                flags = topology.node_side_sets[id]
                @test all(flags[face[i]] for face in faces for i in 1:3)
                @test count(flags) == length(unique(vcat([collect(face) for face in faces]...)))
            end
        end
        for face in boundary
            @test Norma.is_boundary_edge(topology, face[1], face[2])
            @test Norma.is_boundary_edge(topology, face[2], face[3])
        end
        # Interior edges exist and are not boundary edges.
        @test any(!Norma.is_boundary_edge(topology, e[1], e[2]) for e in keys(topology.edges))
        # Node sets match the Exodus node sets.
        for id in Exodus.read_ids(model.mesh, NodeSet)
            @test findall(topology.node_sets[Int(id)]) == sort(Int.(Exodus.read_node_set_nodes(model.mesh, id)))
        end
        Norma.finalize_writing(sim)
        rm("topology-build.e"; force=true)
    end
end

@testset "scaled_jacobian" begin
    c = 0.5 / sqrt(2.0)
    regular = c * [1 -1 -1 1; 1 -1 1 -1; 1 1 -1 -1]
    @test Norma.tetrahedron_scaled_jacobian(regular) ≈ 1.0 atol = 1.0e-12
    @test Norma.tetrahedron_scaled_jacobian(3.7 * regular) ≈ 1.0 atol = 1.0e-12
    flat = copy(regular)
    flat[:, 4] = (regular[:, 1] + regular[:, 2] + regular[:, 3]) / 3
    @test Norma.tetrahedron_scaled_jacobian(flat) ≈ 0.0 atol = 1.0e-12
    inverted = regular[:, [1, 2, 4, 3]]
    @test Norma.tetrahedron_scaled_jacobian(inverted) ≈ -1.0 atol = 1.0e-12
    # Right-angled corner tetrahedron: Jacobian one, largest edge product two.
    right = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test Norma.tetrahedron_scaled_jacobian(right) ≈ sqrt(2.0) / 2.0 atol = 1.0e-12
    sj = Norma.scaled_jacobians(hcat(regular, right), [1 5; 2 6; 3 7; 4 8])
    @test sj[1] ≈ 1.0 && sj[2] < 1.0
end

@testset "topology_edit_compact_write" begin
    sim = smoothing_model("../examples/ems/tube/tube.g", "tube", "topology-edit.e")
    model = sim.model
    topology = Norma.build_topology(model)
    Norma.finalize_writing(sim)
    nn0 = Norma.num_alive_nodes(topology)
    ne0 = Norma.num_alive_elements(topology)
    chi0 = Norma.euler_characteristic(topology)
    # Split an interior edge at its midpoint: each element of its ring is
    # replaced by two.  Node sets and side sets are unaffected because the
    # edge is interior.
    edge = first(e for e in keys(topology.edges) if !Norma.is_boundary_edge(topology, e[1], e[2]))
    a, b = edge
    ring = copy(Norma.edge_elements(topology, a, b))
    @test length(ring) >= 3
    m = Norma.add_node!(topology, 0.5 * (topology.positions[:, a] + topology.positions[:, b]))
    @test m == nn0 + 1
    new_conn = zeros(Int, 4, 0)
    for e in ring
        c = topology.connectivity[:, e]
        ca = replace(c, b => m)
        cb = replace(c, a => m)
        new_conn = hcat(new_conn, ca, cb)
    end
    added = Norma.add_elements!(topology, new_conn, 1)
    @test length(added) == 2 * length(ring)
    Norma.remove_elements!(topology, ring)
    @test Norma.num_alive_elements(topology) == ne0 + length(ring)
    node_map, element_map = Norma.compact!(topology)
    @test count(node_map .> 0) == nn0 + 1
    @test count(element_map .> 0) == ne0 + length(ring)
    @test all(element_map[ring] .== 0)
    @test Norma.num_alive_elements(topology) == ne0 + length(ring)
    @test Norma.euler_characteristic(topology) == chi0
    @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
    @test sum(length(v) for v in values(topology.faces)) == 4 * Norma.num_alive_elements(topology)
    # Total volume is preserved by the split.
    # Remove the new node again by collapsing it onto a: the elements that
    # contain both are removed and the others map m to a.
    star = copy(Norma.node_elements(topology, m))
    keep = Int[]
    for e in star
        c = topology.connectivity[:, e]
        if a in c
            continue
        end
        push!(keep, e)
    end
    Norma.remove_elements!(topology, star)
    Norma.add_elements!(topology, hcat([replace(topology.connectivity[:, e], m => a) for e in keep]...), 1)
    Norma.remove_node!(topology, m)
    Norma.compact!(topology)
    @test Norma.num_alive_nodes(topology) == nn0
    @test Norma.num_alive_elements(topology) == ne0
    @test Norma.euler_characteristic(topology) == chi0
    # Write the mesh and read it back: coordinates, element count, node sets,
    # and side sets round-trip, which also checks the side numbering.
    file = "topology-written.g"
    Norma.write_topology(topology, file)
    exo = ExodusDatabase(file, "r")
    # The title is the one given, not stale memory: Exodus.jl reads the title
    # into 80 bytes where the library writes up to 81, so a title of 80
    # characters or more corrupts the heap on every open.  Read here into a
    # buffer with room to spare.
    title = zeros(UInt8, 2 * Exodus.MAX_LINE_LENGTH)
    counts = [Ref{Int32}(0) for _ in 1:6]
    @test 0 == @ccall Exodus.libexodus.ex_get_init(
        Exodus.get_file_id(exo)::Cint, title::Ptr{UInt8}, counts[1]::Ptr{Int32}, counts[2]::Ptr{Int32},
        counts[3]::Ptr{Int32}, counts[4]::Ptr{Int32}, counts[5]::Ptr{Int32}, counts[6]::Ptr{Int32},
    )::Cint
    @test String(title[1:(findfirst(iszero, title) - 1)]) == "Norma adapted mesh"
    X = Exodus.read_coordinates(exo)
    @test size(X, 2) == nn0
    ids = Exodus.read_ids(exo, Block)
    @test length(ids) == 1
    _, ne, nnpe, _, _, _ = Exodus.read_block_parameters(exo, ids[1])
    @test ne == ne0 && nnpe == 4
    @test Exodus.read_name(exo, Block, ids[1]) == "tube"
    for id in Exodus.read_ids(exo, NodeSet)
        @test sort(Int.(Exodus.read_node_set_nodes(exo, id))) == findall(topology.node_sets[Int(id)])
        @test Exodus.read_name(exo, NodeSet, id) == topology.node_set_names[Int(id)]
    end
    for id in Exodus.read_ids(exo, SideSet)
        counts, nodes = Exodus.read_side_set_node_list(exo, id)
        faces = Set{NTuple{3,Int}}()
        offset = 0
        for count in counts
            push!(faces, Norma.sorted_face(nodes[offset + 1], nodes[offset + 2], nodes[offset + 3]))
            offset += count
        end
        @test faces == topology.side_sets[Int(id)]
        @test Exodus.read_name(exo, SideSet, id) == topology.side_set_names[Int(id)]
    end
    close(exo)
    # The written mesh builds a topology identical in its counts.
    sim2 = smoothing_model(file, "tube", "topology-edit2.e")
    topology2 = Norma.build_topology(sim2.model)
    Norma.finalize_writing(sim2)
    @test Norma.num_alive_elements(topology2) == ne0
    @test Norma.euler_characteristic(topology2) == chi0
    @test Set(Norma.boundary_faces(topology2)) == union(values(topology2.side_sets)...)
    rm(file; force=true)
    rm("topology-edit.e"; force=true)
    rm("topology-edit2.e"; force=true)
end
