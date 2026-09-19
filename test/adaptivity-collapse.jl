# Edge collapses of the adaptivity loop (docs/notes/ems-adaptivity): the
# constraints on the removed node, the operation on a hand-built cavity, and
# a topology phase with swaps and collapses on a mesh, with and without a
# prescribed target.
using LinearAlgebra
using Random
using Test
using Exodus

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function smoothing_model(mesh_file, block_name, output; size_field="0.214", surfaces=false)
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
                Dict{String,Any}("side set" => "outer", "function" => "x^2 + y^2 - 1.0"),
                Dict{String,Any}("side set" => "inner", "function" => "x^2 + y^2 - 0.81"),
                Dict{String,Any}("side set" => "bottom", "function" => "z + 1.0"),
                Dict{String,Any}("side set" => "top", "function" => "z - 1.0"),
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
        "name" => "collapse",
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

options = Norma.AdaptivityOptions(0.05, Inf, 1.0e-8, 2, 10, 2, true, true)

# A closed boundary surface: every boundary edge lies in exactly two
# boundary faces.
function closed_boundary(topology)
    counts = Dict{Tuple{Int,Int},Int}()
    for face in Norma.boundary_faces(topology)
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            edge = Norma.sorted_edge(face[i], face[j])
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    return all(v == 2 for v in values(counts))
end

@testset "collapse_cavity" begin
    sim = smoothing_model("../examples/ems/cube/cube.g", "cube", "collapse-cavity.e"; size_field="0.1")
    model = sim.model
    # An interior node inside a tetrahedron, joined to its four vertices; the
    # star has four elements.  Collapsing the interior node onto a vertex
    # removes the three elements that contain both and maps the fourth onto
    # the outer tetrahedron, whose size is the target.
    h = 0.1
    c = 0.5 / sqrt(2.0)
    outer = h * c * [1 -1 -1 1; 1 -1 1 -1; 1 1 -1 -1]
    inner = 0.15 * outer[:, 1] + [0.01h, -0.005h, 0.0]
    positions = hcat(outer, inner)
    conn = hcat([[5, 2, 3, 4], [5, 1, 4, 3], [5, 1, 2, 4], [5, 1, 3, 2]]...)
    for e in 1:4
        if Norma.tetrahedron_volume(positions[:, conn[:, e]]) < 0.0
            conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
        end
    end
    topology = Norma.build_topology(positions, conn)
    @test Norma.num_alive_elements(topology) == 4
    before = sum(Norma.element_energies(model, 1, conn, positions))
    proposal = Norma.try_edge_collapse(model, topology, 5, 1, options)
    @test proposal !== nothing
    @test proposal.removed_node == 5
    @test proposal.surviving_node == 1
    @test length(proposal.old_elements) == 4
    @test size(proposal.new_connectivity) == (4, 1)
    @test sort(proposal.new_connectivity[:, 1]) == [1, 2, 3, 4]
    @test proposal.energy_before ≈ before
    @test proposal.energy_after < before
    @test proposal.energy_after ≈ 0.0 atol = 1.0e-20
    Norma.apply!(topology, proposal)
    @test topology.node_alive[5] == false
    Norma.compact!(topology)
    @test Norma.num_alive_nodes(topology) == 4
    @test Norma.num_alive_elements(topology) == 1
    # Collapsing a vertex onto the interior node is refused when it does not
    # lower the energy enough; collapsing onto a node that is not a neighbor
    # is not proposed.
    topology = Norma.build_topology(positions, conn)
    strict = Norma.AdaptivityOptions(0.05, Inf, 10.0, 2, 10, 2, true, true)
    @test Norma.try_edge_collapse(model, topology, 5, 1, strict) === nothing
    # Constraints: a node in a node set may only collapse onto a node of the
    # same set, and a boundary node only along a boundary edge onto a node of
    # the same surfaces.
    topology.node_sets[1] = falses(5)
    topology.node_sets[1][5] = true
    @test Norma.may_collapse(topology, 5, 1) == false
    topology.node_sets[1][1] = true
    @test Norma.may_collapse(topology, 5, 1) == true
    delete!(topology.node_sets, 1)
    topology.side_sets[7] = Set([Norma.sorted_face(1, 2, 3)])
    Norma.build_adjacency!(topology)
    @test Norma.may_collapse(topology, 1, 5) == false   # boundary node onto an interior node
    @test Norma.may_collapse(topology, 1, 2) == true    # along the boundary edge (1, 2)
    @test Norma.may_collapse(topology, 1, 4) == false   # onto a node off the surface
    Norma.finalize_writing(sim)
    rm("collapse-cavity.e"; force=true)
end

@testset "collapse_phase_on_meshes" begin
    # The distorted cube with the size of its own mesh: collapses act as the
    # fallback after the swaps; and the tube with a target twice its mesh
    # size: the short edges make collapses the main operation and the element
    # count falls.  In both the boundary stays a closed surface, the sets are
    # carried, and the energy decreases by the accounted amount.
    for (mesh, block, size, surfaces, output) in (
        ("../examples/ems/awful-cube/awful-cube.g", "awful", "0.214", false, "collapse-cube.e"),
        ("../examples/ems/tube/tube.g", "tube", "0.24", true, "collapse-tube.e"),
    )
        sim = smoothing_model(mesh, block, output; size_field=size, surfaces=surfaces)
        Norma.run(sim)
        model = sim.model
        topology = Norma.build_topology(model)
        ne0 = Norma.num_alive_elements(topology)
        nn0 = Norma.num_alive_nodes(topology)
        chi0 = Norma.euler_characteristic(topology)
        node_set_counts = Dict(id => count(flags) for (id, flags) in topology.node_sets)
        densities0 = Norma.energy_densities(model, topology)
        energy0 = sum(densities0 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
        accepted, decrease = Norma.topology_phase!(model, topology, options)
        @test accepted > 0
        @test decrease > 0.0
        @test Norma.euler_characteristic(topology) == chi0
        @test closed_boundary(topology)
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
        @test all(1 <= length(v) <= 2 for v in values(topology.faces))
        densities1 = Norma.energy_densities(model, topology)
        energy1 = sum(densities1 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
        @test energy1 ≈ energy0 - decrease rtol = 1.0e-8
        if surfaces
            # Short edges of the coarse target were collapsed: fewer nodes and
            # elements, and the side sets still cover the whole boundary.
            @test Norma.num_alive_nodes(topology) < nn0
            @test Norma.num_alive_elements(topology) < ne0
            @test Set(Norma.boundary_faces(topology)) == union(values(topology.side_sets)...)
        else
            # Node sets keep their nodes unless a member was collapsed onto
            # another member; corner nodes are never removed.
            for (id, flags) in topology.node_sets
                @test count(flags) <= node_set_counts[id]
            end
        end
        file = replace(output, ".e" => "-adapted.g")
        Norma.write_topology(topology, file)
        sim2 = smoothing_model(file, block, replace(output, ".e" => "-2.e"); size_field=size, surfaces=surfaces)
        topology2 = Norma.build_topology(sim2.model)
        @test Norma.num_alive_elements(topology2) == Norma.num_alive_elements(topology)
        Norma.finalize_writing(sim2)
        rm(file; force=true)
        rm(output; force=true)
        rm(replace(output, ".e" => "-2.e"); force=true)
    end
end
