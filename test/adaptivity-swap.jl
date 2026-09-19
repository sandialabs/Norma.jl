# Edge swaps and the acceptance test of the adaptivity loop
# (docs/notes/ems-adaptivity), on a hand-built cavity and on a mesh.
using LinearAlgebra
using Random
using Test
using Exodus

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function smoothing_model(mesh_file, block_name, output; extra=Dict{String,Any}())
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "swap",
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
            "maximum iterations" => 20,
            "relative tolerance" => 1.0e-12,
            "absolute tolerance" => 1.0e-10,
            "step length" => 1.0e-3,
            "use line search" => true,
            "line search backtrack factor" => 0.5,
            "line search decrease factor" => 1.0e-04,
            "line search maximum iterations" => 16,
        ),
    )
    merge!(params, extra)
    return Norma.create_simulation(params)
end

options = Norma.AdaptivityOptions(0.05, Inf, 1.0e-8, 2, 10, 2, true)

@testset "edge_swap_cavity" begin
    sim = smoothing_model("../examples/ems/cube/cube.g", "cube", "swap-cavity.e")
    model = sim.model
    h = 0.1
    a, b = 1, 2
    # Ring extraction and the triangulations: an edge with a ring of five
    # elements.  The ring is closed, so the edge is interior even though the
    # cavity is open elsewhere; an edge from a to a ring node lies on a
    # boundary face.
    n = 5
    angles = range(0, 2π; length=n + 1)[1:n]
    positions = hcat([0.0, 0.0, -0.9h], [0.0, 0.0, 0.9h], [[cos(t), sin(t), 0.2 * sin(2t)] * h for t in angles]...)
    ring = collect(3:(n + 2))
    conn = hcat([[a, b, ring[i], ring[mod1(i + 1, n)]] for i in 1:n]...)
    for e in 1:n
        if Norma.tetrahedron_volume(positions[:, conn[:, e]]) < 0.0
            conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
        end
    end
    topology = Norma.build_topology(positions, conn)
    @test Norma.edge_ring(topology, a, ring[1]) === nothing
    result = Norma.edge_ring(topology, a, b)
    @test result !== nothing
    elements, ring_nodes = result
    @test length(elements) == n
    @test sort(ring_nodes) == ring
    for i in 1:n
        @test sort(topology.connectivity[:, elements[i]]) == sort([a, b, ring_nodes[i], ring_nodes[mod1(i + 1, n)]])
    end
    @test length(Norma.polygon_triangulations(n)) == 5
    @test length(Norma.polygon_triangulations(7)) == 42
    # Every triangulation gives valid, positively oriented elements here, and
    # the swap picks the one of least energy or none if none decreases it.
    energies = Float64[]
    for triangles in Norma.polygon_triangulations(n)
        c = Norma.swapped_connectivity(topology, a, b, ring_nodes, triangles)
        @test c !== nothing
        @test all(Norma.tetrahedron_volume(positions[:, c[:, e]]) > 0.0 for e in 1:size(c, 2))
        push!(energies, sum(Norma.element_energies(model, 1, c, positions)))
    end
    before = sum(Norma.element_energies(model, 1, conn, positions))
    proposal = Norma.try_edge_swap(model, topology, a, b, options)
    if minimum(energies) < before * (1 - options.minimum_decrease)
        @test proposal !== nothing
        @test proposal.energy_after ≈ minimum(energies)
    else
        @test proposal === nothing
    end
    # A 3-to-2 swap: three elongated elements around a long edge piercing a
    # triangle become two near-regular ones; the swap is accepted.
    tri_angles = (0.0, 2π / 3, 4π / 3)
    tri = hcat([0.0, 0.0, -0.8h], [0.0, 0.0, 0.8h], [[0.6h * cos(t), 0.6h * sin(t), 0.0] for t in tri_angles]...)
    tri_conn = hcat([[a, b, 3, 4], [a, b, 4, 5], [a, b, 5, 3]]...)
    for e in 1:3
        if Norma.tetrahedron_volume(tri[:, tri_conn[:, e]]) < 0.0
            tri_conn[3, e], tri_conn[4, e] = tri_conn[4, e], tri_conn[3, e]
        end
    end
    topology3 = Norma.build_topology(tri, tri_conn)
    before3 = sum(Norma.element_energies(model, 1, tri_conn, tri))
    proposal3 = Norma.try_edge_swap(model, topology3, a, b, options)
    @test proposal3 !== nothing
    @test proposal3.energy_before ≈ before3
    @test proposal3.energy_after < before3
    @test size(proposal3.new_connectivity) == (4, 2)
    @test sum(Norma.element_energies(model, 1, proposal3.new_connectivity, tri)) ≈ proposal3.energy_after
    Norma.apply!(topology3, proposal3)
    Norma.compact!(topology3)
    @test Norma.num_alive_elements(topology3) == 2
    @test isempty(Norma.edge_elements(topology3, a, b))
    @test all(Norma.element_volume(topology3, e) > 0.0 for e in 1:2)
    @test all(1 <= length(v) <= 2 for v in values(topology3.faces))
    @test length(Norma.boundary_faces(topology3)) == 6
    # The same swap is refused when the required decrease is large.
    strict = Norma.AdaptivityOptions(0.05, Inf, 10.0, 2, 10, 2, true)
    @test Norma.try_edge_swap(model, Norma.build_topology(tri, tri_conn), a, b, strict) === nothing
    Norma.finalize_writing(sim)
    rm("swap-cavity.e"; force=true)
end

@testset "topology_phase_on_mesh" begin
    sim = smoothing_model("../examples/ems/awful-cube/awful-cube.g", "awful", "swap-mesh.e")
    Norma.run(sim)
    model = sim.model
    topology = Norma.build_topology(model)
    ne0 = Norma.num_alive_elements(topology)
    nn0 = Norma.num_alive_nodes(topology)
    chi0 = Norma.euler_characteristic(topology)
    boundary0 = Set(Norma.boundary_faces(topology))
    densities0 = Norma.energy_densities(model, topology)
    energy0 = sum(densities0 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
    accepted, decrease = Norma.topology_phase!(model, topology, options)
    @test accepted > 0
    @test decrease > 0.0
    @test Norma.num_alive_nodes(topology) == nn0
    @test Norma.euler_characteristic(topology) == chi0
    @test Set(Norma.boundary_faces(topology)) == boundary0
    @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
    densities1 = Norma.energy_densities(model, topology)
    energy1 = sum(densities1 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
    @test energy1 ≈ energy0 - decrease rtol = 1.0e-8
    @test maximum(densities1) <= maximum(densities0)
    # The adapted mesh is written and read back with its node sets.
    file = "swap-mesh-adapted.g"
    Norma.write_topology(topology, file)
    exo = ExodusDatabase(file, "r")
    @test length(Exodus.read_ids(exo, NodeSet)) == length(Exodus.read_ids(model.mesh, NodeSet))
    close(exo)
    rm(file; force=true)
    rm("swap-mesh.e"; force=true)
end

@testset "adaptive_run" begin
    # The coupled loop from a parameter set: two outer iterations on the
    # distorted cube, each writing the adapted mesh and smoothing it.
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "adaptive",
        "input mesh file" => "../examples/ems/awful-cube/awful-cube.g",
        "output mesh file" => "adaptive.e",
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "adaptivity" => Dict{String,Any}(
            "desired energy density" => 0.05, "adjacency layers" => 2, "maximum passes" => 5, "outer iterations" => 2
        ),
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
                "blocks" => Dict{String,Any}("awful" => "elastic"),
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
            "maximum iterations" => 20,
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
    sim = Norma.run(params)
    @test sim.model.failed == false
    @test isfile("adaptive-adapted-1.g")
    @test isfile("adaptive-adapted-1.e")
    # The final mesh is the last adapted one and its smoothing energy is
    # below that of the smoothed original mesh.
    first = smoothing_model("../examples/ems/awful-cube/awful-cube.g", "awful", "adaptive-first.e")
    Norma.run(first)
    @test sim.model.strain_energy < first.model.strain_energy
    for f in ("adaptive.e", "adaptive-first.e"), k in 1:2
        rm(f; force=true)
        rm("adaptive-adapted-$k.g"; force=true)
        rm("adaptive-adapted-$k.e"; force=true)
    end
end
