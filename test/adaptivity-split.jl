# Edge splits of the adaptivity loop (docs/notes/ems-adaptivity): the
# operation on a hand-built cavity, and topology phases with a target finer
# than the mesh on the tube (analytic surfaces) and on the cube (node sets).
using LinearAlgebra
using Random
using Test
using Exodus
using StaticArrays

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
        "name" => "split",
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

options = Norma.AdaptivityOptions(0.05, Inf, 0.0, 1.0e-8, 2, 10, 2, true, true, true)

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

@testset "split_cavity" begin
    # Four elements around an edge twice as long as the target: the split
    # bisects every element, the new node is relaxed near the midpoint, and
    # the energy decreases.
    sim = smoothing_model("../examples/ems/cube/cube.g", "cube", "split-cavity.e"; size_field="0.1")
    model = sim.model
    h = 0.1
    a, b = 1, 2
    n = 4
    angles = range(0, 2π; length=n + 1)[1:n]
    ring_points = [[0.6h * cos(t), 0.6h * sin(t), 0.02h * sin(2t)] for t in angles]
    positions = hcat([0.0, 0.0, -h], [0.0, 0.0, h], ring_points...)
    ring = collect(3:(n + 2))
    conn = hcat([[a, b, ring[i], ring[mod1(i + 1, n)]] for i in 1:n]...)
    for e in 1:n
        if Norma.tetrahedron_volume(positions[:, conn[:, e]]) < 0.0
            conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
        end
    end
    topology = Norma.build_topology(positions, conn)
    @test Norma.metric_edge_length(model, topology, a, b) ≈ 2.0
    @test (a, b) in Norma.long_edges(model, topology)
    before = sum(Norma.element_energies(model, 1, conn, positions))
    proposal = Norma.try_edge_split(model, topology, a, b, options)
    @test proposal !== nothing
    @test proposal.split !== nothing
    @test proposal.split.edge == (a, b)
    @test isempty(proposal.split.node_sets) && isempty(proposal.split.side_sets)
    @test norm(proposal.split.position - SVector(0.0, 0.0, 0.0)) < 0.2h
    @test size(proposal.new_connectivity) == (4, 2n)
    @test maximum(proposal.new_connectivity) == n + 3
    @test count(==(n + 3), proposal.new_connectivity) == 2n
    @test proposal.energy_before ≈ before
    @test proposal.energy_after < before
    Norma.apply!(topology, proposal)
    Norma.compact!(topology)
    @test Norma.num_alive_nodes(topology) == n + 3
    @test Norma.num_alive_elements(topology) == 2n
    @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:2n)
    @test all(1 <= length(v) <= 2 for v in values(topology.faces))
    @test isempty(Norma.edge_elements(topology, a, b))
    @test length(Norma.edge_elements(topology, a, n + 3)) == n
    energies = Norma.element_energies(model, 1, topology.connectivity, topology.positions)
    @test sum(energies) ≈ proposal.energy_after rtol = 1.0e-12
    # An edge whose split does not lower the energy enough is refused, and so
    # is one whose new elements fall below the geometric floor when the old
    # ones were above it.
    strict = Norma.AdaptivityOptions(0.05, Inf, 0.0, 10.0, 2, 10, 2, true, true, true)
    topology2 = Norma.build_topology(positions, conn)
    @test Norma.try_edge_split(model, topology2, a, b, strict) === nothing
    @test Norma.passes_scaled_jacobian_floor(0.5, 0.4, 0.3)
    @test Norma.passes_scaled_jacobian_floor(0.2, 0.25, 0.3)
    @test !Norma.passes_scaled_jacobian_floor(0.5, 0.25, 0.3)
    @test !Norma.passes_scaled_jacobian_floor(0.2, 0.15, 0.3)
    old_minimum = minimum(Norma.scaled_jacobians(positions, conn))
    split_positions = hcat(positions, Vector(proposal.split.position))
    new_minimum = minimum(Norma.scaled_jacobians(split_positions, proposal.new_connectivity))
    floored = Norma.AdaptivityOptions(0.05, Inf, 0.999, 1.0e-8, 2, 10, 2, true, true, true)
    @test (Norma.try_edge_split(model, topology2, a, b, floored) !== nothing) == (new_minimum >= old_minimum)
    Norma.finalize_writing(sim)
    rm("split-cavity.e"; force=true)
end

@testset "split_phase_on_meshes" begin
    # The tube with a target half its mesh size: the long edges are split,
    # the new boundary nodes lie on the analytic surfaces, the boundary stays
    # closed and covered by the side sets.  The cube with a target finer than
    # its mesh: every node on a face of the cube belongs to the node set of
    # that face, so the Dirichlet conditions hold on the adapted mesh.
    for (mesh, block, size, surfaces, output) in (
        ("../examples/ems/tube/tube.g", "tube", "0.06", true, "split-tube.e"),
        ("../examples/ems/awful-cube/awful-cube.g", "awful", "0.14", false, "split-cube.e"),
    )
        sim = smoothing_model(mesh, block, output; size_field=size, surfaces=surfaces)
        Norma.run(sim)
        model = sim.model
        topology = Norma.build_topology(model)
        ne0 = Norma.num_alive_elements(topology)
        nn0 = Norma.num_alive_nodes(topology)
        chi0 = Norma.euler_characteristic(topology)
        densities0 = Norma.energy_densities(model, topology)
        energy0 = sum(densities0 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
        # Four passes are enough to exercise every operator here and keep the
        # test short.
        phase_options = Norma.AdaptivityOptions(0.05, Inf, 0.0, 1.0e-8, 2, 4, 2, true, true, true)
        accepted, decrease = Norma.topology_phase!(model, topology, phase_options)
        @test accepted > 0
        @test decrease > 0.0
        @test Norma.num_alive_nodes(topology) > nn0
        @test Norma.num_alive_elements(topology) > ne0
        @test Norma.euler_characteristic(topology) == chi0
        @test closed_boundary(topology)
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
        @test all(1 <= length(v) <= 2 for v in values(topology.faces))
        densities1 = Norma.energy_densities(model, topology)
        energy1 = sum(densities1 .* Norma.ideal_element_volumes(model, 1, topology.connectivity, topology.positions))
        @test energy1 ≈ energy0 - decrease rtol = 1.0e-8
        X = topology.positions
        if surfaces
            @test Set(Norma.boundary_faces(topology)) == union(values(topology.side_sets)...)
            for (id, flags) in topology.node_side_sets
                name = topology.side_set_names[id]
                nodes = findall(flags)
                r = sqrt.(X[1, nodes] .^ 2 + X[2, nodes] .^ 2)
                if name == "outer"
                    @test maximum(abs.(r .- 1.0)) < 1.0e-8
                elseif name == "inner"
                    @test maximum(abs.(r .- 0.9)) < 1.0e-8
                elseif name == "bottom"
                    @test maximum(abs.(X[3, nodes] .+ 1.0)) < 1.0e-8
                elseif name == "top"
                    @test maximum(abs.(X[3, nodes] .- 1.0)) < 1.0e-8
                end
            end
        else
            for (id, flags) in topology.node_sets
                name = topology.node_set_names[id]
                axis = Dict('x' => 1, 'y' => 2, 'z' => 3)[name[3]]
                value = name[4] == '-' ? -1.0 : 1.0
                on_face = findall(abs.(X[axis, :] .- value) .< 1.0e-9)
                @test all(flags[on_face])
                @test all(abs.(X[axis, findall(flags)] .- value) .< 1.0e-9)
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
