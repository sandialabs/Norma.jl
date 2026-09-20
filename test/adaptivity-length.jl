# The size operations of the adaptivity loop accepted on the edge length in
# the prescribed target (`size criterion: length`, docs/notes/ems-adaptivity):
# the option, the band guards that keep splits and collapses from undoing each
# other, and a topology phase on the tube that reaches the band.
using LinearAlgebra
using Random
using Test
using StaticArrays

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function length_model(mesh_file, block_name, output; size_field="0.214", surfaces=false)
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
        "name" => "length",
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

function length_closed_boundary(topology)
    counts = Dict{Tuple{Int,Int},Int}()
    for face in Norma.boundary_faces(topology)
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            edge = Norma.sorted_edge(face[i], face[j])
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    return all(v == 2 for v in values(counts))
end


by_length(floor) = Norma.AdaptivityOptions(0.05, Inf, floor, 10.0, 2, 10, 2, true, true, true, true)

@testset "size_criterion_option" begin
    @test Norma.AdaptivityOptions(Dict{String,Any}()).size_by_length == false
    @test Norma.AdaptivityOptions(Dict{String,Any}("size criterion" => "energy")).size_by_length == false
    @test Norma.AdaptivityOptions(Dict{String,Any}("size criterion" => "length")).size_by_length == true
    @test_throws Exception Norma.AdaptivityOptions(Dict{String,Any}("size criterion" => "geometric"))
    @test Norma.AdaptivityOptions(0.05, Inf, 0.0, 1.0e-8, 2, 10, 2, true, true, true).size_by_length == false
    @test Norma.LENGTH_BAND[1] * Norma.LENGTH_BAND[2] ≈ 1.0
end

@testset "length_split_and_guards" begin
    sim = length_model("../examples/ems/cube/cube.g", "cube", "length-cavity.e"; size_field="0.1")
    model = sim.model
    h = 0.1
    a, b = 1, 2
    function cavity(radius, angles)
        n = length(angles)
        ring_points = [[radius * cos(t), radius * sin(t), 0.02h * sin(2t)] for t in angles]
        positions = hcat([0.0, 0.0, -h], [0.0, 0.0, h], ring_points...)
        ring = collect(3:(n + 2))
        conn = hcat([[a, b, ring[i], ring[mod1(i + 1, n)]] for i in 1:n]...)
        for e in 1:n
            if Norma.tetrahedron_volume(positions[:, conn[:, e]]) < 0.0
                conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
            end
        end
        return positions, conn
    end
    square = range(0, 2π; length=5)[1:4]
    # Without the energy test and without the guard (the combination the
    # passes never use) a split is accepted whenever it is valid; with the
    # energy test and an unattainable margin it is refused.
    unguarded = Norma.AdaptivityOptions(0.05, Inf, 0.0, 10.0, 2, 10, 2, true, true, true)
    positions, conn = cavity(0.9h, square)
    topology = Norma.build_topology(positions, conn)
    @test Norma.try_edge_split(model, topology, a, b, unguarded) === nothing
    @test Norma.try_edge_split(model, topology, a, b, unguarded; by_length=true) !== nothing
    # A well-shaped ring: the split of the edge twice as long as the target
    # is accepted on the length, and the new edges from the split node to the
    # ring are inside the band.
    proposal = Norma.try_edge_split(model, topology, a, b, by_length(0.0); by_length=true)
    @test proposal !== nothing
    @test proposal.split.edge == (a, b)
    m = proposal.split.position
    for p in 3:6
        xp = SVector{3,Float64}(positions[:, p])
        @test Norma.metric_length(model, m, xp, (a, b, p, p)) ≥ Norma.LENGTH_BAND[1]
    end
    # A flat ring, its nodes close to the edge: the split would create edges
    # shorter than the band from the new node to the ring, so the guard
    # refuses it although it is valid without the guard.
    flat_positions, flat_conn = cavity(0.3h, square)
    flat = Norma.build_topology(flat_positions, flat_conn)
    @test Norma.try_edge_split(model, flat, a, b, unguarded; by_length=true) !== nothing
    @test Norma.try_edge_split(model, flat, a, b, by_length(0.0); by_length=true) === nothing
    # A short ring edge (3, 4) whose collapse onto node 4 would create the
    # edge (4, 6) longer than the band: refused by the guard, valid without
    # it.
    wide_positions, wide_conn = cavity(0.9h, [-0.25, 0.25, π / 2 + 0.6, π + 0.4])
    wide = Norma.build_topology(wide_positions, wide_conn)
    @test Norma.metric_edge_length(model, wide, 3, 4) < Norma.LENGTH_BAND[1]
    @test Norma.metric_edge_length(model, wide, 4, 6) > Norma.LENGTH_BAND[2]
    @test Norma.try_edge_collapse(model, wide, 3, 4, unguarded; by_length=true) !== nothing
    @test Norma.try_edge_collapse(model, wide, 3, 4, by_length(0.0); by_length=true) === nothing
    Norma.finalize_writing(sim)
    rm("length-cavity.e"; force=true)
end

@testset "length_phase_on_tube" begin
    # The tube with a target of two thirds of its mesh size: under the energy
    # criterion the phase stalls with most edges still long (the shape
    # penalty of a bisection outweighs the size gain at that ratio); under
    # the length criterion the same passes bring most edges inside the band,
    # with the boundary closed, on the surfaces, and covered by the side sets.
    function band_fraction(model, topology)
        inside = 0
        for edge in keys(topology.edges)
            r = Norma.metric_edge_length(model, topology, edge[1], edge[2])
            Norma.LENGTH_BAND[1] ≤ r ≤ Norma.LENGTH_BAND[2] && (inside += 1)
        end
        return inside / length(topology.edges)
    end
    fractions = Dict{Bool,Float64}()
    counts = Dict{Bool,Int}()
    for length_criterion in (false, true)
        output = "length-tube-$length_criterion.e"
        sim = length_model("../examples/ems/tube/tube.g", "tube", output; size_field="0.08", surfaces=true)
        Norma.run(sim)
        model = sim.model
        topology = Norma.build_topology(model)
        chi0 = Norma.euler_characteristic(topology)
        before = band_fraction(model, topology)
        phase_options = Norma.AdaptivityOptions(0.05, Inf, 0.15, 1.0e-8, 2, 6, 2, true, true, true, length_criterion)
        accepted, _ = Norma.topology_phase!(model, topology, phase_options)
        @test accepted > 0
        fractions[length_criterion] = band_fraction(model, topology)
        counts[length_criterion] = Norma.num_alive_elements(topology)
        @test fractions[length_criterion] > before
        @test Norma.euler_characteristic(topology) == chi0
        @test length_closed_boundary(topology)
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
        @test Set(Norma.boundary_faces(topology)) == union(values(topology.side_sets)...)
        X = topology.positions
        for (id, flags) in topology.node_side_sets
            name = topology.side_set_names[id]
            nodes = findall(flags)
            r = sqrt.(X[1, nodes] .^ 2 + X[2, nodes] .^ 2)
            name == "outer" && @test maximum(abs.(r .- 1.0)) < 1.0e-8
            name == "inner" && @test maximum(abs.(r .- 0.9)) < 1.0e-8
            name == "bottom" && @test maximum(abs.(X[3, nodes] .+ 1.0)) < 1.0e-8
            name == "top" && @test maximum(abs.(X[3, nodes] .- 1.0)) < 1.0e-8
        end
        rm(output; force=true)
    end
    @test fractions[true] > 0.7
    @test fractions[true] > fractions[false]
    @test counts[true] > counts[false]
    println("edges inside the band: energy criterion $(fractions[false]), length criterion $(fractions[true]); ",
        "elements $(counts[false]) and $(counts[true])")
end
