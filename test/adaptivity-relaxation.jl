# The cavity relaxation of the adaptivity loop (docs/notes/ems-adaptivity):
# a swap or collapse that the energy test refuses at the fixed positions is
# judged again after the nodes of its cavity are relaxed against the new
# connectivity, interior nodes freely and boundary nodes under their Surface
# and Dirichlet conditions, with the elements around the moved nodes in the
# energy.  Checked on the awful cube: the proposals, a swap pass, the
# agreement of the evaluation ahead of the sequential loop, the constraints
# of the boundary nodes, and the loop with and without the relaxation.
using LinearAlgebra
using Random
using Test
using Exodus

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

function smoothing_model(mesh_file, block_name, output; extra=Dict{String,Any}(), size_field=nothing)
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "swap",
        "input mesh file" => mesh_file,
        "output mesh file" => output,
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => Dict{String,Any}(
            "type" => "mesh smoothing",
            "smooth reference" => size_field === nothing ? "equal volume" : "size field unrestricted",
            "size field" => size_field === nothing ? "1.0" : size_field,
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
                Dict{String,Any}("node set" => ns, "component" => c, "function" => "0.0") for (ns, c) in
                (("nsx-", "x"), ("nsx+", "x"), ("nsy-", "y"), ("nsy+", "y"), ("nsz-", "z"), ("nsz+", "z"))
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
    size_field === nothing && delete!(params["model"], "size field")
    merge!(params, extra)
    return Norma.create_simulation(params)
end

# A mesh adapted by the loop from a parameter set, and the quality of its
# final mesh.
function adapted_quality(params)
    sim = Norma.run(params)
    model = sim.model
    x = model.reference + model.displacement
    conn = Matrix{Int}(model.blocks[1].connectivity)
    sj = Norma.scaled_jacobians(x, conn)
    return (energy=model.strain_energy, elements=size(conn, 2), min_sj=minimum(sj), mean_sj=sum(sj) / length(sj))
end

@testset "cavity_relaxation" begin
    # Under `cavity relaxation` a swap refused at fixed positions is judged
    # again after the interior nodes of its cavity are relaxed against the
    # new connectivity, with the elements around them in the energy.
    sim = smoothing_model("../examples/ems/awful-cube/awful-cube.g", "awful", "swap-relaxed.e"; size_field="0.214")
    Norma.run(sim)
    model = sim.model
    fixed = Norma.AdaptivityOptions(
        Dict{String,Any}("desired energy density" => 0.0, "collapses" => false, "splits" => false)
    )
    relaxed = Norma.AdaptivityOptions(
        Dict{String,Any}(
            "desired energy density" => 0.0, "collapses" => false, "splits" => false, "cavity relaxation" => true
        ),
    )
    @test !fixed.cavity_relaxation
    @test relaxed.cavity_relaxation
    topology = Norma.build_topology(model)
    function interior(n)
        return !any(flags[n] for flags in values(topology.node_sets)) &&
               !any(flags[n] for flags in values(topology.node_side_sets))
    end
    # Interior edges refused at fixed positions and accepted after the
    # relaxation: the moved nodes are interior nodes of the cavity (the
    # boundary nodes of this model are in node sets, so none may move), the
    # energies of the proposal are those of the old and the new elements
    # together with the elements around the moved nodes, before and after
    # the move, and the new elements alone do not lower the cavity energy
    # at the fixed positions.
    found = 0
    for edge in sort(collect(keys(topology.edges)))
        Norma.is_boundary_edge(topology, edge...) && continue
        Norma.try_edge_swap(model, topology, edge..., fixed) === nothing || continue
        proposal = Norma.try_edge_swap(model, topology, edge..., relaxed)
        proposal === nothing && continue
        found += 1
        @test !isempty(proposal.moved_nodes)
        @test all(interior, proposal.moved_nodes)
        cavity_nodes = vec(topology.connectivity[:, proposal.old_elements])
        @test all(n in cavity_nodes for n in proposal.moved_nodes)
        outer = Norma.outer_elements(topology, proposal.old_elements, proposal.moved_nodes)
        @test !isempty(outer)
        @test isempty(intersect(outer, proposal.old_elements))
        x = topology.positions
        old_connectivity = topology.connectivity[:, proposal.old_elements]
        before = sum(Norma.element_energies(model, 1, hcat(old_connectivity, topology.connectivity[:, outer]), x))
        @test proposal.energy_before ≈ before
        moved = Norma.MovedPositions(x, proposal.moved_nodes, proposal.moved_positions)
        after_connectivity = hcat(proposal.new_connectivity, topology.connectivity[:, outer])
        after = sum(Norma.element_energies(model, 1, after_connectivity, moved))
        @test proposal.energy_after ≈ after
        @test after < before * (1 - relaxed.minimum_decrease)
        @test sum(Norma.element_energies(model, 1, proposal.new_connectivity, x)) >
            sum(Norma.element_energies(model, 1, old_connectivity, x)) * (1 - fixed.minimum_decrease)
        @test all(
            Norma.tetrahedron_volume(Matrix(moved[:, after_connectivity[:, e]])) > 0.0 for
            e in 1:size(after_connectivity, 2)
        )
        found ≥ 5 && break
    end
    @test found ≥ 5
    # A swap operator with the relaxation accepts at least as many swaps as
    # one without it and reaches a lower energy; it moves interior nodes
    # only, keeps the boundary and the orientation of every element, its
    # decrease is that of the total energy, and the memory of the phase
    # matches a fresh evaluation after the moves (apply_accepted!).
    results = Dict{Bool,Tuple{Int,Float64}}()
    for opts in (fixed, relaxed)
        topology = Norma.build_topology(model)
        x0 = Matrix{Float64}(topology.positions)
        boundary0 = Set(Norma.boundary_faces(topology))
        total_energy(t) =
            sum(Norma.energy_densities(model, t) .* Norma.ideal_element_volumes(model, 1, t.connectivity, t.positions))
        energy0 = total_energy(topology)
        memory = Norma.PhaseMemory()
        Norma.phase_element_values!(memory, model, topology, opts)
        result = Norma.swap_operator!(memory, model, topology, opts)
        moved_nodes = [n for n in 1:size(x0, 2) if topology.positions[:, n] != x0[:, n]]
        @test isempty(moved_nodes) == !opts.cavity_relaxation
        @test all(interior, moved_nodes)
        @test Set(Norma.boundary_faces(topology)) == boundary0
        @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
        energy1 = total_energy(topology)
        @test energy1 ≈ energy0 - result.decrease rtol = 1.0e-8
        @test memory.density == Norma.energy_densities(model, topology)
        results[opts.cavity_relaxation] = (result.accepted, energy1)
    end
    @test results[true][1] >= results[false][1]
    @test results[true][2] < results[false][2]
    # The evaluation ahead of the sequential loop agrees with the sequential
    # result under the relaxation, which also depends on the positions of
    # the nodes around the cavity (edge_swap_from_evaluation).
    topology = Norma.build_topology(model)
    memory = Norma.PhaseMemory()
    Norma.phase_element_values!(memory, model, topology, relaxed)
    edges = sort([edge for edge in keys(topology.edges) if !Norma.is_boundary_edge(topology, edge...)])
    evaluations = [Norma.evaluate_edge_swap(model, topology, a, b, relaxed) for (a, b) in edges]
    created = Norma.CreatedEntities()
    same(p, q) =
        (p === nothing && q === nothing) || (
            p !== nothing &&
            q !== nothing &&
            p.old_elements == q.old_elements &&
            p.new_connectivity == q.new_connectivity &&
            p.energy_before == q.energy_before &&
            p.energy_after == q.energy_after &&
            p.moved_nodes == q.moved_nodes &&
            p.moved_positions == q.moved_positions
        )
    agree = true
    applied = 0
    for (k, (a, b)) in enumerate(edges)
        expected = Norma.try_edge_swap(model, topology, a, b, relaxed; created)
        taken = Norma.edge_swap_from_evaluation(model, topology, a, b, relaxed, evaluations[k], created)
        agree &= same(expected, taken)
        expected === nothing && continue
        Norma.apply_accepted!(memory, model, topology, relaxed, created, expected)
        applied += 1
    end
    @test applied > 10
    @test agree
    @test !isempty(created.moved)
    rm("swap-relaxed.e"; force=true)
end

@testset "cavity_relaxation_boundary" begin
    # Under Surface conditions the relaxation moves boundary nodes along
    # their surfaces, a node of two surfaces along their line, and leaves the
    # corners; under Dirichlet conditions it keeps the fixed components.
    params = Norma.YAML.load_file("../examples/ems/awful-cube/awful-cube-adaptive.yaml")
    delete!(params, "adaptivity")
    params["input mesh file"] = "../examples/ems/awful-cube/awful-cube.g"
    params["output mesh file"] = "swap-boundary.e"
    params["name"] = "swap-boundary"
    Norma.YAML.write_file("swap-boundary.yaml", params)
    sim = Norma.run("swap-boundary.yaml")
    model = sim.model
    relaxed = Norma.AdaptivityOptions(
        Dict{String,Any}(
            "desired energy density" => 0.0, "collapses" => false, "splits" => false, "cavity relaxation" => true
        ),
    )
    topology = Norma.build_topology(model)
    x0 = Matrix{Float64}(topology.positions)
    side_sets_of(n) = [id for (id, flags) in topology.node_side_sets if flags[n]]
    # Every face of the cube carries a Surface condition, so a node on one or
    # two faces may move and a corner may not; the interior nodes have no
    # constraint.
    for n in 1:size(x0, 2)
        c = Norma.node_constraints(model, topology, n)
        k = length(side_sets_of(n))
        if k ≥ 3
            @test c === nothing
        else
            @test c !== nothing
            @test length(c.surfaces) == k
            @test !any(c.fixed)
        end
    end
    memory = Norma.PhaseMemory()
    Norma.phase_element_values!(memory, model, topology, relaxed)
    result = Norma.swap_operator!(memory, model, topology, relaxed)
    @test result.accepted > 0
    moved = [n for n in 1:size(x0, 2) if topology.positions[:, n] != x0[:, n]]
    boundary_moved = [n for n in moved if !isempty(side_sets_of(n))]
    @test !isempty(boundary_moved)
    @test any(n -> length(side_sets_of(n)) == 2, boundary_moved)
    @test !any(n -> length(side_sets_of(n)) ≥ 3, moved)
    # A moved boundary node stays on the faces of the cube: the Surface
    # conditions of the example are the planes x, y, z = ±1.
    names = Dict(id => topology.side_set_names[id] for id in keys(topology.side_set_names))
    plane(name) = (Dict('x' => 1, 'y' => 2, 'z' => 3)[name[3]], name[4] == '-' ? -1.0 : 1.0)
    for n in boundary_moved, id in side_sets_of(n)
        component, value = plane(names[id])
        @test abs(topology.positions[component, n] - value) < 1.0e-8
    end
    @test Set(Norma.boundary_faces(topology)) == Set(Norma.boundary_faces(Norma.build_topology(model)))
    @test all(Norma.element_volume(topology, e) > 0.0 for e in 1:Norma.num_alive_elements(topology))
    # Under a Dirichlet condition on a component the relaxation keeps it.
    dirichlet = smoothing_model(
        "../examples/ems/awful-cube/awful-cube.g", "awful", "swap-dirichlet.e"; size_field="0.214"
    )
    for n in 1:size(dirichlet.model.reference, 2)
        c = Norma.node_constraints(dirichlet.model, Norma.build_topology(dirichlet.model), n)
        sets = [id for (id, flags) in Norma.build_topology(dirichlet.model).node_sets if flags[n]]
        if length(sets) ≥ 3
            @test c === nothing
        elseif !isempty(sets)
            # The faces carry side sets without Surface conditions in this
            # model, so a boundary node may not move.
            @test c === nothing
        else
            @test c !== nothing && !any(c.fixed)
        end
        n ≥ 200 && break
    end
    Norma.finalize_writing(dirichlet)
    for f in ("swap-boundary.e", "swap-boundary.yaml", "swap-dirichlet.e")
        rm(f; force=true)
    end
end

@testset "relaxation_in_the_loop" begin
    # Under the energy criteria the loop with the relaxation ends at a lower
    # energy than without it; under the length and scaled Jacobian criteria
    # the relaxation does nothing, and the loops agree.
    base = Norma.load_input("../examples/ems/awful-cube/awful-cube-adaptive.yaml")
    base["input mesh file"] = "../examples/ems/awful-cube/awful-cube.g"
    base["name"] = "relaxation-loop"
    base["adaptivity"]["outer iterations"] = 2
    base["adaptivity"]["maximum passes"] = 3
    results = Dict{String,Any}()
    for (tag, criteria, relaxation) in (
        ("energy-fixed", false, false),
        ("energy-relaxed", false, true),
        ("geometric-fixed", true, false),
        ("geometric-relaxed", true, true),
    )
        params = deepcopy(base)
        params["output mesh file"] = "relaxation-loop-$tag.e"
        adaptivity = params["adaptivity"]
        if !criteria
            for key in ("size criterion", "shape criterion", "minimum scaled Jacobian")
                delete!(adaptivity, key)
            end
        end
        adaptivity["cavity relaxation"] = relaxation
        results[tag] = adapted_quality(params)
        for k in 1:2
            rm("relaxation-loop-$tag-adapted-$k.g"; force=true)
            rm(Norma.sequence_file_name("relaxation-loop-$tag.e", k + 1); force=true)
        end
        rm("relaxation-loop-$tag.e"; force=true)
    end
    @test results["energy-relaxed"].energy < results["energy-fixed"].energy
    @test results["energy-relaxed"].mean_sj > results["energy-fixed"].mean_sj
    @test results["geometric-relaxed"].energy == results["geometric-fixed"].energy
    @test results["geometric-relaxed"].elements == results["geometric-fixed"].elements
end
