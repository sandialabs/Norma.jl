# Energy kernel on element subsets for the adaptivity loop
# (docs/notes/ems-adaptivity): the energies of proposed elements, which need
# not exist in the mesh, evaluated exactly as the assembly evaluates them, and
# checked on hand-built cavities before and after an edge swap, an edge
# collapse, and an edge split, with a metric target and with the equal-volume
# fallback.
using LinearAlgebra
using Random
using Test
using StaticArrays

if !isdefined(Main, :Norma)
    include("../src/Norma.jl")
end
Random.seed!(0)

tet_volume(x) = dot(x[:, 2] - x[:, 1], cross(x[:, 3] - x[:, 1], x[:, 4] - x[:, 1])) / 6.0
tet_gradient(X, x) = SMatrix{3,3,Float64,9}((x[:, 2:4] .- x[:, 1]) / (X[:, 2:4] .- X[:, 1]))
# Reorder the nodes of every element so that its volume is positive.
function orient!(conn, positions)
    for e in 1:size(conn, 2)
        if tet_volume(positions[:, conn[:, e]]) < 0.0
            conn[3, e], conn[4, e] = conn[4, e], conn[3, e]
        end
    end
    return conn
end

# The models supply the block data, the material, and the target; the
# positions and connectivities below are the test's own.
function smoothing_model(smooth_reference, extra)
    params = Dict{String,Any}(
        "type" => "single",
        "name" => "cavity",
        "input mesh file" => "../examples/ems/cube/cube.g",
        "output mesh file" => "cavity-" * replace(smooth_reference, " " => "-") * ".e",
        "Exodus output interval" => 0,
        "CSV output interval" => 0,
        "model" => merge(
            Dict{String,Any}(
                "type" => "mesh smoothing",
                "smooth reference" => smooth_reference,
                "material" => Dict{String,Any}(
                    "elastic" => Dict{String,Any}(
                        "model" => "seth-hill",
                        "m" => 2,
                        "n" => 2,
                        "bulk modulus" => 1.0,
                        "shear modulus" => 1.0,
                        "density" => 1.0,
                    ),
                    "blocks" => Dict{String,Any}("cube" => "elastic"),
                ),
            ),
            extra,
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

h = 0.4
θ = 0.3
sizes = (0.5h, 1.2h, 0.9h)
metric_sim = smoothing_model(
    "metric field unrestricted",
    Dict{String,Any}(
        "metric field" => Dict{String,Any}(
            "sizes" => ["$(sizes[1])", "$(sizes[2])", "$(sizes[3])"], "rotation vector" => ["0", "0", "$θ"]
        ),
    ),
)
volume_sim = smoothing_model("equal volume", Dict{String,Any}())
material = metric_sim.model.materials[1]

# Independent evaluation of the energy of one element: the ideal element built
# by hand from the target, one-point quadrature on a linear tetrahedron.
R = Norma.rt_of_rv(SVector(0.0, 0.0, θ))
F_M = Diagonal([1 / sizes[1], 1 / sizes[2], 1 / sizes[3]]) * R'
c = 0.5 / sqrt(2.0)
Y = c * [1 -1 -1 1; 1 -1 1 -1; 1 1 -1 -1]
function reference_energy_metric(x)
    X = R * Diagonal(collect(sizes)) * Y
    F = tet_gradient(X, x)
    return Norma.strain_energy(material, SMatrix{3,3,Float64,9}(F_M * F * inv(F_M))) * abs(tet_volume(X))
end
function reference_energy_volume(x)
    a = cbrt(6.0 * sqrt(2.0) * tet_volume(x))   # edge of the regular tetrahedron of equal volume
    X = a * Y
    return Norma.strain_energy(material, tet_gradient(X, x)) * abs(tet_volume(X))
end
reference_energy(sim, conn, positions) = sum(
    (sim === metric_sim ? reference_energy_metric : reference_energy_volume)(positions[:, conn[:, e]]) for
    e in 1:size(conn, 2)
)

# Cavity of an interior edge: the edge from a to b along z, four ring nodes
# around it, and the four tetrahedra that share the edge.
a_node, b_node = 1, 2
loop_positions = hcat(
    [0.0, 0.0, -0.5h], [0.05h, -0.02h, 0.55h],
    [0.9h, 0.0, 0.1h], [0.0, 0.7h, -0.1h], [-0.8h, 0.1h, 0.0], [0.1h, -1.1h, 0.05h],
)
ring = [3, 4, 5, 6]
loop_conn = orient!(hcat([[a_node, b_node, ring[i], ring[mod1(i + 1, 4)]] for i in 1:4]...), loop_positions)
# The two triangulations of the ring quadrilateral give the two swapped
# configurations, each with four elements.
function swapped(diagonal)
    p, q = diagonal
    others = setdiff(ring, [p, q])
    conn = hcat([[p, q, r, a_node] for r in others]..., [[p, q, r, b_node] for r in others]...)
    return orient!(conn, loop_positions)
end
swap_a = swapped((3, 5))
swap_b = swapped((4, 6))
# Split of the edge at a new node m: every element of the loop is bisected.
split_positions = hcat(loop_positions, 0.5 * (loop_positions[:, a_node] + loop_positions[:, b_node]))
m_node = 7
split_conn = orient!(
    hcat(
        [[a_node, m_node, ring[i], ring[mod1(i + 1, 4)]] for i in 1:4]...,
        [[m_node, b_node, ring[i], ring[mod1(i + 1, 4)]] for i in 1:4]...,
    ),
    split_positions,
)
# Star of an interior node a inside a tetrahedron (p1, p2, p3, p4): four
# elements; collapsing the edge from a to p1 removes the three elements that
# contain it and maps the fourth onto the outer tetrahedron.
star_positions = hcat(0.6h * Y .+ [0.02h; -0.03h; 0.01h], 0.6h * Y[:, 1] * 0.2 + [0.05h; 0.0; -0.02h])
star_a = 5
star_conn = orient!(
    hcat([[star_a, 2, 3, 4], [star_a, 1, 4, 3], [star_a, 1, 2, 4], [star_a, 1, 3, 2]]...), star_positions
)
collapsed_conn = orient!(reshape([1, 2, 3, 4], 4, 1), star_positions)

@testset "cavity_energies" begin
    for sim in (metric_sim, volume_sim)
        model = sim.model
        for (conn, positions) in (
            (loop_conn, loop_positions),
            (swap_a, loop_positions),
            (swap_b, loop_positions),
            (split_conn, split_positions),
            (star_conn, star_positions),
            (collapsed_conn, star_positions),
        )
            energies = Norma.element_energies(model, 1, conn, positions)
            @test length(energies) == size(conn, 2)
            @test all(isfinite, energies)
            @test all(energies .>= 0.0)
            @test sum(energies) ≈ reference_energy(sim, conn, positions) rtol = 1.0e-12 atol = 1.0e-24
        end
        # The acceptance test compares cavity sums; here the collapse restores the regular
        # outer tetrahedron, which the equal-volume rule scores as ideal.
        star = sum(Norma.element_energies(model, 1, star_conn, star_positions))
        collapsed = sum(Norma.element_energies(model, 1, collapsed_conn, star_positions))
        @test collapsed < star
        if sim === volume_sim
            @test collapsed ≈ 0.0 atol = 1.0e-24
        end
        # An inverted element has infinite energy.
        inverted = copy(loop_conn)
        inverted[3, 1], inverted[4, 1] = inverted[4, 1], inverted[3, 1]
        @test Norma.element_energies(model, 1, inverted, loop_positions)[1] == Inf
        # With the target sampled at the current positions the ideal element
        # of the equal-volume rule follows the current volume: scaling the
        # positions leaves the energy density unchanged, so the energy scales
        # with the volume.
        if sim === volume_sim
            scaled = 2.5 * loop_positions
            @test sum(Norma.element_energies(model, 1, loop_conn, scaled)) ≈
                2.5^3 * sum(Norma.element_energies(model, 1, loop_conn, loop_positions)) rtol = 1.0e-12
        end
    end
    # The two swapped configurations differ, so the acceptance test can discriminate.
    ea = sum(Norma.element_energies(metric_sim.model, 1, swap_a, loop_positions))
    eb = sum(Norma.element_energies(metric_sim.model, 1, swap_b, loop_positions))
    @test abs(ea - eb) > 1.0e-6 * max(ea, eb)
end

@testset "kernel_matches_assembly" begin
    # Over the whole mesh, with the target sampled on the original mesh as the
    # smoother does, the kernel reproduces the assembled energy.
    for sim in (metric_sim, volume_sim)
        model = sim.model
        conn = model.blocks[1].connectivity
        model.displacement .= 0.02 * 0.1 * randn(3, size(model.reference, 2))
        Norma.evaluate(model, sim.integrator, sim.solver)
        @test model.failed == false
        positions = model.reference + model.displacement
        energies = Norma.element_energies(model, 1, conn, positions; sample_positions=model.reference)
        @test sum(energies) ≈ model.strain_energy rtol = 1.0e-12
        @test energies ≈ model.stored_energy[1] rtol = 1.0e-12
        volumes = Norma.ideal_volumes(model, 1)
        @test all(volumes .> 0.0)
        if sim === metric_sim
            @test all(volumes .≈ sizes[1] * sizes[2] * sizes[3] / (6.0 * sqrt(2.0)))
        end
    end
    Norma.finalize_writing(metric_sim)
    Norma.finalize_writing(volume_sim)
end
