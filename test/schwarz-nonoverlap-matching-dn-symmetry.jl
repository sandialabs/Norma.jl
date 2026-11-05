# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Matching DN exchange with symmetry planes that cut the interface
# (held_interface_dofs). Two stacked boxes of
# examples/nonoverlap/dynamic-same-step/symmetry-seam, nonconforming at the
# interface, impacted at -35 m/s in z. The quarter model has symmetry
# conditions on nsx- and nsy- of both boxes; its meshes mirror the full ones,
# whose solution is the reference. Without held_interface_dofs the matching
# jump on the plane rows never falls below Π_D v_N there (594 of 600 stops
# unconverged, final jump 3.8e-3), the symmetry reaction is transferred to the
# Neumann side (11% of the force), and the rate invariant drifts by 2.1e-3.

const SEAM_EXAMPLE = joinpath(@__DIR__, "..", "examples", "nonoverlap", "dynamic-same-step", "symmetry-seam")

function seam_run(kind::String)
    dir = mktempdir()
    result = cd(dir) do
        for file in readdir(joinpath(SEAM_EXAMPLE, kind))
            cp(joinpath(SEAM_EXAMPLE, kind, file), file; force=true)
        end
        sim = Norma.run("seam.yaml")
        rows = readlines("seam-energy.csv")
        i = findfirst(==("rate_invariant_total"), split(rows[1], ","))
        invariant = [parse(Float64, split(row, ",")[i]) for row in rows[2:end]]
        (sim=sim, invariant=invariant)
    end
    rm(dir; recursive=true, force=true)
    return result
end

seam_subsim(run, name) = run.sim.subsims[findfirst(s -> s.name == name, run.sim.subsims)]
seam_dn_bc(subsim) =
    only(bc for bc in subsim.model.boundary_conditions if bc isa Norma.SolidMechanicsNonOverlapSchwarzBoundaryCondition)

const SEAM_FULL = seam_run("full")
const SEAM_QUARTER = seam_run("quarter")

@testset "Matching DN Symmetry Seam: Quarter Converges" begin
    # Measured: 8.06 iterations on average, 9 at most, final jump 4.8e-11.
    @test SEAM_QUARTER.sim.failed == false
    @test maximum(SEAM_QUARTER.sim.controller.schwarz_iters) ≤ 16
    @test Norma.dn_interface_residuals(SEAM_QUARTER.sim)[1][2].velocity_jump ≤ 1.0e-9
end

@testset "Matching DN Symmetry Seam: Neumann Force Excludes the Symmetry Reaction" begin
    upper = seam_subsim(SEAM_QUARTER, "upper")
    bc_D = seam_dn_bc(upper)
    bc_N = seam_dn_bc(seam_subsim(SEAM_QUARTER, "lower"))
    # Held: the x rows of the interface nodes on x = 0, the y rows on y = 0.
    X = upper.model.reference[:, bc_D.global_from_local_map]
    held = Norma.held_interface_dofs(bc_D)
    @test held == vcat(abs.(X[1:2, :]) .< 1.0e-12, falses(1, size(X, 2)))
    # The force applied on the Neumann side is Π_Dᵀ r_D with the held rows zeroed.
    m = upper.model
    r_D = reshape(Norma.extract_local_vector(bc_D, -(m.internal_force + Norma.dalembert_inertia_minus_loads(m)), 3), 3, :)
    r_D[held] .= 0.0
    expected = vec(transpose(transpose(bc_D.dirichlet_projector) * transpose(r_D)))
    @test norm(bc_N.transferred_force - expected) ≤ 1.0e-12 * norm(expected)
end

@testset "Matching DN Symmetry Seam: Quarter Reproduces the Full Solution" begin
    # Lateral displacement of the quarter nodes against their twins in the full
    # model: measured 1.5e-3, the difference between the unconstrained and the
    # plane-constrained projection at the seam (6.2e-3 without the fix).
    key(x) = Tuple(round.(x; digits=10) .+ 0.0)
    lateral_q, lateral_f = Float64[], Float64[]
    for name in ("upper", "lower")
        mq, mf = seam_subsim(SEAM_QUARTER, name).model, seam_subsim(SEAM_FULL, name).model
        twin = Dict(key(mf.reference[:, i]) => i for i in axes(mf.reference, 2))
        idx = [twin[key(mq.reference[:, i])] for i in axes(mq.reference, 2)]
        append!(lateral_q, vec(mq.displacement[1:2, :]))
        append!(lateral_f, vec(mf.displacement[1:2, idx]))
    end
    @test norm(lateral_q - lateral_f) ≤ 5.0e-3 * norm(lateral_f)
    # Rate invariant: measured 1.9e-11 (quarter) and 1.2e-11 (full).
    @test maximum(abs.(SEAM_QUARTER.invariant ./ SEAM_QUARTER.invariant[1] .- 1.0)) ≤ 1.0e-9
end
