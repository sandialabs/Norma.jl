# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests a single-domain dynamic simulation that starts in OpInf ROM mode and
# switches to a FOM (solid mechanics) at a specified time threshold.
#
# Problem setup:
#   Mesh: cuboid hex8 (cuboid.g)
#   Material: linear elastic, E = 1 GPa, ν = 0.25, ρ = 1000 kg/m³
#   Time integrator: Newmark (β = 0.25, γ = 0.5 — trapezoidal rule)
#   Applied displacement (nsz+ node set, z-component):
#     u_z = 0.5 - 0.5 * cos(π * t)
#   All other faces pinned in their normal direction.
#
# Phase 1 (ROM):  t ∈ [0.0, 0.75],  Δt = 0.01  (linear OpInf ROM)
# Phase 2 (FOM):  t ∈ [0.75, 1.0],  Δt = 0.1   (solid mechanics FOM)
#
# The TimeSwapCriterion fires when sim.controller.time > 0.75, replacing
# cuboid-rom-t-swap.yaml with cuboid-fom-t-swap-phase2.yaml.
# After the swap the returned simulation is a SingleDomainSimulation whose
# model is a SolidMechanics (FOM) instance.

@testset "Single Dynamic OpInf ROM-to-FOM Time Swap" begin
    example_dir = "../examples/ahead/single/cuboid/dynamic-opinf-rom-fom-swap"

    # ── Stage files ─────────────────────────────────────────────────────────
    # The YAML files reference the mesh as ../cuboid.g (relative to CWD),
    # so we copy it one directory above the test working directory.
    cp("../examples/ahead/single/cuboid/cuboid.g", "../cuboid.g"; force=true)
    cp("$example_dir/cuboid-rom-t-swap.yaml",          "cuboid-rom-t-swap.yaml";          force=true)
    cp("$example_dir/cuboid-fom-t-swap-phase2.yaml",   "cuboid-fom-t-swap-phase2.yaml";   force=true)
    cp("$example_dir/opinf-operator.npz",              "opinf-operator.npz";              force=true)

    sim = Norma.run("cuboid-rom-t-swap.yaml")

    # ── Clean up ────────────────────────────────────────────────────────────
    rm("cuboid-rom-t-swap.yaml";        force=true)
    rm("cuboid-fom-t-swap-phase2.yaml"; force=true)
    rm("opinf-operator.npz";            force=true)
    rm("../cuboid.g";                   force=true)
    rm("cuboid-rom.e";                  force=true)
    rm("cuboid-fom-phase2.e";           force=true)

    # ── Swap fired ──────────────────────────────────────────────────────────
    # After execute the sim should carry the phase-2 name and have no
    # remaining swap plans.
    @test sim.name == "cuboid-fom-t-swap-phase2"
    @test all(p.applied for p in sim.swaps)

    # ── Completion ──────────────────────────────────────────────────────────
    @test sim.failed == false
    # Phase 2 final time is 1.0; check we reached it.
    @test sim.controller.time ≈ 1.0 rtol = 1.0e-9

    # ── Model type after swap ────────────────────────────────────────────────
    # The swapped-in phase-2 simulation uses a plain solid mechanics FOM,
    # not an OpInf ROM.
    @test sim.model isa Norma.SolidMechanics

    # ── Physical sanity checks ───────────────────────────────────────────────
    # At t = 1.0 the top face (nsz+) has u_z = 0.5 - 0.5*cos(π) = 1.0 m.
    # Due to Poisson contraction the x and y displacements are negative.
    # We check the final-state displacement statistics on the FOM model.
    disp = sim.model.displacement   # 3 × num_nodes matrix

    max_disp_z = maximum(disp[3, :])
    min_disp_x = minimum(disp[1, :])
    min_disp_y = minimum(disp[2, :])

    # The top face reaches close to 1.0 m in z.
    @test max_disp_z ≈ 1.0 atol = 1.0e-2

    # Poisson contraction: lateral displacements are non-trivially negative.
    @test min_disp_x < 0.0
    @test min_disp_y < 0.0

    # Average stress: the dominant component is σ_zz ≈ E * ε_zz.
    avg_stress = average_components(sim.model.stress)
    # σ_xx and σ_yy should be small (near zero for uniaxial stretch).
    @test abs(avg_stress[1]) < abs(avg_stress[3]) * 0.5
    @test abs(avg_stress[2]) < abs(avg_stress[3]) * 0.5
    # σ_zz should be the dominant (positive) stress component.
    @test avg_stress[3] > 0.0

    # ── Regression values ────────────────────────────────────────────────────
    # These values were obtained from a reference run and are stored here to
    # guard against unintended changes. Tolerances are intentionally loose
    # (1 %) to accommodate minor floating-point variation across platforms.
    @test max_disp_z  ≈ 1.0        rtol = 1.0e-6
    @test min_disp_x  ≈ -0.16669139545385825 rtol = 1.0e-4
    @test min_disp_y  ≈ -0.16669139545385825 rtol = 1.0e-4
    @test avg_stress[3] ≈ 6.666666666666667e8 rtol = 1.0e-2
end
