# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Tests the AHeaD Schwarz overlap two-domain notched-cylinder problem in
# which domain 1 starts as a linear OpInf reduced-order model (ROM) and is
# swapped mid-run to a full-order model (FOM), while domain 2 runs as a ROM
# throughout.  This is the opposite swap direction from
# schwarz-ahead-overlap-dynamic-clamped-3sd-fom-rom-swap.jl (ROM -> FOM here,
# vs FOM -> ROM there), and on a 3-D tet/hex notched-cylinder mesh rather
# than the 1-D clamped bar, so it exercises a different code path in
# copy_model_state! (the ROM -> FOM overload, which lifts the ROM's shadow
# FOM fields into the destination FOM) and a different swap target slot.
#
# Run exactly as published in
# examples/ahead/overlap/notched-cylinder/dynamic-opinf-swap-rom-fom — no
# parameter overrides; CSV output is already disabled in the example itself
# ("CSV output interval: 0").
#
# Configuration (notched-cylinder.yaml):
#   domains: notched-cylinder-rom-1 (ROM), notched-cylinder-rom-2 (ROM)
#   Material: linear elastic, E = 70 GPa, ν = 0.25, ρ = 1000 kg/m³
#             (yield stress / hardening modulus are present in the YAML but
#             unused, since the material model is "linear elastic", not J2)
#   Time integrator: Newmark (β = 0.25, γ = 0.5), adaptive time stepping
#       (maximum/minimum time step, increase/decrease factor)
#   Domain 1: -X_bottom, -Y_bottom, -Z_bottom pinned to 0 (x, y, z resp.)
#   Domain 2: -X_top, -Y_top pinned to 0 (x, y); +Z_top driven by
#       u_z(t) = 0.0064 * (0.5 - 0.5*cos(π t / 6))
#   initial time = 0.0, final time = 6.0
#
# Swap: domain 1 (notched-cylinder-rom-1, a LinearOpInfRom) is replaced by
# notched-cylinder-fom-1 (a SolidMechanics FOM) once the controller's last
# completed step reaches t_swap = 3.0 — well before the final time of 6.0,
# so the swap has ample margin to fire and settle (note: with adaptive time
# stepping the exact step count is not fixed in advance, but the 3.0 time
# units of margin between t_swap and final time make this robust regardless
# of step-size history).
#
# notched-cylinder-fom-1.yaml declares the same "output mesh file"
# (notched-cylinder-1.e) as notched-cylinder-rom-1.yaml, which has already
# been created by the time the swap fires, so Norma's uniquify_swap_output!
# renames the replacement's output file (and the subsim) with a "-phase2"
# suffix. The post-swap subsim's name is derived from the *replacement
# file's own* stem (stripped_name("notched-cylinder-fom-1.yaml") =
# "notched-cylinder-fom-1"), not from the original domain's base name, so
# the post-swap name is "notched-cylinder-fom-1-phase2" (note: still
# "fom-1", not "1") — while the renamed *output mesh file* is
# "notched-cylinder-1-phase2.e" (base taken from the file's own declared
# "output mesh file" value, which is "notched-cylinder-1.e" regardless of
# the YAML's own filename).
#
# Assertions:
#   1. The swap fired: subsim[1] carries the phase-2 name, is now a plain
#      SolidMechanics FOM, and the swap plan is marked applied.
#   2. handle_by_name resolves both the original and phase-2 names to slot 1
#      (the original name survives as an alias; see the apply_swap! comment
#      in src/swap.jl for why old names are kept rather than removed).
#   3. The simulation completed to final time without failure.
#   4. Physical sanity: the Dirichlet boundary conditions are satisfied to
#      numerical precision on both domains at the final state — in
#      particular, the driven +Z_top face on domain 2 reaches exactly
#      0.0064 * (0.5 - 0.5*cos(π)) = 0.0064 at t = 6.0, and every pinned
#      face (-X/-Y/-Z on domain 1, -X/-Y on domain 2) is at 0. This is an
#      end-to-end check that the swapped-in FOM is correctly coupled and
#      satisfying its boundary data, not a comparison against a saved
#      reference solution.

@testset "Schwarz AHeaD Overlap Dynamic Notched Cylinder OpInf ROM-to-FOM Swap" begin
    example_dir = "../examples/ahead/overlap/notched-cylinder/dynamic-opinf-swap-rom-fom"
    mesh_dir    = "../examples/ahead/overlap/notched-cylinder"

    # ── Stage input files ──────────────────────────────────────────────────
    # YAML files are copied to the test working directory; mesh files are
    # expected one level up (as specified by "input mesh file: ../…" in the
    # sub-domain YAMLs), matching the convention used by
    # schwarz-ahead-overlap-dynamic-notched-cylinder-j2-swap-fom.jl.
    cp("$example_dir/notched-cylinder.yaml",         "notched-cylinder.yaml";         force=true)
    cp("$example_dir/notched-cylinder-rom-1.yaml",    "notched-cylinder-rom-1.yaml";   force=true)
    cp("$example_dir/notched-cylinder-fom-1.yaml",    "notched-cylinder-fom-1.yaml";   force=true)
    cp("$example_dir/notched-cylinder-rom-2.yaml",    "notched-cylinder-rom-2.yaml";   force=true)
    cp("$example_dir/linear-opinf-operator-M3-1.npz", "linear-opinf-operator-M3-1.npz"; force=true)
    cp("$example_dir/linear-opinf-operator-M3-2.npz", "linear-opinf-operator-M3-2.npz"; force=true)
    cp("$mesh_dir/notched-cylinder-1-coarse.g",       "../notched-cylinder-1-coarse.g"; force=true)
    cp("$mesh_dir/notched-cylinder-2-coarse.g",       "../notched-cylinder-2-coarse.g"; force=true)

    # ── Run ────────────────────────────────────────────────────────────────
    # Run the published example exactly as-is: no time-shortening or CSV
    # overrides (CSV output is already off in the example's own YAML).
    sim = Norma.run("notched-cylinder.yaml")

    # ── Cleanup ────────────────────────────────────────────────────────────
    rm("notched-cylinder.yaml";          force=true)
    rm("notched-cylinder-rom-1.yaml";    force=true)
    rm("notched-cylinder-fom-1.yaml";    force=true)
    rm("notched-cylinder-rom-2.yaml";    force=true)
    rm("linear-opinf-operator-M3-1.npz"; force=true)
    rm("linear-opinf-operator-M3-2.npz"; force=true)
    rm("../notched-cylinder-1-coarse.g"; force=true)
    rm("../notched-cylinder-2-coarse.g"; force=true)
    rm("notched-cylinder-1.e";           force=true)
    # The swap's replacement (notched-cylinder-fom-1.yaml) reuses domain 1's
    # original "output mesh file" (notched-cylinder-1.e), already written
    # before the swap; uniquify_swap_output! detects the collision and
    # renames the output mesh file (base taken from the declared "output
    # mesh file" value, "notched-cylinder-1.e") to
    # "notched-cylinder-1-phase2.e" — note this is independent of the
    # replacement YAML's own filename, which is reflected in the subsim's
    # *name* instead (see the assertion below).
    rm("notched-cylinder-1-phase2.e";    force=true)
    rm("notched-cylinder-2.e";           force=true)

    subsims = sim.subsims
    subsim_1 = subsims[1]
    subsim_2 = subsims[2]
    model_1  = subsim_1.model
    model_2  = subsim_2.model

    # ── Swap fired ─────────────────────────────────────────────────────────
    @test length(sim.swaps) == 1
    @test sim.swaps[1].applied == true
    @test sim.swaps[1].subsim_name == "notched-cylinder-rom-1"

    # After the swap, slot 1 holds the phase-2 FOM replacement; slot 2 (never
    # swapped) keeps its original ROM model throughout.
    @test subsim_1.name == "notched-cylinder-fom-1-phase2"
    @test model_1 isa Norma.SolidMechanics
    @test model_2 isa Norma.LinearOpInfRom

    # handle_by_name resolves both the original domain name and the new
    # phase-2 name to slot 1.
    @test subsim_1.handle.id == 1
    @test sim.handle_by_name["notched-cylinder-rom-1"].id == 1
    @test sim.handle_by_name["notched-cylinder-fom-1-phase2"].id == 1

    # ── Completion ─────────────────────────────────────────────────────────
    @test sim.failed == false
    @test sim.controller.time ≈ 6.0 rtol = 1.0e-9

    # ── Physical sanity: Dirichlet BCs satisfied on the post-swap state ─────
    # Domain 1 is a plain SolidMechanics FOM after the swap, so its
    # displacement field is read directly.  Domain 2 is still a ROM: its
    # *boundary condition definitions* live on the ROM wrapper itself
    # (model_2.boundary_conditions — create_bcs assigns to sim.model, which
    # for a ROM subsim is the LinearOpInfRom, not its shadow FOM), but the
    # actual *displacement values* at those constrained nodes are written
    # onto the shadow FOM (apply_bc(model::OpInfModel, bc) calls
    # apply_bc(model.fom_model, bc) — see src/opinf/opinf_ics_bcs.jl), so the
    # shadow FOM's boundary_conditions list itself is always empty and must
    # not be searched directly.
    fom_2 = Norma.get_fom_model(subsim_2)

    # Helper: find the Dirichlet BC on `node_set_name` (component `offset`,
    # 1=x/2=y/3=z) on `bc_model`, and return the displacement values at its
    # constrained nodes from `disp_model` (defaults to `bc_model` itself,
    # the common case where both the BC list and the live displacement data
    # live on the same object, i.e. a plain SolidMechanics FOM).
    function dirichlet_disp(bc_model, node_set_name::String, offset::Int64, disp_model=bc_model)
        for bc in bc_model.boundary_conditions
            if bc isa Norma.SolidMechanicsDirichletBoundaryCondition && bc.name == node_set_name && bc.offset == offset
                return disp_model.displacement[offset, bc.node_set_node_indices]
            end
        end
        error("Dirichlet BC on node set '$node_set_name' (offset $offset) not found")
    end

    # Domain 1 (post-swap FOM): -X_bottom/-Y_bottom/-Z_bottom all pinned to 0.
    @test maximum(abs.(dirichlet_disp(model_1, "-X_bottom", 1))) ≈ 0.0 atol = 1.0e-12
    @test maximum(abs.(dirichlet_disp(model_1, "-Y_bottom", 2))) ≈ 0.0 atol = 1.0e-12
    @test maximum(abs.(dirichlet_disp(model_1, "-Z_bottom", 3))) ≈ 0.0 atol = 1.0e-12

    # Domain 2 (ROM): BC definitions come from model_2 (the ROM wrapper),
    # displacement values come from fom_2 (its shadow FOM).  -X_top/-Y_top
    # pinned to 0; +Z_top driven to
    # u_z(6.0) = 0.0064 * (0.5 - 0.5*cos(π*6/6)) = 0.0064 exactly.
    @test maximum(abs.(dirichlet_disp(model_2, "-X_top", 1, fom_2))) ≈ 0.0 atol = 1.0e-12
    @test maximum(abs.(dirichlet_disp(model_2, "-Y_top", 2, fom_2))) ≈ 0.0 atol = 1.0e-12
    u_z_top_exact = 0.0064 * (0.5 - 0.5 * cos(π * 6.0 / 6.0))
    disp_z_top = dirichlet_disp(model_2, "+Z_top", 3, fom_2)
    @test all(d -> isapprox(d, u_z_top_exact; atol=1.0e-8), disp_z_top)

    # Schwarz iteration counts should be non-trivial (coupling active).
    @test all(sim.controller.schwarz_iters .>= 0)
    @test any(sim.controller.schwarz_iters .> 0)
end
