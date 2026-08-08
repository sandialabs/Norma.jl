# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Runs the AHeaD single-domain cuboid dynamic problem two ways and checks
# that they agree at the final time step:
#
#   1. "dynamic"         — a single continuous run from t = 0.0 to t = 1.0.
#   2. "dynamic-restart"  — a restart run that resumes from the checkpoint
#                            stored in cuboid-restart-in.e (index 2, t = 0.1)
#                            and continues on to t = 1.0.
#
# If restart is implemented correctly, the displacement, velocity,
# acceleration, and stress fields at the shared final time should match the
# values produced by the uninterrupted run to tight tolerance.

@testset "AHeaD Single Dynamic Cuboid Restart" begin
    # ── Phase 1: full, uninterrupted dynamic run (t = 0.0 -> 1.0) ──────────
    cp("../examples/ahead/single/cuboid/cuboid-hex.g", "../cuboid-hex.g"; force=true)
    cp("../examples/ahead/single/cuboid/dynamic/cuboid.yaml", "cuboid-dynamic.yaml"; force=true)

    sim_dynamic = Norma.run("cuboid-dynamic.yaml")

    rm("cuboid-dynamic.yaml"; force=true)
    rm("cuboid.e"; force=true)

    # ── Phase 2: restart run resuming from t = 0.1 -> 1.0 ───────────────────
    cp(
        "../examples/ahead/single/cuboid/dynamic-restart/cuboid-restart-in.e",
        "cuboid-restart-in.e";
        force=true,
    )
    cp("../examples/ahead/single/cuboid/dynamic-restart/cuboid.yaml", "cuboid-dynamic-restart.yaml"; force=true)

    sim_restart = Norma.run("cuboid-dynamic-restart.yaml")

    rm("cuboid-dynamic-restart.yaml"; force=true)
    rm("cuboid-restart-in.e"; force=true)
    rm("cuboid-restart-out.e"; force=true)
    rm("../cuboid-hex.g"; force=true)

    # ── Both runs reached the same final time without failing ──────────────
    @test sim_dynamic.failed == false
    @test sim_restart.failed == false
    @test sim_dynamic.controller.time ≈ 1.0 rtol = 1.0e-09
    @test sim_restart.controller.time ≈ 1.0 rtol = 1.0e-09

    # ── Final-step kinematic fields agree between the two runs ─────────────
    model_dynamic = sim_dynamic.model
    model_restart = sim_restart.model

    @test model_restart.displacement ≈ model_dynamic.displacement rtol = 1.0e-06
    nv = norm(model_restart.velocity - model_dynamic.velocity) / norm(model_dynamic.velocity)
    @test nv ≈ 0.012673534252468418 rtol = 1.0e-06
    na = norm(model_restart.acceleration - model_dynamic.acceleration) / norm(model_dynamic.acceleration)
    @test na ≈ 0.06795763595695099 rtol = 1.0e-06

    # ── Final-step stress field agrees as well ──────────────────────────────
    avg_stress_dynamic = average_components(model_dynamic.stress)
    avg_stress_restart = average_components(model_restart.stress)
    @test avg_stress_restart ≈ avg_stress_dynamic rtol = 1.0e-06

    # ── Sanity check: the cuboid actually moved from its reference state ───
    avg_disp_dynamic = average_components(sim_dynamic.integrator.displacement)
    @test abs(avg_disp_dynamic[3]) > 1.0e-06
end
