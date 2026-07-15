# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# examples/overlap/static-same-step/cuboids-restart is a Schwarz-overlap,
# quasistatic (no velocity field) restart problem: its checkpoints
# (cuboid-1-in.e / cuboid-2-in.e) hold two snapshots, t = 0.0 and t = 0.2,
# and `restart: index: -1` resumes from the last one (t = 0.2), continuing
# with `time step: 0.8` to `final time: 1.0` -- landing on the same final
# time as the non-restart `cuboids` example it was derived from (same mesh,
# material, and boundary conditions, just a single step from t = 0 instead
# of a restarted step from t = 0.2).
#
# Since this is a quasistatic linear-elastic problem, both runs solve for
# equilibrium under the same applied boundary displacement at t = 1.0
# (`"1.0 * t"` on cuboid-2's driven face), so the final solution should
# agree regardless of the time-stepping path taken to get there. This test
# runs the full `cuboids` problem, then the restarted `cuboids-restart`
# problem, and compares their final nodal displacement fields.

# Run a Schwarz simulation from the given example directory and set of
# input files, then clean up. `output_files` are the per-subdomain output
# meshes Norma.run() produces, removed afterward along with the copied
# inputs.
function _run_and_cleanup(example_dir, top_level_yaml, yaml_files, output_files)
    for f in yaml_files
        cp("$example_dir/$f", f; force=true)
    end
    params = YAML.load_file(top_level_yaml; dicttype=Norma.Parameters)
    params["name"] = top_level_yaml
    sim = Norma.run(params)
    for f in vcat(yaml_files, output_files)
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Overlap Static Restart Matches Full Run" begin
    # --- Full (non-restart) run of the original problem ---------------------
    sim_full = _run_and_cleanup(
        "../examples/overlap/static-same-step/cuboids",
        "cuboids.yaml",
        ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.g", "cuboid-2.g"],
        ["cuboid-1.e", "cuboid-2.e"],
    )
    @test sim_full.failed == false

    # --- Restarted run, resuming from t = 0.2 to the same final time = 1.0 --
    sim_restart = _run_and_cleanup(
        "../examples/overlap/static-same-step/cuboids-restart",
        "cuboids.yaml",
        ["cuboids.yaml", "cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1-in.e", "cuboid-2-in.e"],
        ["cuboid-1-out.e", "cuboid-2-out.e"],
    )
    @test sim_restart.failed == false
    # Checkpoint's last snapshot is at t = 0.2 (see comment above).
    @test sim_restart.controller.initial_time ≈ 0.2 rtol = 1.0e-09
    # Both runs reach the same final time.
    @test sim_full.controller.time ≈ sim_restart.controller.time rtol = 1.0e-09

    # --- Compare final solutions ---------------------------------------------
    @test length(sim_full.subsims) == length(sim_restart.subsims) == 2
    for i in 1:2
        d_full = sim_full.subsims[i].model.displacement
        d_restart = sim_restart.subsims[i].model.displacement
        @test size(d_full) == size(d_restart)
        @test d_restart ≈ d_full atol = 1.0e-06 rtol = 1.0e-06
    end
end
