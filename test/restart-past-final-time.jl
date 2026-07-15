# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus
using YAML

# A restart at or past `final time` leaves nothing to integrate. Before the
# fix, num_stops = max(round((final_time - restart_time) / time_step) + 1, 2)
# was forced to at least 2 regardless, and evolve() always advances one
# control step before checking stop_evolve(), so the run silently took one
# step past final_time and wrote a bogus extra stop there instead of
# stopping cleanly. This test checks that both the single-domain
# (process_restart!) and multi-domain (process_multidomain_restart!) restart
# paths now abort instead.
@testset "Restart At Or Past Final Time" begin
    # cuboid-restart-in.e's index 2 snapshot is at t = 0.1 (see
    # single-ahead-cuboid-dynamic-restart.jl).
    @testset "Single-domain restart" begin
        cp("../examples/ahead/single/cuboid/cuboid-hex.g", "../cuboid-hex.g"; force=true)
        cp(
            "../examples/ahead/single/cuboid/dynamic-restart/cuboid-restart-in.e",
            "cuboid-restart-in.e";
            force=true,
        )
        yaml_path = "../examples/ahead/single/cuboid/dynamic-restart/cuboid.yaml"

        try
            # Restart time == final time.
            params = YAML.load_file(yaml_path; dicttype=Norma.Parameters)
            params["name"] = "restart-past-final-time-equal"
            params["time integrator"]["final time"] = 0.1
            Norma.NORMA_TEST_MODE[] = true
            @test_throws Norma.NormaAbortException Norma.run(params)

            # Restart time > final time.
            params = YAML.load_file(yaml_path; dicttype=Norma.Parameters)
            params["name"] = "restart-past-final-time-past"
            params["time integrator"]["final time"] = 0.05
            @test_throws Norma.NormaAbortException Norma.run(params)

            # Sanity check: restart time comfortably before final time still
            # passes this particular validation (the original `final time:
            # 1.0` from the example file goes through create_simulation()
            # without the new check firing; this does not run the full
            # simulation, just confirms construction gets past process_restart!).
            params = YAML.load_file(yaml_path; dicttype=Norma.Parameters)
            params["name"] = "restart-past-final-time-ok"
            sim = Norma.create_simulation(params)
            @test sim.controller.initial_time ≈ 0.1 rtol = 1.0e-09
            Exodus.close(sim.params["input_mesh"])
            Exodus.close(sim.params["output_mesh"])
        finally
            Norma.NORMA_TEST_MODE[] = false
            rm("cuboid-restart-in.e"; force=true)
            rm("cuboid-restart-out.e"; force=true)
            rm("../cuboid-hex.g"; force=true)
        end
    end

    # process_multidomain_restart!() only reads each subdomain's `input mesh
    # file` (for the checkpoint time at the shared restart index) and
    # `model:` (for the j2-plasticity check), so a single minimal synthetic
    # subdomain file pointing at the same checkpoint mesh, used twice, is
    # enough to exercise the shared-restart-time validation without having
    # to stand up a full two-subdomain Schwarz problem.
    @testset "Multi-domain restart" begin
        cp(
            "../examples/ahead/single/cuboid/dynamic-restart/cuboid-restart-in.e",
            "restart-past-final-domain.e";
            force=true,
        )
        write(
            "restart-past-final-domain.yaml",
            "input mesh file: restart-past-final-domain.e\nmodel:\n  type: solid mechanics\n",
        )
        try
            params = Norma.Parameters(
                "domains" => ["restart-past-final-domain.yaml", "restart-past-final-domain.yaml"],
                "restart" => Norma.Parameters("index" => 2),
                "final time" => 0.1,
            )
            Norma.NORMA_TEST_MODE[] = true
            @test_throws Norma.NormaAbortException Norma.process_multidomain_restart!(params)

            params = Norma.Parameters(
                "domains" => ["restart-past-final-domain.yaml", "restart-past-final-domain.yaml"],
                "restart" => Norma.Parameters("index" => 2),
                "final time" => 0.05,
            )
            @test_throws Norma.NormaAbortException Norma.process_multidomain_restart!(params)

            # Restart time (0.1) comfortably before final time still passes.
            params = Norma.Parameters(
                "domains" => ["restart-past-final-domain.yaml", "restart-past-final-domain.yaml"],
                "restart" => Norma.Parameters("index" => 2),
                "final time" => 1.0,
            )
            @test Norma.process_multidomain_restart!(params) === nothing
            @test params["initial time"] ≈ 0.1 rtol = 1.0e-09
        finally
            Norma.NORMA_TEST_MODE[] = false
            rm("restart-past-final-domain.e"; force=true)
            rm("restart-past-final-domain.yaml"; force=true)
        end
    end
end
