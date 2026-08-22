# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Regression test for #221: restore_prev_state read sim.model.prev_state_old
# unconditionally, but that field only exists on SolidMechanics. Any ROM
# subdomain that failed a step and was retried crashed with a FieldError
# instead of rolling back. This example is the case in the repository that
# actually exercises that path: adaptive time stepping (decrease/increase
# factor) is enabled on both ROM subdomains, and the quadratic OpInf ROM
# fails at least one step here and gets retried with a smaller time step --
# unlike the linear/cubic ROM-ROM cuboid cases, which converge without ever
# rolling back. Before the fix, this test would crash instead of completing.
@testset "Schwarz AHeaD Non-Overlap Dynamic Cuboid HEX8-HEX8 Quadratic ROM-ROM" begin
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/cuboids.yaml", "cuboids.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/cuboid-1.yaml", "cuboid-1.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/cuboid-2.yaml", "cuboid-2.yaml"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-1.g", "../cuboid-1.g"; force=true)
    cp("../examples/ahead/nonoverlap/cuboid/cuboid-2.g", "../cuboid-2.g"; force=true)
    cp(
        "../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/qopinf-operator-1.npz",
        "qopinf-operator-1.npz";
        force=true,
    )
    cp(
        "../examples/ahead/nonoverlap/cuboid/dynamic-neohookean-quadratic-rom-rom/qopinf-operator-2.npz",
        "qopinf-operator-2.npz";
        force=true,
    )
    input_file = "cuboids.yaml"
    params = YAML.load_file(input_file; dicttype=Norma.Parameters)
    params["initial time"] = 0.0
    params["time step"] = 0.01
    params["final time"] = 1.0
    params["name"] = input_file

    # Capture stdout so we can confirm a step actually failed and was retried
    # through the rollback path this issue is about, rather than merely
    # confirming the run finishes (which would also pass if this case never
    # happened to exercise restore_prev_state at all). redirect_stdout needs
    # a real file-backed stream (not a bare IOBuffer), so route it through a
    # temp file.
    log_path, log_io = mktemp()
    sim = redirect_stdout(log_io) do
        Norma.run(params)
    end
    close(log_io)
    log_output = read(log_path, String)
    rm(log_path; force=true)

    subsims = sim.subsims
    model_cuboid1 = subsims[1].model
    model_cuboid2 = subsims[2].model

    rm("cuboid.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("qopinf-operator-1.npz"; force=true)
    rm("qopinf-operator-2.npz"; force=true)

    # The regression itself: before the fix, restore_prev_state threw a
    # FieldError as soon as either ROM subdomain failed a step, so Norma.run
    # never returned. Reaching the final time at all is the primary assertion.
    @test sim.controller.time ≈ params["final time"] atol = 1e-9

    min_disp_x_cuboid1 = minimum(model_cuboid1.fom_model.displacement[1, :])
    min_disp_y_cuboid1 = minimum(model_cuboid1.fom_model.displacement[2, :])
    max_disp_z_cuboid1 = maximum(model_cuboid1.fom_model.displacement[3, :])
    min_disp_x_cuboid2 = minimum(model_cuboid2.fom_model.displacement[1, :])
    min_disp_y_cuboid2 = minimum(model_cuboid2.fom_model.displacement[2, :])
    max_disp_z_cuboid2 = maximum(model_cuboid2.fom_model.displacement[3, :])

    @test min_disp_x_cuboid1 ≈ -0.11889695523833804 atol = 1.0e-6   
    @test min_disp_y_cuboid1 ≈ -0.11889695523833796 atol = 1.0e-6   
    @test max_disp_z_cuboid1 ≈ 0.5000014637341199 atol = 1.0e-8
    @test min_disp_x_cuboid2 ≈ -0.11889368674381376 atol = 1.0e-6   
    @test min_disp_y_cuboid2 ≈ -0.11889368674381375 atol = 1.0e-6   
    @test max_disp_z_cuboid2 ≈ 1.0 atol = 0.0
end
