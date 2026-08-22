# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

# Regression test for #221: restore_prev_state read sim.model.prev_state_old
# unconditionally, but that field only exists on SolidMechanics. Any ROM
# subdomain that failed a step and was retried crashed with a FieldError
# instead of rolling back.
#
# The regression assertion is the direct save/restore rollback check below:
# it exercises restore_prev_state on both quadratic OpInf ROM subdomains
# deterministically, and before the fix it threw
# `FieldError: type QuadraticOpInfRom has no field prev_state_old`
# regardless of whether the physics happens to fail a step. The full run
# that follows is integration coverage for the example itself: with the
# current operators it does not necessarily fail a step (measured: no
# rollback on Linux/Julia 1.12), so run completion alone must not be relied
# on to cover the rollback path.
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

    # Direct regression check for #221: drive both ROM subdomains through the
    # save_curr_state / restore_prev_state rollback pair that advance_one_step
    # uses on a step failure. Before the fix, restore_prev_state threw a
    # FieldError on the first ROM subdomain. finalize_writing releases the
    # Exodus output handles so the full run below can recreate its outputs.
    probe_sim = Norma.create_simulation(deepcopy(params))
    for subsim in probe_sim.subsims
        @test subsim.model isa Norma.QuadraticOpInfRom
        Norma.save_curr_state(subsim)
        Norma.restore_prev_state(subsim)
    end
    Norma.finalize_writing(probe_sim)

    sim = Norma.run(params)

    subsims = sim.subsims
    model_cuboid1 = subsims[1].model
    model_cuboid2 = subsims[2].model

    rm("cuboids.yaml"; force=true)
    rm("cuboid-1.yaml"; force=true)
    rm("cuboid-2.yaml"; force=true)
    rm("../cuboid-1.g"; force=true)
    rm("../cuboid-2.g"; force=true)
    rm("cuboid-1.e"; force=true)
    rm("cuboid-2.e"; force=true)
    rm("qopinf-operator-1.npz"; force=true)
    rm("qopinf-operator-2.npz"; force=true)

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
