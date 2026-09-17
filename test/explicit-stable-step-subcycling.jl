# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# A central difference run whose requested time step exceeds the stable step
# must reach each stop by substeps of the stable step. Before the cap moved
# into advance_time, the predictor capped the step used by the kinematic
# update while the time had already advanced by the requested step, so the
# whole response was compressed in time by the ratio of the two steps. The
# clamped explicit example has a d'Alembert reference solution (see
# single-explicit-dynamic-solid-clamped.jl); its stable step at CFL 0.2 is
# about 3.2e-7 s, so a requested step of 1.0e-6 s is capped on every stop.

using YAML

@testset "Explicit Stable Step Subcycling" begin
    cp("../examples/single/explicit-dynamic-solid/clamped/clamped.g", "clamped.g"; force=true)
    input = YAML.load_file("../examples/single/explicit-dynamic-solid/clamped/clamped.yaml"; dicttype=Dict{String,Any})
    input["time integrator"]["time step"] = 1.0e-06
    YAML.write_file("clamped-subcycled.yaml", input)
    simulation = Norma.run("clamped-subcycled.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("clamped-subcycled.yaml"; force=true)
    rm("clamped.g"; force=true)
    rm("clamped.e"; force=true)

    @test integrator.time ≈ 1.0e-05 rtol = 1.0e-10
    @test integrator.time_step < 1.0e-06

    α = 0.01
    s = 0.02
    c = sqrt(1.0e9 / 1000.0)
    t = 1.0e-5
    ct = c * t
    g(ξ) = exp(-ξ^2 / (2 * s^2))
    u_ref(z) = (α / 2) * (g(z - ct) + g(z + ct))
    v_ref(z) = (α * c / (2 * s^2)) * ((z - ct) * g(z - ct) - (z + ct) * g(z + ct))
    z_nodes = model.reference[3, :]
    u_field = u_ref.(z_nodes)
    v_field = v_ref.(z_nodes)

    max_disp = maximum_components(integrator.displacement)
    max_velo = maximum_components(integrator.velocity)
    min_velo = minimum_components(integrator.velocity)
    # The substeps are about three times the 1.0e-7 s step of the reference
    # test, so the tolerances are looser than there. With the pre-fix time
    # compression the maximum displacement was off by about 1.4 percent and
    # the velocity extrema by more than 10 percent.
    @test max_disp[3] ≈ maximum(u_field) rtol = 5.0e-04
    @test max_velo[3] ≈ maximum(v_field) rtol = 5.0e-03
    @test min_velo[3] ≈ minimum(v_field) rtol = 5.0e-03
end
