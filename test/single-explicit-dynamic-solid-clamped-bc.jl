# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Same physical setup as single-implicit-dynamic-solid-clamped-bc.jl, driven with
# central-difference explicit time integration (Δt = 5e-7, CFL = 0.5).
@testset "Single Explicit Dynamic Solid Clamped BC Pulse" begin
    cp("../examples/single/explicit-dynamic-solid/clamped-bc/clamped-bc.yaml", "clamped-bc.yaml"; force=true)
    cp("../examples/single/explicit-dynamic-solid/clamped-bc/clamped-bc.g", "clamped-bc.g"; force=true)
    simulation = Norma.run("clamped-bc.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("clamped-bc.yaml"; force=true)
    rm("clamped-bc.g"; force=true)
    rm("clamped-bc.e"; force=true)

    a = 1.0e-3
    τ = 5.0e-5
    sqrt_e = sqrt(exp(1))
    u_amp = a
    v_amp = a / (τ * sqrt_e)
    a_amp = a / τ^2

    u_eta1 = a / sqrt_e

    z_coords = model.reference[3, :]
    function node_at(z_target)
        idx = argmin(abs.(z_coords .- z_target))
        @assert abs(z_coords[idx] - z_target) < 1.0e-9 "no node at z = $z_target (closest $(z_coords[idx]))"
        return idx
    end
    i_peak = node_at(0.0)
    i_lead = node_at(0.05)
    i_trail = node_at(-0.05)

    uz(i) = integrator.displacement[3 * i]
    vz(i) = integrator.velocity[3 * i]
    az(i) = integrator.acceleration[3 * i]

    # τ/Δt = 100, CFL = 0.5: errors ≲ 0.03% throughout the bar after 1 m of travel.
    @test uz(i_peak)  ≈ u_amp       rtol = 1.0e-06
    @test abs(vz(i_peak))  < 1.0e-03 * v_amp
    @test az(i_peak)  ≈ -a_amp      rtol = 1.0e-03

    @test uz(i_lead)  ≈ u_eta1      rtol = 1.0e-03
    @test vz(i_lead)  ≈ +v_amp      rtol = 1.0e-03
    @test abs(az(i_lead))  < 1.0e-03 * a_amp

    @test uz(i_trail) ≈ u_eta1      rtol = 1.0e-03
    @test vz(i_trail) ≈ -v_amp      rtol = 1.0e-03
    @test abs(az(i_trail)) < 1.0e-03 * a_amp
end
