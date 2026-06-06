# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Drive z = -L/2 with a Gaussian-in-time displacement and verify the right-
# traveling pulse against the 1D wave analytical solution u(z, t) = g(t - (z + L/2)/c).
# Bar: 1 mm × 1 mm × 1 m, E = 1e9, ν = 0, ρ = 1000, so c = sqrt(E/ρ) = 1000 m/s.
# Driver: g(t) = a*exp(-((t-t_c)/τ)^2/2) with a = 1e-3, t_c = 2.5e-4, τ = 5e-5.
# Final time t_f = 7.5e-4 places the pulse peak (η = 0) at z = c·(t_f - t_c) - L/2 = 0,
# strictly before reflection (first arrival at z = +0.5 is at t = t_c + L/c = 1.25e-3).
# At t_f the references at z ∈ {-0.05, 0, +0.05} sit on exact nodes (h = 1 mm) and
# correspond to η ∈ {+1, 0, -1}, sampling both the displacement peak and the velocity
# peaks on either side of it.
@testset "Single Implicit Dynamic Solid Clamped BC Pulse" begin
    cp("../examples/single/implicit-dynamic-solid/clamped-bc/clamped-bc.yaml", "clamped-bc.yaml"; force=true)
    cp("../examples/single/implicit-dynamic-solid/clamped-bc/clamped-bc.g", "clamped-bc.g"; force=true)
    simulation = Norma.run("clamped-bc.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("clamped-bc.yaml"; force=true)
    rm("clamped-bc.g"; force=true)
    rm("clamped-bc.e"; force=true)

    a = 1.0e-3
    τ = 5.0e-5
    sqrt_e = sqrt(exp(1))
    u_amp = a                       # peak displacement
    v_amp = a / (τ * sqrt_e)        # peak |∂u/∂t|, attained at η = ±1
    a_amp = a / τ^2                 # peak |∂²u/∂t²|, attained at η = 0

    u_eta1 = a / sqrt_e             # u at η = ±1

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

    # τ/Δt = 25 with Newmark (γ=0.5, β=0.25) gives ≲ 0.3% error in u, v, a after
    # the pulse traverses 1 m.  Tolerances below are 2-3× the observed errors.
    @test uz(i_peak)  ≈ u_amp       rtol = 1.0e-04
    @test abs(vz(i_peak))  < 0.01 * v_amp     # 1% of velocity amplitude
    @test az(i_peak)  ≈ -a_amp      rtol = 2.0e-03

    @test uz(i_lead)  ≈ u_eta1      rtol = 5.0e-03
    @test vz(i_lead)  ≈ +v_amp      rtol = 5.0e-03
    @test abs(az(i_lead))  < 0.01 * a_amp     # 1% of acceleration amplitude

    @test uz(i_trail) ≈ u_eta1      rtol = 5.0e-03
    @test vz(i_trail) ≈ -v_amp      rtol = 5.0e-03
    @test abs(az(i_trail)) < 0.01 * a_amp
end
