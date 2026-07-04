# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# HHT-α numerical dissipation on the Newmark integrator: `HHT alpha: ᾱ`
# blends the internal force to the state (1-ᾱ) u_{n+1} + ᾱ u_n and overrides
# γ = 1/2 + ᾱ, β = (1+ᾱ)²/4, giving high-frequency dissipation
# ρ∞ = (1-ᾱ)/(1+ᾱ) at second-order accuracy. The test runs the conforming
# overlap cantilever with and without dissipation and asserts the derived
# parameters and a strict energy ordering at the same horizon.

@testset "Newmark HHT Alpha" begin
    example = "../examples/overlap/dynamic-same-step/cantilever-conforming"
    files = [
        "cantilever-impedance.yaml",
        "cantilever-clamped-impedance.yaml",
        "cantilever-free-impedance.yaml",
        "cantilever-clamped.g",
        "cantilever-free.g",
    ]

    function total_energy(sim)
        return sum(s.integrator.stored_energy + s.integrator.kinetic_energy for s in sim.subsims)
    end

    function run_case(hht_alpha)
        for f in files
            cp(joinpath(example, f), f; force=true)
        end
        if hht_alpha > 0.0
            for f in ("cantilever-clamped-impedance.yaml", "cantilever-free-impedance.yaml")
                y = read(f, String)
                y = replace(y, "γ: 0.5" => "γ: 0.5\n  HHT alpha: $hht_alpha")
                write(f, y)
            end
        end
        multi = read("cantilever-impedance.yaml", String)
        # Long enough that the HHT dissipation (tens of percent) dominates
        # the ±15% overlap double-count artifact of this energy sum.
        multi = replace(multi, "final time: 3.0e-4" => "final time: 1.5e-4")
        write("cantilever-impedance.yaml", multi)
        sim = Norma.run("cantilever-impedance.yaml")
        for f in vcat(files, ["cantilever-clamped.e", "cantilever-free.e"])
            rm(f; force=true)
        end
        return sim
    end

    sim_plain = run_case(0.0)
    sim_hht = run_case(0.1)

    @test sim_plain.failed == false
    @test sim_hht.failed == false

    for subsim in sim_hht.subsims
        @test subsim.integrator.hht_alpha ≈ 0.1
        @test subsim.integrator.γ ≈ 0.6
        @test subsim.integrator.β ≈ 0.3025
    end
    for subsim in sim_plain.subsims
        @test subsim.integrator.hht_alpha == 0.0
        @test subsim.integrator.γ ≈ 0.5
    end

    # Numerical dissipation: strictly less energy at the same horizon. The
    # release is broadband, so ρ∞ = 0.818 removes a large fraction of the
    # (mostly high-frequency) energy quickly; the floor only guards against
    # NaN or collapse.
    E_plain = total_energy(sim_plain)
    E_hht = total_energy(sim_hht)
    @test E_hht < E_plain
    @test E_hht > 0.05 * E_plain
end
