# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

const cantilever_dn_secant_example = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"

# Monolithic single-domain reference run.
function run_cantilever_single_secant(num_steps)
    for f in ["cantilever.yaml", "cantilever.g"]
        cp("$cantilever_dn_secant_example/$f", f; force=true)
    end
    params = YAML.load_file("cantilever.yaml"; dicttype=Norma.Parameters)
    params["time integrator"]["final time"] = num_steps * 5.0e-07
    params["name"] = "cantilever.yaml"
    sim = Norma.run(params)
    for f in ["cantilever.yaml", "cantilever.g", "cantilever.e"]
        rm(f; force=true)
    end
    return sim
end

# DN Schwarz cantilever (Newmark, displacement-form) with the requested
# relaxation method.
function run_cantilever_dn_secant(num_steps, relaxation)
    driver = "cantilever-multi-aitken.yaml"
    for f in [driver, "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_dn_secant_example/$f", f; force=true)
    end
    params = YAML.load_file(driver; dicttype=Norma.Parameters)
    params["final time"] = num_steps * 5.0e-07
    params["relaxation"] = relaxation
    params["name"] = driver
    sim = Norma.run(params)
    for f in [driver, "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g",
              "cantilever-clamped.e", "cantilever-free.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cantilever DN Aitken Secant" begin
    num_steps = 10
    sim_ref = run_cantilever_single_secant(num_steps)
    sim_secant = run_cantilever_dn_secant(num_steps, "aitken secant")
    sim_aitken = run_cantilever_dn_secant(num_steps, "aitken recursive")

    @test sim_secant.controller.relaxation_method == :aitken_secant
    @test sim_secant.failed == false

    # Physical correctness: the Aitken-secant DN tip displacement matches the
    # monolithic single-domain reference. This validates the d-form recovery path
    # (recover_interface_kinematics! Newmark branch), which the static test does
    # not reach. (Note: Aitken secant and Aitken recursive legitimately differ in
    # the transient interior because they impose different interface velocity /
    # acceleration — recovered-consistent vs independently relaxed — so they are
    # compared to the monolithic reference, not to each other.)
    m_ref = sim_ref.model
    uy_ref = m_ref.displacement[2, argmax(m_ref.reference[1, :])]
    m_free = sim_secant.subsims[1].model
    uy_secant = m_free.displacement[2, argmax(m_free.reference[1, :])]
    @test uy_secant ≈ uy_ref rtol = 1.0e-02

    # Interface displacement continuity between the free and clamped subdomains
    # (interface at x = L/2 = 0.127 m).
    m_clamped = sim_secant.subsims[2].model
    iface_free = findall(x -> abs(x - 0.127) < 1e-4, m_free.reference[1, :])
    iface_clamped = findall(x -> abs(x - 0.127) < 1e-4, m_clamped.reference[1, :])
    avg_free = sum(m_free.displacement[2, iface_free]) / length(iface_free)
    avg_clamped = sum(m_clamped.displacement[2, iface_clamped]) / length(iface_clamped)
    @test avg_free ≈ avg_clamped rtol = 1.0e-02

    # Sweep counts. The secant and the Aitken-recursive Irons-Tuck form cost
    # essentially the same here, 41 Schwarz iterations against 40 over these 10 steps, the
    # difference being one extra iteration in the first step. The secant used to
    # take 20, but that margin came from the degenerate first pair fixed for
    # issue #219: its zero factor discarded the first interface iterate, which
    # happened to help this problem while silently returning a wrong answer on
    # subdomain chains. The bound is the recursive count plus that one iteration, so
    # a genuine regression in the secant factor still fails it.
    iters_secant = sim_secant.controller.schwarz_iters[1:num_steps]
    iters_aitken = sim_aitken.controller.schwarz_iters[1:num_steps]
    @test all(iters_secant .≤ 10)
    @test sum(iters_secant) ≤ sum(iters_aitken) + 1
end
