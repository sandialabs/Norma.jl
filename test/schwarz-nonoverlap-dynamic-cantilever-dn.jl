# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using YAML

const cantilever_dn_example = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"

# Run the single-domain cantilever reference for a given number of steps.
function run_cantilever_single(num_steps)
    for f in ["cantilever.yaml", "cantilever.g"]
        cp("$cantilever_dn_example/$f", f; force=true)
    end
    params = YAML.load_file("cantilever.yaml"; dicttype=Norma.Parameters)
    dt = 5.0e-07
    params["time integrator"]["final time"] = num_steps * dt
    params["name"] = "cantilever.yaml"
    sim = Norma.run(params)
    for f in ["cantilever.yaml", "cantilever.g", "cantilever.e"]
        rm(f; force=true)
    end
    return sim
end

# Run the DN Schwarz cantilever for a given number of steps using the
# specified multi-domain driver YAML.
function run_cantilever_dn(num_steps, driver)
    for f in [driver, "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_dn_example/$f", f; force=true)
    end
    params = YAML.load_file(driver; dicttype=Norma.Parameters)
    dt = 5.0e-07
    params["final time"] = num_steps * dt
    params["name"] = driver
    sim = Norma.run(params)
    for f in [driver, "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g",
              "cantilever-clamped.e", "cantilever-free.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cantilever DN" begin
    num_steps = 10
    sim_ref = run_cantilever_single(num_steps)
    sim_dn = run_cantilever_dn(num_steps, "cantilever-multi.yaml")
    sim_aitken = run_cantilever_dn(num_steps, "cantilever-multi-aitken.yaml")
    sim_pred = run_cantilever_dn(num_steps, "cantilever-multi-predictor.yaml")

    # Reference: get tip displacement (max x node, y-component)
    m_ref = sim_ref.model
    ref_x = m_ref.reference[1, :]
    tip_ref = argmax(ref_x)
    uy_ref = m_ref.displacement[2, tip_ref]

    # DN: get tip displacement from the free subdomain
    m_free = sim_dn.subsims[1].model  # cantilever-free is listed first in domains
    free_x = m_free.reference[1, :]
    tip_free = argmax(free_x)
    uy_dn = m_free.displacement[2, tip_free]

    # Tip displacement should match single domain within 1%
    @test uy_dn ≈ uy_ref rtol = 1.0e-02

    # Interface displacement continuity
    m_clamped = sim_dn.subsims[2].model  # cantilever-clamped is listed second
    u1y = m_clamped.displacement[2, :]
    u2y = m_free.displacement[2, :]
    clamped_x = m_clamped.reference[1, :]
    free_x2 = m_free.reference[1, :]
    # Interface at x = L/2 = 0.127 m
    iface_clamped = findall(x -> abs(x - 0.127) < 1e-4, clamped_x)
    iface_free = findall(x -> abs(x - 0.127) < 1e-4, free_x2)
    avg_clamped = sum(u1y[iface_clamped]) / length(iface_clamped)
    avg_free = sum(u2y[iface_free]) / length(iface_free)
    @test avg_clamped ≈ avg_free rtol = 1.0e-02

    # Schwarz convergence bounded
    iters = sim_dn.controller.schwarz_iters[1:num_steps]
    @test all(iters .≤ 30)

    # Aitken: same physics (reuse tip node index from free subdomain).
    m_free_aitken = sim_aitken.subsims[1].model
    uy_aitken = m_free_aitken.displacement[2, tip_free]
    @test uy_aitken ≈ uy_ref rtol = 1.0e-02
    @test sim_aitken.controller.relaxation_method == :aitken

    # Aitken should reduce Schwarz iterations vs classical θ=0.5 on this problem.
    iters_aitken = sim_aitken.controller.schwarz_iters[1:num_steps]
    @test sum(iters_aitken) < sum(iters)
    @test all(iters_aitken .≤ 10)

    # Interface predictor: same physics, exercises compute_interface_predictor!
    m_free_pred = sim_pred.subsims[1].model
    uy_pred = m_free_pred.displacement[2, tip_free]
    @test uy_pred ≈ uy_ref rtol = 1.0e-02
    @test sim_pred.controller.use_interface_predictor == true
    # Predictor should not make convergence worse than classical relaxation.
    iters_pred = sim_pred.controller.schwarz_iters[1:num_steps]
    @test all(iters_pred .≤ 30)
end
