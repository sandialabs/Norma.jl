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

# Run the DN Schwarz cantilever for a given number of steps.
function run_cantilever_dn(num_steps)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g"]
        cp("$cantilever_dn_example/$f", f; force=true)
    end
    params = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
    dt = 5.0e-07
    params["final time"] = num_steps * dt
    params["name"] = "cantilever-multi.yaml"
    sim = Norma.run(params)
    for f in ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
              "cantilever-clamped.g", "cantilever-free.g",
              "cantilever-clamped.e", "cantilever-free.e"]
        rm(f; force=true)
    end
    return sim
end

@testset "Schwarz Nonoverlap Dynamic Cantilever DN" begin
    num_steps = 10
    sim_ref = run_cantilever_single(num_steps)
    sim_dn = run_cantilever_dn(num_steps)

    # Reference: get tip displacement (max x node, y-component)
    m_ref = sim_ref.model
    ref_x = m_ref.reference[1, :]
    tip_ref = argmax(ref_x)
    uy_ref = m_ref.current[2, tip_ref] - m_ref.reference[2, tip_ref]

    # DN: get tip displacement from the free subdomain
    m_free = sim_dn.subsims[1].model  # cantilever-free is listed first in domains
    free_x = m_free.reference[1, :]
    tip_free = argmax(free_x)
    uy_dn = m_free.current[2, tip_free] - m_free.reference[2, tip_free]

    # Tip displacement should match single domain within 1%
    @test uy_dn ≈ uy_ref rtol = 1.0e-02

    # Interface displacement continuity
    m_clamped = sim_dn.subsims[2].model  # cantilever-clamped is listed second
    u1y = m_clamped.current[2, :] .- m_clamped.reference[2, :]
    u2y = m_free.current[2, :] .- m_free.reference[2, :]
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
end
