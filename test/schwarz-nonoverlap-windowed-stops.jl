# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Windowed controller stops on nonoverlapping Schwarz couplings. A stop is
# windowed when the controller `time step` is a multiple of the subdomain time
# steps, so that each subdomain takes several substeps per Schwarz iteration.
# The relaxation state of the Schwarz iteration is kept per interface and per
# substep time (relaxation_slot! in schwarz.jl), so a windowed run must
# converge to the trajectory of the run whose controller step equals the
# subdomain step. Before the relaxation state was keyed by substep time, one
# vector per interface was blended across the substeps of a stop, which
# shifted the fixed point of the windowed iteration. With Aitken relaxation
# configured, a windowed stop uses the fixed relaxation parameter (Aitken
# applies only to stops with one substep, see aitken_applies) and must reach
# the same trajectory.
#
# Both couplings are run on the conforming cantilever (Newmark, subdomain time
# step 5.0e-7 s) for 10 substeps: the Robin-Robin condition, whose relaxed
# quantity is the Robin datum of the side listed later, and the
# Dirichlet-Neumann condition, whose relaxed quantity is the interface
# displacement of the Dirichlet side. The windowed runs take 5 substeps per
# stop. The Schwarz tolerances are tightened so that the difference between the
# runs measures the windowing, not the Schwarz truncation error.

using YAML
using LinearAlgebra

function run_windowed_cantilever(example; controller_dt=5.0e-7, aitken=false)
    files = ["cantilever-multi.yaml", "cantilever-clamped.yaml", "cantilever-free.yaml",
             "cantilever-clamped.g", "cantilever-free.g"]
    for f in files
        cp("$example/$f", f; force=true)
    end
    params = YAML.load_file("cantilever-multi.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cantilever-multi.yaml"
    params["final time"] = 10 * 5.0e-7
    params["time step"] = controller_dt
    params["Exodus output interval"] = 0
    params["CSV output interval"] = 0
    params["relative tolerance"] = 1.0e-10
    params["absolute tolerance"] = 1.0e-12
    params["relaxation parameter"] = 0.5
    if aitken
        params["relaxation"] = "aitken recursive"
    end
    sim = Norma.run(params)
    for f in vcat(files, ["cantilever-clamped.e", "cantilever-free.e"])
        rm(f; force=true)
    end
    return sim
end

for (label, example) in (
    ("Robin-Robin", "../examples/nonoverlap/dynamic-same-step/cantilever-rr"),
    ("Dirichlet-Neumann", "../examples/nonoverlap/dynamic-same-step/cantilever-dn"),
)
    @testset "Schwarz Nonoverlap Windowed Stops: $label" begin
        sim_ss = run_windowed_cantilever(example)
        sim_win = run_windowed_cantilever(example; controller_dt=2.5e-6)
        sim_win_ait = run_windowed_cantilever(example; controller_dt=2.5e-6, aitken=true)
        for sim in (sim_ss, sim_win, sim_win_ait)
            @test sim.failed == false
        end
        for i in 1:2
            u_ss = sim_ss.subsims[i].model.displacement
            u_win = sim_win.subsims[i].model.displacement
            u_win_ait = sim_win_ait.subsims[i].model.displacement
            @test norm(u_win - u_ss) / norm(u_ss) < 1.0e-6
            @test norm(u_win_ait - u_ss) / norm(u_ss) < 1.0e-6
        end
    end
end
