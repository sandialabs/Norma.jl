# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Hash of the final displacement, velocity, and acceleration of a case, with
# its energies, to compare the results of two versions of the code or of two
# thread counts bit for bit:
#   julia -t N --project=.. fingerprint.jl case.yaml
using Norma, SHA
input = ARGS[1]
sim = Norma.run(input)
sims = sim isa Norma.MultiDomainSimulation ? sim.subsims : [sim]
for (i, s) in enumerate(sims)
    ig = s.integrator
    h = bytes2hex(
        sha256(
            vcat(
                reinterpret(UInt8, ig.displacement),
                reinterpret(UInt8, ig.velocity),
                reinterpret(UInt8, ig.acceleration),
            ),
        ),
    )
    println(
        "FP $(basename(input)) sub=$i fields=$(h[1:16]) stored=$(ig.stored_energy) kinetic=$(hasproperty(ig, :kinetic_energy) ? ig.kinetic_energy : 0.0) t=$(ig.time)",
    )
end
