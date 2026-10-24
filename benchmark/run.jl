# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Wall time of one benchmark case:
#   julia -t N --project=.. run.jl case.yaml
# The case runs twice, the first time to compile; the second run is timed,
# with the time spent in garbage collection and the bytes allocated. One line
# per run, machine readable:
#   BENCH threads=N case=... seconds=... gctime=... GB=...
using Norma
input = ARGS[1]
Norma.run(input)
r = @timed Norma.run(input)
println(
    "BENCH threads=",
    Threads.nthreads(),
    " case=",
    input,
    " seconds=",
    round(r.time; digits=3),
    " gctime=",
    round(r.gctime; digits=3),
    " GB=",
    round(r.bytes / 1e9; digits=2),
)
