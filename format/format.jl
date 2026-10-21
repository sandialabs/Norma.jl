# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Formats the Julia sources of Norma with the version of JuliaFormatter pinned
# in format/Manifest.toml and the style in .JuliaFormatter.toml:
#
#     julia --project=format format/format.jl          # rewrite the files
#     julia --project=format format/format.jl --check  # exit 1 if any file would change
#
# The check is what CI runs, so a change is accepted only when its files are
# formatted with this version.

using JuliaFormatter

const ROOT = dirname(@__DIR__)
const PATHS = [joinpath(ROOT, p) for p in ("src", "test", "ext", "examples", "studies", "docs", "format")]

check = "--check" in ARGS
unformatted = String[]
for path in PATHS
    isdir(path) || continue
    for (dir, _, files) in walkdir(path), file in files
        endswith(file, ".jl") || continue
        full = joinpath(dir, file)
        format(full; overwrite=(!check)) || push!(unformatted, relpath(full, ROOT))
    end
end
if check && !isempty(unformatted)
    println(stderr, "These files are not formatted; run `julia --project=format format/format.jl`:")
    foreach(f -> println(stderr, "  ", f), unformatted)
    exit(1)
end
