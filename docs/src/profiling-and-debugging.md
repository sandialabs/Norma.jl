# Profiling and debugging

## Debug logging

Enable debug-level logging with the `JULIA_DEBUG` environment variable:

```bash
JULIA_DEBUG=Norma julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

Add debug messages in code with:

```julia
@debug "Starting simulation with input file: input.yaml"
```

Disable debug printing by unsetting the variable, or suppress it at launch:

```bash
unset JULIA_DEBUG
# or, per invocation:
JULIA_DEBUG= julia --project=. src/Norma.jl input.yaml
```

## Profiling

Norma works with Julia's built-in `Profile` module; results can be visualized
with [`ProfileCanvas.jl`](https://github.com/SciML/ProfileCanvas.jl).

Collect a profile:

```julia
using Profile
include("/path/to/Norma.jl/src/Norma.jl")
cd("/path/to/Norma.jl/examples/contact/implicit-dynamic/2-bars")
@profile Norma.run("bars.yaml")
```

Visualize it:

```julia
using Pkg; Pkg.add("ProfileCanvas")
using ProfileCanvas
ProfileCanvas.canvas()
```

Or as a single command line:

```bash
julia --project=@/path/to/Norma.jl \
  -e 'using Profile; using Norma; cd("examples/contact/implicit-dynamic/2-bars"); @profile Norma.run("bars.yaml")' \
  -E 'using ProfileCanvas; ProfileCanvas.canvas()'
```
