# Running Norma

Norma runs a simulation from a single YAML [input file](reference/index.md).

## Command line

The self-activating wrapper is the simplest way to run:

```bash
bin/norma input.yaml
```

Equivalently, invoke Julia directly:

```bash
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

## Interactive session

```julia
using Pkg
Pkg.activate("/path/to/Norma.jl")
using Norma
Norma.run("input.yaml")
```

A typical workflow navigates to an example directory first:

```julia
cd("examples/contact/implicit-dynamic/2-bars")
Norma.run("bars.yaml")
```

If you change Norma's source, reload the module (`using Norma`) for the changes
to take effect.

## Multiple threads

```bash
bin/norma input.yaml --threads 8
```

or via the environment variable:

```bash
JULIA_NUM_THREADS=8 julia --project=. src/Norma.jl input.yaml
```

The count applies to the whole run, not only to the Julia code. Norma binds the
BLAS thread pool to the same number at startup, so the dense kernels reached
through the linear solver, the consistent stress recovery and the Schwarz
projectors stay within the requested budget instead of defaulting to half the
hardware threads. Without a thread flag the run is serial in both pools.

To size the two pools independently, set `OPENBLAS_NUM_THREADS` (or
`OMP_NUM_THREADS`) in the environment. Either variable disables the binding and
Norma leaves the library at the requested value:

```bash
OPENBLAS_NUM_THREADS=4 bin/norma input.yaml --threads 8
```

The thread configuration in effect is reported at the start of every run.

## Reduced-order models

Running Operator-Inference reduced-order models is a separate three-step
workflow (generate snapshots with a full-order run, build the operator with
`norma-opinf`, then run in reduced-order mode). It is documented with the
[`norma-opinf`](https://github.com/sandialabs/norma-opinf) companion package and
the example READMEs under `examples/ahead/`, and is outside the scope of this
guide.

## Next steps

- Learn the input file structure in the [Input File Reference](reference/index.md).
- Browse runnable inputs in [Examples](examples.md).
- Enable debug logging or profile a run in [Profiling and debugging](profiling-and-debugging.md).
