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

Norma runs with one thread unless a count is requested:

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

The element loops of the assembly, of the stable time step, and of the blended
energy run on all threads, in chunks of fixed size with one task per thread,
and so do the reductions of the element values into the global vectors and
matrices and the sparse product of the conjugate gradient. The vector updates
of the time integrators, the rest of the conjugate gradient, and the direct
solver (through its BLAS calls) run on one thread. On a 12-core Ryzen 9 9900X
(October 24, 2026), an explicit run of 110 592 hexahedra over 100 steps took
11.6 s on one thread and 2.25 s on twelve pinned threads, and a Newmark run
of 45 000 free degrees of freedom with the conjugate gradient took 4.9 s and
1.34 s; `benchmark/README.md` has the cases and the measurements.
The thread count does not change the results: the chunks and the order of
every sum are fixed, so the fields and the energies are the same on any
number of threads.

Julia's worker threads sleep 100 µs after a parallel loop, and the serial
parts of a step last longer, so every loop that follows one waits for the
workers to wake, and the operating system wakes them on other cores than the
ones that hold their data. The `bin/norma` launcher therefore exports
`JULIA_EXCLUSIVE=1`, which pins the threads to cores, when more than one
thread is requested and the variable is not in the environment. On 12
threads of the Ryzen 9 9900X a coupled run of two subdomains of 45 000 and
14 700 unknowns took 22.3 s pinned against 34.8 s, the Newmark run above
1.62 against 2.02 s at that commit, and the explicit run 2.25 against
2.46 s; the direct
solver lost nothing (`benchmark/README.md`). A run started with `julia`
directly needs the variable set by hand:

```bash
JULIA_EXCLUSIVE=1 julia -t 8 --project=. src/Norma.jl input.yaml
```

The threads are pinned to the first cores, so two threaded runs on one
machine would share them: start the second with `JULIA_EXCLUSIVE=0`, which
the launcher leaves alone. The same applies to a run with fewer threads than
cores on a processor with several core complexes: the pinned threads fill
the first complex and share its cache, and the coupled run above took 28.4 s
on 6 pinned threads against 22.3 s on 6 unpinned ones (and 17.8 s on 12
pinned, 30.5 s on 12 unpinned). Do not confine the process with `taskset` to
as many CPUs as it has threads: the collector and the BLAS threads then
displace the pinned chunks of the loops, and the same run took 116 s.

An automatic count is refused: `-t auto`, `--threads auto` and
`JULIA_NUM_THREADS=auto` all abort the run at startup. On a large node the
automatic count sizes both pools to the hardware thread count, and the level 1
BLAS calls of the linear solver then spend their time synchronizing hundreds of
threads instead of computing. Ask for the number of cores you intend to use.

The binding resizes the pool after OpenBLAS has created it, so a run started
with `julia` directly still owns the idle default threads, and a process
monitor lists them. The `bin/norma` launcher avoids this by exporting
`OPENBLAS_NUM_THREADS` to match the requested count before Julia starts, so
the pool is created at the right size, and always passes `-t` so that a
`JULIA_NUM_THREADS` setting in the environment cannot change the count.

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
