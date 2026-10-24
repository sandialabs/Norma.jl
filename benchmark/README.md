# Benchmarks

Cases and scripts that measure the time of a Norma run against the number of
threads, with the components of a step and the bytes allocated. They were
written for the thread-scaling study of October 2026 and are rerun after every
change of the assembly, the reductions, or the linear solvers.

## Cases

All cases are neo-Hookean cubes of side 0.1 m (E = 1 GPa, ν = 0.3,
ρ = 1000 kg/m³), fixed on one face, meshed with hexahedra by `cube.jou`:
`./meshes.sh [cubit]` writes `cube24.g` (13 824 elements, 15 625 nodes) and
`cube48.g` (110 592 elements, 117 649 nodes).

| Input | Integrator | Elements | Linear solver |
|---|---|---:|---|
| `explicit48.yaml` | central difference, 100 steps | 110 592 | none |
| `implicit24-cg.yaml` | Newmark, 10 steps, 30 Newton iterations | 13 824 | conjugate gradient |
| `implicit24-cg-jacobi.yaml` | same | 13 824 | conjugate gradient, Jacobi |
| `implicit24-cg-cholesky.yaml` | same | 13 824 | conjugate gradient, reused Cholesky factor |
| `implicit24-direct.yaml` | same | 13 824 | sparse Cholesky |
| `quasistatic24-cg.yaml` | quasi-static, 10 steps, 50 Newton iterations | 13 824 | conjugate gradient |
| `quasistatic24-jacobi.yaml` | same | 13 824 | conjugate gradient, Jacobi |
| `quasistatic24-cholesky.yaml` | same | 13 824 | conjugate gradient, reused Cholesky factor |
| `quasistatic24-direct.yaml` | same | 13 824 | sparse Cholesky |

The implicit cases have 45 000 free degrees of freedom and a stiffness of
3.36 million nonzeros.

## Scripts

- `run.jl case.yaml`: runs the case twice, the first time to compile, and
  prints the wall time of the second run, its time in garbage collection, and
  the bytes it allocated.
- `sweep.sh ["1 2 4 8 12"] [cases]`: one Julia process per case and thread
  count, appended to `results.txt` with the date, the commit, and the host.
- `components.jl case.yaml`: the time of each component of one step on the
  initialized model (the evaluation, the reductions, the stable step, the
  kinematic updates, the reduced system, the linear solve), as the minimum of
  five calls after one call to compile.
- `fingerprint.jl case.yaml`: a hash of the final kinematic fields with the
  energies, to compare two versions of the code or two thread counts bit for
  bit. The fields and the energies do not depend on the thread count; the
  energies and the conjugate gradient iterates differ at roundoff level from
  the versions before October 2026.

Run every script from this directory with `--project=..`, and the sweep with
the machine otherwise idle: a busy core delays every static chunk of the
loops.

## Measurements

Ryzen 9 9900X (12 cores, 2 × 32 MB L3, dual-channel DDR5), Julia 1.13.0,
wall seconds of the timed run.

| Case, threads | 2026-10-23, before | 2026-10-24, after |
|---|---:|---:|
| explicit48, 1 | 41.8 | 12.5 |
| explicit48, 12 | 19.1 | 2.4 |
| implicit24-cg, 1 | 6.1 | 4.9 |
| implicit24-cg, 12 | 2.7 | 2.0 |
| implicit24-cg-jacobi, 12 | | 1.8 |
| implicit24-direct, 12 | 15.9 | 15.5 |
| quasistatic24-cg, 12 | | 2.45 |
| quasistatic24-direct, 12 | | 24.1 |

Before: garbage collection took 52% of the explicit run on 12 threads (549 MB
allocated per step), the stable time step read the untyped shape function
tables inside its loop, the per-thread accumulators shared cache lines, and
the reductions, the reduced system, and the vector operations of the
conjugate gradient ran on one thread. After: the explicit step allocates
29 MB, the element loops run in fixed chunks pinned to the threads, the
reductions gather on all threads, the reduced system is a cached gather, and
the conjugate gradient runs on a reused workspace with Jacobi or Cholesky
preconditioning. The element kernel of the explicit case takes 93 ms on one
thread and 10.7 ms on twelve.
