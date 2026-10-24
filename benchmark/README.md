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
| `incompressible24-cg.yaml`, `-direct.yaml` | quasi-static, ν = 0.49, tangent predictor | 13 824 | conjugate gradient, sparse Cholesky |
| `beam-cg.yaml`, `-jacobi.yaml`, `-cholesky.yaml`, `-direct.yaml` | quasi-static bending of a beam of 0.02 × 0.02 × 1 m, tangent predictor | 25 600 | the four options |

The implicit cube cases have 45 000 free degrees of freedom and a stiffness
of 3.36 million nonzeros; the beam 97 119 and 6.73 million. The beam is the
ill-conditioned case: the plain conjugate gradient takes 8056 iterations on
its first system.

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

Coupled cases, one thread, wall seconds by the linear solver of the
subdomains (the examples `nonoverlap/dynamic-same-step/cantilever-dn-nonconforming`
over 200 steps and `nonoverlap/static-same-step/matrix-inclusion-split`):

| Case | `cg` | `cg`, `jacobi` | `cg`, `cholesky` | `direct` |
|---|---:|---:|---:|---:|
| cantilever, 315 and 66 nodes, 1992 solves | 2.45 | 2.02 | 0.97 | 1.72 |
| inclusion split, 703 and 262 nodes, 48 solves | 1.19 | 1.22 | 1.02 | 0.81 |

The coupling layer took 1.5% and 7.6% of these runs; the rest is the
subdomain solves.

Ill-conditioned cases, 12 threads, wall seconds by linear solver:

| Case | `cg` | `cg`, `jacobi` | `cg`, `cholesky` | `direct` |
|---|---:|---:|---:|---:|
| incompressible24, 20 Newton iterations | 3.0 | | | 14.3 |
| beam, 12 Newton iterations | 119.3 | 108.4 | 7.7 | 9.2 |

Algebraic multigrid (AlgebraicMultigrid.jl 2.0.2, tried on October 24, 2026
and not adopted), as the preconditioner of the library conjugate gradient on
the first system of each case: on the beam the Ruge-Stüben and smoothed
aggregation hierarchies cut 8056 iterations to 1373 and 1256 but their cycles
took 63 and 43 s, against 0.44 s for the Cholesky factor and 0.07 s for its
solve; on the quasi-static cube 254 iterations to 52 and 53, with cycles of
1.19 and 0.82 s against 0.24 s for the plain conjugate gradient. The scalar
hierarchies of the library do not use the three-component structure or the
rigid modes of elasticity, and their cycles run on one thread.

Sparse Cholesky (CHOLMOD) against the BLAS thread count, which Norma binds to
the Julia thread count; seconds of the full factorization, of the numeric
phase with the symbolic analysis reused, and of one solve:

| Case | threads | nnz(L) | factorization | numeric | solve |
|---|---:|---:|---:|---:|---:|
| quasistatic24, 44 375 unknowns | 1 | 50.1 M | 1.46 | 1.37 | 0.021 |
| | 6 | | 0.57 | 0.45 | 0.020 |
| | 12 | | 0.53 | 0.42 | 0.020 |
| beam, 97 119 unknowns | 1 | 37.9 M | 0.53 | 0.39 | 0.021 |
| | 6 | | 0.37 | 0.24 | 0.021 |
| | 12 | | 0.44 | 0.32 | 0.021 |

The factorization gains nothing beyond 6 threads (one core complex of the
processor), and the fill of the cube is 15 times the nonzeros of the matrix
under the default ordering. The solve is cheap, which is why the factor
reused as a preconditioner pays on the beam.

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
