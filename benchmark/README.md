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
| `coupled.yaml` (`coupled-clamped.yaml`, `coupled-free.yaml`) | Newmark, 20 steps, matching Dirichlet-Neumann exchange with Anderson acceleration between a clamped block of 24³ elements and a free block of 16³ elements, a nonconforming 3:2 interface; linear elastic, E = 6.895 GPa | 13 824 + 4096 | conjugate gradient |
| `coupled-mono.yaml` | the same bar in one mesh of 48 × 24 × 24 elements | 27 648 | conjugate gradient |

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
- `coupling-profile.jl case.yaml`: the samples of the thread that drives the
  Schwarz iteration attributed inclusively to the source files of Norma and
  to the functions of the coupling layer, and the innermost functions by
  exclusive count, after one run to compile.
- `product.jl case.yaml`: the sparse product of the first reduced tangent
  with 64-bit and 32-bit column indices on all threads, and with the full
  rows and with one stored triangle on one thread.
- `ordering.jl case.yaml`: the sparse Cholesky of the first reduced tangent
  under the default ordering of CHOLMOD, the nested dissection of METIS that
  Norma uses, and the natural ordering: nonzeros of the factor, time of the
  factorization, of the numeric phase with the analysis reused, and of one
  solve.
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

The direct solver and the Cholesky preconditioner under the METIS ordering
(October 24, 2026, 12 threads, the ordering below): implicit24-direct 15.5 →
9.5, implicit24-cg-cholesky 6.5 → 4.9, quasistatic24-direct 24.1 → 13.4,
quasistatic24-cholesky 16.5 → 11.2.

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
| incompressible24, METIS ordering | | | | 8.3 |
| beam, METIS ordering | | | 7.4 | 8.2 |

### Coupled case

`coupled.yaml` on 12 threads (October 24, 2026) runs 20 steps with 15
Schwarz iterations per step and 2 Newton iterations per subdomain solve, 600
subdomain solves in all, against 20 solves of the monolithic bar. Wall
seconds of the timed run, before the changes below:

| Case | 1 thread | 6 threads | 12 threads |
|---|---:|---:|---:|
| coupled | 60.6 | 33.1 | 32.8 |
| coupled-mono | 8.1 | 4.8 | 4.6 |

The coupling layer (`schwarz.jl`, by `coupling-profile.jl`) took 5.0% of the
run on one thread and 12.3% on 12, of which the construction of the
projectors of the nonconforming interface 1.6% and 4.6% (once, 1.5 s on 12
threads), the Neumann datum and the matching criterion 3% and 7% (both
compute the d'Alembert reaction `M a - f_body - f_boundary` with the column
product of the sparse library, 1.1 ms on the 46 875 unknowns of the clamped
block, and the criterion recomputes the row sums of the mass, 2.5 ms), and
the time interpolation of the partner histories 0.7% and 2.1%. The rest of
the run is the subdomain solves; the coupling layer does not limit the
scaling.

Two limits of the solves were found on this case. The row product of the
conjugate gradient ran on one thread below 20 000 unknowns, so the 600
solves of the free block (14 739 unknowns) took 10.1 ms each on 12 threads,
as on one; with the threshold at 4000 rows and the rows divided evenly over
the threads they take 2.4 ms. And Julia's worker threads sleep 100 µs after
a parallel loop, so after each serial part of the Schwarz iteration the
first loop waits for them to wake, and the operating system wakes them on
other cores than the ones that hold their chunks: inside the run the solves
of the clamped block took 18 ms each against 7 ms in isolation, and with
the workers kept spinning (`JULIA_THREAD_SLEEP_THRESHOLD=infinite`) 7.8 ms.
Pinning the threads to cores (`JULIA_EXCLUSIVE=1`) has the same effect on
the coupled run without the spinning, which slows the direct solver (its
BLAS threads compete with the spinning workers). With the even chunks and
the default thread settings the coupled run took 34.8 s, more than the
32.8 s before them: the threaded solves of the free block then wait for the
workers at every pass. Wall seconds on 12 threads with the even chunks:

| Case | default | spinning 0.1 s | spinning | pinned |
|---|---:|---:|---:|---:|
| coupled | 34.8 | 24.1 | 22.8 | 22.3 |
| coupled-mono | 4.6 | | 4.3 | 4.3 |
| implicit24-cg | 2.02 | | 1.65 | 1.62 |
| quasistatic24-cg | 2.77 | | 2.64 | 2.26 |
| quasistatic24-direct | 13.5 | 14.6 | | 13.6 |
| beam-cholesky | 7.5 | 8.1 | | 7.3 |
| explicit48 | 2.46 | 2.56 | | 2.25 |

The coupled run on 6 pinned threads took 32.8 s, as without pinning, and on
one thread 62.7 s. `bin/norma` exports `JULIA_EXCLUSIVE=1` when more than
one thread is requested and the variable is not in the environment, and
`sweep.sh` exports the same value; the measurements earlier in this file
were made without it.

The serial items of the coupling layer were then reduced to the interface
rows (October 24, 2026). The d'Alembert reaction `M a - f_body - f_boundary`
of the Dirichlet side was computed on all its degrees of freedom with the
column product of the sparse library (1.1 ms on 46 875) twice per Schwarz
iteration, for the Neumann datum and for the criterion, and the momentum
rate of the criterion computed the whole momentum and the row sums of the
mass (2.4 ms); every caller kept the interface rows only. They are now
computed on those rows, a few hundred dot products with the rows of the
symmetric mass (`mass_rows_product`, `dalembert_inertia_minus_loads` with
the rows), and the criterion interpolates the partner histories on the
interface rows instead of the whole vectors. The coupled run on 12 pinned
threads: 21.2 to 17.9 s. The roundoff of the reaction changes
with the order of the sums, so the fields differ from the previous commit
at that level; the Schwarz iteration counts are unchanged (15 per step).

The stiffness and the mass are reduced into the model's existing matrices
as well (October 24, 2026, `reduce_element_matrix!`): the values of the
stiffness, 27 MB on the 24³ cube, were allocated at every evaluation. On 12
pinned threads the Newmark cube takes 1.34 s against 1.62 (1.1 GB allocated
against 2.2), the quasi-static cube 1.99 s against 2.26, and the coupled
case 17.7 s against 17.9; the fields of the three cases are bitwise
identical. The tangent predictor copies the converged stiffness of a step,
which it held by reference while every evaluation built a new matrix.

The allocation of the coupled run, 45 GB over 1200 Newton iterations, was
the Newmark Hessian: the values of K + M/(βΔt²) were allocated at every
iteration (28 MB on the clamped block, 8.5 MB on the free one), and the
garbage collection took 14% of the run on 12 threads. Since October 24, 2026
(`newmark_hessian!`) they are written into the Hessian of the previous
iteration when the pattern is unchanged: the Newmark wrapper allocates
2.6 MB per iteration on the clamped block, the coupled run 12 GB with 1.5 s
of collection, and its wall time on 12 pinned threads 21.2 s against 22.3;
the monolithic bar, timed alternately with the previous commit under a
loaded machine, 5.0 against 5.5 s. The fields of both cases are bitwise
identical. The assembly itself allocates 0.3 MB per evaluation and the
conjugate gradient 4 MB per solve in the tasks of its passes.

### Memory traffic of the sparse product

The product of the reduced tangent of the 24³ cube (45 000 unknowns, 3.36
million nonzeros, 27 MB of values and 27 MB of 64-bit column indices) was
timed with `product.jl` on October 24, 2026, threads pinned:

| Variant | threads | ms |
|---|---:|---:|
| full rows, 64-bit indices | 1 | 0.866 |
| full rows, 32-bit indices | 1 | 0.818 |
| one stored triangle | 1 | 0.62 |
| full rows, 64-bit indices | 12 | 0.115 |
| full rows, 32-bit indices | 12 | 0.105 |

The product scales 7.5 times on 12 threads: the matrix fits the two 32 MB
level 3 caches when the threads are pinned, so it is not bound by the
memory bandwidth at this size. The 32-bit indices save 9% of the product,
about 2% of an implicit run, and were not adopted; one stored triangle saves
28% on one thread but its transposed part scatters into every row, so it
does not run on several threads without a second pass, and the full rows on
12 threads are 5 times faster than the triangle on one.

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
under the default ordering. The solve takes 0.02 s, a twentieth of the
factorization, which is why the factor reused as a preconditioner is the
fastest option on the beam.

Ordering of the factor (`ordering.jl`): CHOLMOD's default ordering (the
approximate minimum degree, with its own nested dissection only when the
fill exceeds its threshold), the nested dissection of METIS 5.1 through
Metis.jl, and the natural ordering of the unknowns; nonzeros of the factor
and seconds of the factorization, of the numeric phase, and of one solve.
The METIS ordering takes 0.13 s on the cube and 0.26 s on the beam, once per
pattern.

| Case | ordering | nnz(L) | threads | factorization | numeric | solve |
|---|---|---:|---:|---:|---:|---:|
| quasistatic24 | default | 50.1 M | 1 | 1.47 | 1.35 | 0.021 |
| | | | 6 | 0.55 | 0.46 | 0.020 |
| | | | 12 | 0.53 | 0.41 | 0.019 |
| | METIS | 28.9 M | 1 | 0.50 | 0.44 | 0.014 |
| | | | 6 | 0.27 | 0.20 | 0.013 |
| | | | 12 | 0.28 | 0.22 | 0.013 |
| | natural | 81.6 M | 1 | 1.73 | 1.63 | 0.031 |
| | | | 12 | 0.78 | 0.68 | 0.032 |
| beam | default | 37.9 M | 1 | 0.51 | 0.38 | 0.021 |
| | | | 6 | 0.37 | 0.24 | 0.021 |
| | | | 12 | 0.44 | 0.31 | 0.021 |
| | METIS | 30.5 M | 1 | 0.39 | 0.29 | 0.020 |
| | | | 6 | 0.30 | 0.21 | 0.020 |
| | | | 12 | 0.35 | 0.25 | 0.019 |
| | natural | 1056 M | 1 | 101.6 | 100.2 | 0.380 |
| | | | 12 | 17.0 | 15.6 | 0.330 |

The direct solver and the Cholesky preconditioner use the METIS ordering
since October 24, 2026, computed when the pattern of the reduced matrix
changes, as the symbolic analysis is. Wall
seconds of the runs at 12 threads before and after: quasistatic24-direct
24.1 and 13.4, quasistatic24-cholesky 16.5 and 11.2, implicit24-direct 15.5
and 9.5, implicit24-cg-cholesky 6.5 and 4.9, incompressible24-direct 14.3
and 8.3, beam-direct 9.2 and 8.2, beam-cholesky 7.7 and 7.4. The numeric
factorizations were 85% of the quasi-static direct run on the cube (50
Newton iterations at 0.41 s) and 40% on the beam (12 iterations at 0.31 s),
which bounds the gain of the ordering on each.

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
