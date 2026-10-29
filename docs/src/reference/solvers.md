# Solvers

The `solver` block selects the algorithm that solves each time step's system.
The `type` is required and case-sensitive, and must match the time integrator:

| `type` value | Solver | Pair with |
|---|---|---|
| `Hessian minimizer` | Newton's method (Hessian minimizer) | `quasi static`, `Newmark` |
| `explicit solver` | Explicit lumped-mass update | `central difference` |
| `steepest descent` | Matrix-free steepest descent | implicit schemes |

## `Hessian minimizer`

Newton's method with a configurable step and optional backtracking line search.

```yaml
solver:
  type: Hessian minimizer
  step: full Newton
  minimum iterations: 1
  maximum iterations: 16
  relative tolerance: 1.0e-10
  absolute tolerance: 1.0e-06
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| `step` | yes | — | Newton step type (see [Steps](#steps)) |
| `minimum iterations` | yes | — | minimum nonlinear iterations per step |
| `maximum iterations` | yes | — | maximum nonlinear iterations per step |
| `absolute tolerance` | yes | — | absolute residual-norm convergence tolerance |
| `relative tolerance` | yes | — | relative residual-norm convergence tolerance, relative to the residual of the first iterate; with a quasi-static [predictor](time-integrators.md#Predictor), relative to the larger of that residual and the residual the iterate would have without the predictor |
| `unconverged solve action` | no | `fail` (`warn` for matrix-free solvers) | what to do with a solve that exhausts `maximum iterations` without meeting either tolerance: `fail` marks the step failed, so adaptive time stepping retries it with a smaller step or the run aborts; `warn` reports it and continues |
| `linear solver` | no | `cg` | `cg` (conjugate gradient, for symmetric positive definite matrices) or `direct` (sparse Cholesky, CHOLMOD, under the nested-dissection ordering of METIS, with the ordering and the symbolic analysis computed once per pattern; sparse LU, UMFPACK, under its own METIS ordering, when the matrix is not symmetric, as the algorithmic tangent of crystal plasticity is, or not positive definite). Also the solver of the mass for the initial acceleration of `Newmark` |
| `linear solver absolute tolerance` | no | `0.0` | absolute tolerance for the inner CG linear solve |
| `linear solver relative tolerance` | no | `sqrt(eps)` | relative tolerance for the inner CG linear solve. On a linear material Newton converges in one iteration up to this tolerance, and a Newton `relative tolerance` tighter than it (1e-10 against the default 1.5e-8) makes a second iteration solve again: `1.0e-12` then halves the solves (the coupled benchmark took 13.3 s against 17.8 on 12 threads, `benchmark/README.md`) |
| `preconditioner` | no | `none` | preconditioner of the conjugate gradient: `none`, `jacobi` (the inverse diagonal of the reduced matrix), or `cholesky` (the sparse Cholesky factor of the reduced matrix, computed at the first solve and reused across the Newton iterations and steps that follow while the pattern is unchanged and the solve stays within `preconditioner refresh iterations`; a factorization is then computed once per several solves instead of at every one, as with `linear solver: direct`). The conjugate gradient runs on a workspace reused across solves, with every pass over its vectors in fixed chunks on all threads, so its iterates do not depend on the thread count |
| `preconditioner refresh iterations` | no | `10` | with `preconditioner: cholesky`, the iteration count of a solve above which the factor is recomputed at the next solve |

Choice of the linear solver, from the measurements in `benchmark/README.md`
(12 threads of a Ryzen 9 9900X, October 2026): on well-conditioned systems,
a cube of 45 000 unknowns under Newmark or quasi-static loading, the
conjugate gradient is 5 times faster than the direct solver, and the
Jacobi preconditioner saves up to a fifth of its time; on an ill-conditioned
system, a slender beam of 97 000 unknowns in bending, the conjugate gradient
takes thousands of iterations and the direct solver or the reused Cholesky
factor as preconditioner is 14 to 16 times faster; on small subdomains solved many
times, as in the Schwarz couplings, the reused factor is the fastest option.
Algebraic multigrid was measured and not adopted (`benchmark/README.md`). The
sparse Cholesky gains nothing beyond 6 BLAS threads on that processor. Its
factor is ordered by the nested dissection of METIS, which on the cube has
28.9 million nonzeros against 50.1 million under the default ordering of
CHOLMOD, and halves the factorization time on 12 threads; the ordering itself
takes 0.13 s on the cube and 0.26 s on the beam, once per pattern. The
sparse LU takes the METIS ordering of UMFPACK instead of its default (AMD):
on a crystal-plasticity polycrystal of 388 661 unknowns, its factor has 1.79
billion nonzeros against 4.43 billion and takes 148 s against 593 to 709 s
on a 64-core EPYC (issue #251).
| `use line search` | no | `false` | enable backtracking line search; a trial step whose evaluation fails (a non-positive Jacobian, or a constitutive failure such as the slip-rate cap of crystal plasticity) fails the solve, so adaptive time stepping cuts the step |
| `line search backtrack factor` | no | `0.5` | step-reduction factor per backtrack |
| `line search decrease factor` | no | `1.0e-04` | Armijo sufficient-decrease coefficient |
| `line search maximum iterations` | no | `16` | maximum backtracking iterations |

## `explicit solver`

Single-shot explicit update for central-difference dynamics; no nonlinear
iteration or tolerances. Requires only a `step`, which is `explicit`.

```yaml
solver:
  type: explicit solver
  step: explicit
```

## `steepest descent`

Matrix-free steepest descent with a mandatory line search.

```yaml
solver:
  type: steepest descent
  step: steepest descent
  minimum iterations: 1
  maximum iterations: 32
  relative tolerance: 1.0e-10
  absolute tolerance: 1.0e-06
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| `step` | yes | — | `steepest descent` |
| `minimum iterations` | yes | — | minimum iterations per step |
| `maximum iterations` | yes | — | maximum iterations per step |
| `absolute tolerance` | yes | — | absolute convergence tolerance |
| `relative tolerance` | yes | — | relative convergence tolerance |
| `line search backtrack factor` | no | `0.5` | step-reduction factor per backtrack |
| `line search decrease factor` | no | `1.0e-04` | Armijo sufficient-decrease coefficient |
| `line search maximum iterations` | no | `16` | maximum backtracking iterations |
| `energy stagnation window` | no | `0` (disabled) | mesh smoothing only: window length for the energy stagnation exit |
| `energy stagnation tolerance` | no | `1.0e-06` | mesh smoothing only: relative energy decrease below which a window counts as stalled |

### Energy stagnation exit (mesh smoothing)

A mesh being smoothed can reach the energy floor its topology permits long
before the residual tolerances are met: near that floor the solver keeps
reducing the gradient without reducing the energy, since the remedy is
topology modification, not more iterations. The `adaptivity` block of the
[model reference](model.md) supplies it, and ends every smoothing phase of its
loop on this exit, so its examples set the window. With
`energy stagnation window` set to ``W > 0``, the solve stops once the energy
decrease over the last ``W`` iterations, relative to the current energy, stays
below `energy stagnation tolerance` for ``W`` consecutive iterations. The
energy still recoverable past that point is of the order of the tolerance
times the energy itself. When the criterion is met the step is accepted as converged, a
message advising topology modification is logged, and the solver's
`stagnated` flag is raised. The criterion applies to `mesh smoothing` models
with the `steepest descent` solver only; setting the keys anywhere else
aborts. See `examples/ems/awful-cube/awful-cube-lbfgs.yaml`, and for the
exit inside the adaptivity loop `examples/ems/awful-cube/awful-cube-adaptive.yaml`.

## Steps

The `step` key inside `solver` names the update direction. It is required and
case-sensitive:

| `step` value | Step | Notes |
|---|---|---|
| `full Newton` | Newton step from the tangent (Hessian) | for `Hessian minimizer` |
| `explicit` | Explicit lumped-mass update | for `explicit solver` |
| `steepest descent` | Negative-gradient step | for `steepest descent` |
| `lbfgs` | Limited-memory BFGS quasi-Newton step | for `Hessian minimizer` |

| Key | Required | Default | Applies to | Meaning |
|---|---|---|---|---|
| `step length` | no | `1.0` | all steps | scaling applied to the computed step |
| `memory` | no | `10` | `lbfgs` | number of stored curvature pairs |

## Canonical examples

- Newton (implicit): `examples/single/implicit-dynamic-solid`
- Explicit: `examples/single/explicit-dynamic-solid`
- Quasi-static Newton: `examples/single/static-solid`
