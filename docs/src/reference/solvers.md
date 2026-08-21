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
| `relative tolerance` | yes | — | relative residual-norm convergence tolerance |
| `linear solver absolute tolerance` | no | `0.0` | absolute tolerance for the inner CG linear solve |
| `linear solver relative tolerance` | no | `sqrt(eps)` | relative tolerance for the inner CG linear solve |
| `use line search` | no | `false` | enable backtracking line search |
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
topology modification with an external tool, not more iterations. With
`energy stagnation window` set to ``W > 0``, the solve stops once the energy
decrease over the last ``W`` iterations, relative to the current energy, stays
below `energy stagnation tolerance` for ``W`` consecutive iterations. The
energy still recoverable past that point is of the order of the tolerance
times the energy itself. On trigger the step is accepted as converged, a
message advising topology modification is logged, and the solver's
`stagnated` flag is raised. The criterion applies to `mesh smoothing` models
with the `steepest descent` solver only; setting the keys anywhere else
aborts. See `examples/ems/awful-cube/awful-cube-lbfgs.yaml`.

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
