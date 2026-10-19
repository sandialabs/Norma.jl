# Time integrators

The `time integrator` block selects the time-stepping scheme and the time
window. The `type` is required and case-sensitive:

| `type` value | Scheme | Use |
|---|---|---|
| `quasi static` | Quasi-static (no inertia) | static/incremental loading |
| `Newmark` | Newmark-β implicit dynamics | implicit transient |
| `central difference` | Central-difference explicit dynamics | explicit transient |

## Time window (all schemes)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `initial time` | yes | — | start time |
| `final time` | yes | — | end time |
| `time step` | yes | — | base (target) time increment Δt |

In a [multidomain](multidomain.md) simulation the time window is set on the
top-level controller and passed to each subdomain, so subdomain input files
may omit `initial time` and `final time`.

## Adaptive time stepping (all schemes)

Adaptive stepping is all-or-nothing: either supply none of the four keys below
(fixed step) or supply all four. The step is increased after successful steps
and reduced after failed ones, within the minimum and maximum bounds.

| Key | Required | Default | Constraint |
|---|---|---|---|
| `minimum time step` | with the group | `time step` | ≤ `maximum time step` |
| `maximum time step` | with the group | `time step` | — |
| `decrease factor` | with the group | `1.0` | ≤ 1.0 |
| `increase factor` | with the group | `1.0` | ≥ 1.0 |

See `examples/adaptive-time-stepping/`.

## `quasi static`

Drops inertia; solves a sequence of static equilibria over the time window.
Pair with the [`Hessian minimizer`](solvers.md) solver.

```yaml
time integrator:
  type: quasi static
  initial time: 0.0
  final time: 1.0
  time step: 0.1
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| `initial equilibrium` | no | `false` | solve for static equilibrium before the first step |
| `predictor` | no | `constant` | the iterate each step's Newton solve starts from: `constant`, `linear` or `tangent` |

### Predictor

Each step starts from the previous solution with the new values of the
prescribed displacements. With the default `constant` predictor, the first
Newton iterate therefore puts the whole increment of a loaded face on the
layer of elements next to it, a strain increment of Δu/h for elements of size
h, however small the strain increment of the body. A path-dependent material
evaluated at that strain can fail where the converged step would not: the flow
rule of [crystal plasticity](materials.md) has a bounded range of slip rates,
and the step is then cut by adaptive time stepping. The other two predictors
start Newton from an iterate that spreads the increment over the body:

- `linear` adds to the free degrees of freedom the increment of the last
  accepted step, scaled by the ratio of the time steps. It costs nothing with
  the Hessian minimizer, and one more evaluation of the residual per step with
  a matrix-free solver (see below). The first step has no increment to
  extrapolate and starts as with `constant`.
- `tangent` adds the solution Δu_f of K_ff Δu_f = −K_fc Δu_c, the linear
  response of the free degrees of freedom to the increment Δu_c of the
  prescribed ones, with the stiffness K of the last accepted step, or, before
  the first step, the stiffness of the initial state. It costs a linear solve
  per step and an assembly on the first step, and requires the
  [`Hessian minimizer`](solvers.md). It responds to the increment of the
  prescribed displacements, not to that of the applied forces.

```yaml
time integrator:
  type: quasi static
  initial time: 0.0
  final time: 5.0e-3
  time step: 5.0e-4
  predictor: tangent
```

The polycrystals of `examples/crystal-plasticity` take LCM's steps of 5 × 10⁻⁴
with `tangent`; their README lists the steps each predictor cuts.

A predictor changes the iterate Newton starts from, not the solution it
converges to. Since it lowers the residual of that iterate, the `relative
tolerance` is measured from the residual the iterate would have without the
predictor, when that residual is the larger: for a linear problem the
predicted iterate is the solution, and a tolerance relative to its residual
would ask for less than the roundoff. The Hessian minimizer estimates that
residual with the stiffness. A matrix-free solver has no stiffness, so it
evaluates the residual of the iterate without the prediction, one more
evaluation per step; if that evaluation fails, the tolerance is measured from
the predicted iterate. A failed step that is retried with a smaller step is
predicted again from the last accepted step. The history starts with the run:
after a restart, a swap, or a switch of an in-memory restart, the first step
is predicted as the first step of a run.

In a [multidomain](multidomain.md) simulation each subdomain has its own
predictor, and the values a Schwarz coupling prescribes on the subdomain's
boundary count as prescribed displacements. Every Schwarz iteration restarts
the subdomain at the start of the stop, and its first step is predicted from
the last step the subdomain took, the one that ended the stop in the previous
iteration. With one step per stop, `linear` therefore starts each Schwarz
iteration from the subdomain's solution of the previous one.

## `Newmark`

Implicit dynamics with the Newmark-β family. Pair with the
[`Hessian minimizer`](solvers.md) solver.

```yaml
time integrator:
  type: Newmark
  initial time: 0.0
  final time: 3.0e-6
  time step: 1.0e-6
  β: 0.25
  γ: 0.5
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| `β` | yes | — | Newmark β (use the Unicode key, not `beta`) |
| `γ` | yes | — | Newmark γ (use the Unicode key, not `gamma`) |
| `HHT α` (or `HHT alpha`) | no | `0.0` | Hilber–Hughes–Taylor α, in `[0, 1/3]`; when positive it overrides γ = 0.5 + α and β = 0.25(1 + α)² for controlled numerical damping |

The common choice `β = 0.25`, `γ = 0.5` is the unconditionally stable,
non-dissipative average-acceleration (trapezoidal) rule.

## `central difference`

Explicit dynamics. Pair with the [`explicit solver`](solvers.md). The stable
time step is estimated from the mesh and wave speed scaled by the Courant
number `CFL`. When the requested `time step` (or, in a multi-domain run, the
controller `time step`) exceeds the stable step, each stop is reached by
substeps no longer than the stable step, and a warning is printed at every
substep.

```yaml
time integrator:
  type: central difference
  initial time: 0.0
  final time: 3.0e-6
  CFL: 0.2
  γ: 0.5
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| `CFL` | yes | — | Courant number used to bound the stable step from the element size and wave speed |
| `γ` | yes | — | central-difference γ (0.5 is the standard second-order value) |

The `CFL` value scales an element-size-and-wave-speed estimate of the stable
step; it is a safety factor rather than the exact operating Courant ratio. The
element-size estimate is a characteristic diameter, so a value near 1.0 can
exceed the true stability limit on hexahedral meshes. Values of about 0.5 or
below are a conservative choice. The wave speed is that of the material's
stiffness in the undeformed state (for the anisotropic models, the largest
over all directions), so a material that stiffens strongly with strain, such
as the exponential fibers of `anisotropic hyperelastic damage`, can exceed it
by orders of magnitude, which the CFL must allow for.

## Canonical examples

- Quasi-static: `examples/single/static-solid`
- Quasi-static with the tangent predictor: `examples/crystal-plasticity/rubiks-cube/rubiks-cube.yaml`
- Newmark: `examples/single/implicit-dynamic-solid`
- Central difference: `examples/single/explicit-dynamic-solid`
- Adaptive stepping: `examples/adaptive-time-stepping/`
