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
top-level controller and injected into each subdomain, so subdomain input files
may omit `initial time` and `final time`.

## Adaptive time stepping (all schemes)

Adaptive stepping is all-or-nothing: either supply none of the four keys below
(fixed step) or supply all four. The step is grown after successful steps and
shrunk after failed ones, staying within the min/max bounds.

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
number `CFL`.

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
below are a conservative choice.

## Canonical examples

- Quasi-static: `examples/single/static-solid`
- Newmark: `examples/single/implicit-dynamic-solid`
- Central difference: `examples/single/explicit-dynamic-solid`
- Adaptive stepping: `examples/adaptive-time-stepping/`
