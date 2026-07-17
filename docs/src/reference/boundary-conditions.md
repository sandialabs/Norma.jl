# Boundary conditions

The `boundary conditions` block groups conditions by type. Each type key maps to
a **list** of entries:

```yaml
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsx-
      component: y
      function: "0.0"
  Neumann:
    - side set: ssz+
      component: z
      function: "1.0e+09 * t"
```

The block is optional; a model with no `boundary conditions` has none. This page
covers the single-domain condition types. The Schwarz coupling condition types
(`Schwarz overlap`, `Schwarz DN nonoverlap`, `Schwarz impedance nonoverlap`,
`Schwarz impedance overlap`, `Schwarz contact`) are documented under
[Multidomain and Schwarz coupling](multidomain.md).

The `component` value is always one of `x`, `y`, or `z`. The `function` value is
a [function expression](functions.md) in `t, x, y, z`. Geometry is selected by a
`node set` or `side set` name from the Exodus mesh, as noted per type.

## `Dirichlet`

Prescribes a displacement component. Velocity and acceleration follow by
automatic time differentiation (see [Function expressions](functions.md)).

| Key | Required | Meaning |
|---|---|---|
| `node set` or `side set` | one required | geometry the condition acts on |
| `component` | yes | `x`, `y`, or `z` |
| `function` | yes | prescribed displacement |

```yaml
Dirichlet:
  - node set: nsz+
    component: z
    function: "0.01 * t"
```

## `Neumann`

Applies a traction component over a side set.

| Key | Required | Meaning |
|---|---|---|
| `side set` | yes | surface the traction acts on |
| `component` | yes | traction direction `x`, `y`, or `z` |
| `function` | yes | traction value |

```yaml
Neumann:
  - side set: ssz+
    component: z
    function: "1.0e+09 * t"
```

## `Neumann pressure`

Applies a scalar follower pressure normal to a side set. There is no
`component` — the direction is the surface normal.

| Key | Required | Meaning |
|---|---|---|
| `side set` | yes | surface the pressure acts on |
| `function` | yes | scalar pressure (positive acts along the outward normal) |

```yaml
Neumann pressure:
  - side set: interior_ss_bottom
    function: "-937.5e3 * 1500 * (0.5 - 0.5 * cos(pi * t / 1500.0))"
```

## `Robin`

A mixed condition `traction + (robin parameter)·displacement = function` on a
side-set component.

| Key | Required | Meaning |
|---|---|---|
| `side set` | yes | surface the condition acts on |
| `component` | yes | `x`, `y`, or `z` |
| `function` | yes | right-hand side |
| `robin parameter` | yes | coefficient on the displacement term (nonzero) |

## `Surface`

An inclined-support (roller) condition: the constrained nodes are held on the
level surface `g(x, y, z) = 0` given by `function`, free to slide within it. The
constraint direction is the surface gradient ∇g, computed automatically; a
constant `function` (no gradient) is rejected.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `side set` | yes | — | surface whose nodes are constrained |
| `function` | yes | — | level-set expression g in `x, y, z` |
| `enforcement` | no | `exact` | `exact` or `penalty` |
| `penalty` | no | `1.0e6` | penalty coefficient (positive; used when `enforcement: penalty`) |

## Canonical examples

- Fixed and time-dependent Dirichlet: `examples/single/static-solid`,
  `examples/single/implicit-dynamic-solid`
- Neumann traction: single-domain examples under `examples/single/`
- Neumann pressure: pressurized cases under `examples/ahead/`
