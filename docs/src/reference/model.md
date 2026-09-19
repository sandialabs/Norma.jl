# Model

The `model` block defines the physics, the material assignment, and optional
field-recovery and mesh-quality settings. For a solid-mechanics simulation the
type is `solid mechanics`:

```yaml
model:
  type: solid mechanics
  material:
    blocks:
      body: steel
    steel:
      model: linear elastic
      elastic modulus: 6.895e+09
      Poisson's ratio: 0.25
      density: 2768.0
```

## `type`

Required string. Selects the model. For full-order solid mechanics use
`solid mechanics`. The value `mesh smoothing` runs the mesh-smoothing model
(and implicitly enables the `mesh smoothing` option below). The reduced-order
model types are documented separately with the `norma-opinf` companion package
and are outside the scope of this guide.

## `material`

Required. Assigns a constitutive model to each mesh element block. See
[Materials](materials.md) for the full specification of `blocks` and the
per-material parameters.

## Integration points

| Key | Required | Default | Meaning |
|---|---|---|---|
| `num integration points` | no | element default | map of element-block name → integer quadrature-point count, overriding the default rule for that block |

```yaml
model:
  type: solid mechanics
  num integration points:
    body: 8
  material:
    ...
```

## Nodal recovery

Optional. Enables L2-projection recovery of element quantities to the nodes for
output. Omitting the block disables recovery.

| Key | Required | Default | Values / meaning |
|---|---|---|---|
| `method` | yes (if block present) | — | `lumped`, `consistent`, or `both` (projection mass matrix) |
| `stress` | no | `true` | recover the stress tensor |
| `von mises stress` | no | `false` | recover the von Mises stress |
| `internal variables` | no | `false` | recover material internal variables |
| `deformation gradient` | no | `false` | recover the deformation gradient |

```yaml
model:
  type: solid mechanics
  nodal recovery:
    method: lumped
    stress: true
    von mises stress: true
  material:
    ...
```

The legacy keys `stress recovery` and `recover internal variables` are no longer
accepted and cause an error; use the `nodal recovery` block above.

## Mesh smoothing

| Key | Required | Default | Meaning |
|---|---|---|---|
| `mesh smoothing` | no | `false` | enable mesh smoothing (set automatically when `type: mesh smoothing`); a top-level key |
| `smooth reference` | no | `""` | reference-metric rule for TETRA4 smoothing: `equal volume`, `average edge length`, `max`, `size field`, `size field unrestricted`, `metric field`, or `metric field unrestricted` |
| `size field` | required if `smooth reference: size field`/`size field unrestricted` | `nothing` | expression in `t, x, y, z` giving the target element size |
| `metric field` | required if `smooth reference: metric field`/`metric field unrestricted` | `nothing` | block with `sizes` (three expressions in `t, x, y, z`: the principal sizes h₁, h₂, h₃) and optional `rotation vector` (three expressions: the rotation vector whose exponential carries the global axes onto the principal directions); see below |

### Anisotropic smoothing with a metric field

`smooth reference: metric field` replaces the scalar target size by a metric
tensor `M = R diag(1/h₁², 1/h₂², 1/h₃²) Rᵀ`, where the principal sizes
`h_i` and the rotation `R = exp(hat(v))` come from the `sizes` and
`rotation vector` expressions (radians; identity when omitted). The ideal
element of each tetrahedron is the unit regular tetrahedron scaled by `h_i`
along the global axes and rotated by `R`, and the smoothing energy is evaluated
on the deformation gradient measured in the metric, `F_U = F_M F F_M⁻¹` with
`F_M = diag(1/h_i) Rᵀ`, so the orientation of an elongated element is part of
the objective, which an isotropic energy on `F` alone cannot see. With equal
sizes the rule reduces exactly to `size field`. The field is evaluated at the
centroid of each element in the original mesh and held fixed during a solve,
as for `size field`. `metric field` scales the three sizes uniformly so the
ideal volume is never below the volume of the original element (the
anisotropic form of the `size field` floor); `metric field unrestricted` does
not. The sampled sizes and rotation vector are written as the nodal variables
`size_1..3` and `rotation_1..3`. Formulation and tests:
`docs/notes/ems-anisotropic`.

```yaml
model:
  type: mesh smoothing
  smooth reference: metric field
  metric field:
    sizes: ["0.025", "0.1", "0.1"]
    rotation vector: ["0.0", "0.0", "atan(y, x)"]   # first principal direction radial
```

Mesh smoothing is a specialized capability; most simulations omit these keys.
See `examples/ems/` for smoothing cases.

## Canonical examples

- Basic model/material block: `examples/single/static-solid`
- Mesh smoothing: `examples/ems/cube/cube.yaml`
