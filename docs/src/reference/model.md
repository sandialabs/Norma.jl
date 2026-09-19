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
| `smooth reference` | no | `""` | rule for the ideal element of TETRA4 smoothing: `size field`, `size field unrestricted`, `metric field`, or `metric field unrestricted`, which prescribe the target; `equal volume`, `average edge length`, and `max`, which take the target from the original element, are legacy rules kept for older inputs |
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

### Adaptivity: smoothing alternated with topological operations

A top-level `adaptivity` block turns a mesh smoothing run into the coupled
loop of `docs/notes/ems-adaptivity`: the mesh is smoothed, the interior edges
of the elements of highest energy density are swapped where a swap lowers the
energy of the elements around the edge, edges are collapsed where a collapse
does (a node is removed only onto a node that carries its node sets and lies
on its surfaces, along a boundary edge if it is on the boundary), edges are
split where a split does, the new mesh is written as
`<output name>-adapted-<k>.g` and smoothed again, and so on until a topology
phase accepts no operation or the outer iterations are exhausted. Every
operation is accepted by one test: the energy of the new elements must be
below that of the old ones by the relative margin, no new element may exceed
the allowed density, and the worst new element must respect the geometric
floor. The energy sums the elements of the cavity, so without the floor an
operation can improve the sum while creating one flat element; the floor
is stated in the scaled Jacobian that the analysis codes require. The mesh must consist of four-node tetrahedra.
Node sets and side sets are carried over to the written meshes.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `desired energy density` | no | `0.1` | elements above this energy per unit ideal volume are candidates |
| `allowed energy density` | no | `Inf` | no accepted operation may create an element above this |
| `minimum scaled Jacobian` | no | `0` | geometric floor: no accepted operation may create an element with a scaled Jacobian below this, unless the worst element of the cavity was already below it and the new worst is no worse |
| `minimum decrease` | no | `1.0e-8` | relative decrease of the cavity energy an operation must achieve |
| `adjacency layers` | no | `4` | rings of adjacent elements added to the candidate set |
| `maximum passes` | no | `20` | passes of operations per topology phase |
| `outer iterations` | no | `5` | alternations of smoothing and topology |
| `swaps` | no | `true` | try edge swaps |
| `collapses` | no | `true` | try edge collapses: with a prescribed target, the edges shorter than 1/√2 of it in every pass, and the edges of the elements above the desired density in a pass where no swap was accepted |
| `splits` | no | `true` | try edge splits: with a prescribed target, the edges longer than √2 of it in every pass, and the edges of the elements above the desired density in a pass where no swap was accepted; the new node starts at the midpoint, is returned to the surfaces of the boundary faces of the edge, inherits their side sets and the node sets common to the ends of a boundary edge, and is relaxed locally before the test |

```yaml
adaptivity:
  desired energy density: 0.05
  adjacency layers: 2
  maximum passes: 10
  outer iterations: 3
```

Examples: `examples/ems/awful-cube/awful-cube-adaptive.yaml` (a distorted
cube improved at its own mesh size) and `examples/ems/tube/tube-adaptive.yaml`
(a tube refined toward a finer target with its nodes on analytic surfaces).
Formulation, design, and measurements: `docs/notes/ems-adaptivity`.

Mesh smoothing is a specialized capability; most simulations omit these keys.
See `examples/ems/` for smoothing cases.

## Canonical examples

- Basic model/material block: `examples/single/static-solid`
- Mesh smoothing: `examples/ems/cube/cube.yaml`
- Anisotropic smoothing: `examples/ems/tube/tube-metric.yaml`
- Smoothing with topological operations: `examples/ems/awful-cube/awful-cube-adaptive.yaml`
