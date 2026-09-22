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
| `metric field` | required if `smooth reference: metric field`/`metric field unrestricted` | `nothing` | block giving the metric in one of four forms: `sizes` (three expressions in `t, x, y, z`: the principal sizes h₁, h₂, h₃) with optional `rotation vector` (three expressions: the rotation vector whose exponential carries the global axes onto the principal directions); `tensor` (six expressions: M_xx, M_yy, M_zz, M_xy, M_yz, M_zx); `nodal sizes` (three nodal variables of the input mesh) with optional `nodal rotation vector` (three nodal variables); or `nodal tensor` (six nodal variables), read at `time index` (default: the last step) and interpolated by `interpolation: principal` (default) or `log-Euclidean`; see below |

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
not. The output carries the metric at every node in three forms: the
principal sizes and the rotation vector as `size_1..3` and `rotation_1..3`;
the metric tensor `M = R diag(1/h²) Rᵀ` as `metric_xx`, `metric_yy`,
`metric_zz`, `metric_xy`, `metric_yz`, `metric_zx`, whose quadratic form
gives the squared length of a vector in the target; and the target tensor
`T = R diag(h) Rᵀ` under the same six suffixes with the prefix `target`,
whose eigenvectors are the principal directions and whose eigenvalues are
the target edge lengths, which is the form a tool that reads scaled
eigenvector tensors expects. Formulation and tests:
`docs/notes/ems-anisotropic`.

A metric is a target, not a promise: with the boundary held on its
surfaces and the connectivity fixed, smoothing cannot change the number of
elements across the domain, so a target whose sizes ask for more elements
across a direction than the mesh has is unreachable. On the unit cube
meshed at 0.1, an axis-aligned target with the first size 0.025 leaves the
mesh isotropic: the solve reaches a strict local minimum (gradient below
10⁻⁸, Hessian positive definite), the same one from a randomly displaced
initial mesh, because the faces are fixed planes, so the extent across x is
1 and the elements along a path across it are as many as the connectivity
has. The same target rotated by 45 degrees is partly reachable, since a
path can zigzag in the plane normal to the thin direction: the mean edge
component along it goes from 0.052 to 0.029 against the 0.025 asked, with
the energy five times lower. Ramping the sizes over the pseudo-time reaches
the same minimum. The remedy is the topology of the `adaptivity` block: on
that case it takes the energy from 7.3 × 10⁴ to 8.7 × 10³ in two
iterations, at about 11000 elements from 7677, and the mean edge component
along the thin direction from 0.052 to 0.013 against the 0.0125 that the
target implies, with 97 percent of the edges inside the band in the metric
(`examples/ems/cube/cube-metric-adaptive.yaml`). The two orientations then
agree, as they must: with topology allowed, the axis-aligned target and the
same target at 45 degrees end at the same edge components and the same
fraction inside the band.

```yaml
model:
  type: mesh smoothing
  smooth reference: metric field
  metric field:
    sizes: ["0.025", "0.1", "0.1"]
    rotation vector: ["0.0", "0.0", "atan(y, x)"]   # first principal direction radial
```

The metric can also be given as a tensor, by its six components
`M_xx, M_yy, M_zz, M_xy, M_yz, M_zx` as expressions in `t, x, y, z`. The
tensor is factored into sizes and a rotation by eigendecomposition; the energy
does not depend on the frame chosen, so no convention is needed in the input.

```yaml
  metric field:
    tensor: ["225.0", "225.0", "50.0", "175.0", "0.0", "0.0"]   # sizes 0.05, 0.1414, 0.1414 at 45 degrees
```

A metric computed elsewhere, by an error estimator or a previous run, is
carried by the nodes of the input mesh as nodal variables, either as sizes and
an optional rotation vector or as the six tensor components, read at
`time index` (default: the last step). The value on an element is the mean of
its nodes: the sizes geometrically and the rotation through its vector, which
keeps the sizes exact where the frame turns. For nodal tensors the sizes and
rotation vectors are recovered at the nodes with a frame that follows from
node to node through the mesh (`interpolation: principal`, the default); where
the frame turns by more than 45 degrees between adjacent nodes a warning
counts the turns, and `interpolation: log-Euclidean` interpolates the
logarithm of the tensor instead, which needs no frame. Under `adaptivity`
the nodal metric follows the operations (a split node takes the mean of the
ends of its edge) and the adapted meshes carry it under the same variable
names, so the input applies to them unchanged. The output of a metric run
carries `size_1..3` and `rotation_1..3`, so it can serve as the input mesh
of a nodal run.

```yaml
input mesh file: cube-metric.e
model:
  type: mesh smoothing
  smooth reference: metric field unrestricted
  metric field:
    nodal sizes: [size_1, size_2, size_3]
    nodal rotation vector: [rotation_1, rotation_2, rotation_3]
    time index: 1
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

The energy test alone cannot refine or coarsen a smoothed mesh toward a
prescribed target by more than a factor of about 1.5 in edge length: a split
bisects the elements around the edge and roughly halves their scaled
Jacobian, and with equal bulk and shear moduli that shape penalty outweighs
the size gain until the edges are far longer than the target (on a plate
refined toward a sinusoid, every split of the 920 edges beyond the band
raised the cavity energy, by 52 percent at the median). A larger bulk
modulus makes the splits favorable but also makes the smoothing sacrifice
shape for size. `size criterion: length` therefore accepts the size
operations, the splits of edges longer than √2 and the collapses of edges
shorter than 1/√2 in the target, on the edge length alone, subject to the
geometric floor and, for collapses, the set and surface constraints, with
the energy as a validity check only; the swaps and the shape-driven
operations keep the energy test. Under refinement the floor throttles the
splits (a floor of 0.35 let 981 of 29797 splits through on that plate,
none 29222), so a low floor such as 0.15 suits a refinement stage.

The cavity energy sum can also accept a swap that leaves the worst element
of the cavity as it was, and refuse one that raises it while raising the
sum. `shape criterion: scaled Jacobian` accepts the swaps and the
shape-driven collapses and splits when the minimum scaled Jacobian of the
cavity rises, and chooses among the configurations of a swap the one with
the best worst element; the energy remains the validity check. On the
distorted cube this criterion together with `boundary swaps` brings the
coupled loop to the quality of the staggered workflow of Norma smoothing
alternated with the Sierra tool improve_mesh, which accepts on the same
measure (minimum 0.41 against 0.38, mean 0.692 against 0.703 with the same
Surface conditions, target, and smoothing); with the energy criterion alone
the minimum was 0.36 and the mean 0.668, the deficit being on the boundary,
where bisected boundary triangles were never repaired since only interior
edges were swapped. `boundary swaps: true` swaps a boundary edge whose two
boundary faces belong to one side set (or to none) and make an angle below
`boundary swap angle`: the edge is replaced by the edge between the far
nodes of the two faces, the faces by the two that contain the new edge, in
the side set as well, and the chain of elements behind them by a
triangulation of the resulting polygon. `face swaps: true` adds the swap of
an interior face shared by two elements into the three elements around the
edge between their apexes; measured on the distorted cube it lowers the
quality under both criteria and is off by default. The order of the
operators within a pass matters as well: with `size operations first: true`
the edges outside the band are split and collapsed before the swaps, as
improve_mesh does, and on the distorted cube this raises the final mean from
0.692 to 0.698 and the minimum from 0.39 to 0.41; the loop applied to the
mesh that the staggered workflow produced improves it further (minimum 0.38
to 0.44, mean 0.703 to 0.707), so the remaining difference between the two
workflows from the distorted input is the path each takes, not the
acceptance test. Ranking the cavities worst first, more alternation with
fewer passes, and shape operations in every pass were measured and gain
nothing.

Under a metric field every shape measure of the loop, the criterion, the
floor, the candidates, and the ranking, is the scaled Jacobian of the
element mapped by the metric factor, so that an element matching an
anisotropic target is regular and measures one. Measured in physical
space instead, the criterion penalizes the elongation the target asks for,
the swaps undo what the smoother does, and the loop does not converge (on
the tube under its graded metric the operations per iteration stayed at
220 to 440 and the energy tripled; in the metric space they fall 167, 72,
49, 37, 20, 13 and the mesh reaches 99.6 percent of its edges inside the
band). The splits of edges beyond the band are tried longest first in the
target, as in longest-edge bisection, which bounds the shape of the
children; on the plate this raised the mean scaled Jacobian of the final
stage from 0.696 to 0.702 and the minimum of the third from 0.22 to 0.33,
within 0.011 of the staggered workflow.

Within a pass the operations of one operator are independent, since every
accepted operation marks the elements of its cavity dead and a later
proposal touching a dead element is refused, but they may not create an
edge or a face that exists in the mesh or was created earlier in the pass,
so that no face is shared by more than two elements and the link of every
edge stays a single ring. The topology is compacted after every operator,
so that the collapses and splits see the elements the swaps made. Nodes do
not move within a phase, so a refused proposal is remembered for the phase
and retried only when an element around its edge or face has changed; on a
converged mesh of half a million elements this takes a pass from about
200 seconds to 7.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `size criterion` | no | `energy` | `energy`: every operation must lower the cavity energy; `length`: the splits and collapses of edges outside the length band of the prescribed target are accepted on the length alone, subject to the floor and the constraints |
| `shape criterion` | no | `energy` | `energy`: a swap or a shape-driven collapse or split must lower the cavity energy; `scaled Jacobian`: it must raise the minimum scaled Jacobian of the cavity, measured in the metric space under a metric field |
| `boundary swaps` | no | `false` | try the swap of boundary edges whose two boundary faces lie in one side set (or in none) and form a flat patch |
| `boundary swap angle` | no | `20` | largest angle in degrees between the two boundary faces of an edge for its swap to be tried |
| `size operations first` | no | `false` | collapse and split the edges outside the band before the swaps in every pass, rather than after; recommended, since the swaps then act on elements already near their target size |
| `shape operations every pass` | no | `false` | try the shape-driven collapses and splits in every pass, not only in a pass that accepted no swap |
| `desired scaled Jacobian` | no | `0.9` | under `shape criterion: scaled Jacobian`, the elements below this are the candidates of the swaps and the shape-driven operations, without dilation; lowering it focuses the phase on the worst elements and shortens it |
| `face swaps` | no | `false` | try the swap of an interior face into the three elements around the edge between the apexes |
| `desired energy density` | no | `0.1` | elements above this energy per unit ideal volume are candidates |
| `allowed energy density` | no | `Inf` | no accepted operation may create an element above this |
| `minimum scaled Jacobian` | no | `0` | geometric floor: no accepted operation may create an element with a scaled Jacobian below this, unless the worst element of the cavity was already below it and the new worst is no worse |
| `minimum decrease` | no | `1.0e-8` | relative decrease of the cavity energy an operation must achieve |
| `adjacency layers` | no | `4` | rings of adjacent elements added to the candidate set |
| `maximum passes` | no | `20` | passes of operations per topology phase |
| `outer iterations` | no | `5` | alternations of smoothing and topology |
| `swaps` | no | `true` | try the swaps of interior edges |
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
- Metric given as a tensor: `examples/ems/cube/cube-tensor.yaml`
- Metric carried by the nodes of the input mesh: `examples/ems/cube/cube-metric-nodal.yaml`
- Smoothing with topological operations: `examples/ems/awful-cube/awful-cube-adaptive.yaml`
- Refinement toward a prescribed size field in stages: `examples/ems/plate/plate-sinusoid.yaml` with `refine.jl`
- Refinement and coarsening toward a metric: `examples/ems/tube/tube-metric-adaptive.yaml`
- A metric that smoothing alone cannot reach: `examples/ems/cube/cube-metric-adaptive.yaml`
