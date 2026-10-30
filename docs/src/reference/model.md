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

## Composite tet 10

| Key | Required | Default | Meaning |
|---|---|---|---|
| `composite tet 10` | no | `false` | `true` for every TETRA10 block, or a list of block names (or a map of block names to booleans): those blocks use the composite 10-node tetrahedron |

The composite tetrahedron is LCM's element for `Use Composite Tet 10: true`,
as Intrepid2's `Basis_HGRAD_TET_COMP12_FEM` implements it. That is the base
formulation of J. T. Ostien, J. W. Foulk, A. Mota and M. G. Veilleux, "A
10-node composite tetrahedral finite element for solid mechanics", IJNME
107:1145–1170, 2016 (doi:10.1002/nme.5218). J. W. Foulk III, J. T. Ostien,
B. Talamini, M. R. Tupek, N. K. Crane, A. Mota and M. G. Veilleux extend it in
"Extending a 10-node composite tetrahedral finite element for solid
mechanics", IJNME 122:3845–3875, 2021 (doi:10.1002/nme.6684).
[Differences from the papers](#Differences-from-the-papers) lists where
Norma departs from them, mostly as LCM does. It runs on a TETRA10 mesh. `true` on a mesh without TETRA10 blocks switches nothing, with a
warning, and the values of a map of block names must be `true` or `false`. A
composite block and a standard TETRA10 block that share a face do not
conform: on a face the composite functions are piecewise linear on its four
sub-triangles and the standard ones quadratic, so only a linear displacement
is continuous across it. A connected mesh should use one of the two.

- **Shape functions.** The tetrahedron is split into 12 sub-tetrahedra:
  - the 4 corner tetrahedra of its regular refinement;
  - the 8 cones from its centroid over the faces of the inner octahedron.

  The shape functions are piecewise linear on them. At the centroid, those of
  the edge nodes are 1/6 and those of the corner nodes are 0.
- **Gradients.** The deformation gradient, the internal force and the tangent
  use the gradients' L2 projection onto linear polynomials over the parent
  tetrahedron, as LCM does: the closed form of the 2016 paper (§3.2.1),
  which assumes edge nodes at the midpoints. That is an assumed deformation
  gradient, linear over a straight-sided element.
- **Integration points.** The block uses the integration points of the
  parent tetrahedron. The default is the 4-point rule, which is the rule of
  the paper and LCM's `Cubature Degree: 2`. `num integration points` may
  choose 5 (LCM's default, degree 3, with a negative centroid weight) or 14.
- **Mass.** The mass is integrated exactly on the sub-tetrahedra, as in the
  2021 paper (eq. 50).
  - On a straight-sided tetrahedron, the consistent mass is LCM's analytic
    composite tet mass (`Use Analytic Mass`): V ρ / 1440 times a sparse
    matrix. On a curved one, the two differ.
  - The lumped mass is its row sums: on a straight-sided tetrahedron,
    V ρ (3, 14) / 96 at the corner and edge nodes, LCM's `Lump Analytic Mass`.
    It is positive, unlike the row sums of the quadratic TETRA10, so explicit
    dynamics can use it.
  - The explicit stable step uses the estimates of the other elements
    (see the central difference integrator): the element eigenvalue with
    the composite lumped mass, or the altitude of the four vertices scaled
    by the critical step of the regular composite tetrahedron.
  - For a single element, the critical step is 0.13 to 0.38 h/c over the
    shapes tried, against 0.16 to 0.45 for a TETRA4 of the same shape.
- **Interpolation.** Nodal recovery, field transfer, Schwarz coupling and the
  search for a point in the mesh all use the composite shape functions.
  - A transfer onto composite blocks integrates the nodal source with the
    exact mass rule. With consistent recovery, it is the L2 projection onto
    the composite space, so a field of that space comes back unchanged.
  - A point lies in a composite element only when it lies in one of its
    sub-tetrahedra. The standard tetrahedra also accept points beyond some of
    their faces and extrapolate to them; the composite tet does not, so a
    point outside a composite mesh, such as a Schwarz boundary node beyond its
    surface, is not found.
  - The overlap matching coupling rejects composite blocks: it weights the
    volume element at the integration points, which the exact mass rule does
    not use.
- **Faces.** Neumann and Robin conditions, and the non-overlapping Schwarz
  coupling, integrate over the faces with the quadratic TRI6 shape functions,
  as LCM does for its Neumann conditions. These are not the composite
  element's face functions, so the nodal loads are not consistent with the
  element: under a uniform traction, the stress in the elements at the loaded
  face is off by about half of the traction, at every mesh size. The total
  force is right on straight faces, and prescribed displacements are not
  affected.
- **Patch test.** An affine displacement is reproduced exactly on distorted
  straight-sided tetrahedra.
  - On curved tetrahedra (edge nodes off the chords), the projected
    gradients are not those of the geometry, so the patch test holds only in
    the limit. On a warped cube, the stress error halves with the element
    size. The operator of the 2021 paper passes it; Norma, like LCM, keeps
    that of the 2016 paper.
  - The standard TETRA10 passes it on the same mesh.
- **Accuracy.** On a thick-walled cylinder, the displacement converges at
  second order.
  - As ν → 1/2 on a fixed mesh, the element partially locks, as both papers
    report for this formulation: with 4 elements through the wall, its error
    grows from 0.008 of the inner displacement at ν = 0.3 to 0.26 at
    ν = 0.499999. Norma has neither paper's remedy.
  - The TETRA4's error grows without bound on the same mesh (9.7 at
    ν = 0.499999).
  - The error of the quadratic TETRA10 grows too, to 0.064, about 4 times
    less, on this displacement-driven problem.

`test/composite-tet10-lcm.jl` and `test/lcm-rubiks-cube.jl` compare it with
LCM on these cases:
- LCM's ClampedSDBC1DCt test and a Neohookean cylinder of curved elements,
  which agree to round-off;
- LCM's ProjectTet10 gold, which agrees to LCM's residual tolerance (1e-10),
  and its nodal projection of a linear field;
- the crystal-plasticity RubiksCube cut into TETRA10, whose stresses, slips
  and Fp agree to about 2e-8, as on the hexahedral RubiksCube.

LCM's implicit Newmark (Tempus d-form) starts from a zero acceleration,
whereas Norma solves for the initial acceleration. Implicit dynamic runs
therefore differ by that start, unless the initial state is in equilibrium.

```yaml
model:
  type: solid mechanics
  composite tet 10: [grain_1, grain_2]
  material:
    ...
```

### Differences from the papers

!!! warning "Loads on composite faces"
    Tractions, pressures, Robin conditions and the non-overlapping Schwarz
    coupling load composite blocks through TRI6 face functions, as LCM does.
    The stress in the elements at a loaded face is then wrong by a fraction of
    order one (about half of a uniform traction), and refining the mesh does
    not reduce it. Conditions on the displacement do not have this error.

Norma's composite tet 10 is LCM's. It departs from the papers where LCM
does, and in one place from LCM:

| | The papers | Norma | LCM |
|---|---|---|---|
| Face loads | The element's own face functions, piecewise linear on four sub-triangles (2016, eq. 17) | TRI6 face functions (see above) | As Norma |
| Curved elements | The projection of the true gradients, weighted by the sub-tetrahedra Jacobians, and integration weights from a composite Jacobian (2021, §3.2 and §3.4), which pass the patch test on curved elements | The closed form of 2016 (§3.2.1), which assumes edge nodes at the midpoints: on a curved element, F̄ is not linear and the patch test holds only in the limit | As Norma |
| Near incompressibility | The volume-averaged F̄* and p̄* (2016, §3.4); a five-field functional with a stabilization (2021, §3.6–3.7) | Neither: the element partially locks as ν → 1/2 | F̄* and p̄* as the block options `Weighted Volume Average J` and `Volume Average Pressure` |
| Stress with 5 or 14 points | The projected stress P̄ (2016, eqs. 37–40, for 4 and 5 points) | The stress of F̄ at each point. The nodal forces are the same, and with 4 points so is the stress | As Norma |
| Mass on curved elements | Exact on the sub-tetrahedra (2021, eq. 50) | As the 2021 paper | The straight-sided matrix scaled by the element volume, a few percent off on a curved element |

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

The keys `stress recovery` and `recover internal variables` of earlier versions are no longer
accepted and cause an error; use the `nodal recovery` block above.

## Mesh smoothing

| Key | Required | Default | Meaning |
|---|---|---|---|
| `mesh smoothing` | no | `false` | enable mesh smoothing (set automatically when `type: mesh smoothing`); a top-level key |
| `smooth reference` | no | `""` | rule for the ideal element of TETRA4 smoothing: `size field`, `size field unrestricted`, `metric field`, or `metric field unrestricted`, which prescribe the target; `equal volume`, `average edge length`, and `max`, which take the target from the original element, are rules of earlier versions kept for older inputs |
| `size field` | required if `smooth reference: size field`/`size field unrestricted` | `nothing` | expression in `t, x, y, z` giving the target element size (edge length of the ideal regular tetrahedron); under `size field` the target is at least the edge of the regular tetrahedron with the volume of the original element, so the ideal never asks an element to shrink, and under `size field unrestricted` it is the field as given |
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
the objective, on which an isotropic energy on `F` alone does not depend. With equal
sizes the rule reduces exactly to `size field`. As for `size field`, the
target is evaluated at the centroid of each element in the original mesh, so
it does not follow the nodes as they move, and at the current time, so an
expression in `t` changes the target from one step to the next. `metric field` scales the three sizes uniformly so the
ideal volume is never below the volume of the original element (the
anisotropic form of the `size field` floor); `metric field unrestricted` does
not. The output carries the metric at every node in three forms: the
principal sizes and the rotation vector as `size_1..3` and `rotation_1..3`
(the rotation omitted when the target is given by `sizes` without a
`rotation vector`);
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
meshed at 0.1, an axis-aligned target of sizes (0.025, 0.2, 0.2) leaves the
mesh isotropic: the solve reaches a strict local minimum (gradient below
10⁻⁸, Hessian positive definite), the same one from a randomly displaced
initial mesh, because the faces are fixed planes, so the extent across x is
1 and the elements along a path across it are as many as the connectivity
has. The same target rotated by 45 degrees is partly reachable, since a
path can zigzag in the plane normal to the thin direction: the mean edge
component along it goes from 0.052 to 0.029 against the 0.025 asked, with
the energy five times lower. Ramping the sizes over the pseudo-time reaches
the same minimum. The remedy is the topology of the `adaptivity` block: on
that case it takes the energy from 7.3 × 10⁴ to 9.4 × 10² over twelve
iterations, at 10801 elements from 7677, and the mean edge components from
(0.052, 0.052, 0.052) to (0.0134, 0.091, 0.091) against the
(0.0125, 0.1, 0.1) that the target implies, with 96.8 percent of the edges
inside the band in the metric and every boundary node on its face
(`examples/ems/cube/cube-metric-adaptive.yaml`). The two orientations then
agree, as they must: the same target at 45 degrees
(`examples/ems/cube/cube-metric-rotated-adaptive.yaml`) ends at 10768
elements, edge components (0.0134, 0.094, 0.090), and 96.2 percent inside
the band.

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
of a nodal run, and `metric_*` and `target_*`, so either tensor form can.

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
loop of `docs/notes/ems-adaptivity`: the mesh is smoothed, then a topology
phase swaps edges (interior ones, and with `boundary swaps` boundary ones),
collapses edges (a node is removed only onto a node that carries its node
sets and lies on its surfaces, along a boundary edge if it is on the
boundary), and splits edges, the new mesh is written as
`<output name>-adapted-<k>.g` and smoothed again, and so on until a topology
phase accepts no operation or the outer iterations are exhausted. The
output of the first smoothing phase is the `output mesh file`, and that of
each later phase appends `-s_0002`, `-s_0003`, and so on
(`cube.e`, `cube.e-s_0002`, ...), which ParaView opens as one sequence; the
times written in each file continue one time step after the last time of
the file before it, since ParaView orders the sequence by time and shows
only one file at a repeated time. The mesh
must consist of four-node tetrahedra. Node sets and side sets are carried
over to the written meshes, and a split node on the boundary joins the side
sets of the faces it lies on.

Each operation is accepted by the test its criterion names. With the
defaults every operation must lower the energy of the elements around it
(its cavity) by the relative margin; `size criterion: length` and
`shape criterion: scaled Jacobian`, described below, replace that test for
the size operations and the shape operations. Under every criterion, no
new element may exceed the allowed density, and the worst new element must
respect the geometric floor. The energy sums the elements of the cavity, so
without the floor an operation can improve the sum while creating one flat
element; the floor is stated in the scaled Jacobian that the analysis codes
require. The combination recommended for refinement, coarsening, and
repair alike is the length criterion, the scaled Jacobian criterion,
boundary swaps, and the size operations first, as in the sample below and
in every adaptive example except `tube-adaptive.yaml`, which keeps the
energy test for the swaps.

The energy test alone cannot refine or coarsen a smoothed mesh toward a
prescribed target by more than a factor of about 1.5 in edge length: a split
bisects the elements around the edge and approximately halves their scaled
Jacobian, and with equal bulk and shear moduli the increase of the shape
energy exceeds the decrease of the size energy until the edges are far
longer than the target (on a plate
refined toward a sinusoid, every split of the 920 edges beyond the band
raised the cavity energy, by 52 percent at the median). A larger bulk
modulus makes the splits favorable but also makes the smoothing lower the
shape quality to reach the size. `size criterion: length` therefore accepts the size
operations, the splits of edges longer than √2 and the collapses of edges
shorter than 1/√2 in the target, on the edge length alone, subject to the
geometric floor and, for collapses, the set and surface constraints, with
the energy as a validity check only; the swaps and the shape-driven
operations keep the energy test. Under refinement the floor limits the
number of splits (a floor of 0.35 accepted 981 of 29797 splits on that
plate, none 29222), so a low floor such as 0.15 suits a refinement stage.

The cavity energy sum can also accept a swap that leaves the worst element
of the cavity as it was, and refuse one that raises it while raising the
sum. `shape criterion: scaled Jacobian` accepts the swaps and the
shape-driven collapses and splits when the minimum scaled Jacobian of the
cavity rises, and chooses among the configurations of a swap the one with
the largest minimum scaled Jacobian; the energy remains the validity check. On the
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
workflows from the distorted input is the sequence of operations each
takes, not the acceptance test. Ranking the cavities worst first, more alternation with
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
accepted operation marks the elements of its cavity as removed and a later
proposal that includes a removed element is refused, but they may not create an
edge or a face that exists in the mesh or was created earlier in the pass,
so that no face is shared by more than two elements and the link of every
edge stays a single ring. The topology is compacted after every operator,
so that the collapses and splits operate on the elements the swaps made. Nodes do
not move within a phase, so a refused proposal is remembered for the phase
and retried only when an element around its edge or face has changed; on a
converged mesh of half a million elements this takes a pass from about
200 seconds to 7.

With `cavity relaxation: true` a swap or collapse that the energy test
refuses at the fixed positions is judged again after a relaxation of the
nodes of its cavity that may move, under the constraints of the smoothing:
a node on side sets with Surface conditions slides on their surfaces (along
their line for two), a node in a node set with Dirichlet conditions keeps
the components they fix, and a corner, a node with all components fixed, a
node set with a condition of another kind, and a node on a side set
without a Surface condition stay:
each node in turn takes up to five descent steps on the energy of the
elements around it under the new connectivity, with the gradient of the
energy with respect to the node (the element internal force of the
smoothing material, one evaluation of the elements around the node) and the
step halved until the energy decreases, as the node of a split does. When
after the first step over the nodes the energy of the new and surrounding
elements is still above 1.5 times that of the old and surrounding ones, the
relaxation stops and the operation is refused: on the awful cube 98% of the
operations the full relaxation accepts are below that ratio after the first
step, and 71% of the relaxations tried are above it. The energy test and
the geometric floor then compare the old elements with the new ones at the
relaxed positions, both together with the elements around the moved nodes
that the operation keeps, since those change with the positions; an accepted operation moves the nodes, the
stored values of the elements around them are evaluated again, and the
refusals of their edges and faces are forgotten. The relaxation applies to
the energy test only, not to the tests on the length or on the scaled
Jacobian. A swap that is accepted at the fixed positions is not relaxed,
since the smoothing phase that follows moves the nodes anyway. Measured on
the awful cube example under the energy criteria: the final energy 6.24
against 9.51, the minimum scaled Jacobian 0.42 against 0.36, and the mean
0.73 against 0.68 at the same element count (the shipped configuration of
the example reaches 0.40 and 0.70 with half again as many elements), with
the topology phases taking 40 s against 2.4 s, which is why the option is
off by default.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `size criterion` | no | `energy` | `energy`: every operation must lower the cavity energy; `length`: the splits and collapses of edges outside the length band of the prescribed target are accepted on the length alone, subject to the floor and the constraints |
| `shape criterion` | no | `energy` | `energy`: a swap or a shape-driven collapse or split must lower the cavity energy; `scaled Jacobian`: it must raise the minimum scaled Jacobian of the cavity, measured in the metric space under a metric field |
| `boundary swaps` | no | `false` | try the swap of boundary edges whose two boundary faces lie in one side set (or in none) and form a flat patch |
| `boundary swap angle` | no | `20` | largest angle in degrees between the two boundary faces of an edge for its swap to be tried |
| `size operations first` | no | `false` | collapse and split the edges outside the band before the swaps in every pass, rather than after; recommended, since the swaps then act on elements already near their target size |
| `shape operations every pass` | no | `false` | try the shape-driven collapses and splits in every pass, not only in a pass that accepted no swap |
| `cavity relaxation` | no | `false` | judge a swap or collapse that the energy test refuses at the fixed positions again after relaxing the interior nodes of its cavity against the new connectivity, with the elements around those nodes included in the energy and the floor |
| `desired scaled Jacobian` | no | `0.9` | under `shape criterion: scaled Jacobian`, the elements below this are the candidates of the swaps and the shape-driven operations, without dilation; lowering it restricts the phase to the elements of lowest quality and shortens it |
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
  size criterion: length
  shape criterion: scaled Jacobian
  boundary swaps: true
  size operations first: true
  minimum scaled Jacobian: 0.15
  desired energy density: 0.05
  adjacency layers: 2
  maximum passes: 10
  outer iterations: 3
```

Examples: `examples/ems/awful-cube/awful-cube-adaptive.yaml` (a distorted
cube improved at its own mesh size), `examples/ems/tube/tube-adaptive.yaml`
(a tube refined toward a finer target with its nodes on analytic surfaces),
`examples/ems/plate/plate-sinusoid.yaml` (a plate refined in stages toward a
sinusoidal size field), `examples/ems/tube/tube-metric-adaptive.yaml`
(refinement and coarsening toward a graded metric), and
`examples/ems/cube/cube-metric-adaptive.yaml` with
`cube-metric-rotated-adaptive.yaml` (an anisotropic target in two
orientations).
Formulation, design, and measurements: `docs/notes/ems-adaptivity`.

Mesh smoothing is a specialized capability; most simulations omit these keys.
See `examples/ems/` for smoothing cases.

## Canonical examples

- Basic model/material block: `examples/single/static-solid`
- Composite tet 10: `examples/element-types/composite-tet10/cube/cube.yaml`
- Mesh smoothing: `examples/ems/cube/cube.yaml`
- Anisotropic smoothing: `examples/ems/tube/tube-metric.yaml`
- Metric given as a tensor: `examples/ems/cube/cube-tensor.yaml`
- Metric carried by the nodes of the input mesh: `examples/ems/cube/cube-metric-nodal.yaml`
- Smoothing with topological operations: `examples/ems/awful-cube/awful-cube-adaptive.yaml`
- Swaps only, in stages with a size field that follows the element count: `examples/ems/awful-cube/awful-cube-swaps.yaml` with `size-adaptive.jl`
- Refinement toward a prescribed size field in stages: `examples/ems/plate/plate-sinusoid.yaml` with `refine.jl`
- Refinement and coarsening toward a metric: `examples/ems/tube/tube-metric-adaptive.yaml`
- A metric that smoothing alone cannot reach: `examples/ems/cube/cube-metric-adaptive.yaml`, and the same target turned by 45 degrees: `examples/ems/cube/cube-metric-rotated-adaptive.yaml`
