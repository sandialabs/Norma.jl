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
(`Schwarz overlap`, `Schwarz DN nonoverlap`, `Schwarz RR nonoverlap`,
`Schwarz contact`) are documented under
[Multidomain and Schwarz coupling](multidomain.md).

Conditions are applied in a fixed order that does not depend on the order of
the type keys in the file: all other types first, in alphabetical order of the
type name, and the prescribed types (`Dirichlet`, `OpInf Dirichlet`) last.
Entries of one type keep their order in the list. A node that belongs to more
than one condition, such as an interface edge node that also lies on a
symmetry or clamped face, therefore ends up with its prescribed value.

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
constant `function` (no gradient) is rejected. A node in the side sets of two
`Surface` conditions slides along their intersection, and a node in three is
fixed at their common point; surfaces that meet tangentially at a node are
rejected. Under `exact` enforcement the descent direction is projected onto
the tangent space and every node is returned to its surfaces by a
closest-point projection, which requires the `steepest descent` solver;
`penalty` enforcement adds the force of the penalty ½ κ g² and works with
every solver. This is the condition that keeps the boundary nodes of
[mesh smoothing](model.md) on curved surfaces, and the adaptivity loop places
the nodes that its splits add on the boundary on the same surfaces.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `side set` | yes | — | surface whose nodes are constrained |
| `function` | yes | — | level-set expression g in `t, x, y, z` |
| `enforcement` | no | `exact` | `exact` or `penalty` |
| `penalty` | no | `1.0e6` | penalty coefficient (positive; used when `enforcement: penalty`) |

## `Periodic`

Periodic conditions on a box cell, the boundary conditions of a
representative volume element. Each pair of node sets lists two opposite
faces of the box, the second the image of the first under a translation
along the positive x, y or z axis. For every node X⁻ of the first set and its
image X⁺ on the second,

u(X⁺) − u(X⁻) = (F̄ − I)(X⁺ − X⁻),

where F̄ is the macroscopic deformation gradient, which for a periodic field
is the volume average of F over the cell. Every component F̄ᵢⱼ whose column j
is the axis of a pair is controlled in one of two ways: prescribed as a
function of t under `deformation gradient`, or free under `stress`, with the
same component of P̄ prescribed instead. P̄ is the volume average of the
first Piola–Kirchhoff stress, the work conjugate of F̄. A component listed in
neither map is held at its value in the identity.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `node set pairs` | yes | — | one to three pairs `[first, second]` of node sets on opposite faces |
| `deformation gradient` | no | identity | prescribed components of F̄, as functions of `t` |
| `stress` | no | none | free components of F̄, with the prescribed components of P̄ as functions of `t` |

Components are named by row and column: `xy` is F̄ₓᵧ = ∂x/∂Y, and the names
are `xx`, `xy`, `xz`, `yx`, `yy`, `yz`, `zx`, `zy` and `zz`. A component may
appear in only one of the two maps. The rigid rotation has to stay fixed, so
of `xy` and `yx` (and of `xz` and `zx`, and of `yz` and `zy`) at most one may
be under `stress`. With one or two pairs, a plane that contains an axis that
is not periodic has one component of F̄, or none: a cell periodic along x and
y has `zx` but not `xz`, and a cell periodic along x only has no component in
the y–z plane. The rotation in such a plane is held when its component is
prescribed; otherwise another condition has to hold it, such as a Dirichlet
condition on a node away from the origin, and the run aborts if none does. A
cell periodic along x only, with free lateral faces, needs one for the
rotation about x, for example on component z of the corner node set
`nsx-y+z-`, which a homogeneous deformation leaves at zero.

Uniaxial stress along x prescribes F̄ₓₓ and leaves every other component of P̄
at zero. The lower shears `yx`, `zx` and `zy` stay at zero and fix the
rotation. The free upper shears let an anisotropic cell shear under the
axial load, and the average stress is then uniaxial:

```yaml
Periodic:
  - node set pairs: [[nsx-, nsx+], [nsy-, nsy+], [nsz-, nsz+]]
    deformation gradient:
      xx: "1.0 + 1.0e-3 * t"
    stress:
      yy: "0.0"
      zz: "0.0"
      xy: "0.0"
      xz: "0.0"
      yz: "0.0"
```

The condition eliminates the displacements of the nodes of the second sets:

- The nodes of the two sets of a pair are matched by position, within 10⁻⁶
  of the diameter of the mesh. A node without an image aborts the run, so the
  mesh must conform across the periodic faces.
- The node common to the first sets of all pairs, the lower corner of the
  cell, is the origin O and is held fixed, which removes the rigid
  translation. With one or two pairs these nodes form a face or an edge, and
  the origin is the lowest of them by x, then y, then z, with coordinates
  that differ by less than the matching tolerance taken as equal. It keeps
  zero displacement unless a Dirichlet condition prescribes it. Its image C
  across each pair, at a distance L along the axis of the pair, carries the
  column of F̄ along that axis:
  u(C) − u(O) = (F̄ − I) L eₖ. A prescribed component of F̄ prescribes the
  displacement of C. A stress-controlled component is a force V P̄ᵢⱼ / L on C,
  with V the volume of the cell.
- Every other node of the second sets follows its image on the first sets.
  A node at an edge or a corner of the cell follows the image that is on the
  first sets of all its pairs. The Newton residual and tangent are condensed
  onto the remaining unknowns.
- At every step the increment of the prescribed components of F̄ is applied
  to the whole cell, with the increment Δu_O of the displacement of the
  origin when a Dirichlet condition moves it: u ← u + Δu_O + ΔF̄ (X − X_O).
  The first Newton iterate is the previous solution plus a translation and a
  homogeneous deformation, which is exact for a homogeneous material, so the
  increment does not fall on the elements next to the faces as it would
  with Dirichlet conditions. Without `initial equilibrium`, the cell starts
  in the same way from the affine field of F̄ at the initial time. The
  increment is the prediction of the step, so the relative tolerance of the
  solve is measured from the residual without it, as for a
  [predictor](time-integrators.md#Predictor) of the time integrator. A step
  that a `linear` or `tangent` predictor predicts takes the predicted field
  instead: on the condensed system, the predictor moves the independent nodes
  with the corners. A step it does not predict, such as the first step of
  `linear`, takes the increment of F̄.

F̄ and P̄ are written to the Exodus output as global variables `F_xx`, …,
`F_zz` and `P_xx`, …, `P_zz`. In Julia they are
`Norma.periodic_deformation_gradient(model)` and
`Norma.periodic_stress(model)`. The columns of F̄ along the periodic axes are
those the corners carry. With one or two pairs, the other columns, which the
condition does not control, are the volume average of F computed over the
elements. P̄ is computed as (1/V) Σₐ fₐ ⊗ Xₐ, summed over the internal forces
of the material at all nodes, without the forces of `Robin` and penalty
`Surface` conditions. Because the shape functions interpolate X exactly, this
sum is (1/V) ∫ P dV of the finite element field.

The condition requires the `quasi static` time integrator with the
`Hessian minimizer` solver, in a single-domain simulation, and does not
support `matrix output` (the solve uses the stiffness condensed onto the
independent degrees of freedom). A model has at
most one `Periodic` entry, with up to three pairs. Dirichlet conditions may
hold the origin and the nodes of the first sets that are on no second set,
but not the nodes of the second sets or the corners C, which the periodic
condition determines. With two or three pairs, a first set also holds its
edges on the second sets of the other pairs (`nsx-` of a cell periodic along
x, y and z holds nodes of `nsy+` and `nsz+`), so a Dirichlet condition on a
whole first set aborts, and the message lists the nodes it would prescribe.
A node set of the origin, such as `nsx-y-z-`, or of nodes inside a face, can
be held. With one or two pairs, the faces normal to the other axes take
other conditions, and only the columns of F̄ along the periodic axes are
controlled.

The MicroStructPy converter
`examples/crystal-plasticity/polycrystal/microstructpy_to_exodus.jl` checks
that a periodic mesh conforms across its periodic faces. It writes the face
node sets `nsx-`, …, `nsz+` (for a periodic axis, the i-th node of `nsx+` is
the image of the i-th node of `nsx-`) and a node set for each corner of the
box, `nsx-y-z-` to `nsx+y+z+`.

## Canonical examples

- Fixed and time-dependent Dirichlet: `examples/single/static-solid`,
  `examples/single/implicit-dynamic-solid`
- Neumann traction: single-domain examples under `examples/single/`
- Neumann pressure: pressurized cases under `examples/ahead/`
- Sliding on analytic surfaces: `examples/ems/cube/cube-surface.yaml` (planes)
  and `examples/ems/tube/tube.yaml` (cylinders and caps)
- Periodic representative volume element:
  `examples/crystal-plasticity/polycrystal/polycrystal-periodic.yaml`
