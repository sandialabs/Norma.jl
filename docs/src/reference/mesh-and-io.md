# Mesh and output

## Mesh files

Every single-domain input file names its input and output meshes, both in
Exodus II format:

| Key | Required | Meaning |
|---|---|---|
| `input mesh file` | yes | Exodus mesh read as input (opened read-only) |
| `output mesh file` | yes | Exodus results file written by the run |

```yaml
input mesh file: cantilever.g
output mesh file: cantilever.e
```

The input and output paths must differ. The output file is (re)created at the
start of the run as a copy of the input mesh and then populated with results.
Meshes carry the element blocks referenced by [`material: blocks`](materials.md)
and the node sets and side sets referenced by boundary and initial conditions.

## Output cadence

| Key | Required | Default | Meaning |
|---|---|---|---|
| `Exodus output interval` | no | the time step | simulation-time interval between Exodus writes; a non-positive value disables Exodus output |
| `CSV output interval` | no | `0.0` (disabled) | simulation-time interval between CSV writes of nodal quantities |
| `CSV write sidesets` | no | off | when present, also write per-side-set CSV files at each CSV write |

```yaml
Exodus output interval: 1.0e-6
CSV output interval: 1.0e-7
CSV write sidesets: true
```

Intervals are in simulation time, not step counts. In a
[multidomain](multidomain.md) run the output intervals are set on the top-level
controller and applied uniformly to all subdomains.

## Output fields

Exodus output always includes the reference coordinates (`refe_x/y/z`) and
displacements (`disp_x/y/z`). Dynamic time integrators add velocities
(`velo_x/y/z`) and accelerations (`acce_x/y/z`). Enabling
[`nodal recovery`](model.md) adds recovered nodal fields such as stress
components, von Mises stress, the deformation gradient, and internal variables.
Element output includes per-quadrature stresses and the element stored energy.
A [mesh smoothing](model.md) run adds the element variable `energy_density`,
the smoothing energy per unit ideal volume, and the target at every node:
`size` under a size field; under a metric field the principal sizes
`size_1..3`, the rotation vector `rotation_1..3` (when the target has one),
and the metric and target tensors `metric_xx`, `metric_yy`, `metric_zz`,
`metric_xy`, `metric_yz`, `metric_zx` and `target_xx` through `target_zx`.
A smoothing run with an `adaptivity` block writes each adapted mesh as
`<output name>-adapted-<k>.g` and the output of each smoothing phase after
the first as `<output mesh file>-s_0002`, `-s_0003`, and so on, which
ParaView opens as one sequence.
These field names are produced by the output writer and are not input keys.

## Matrix output

`matrix output` writes the assembled global operators of a monolithic
(`type: single`) solid mechanics run, so the discrete problem can be
reproduced outside Norma. It aborts in a Schwarz (`type: multi`) run and with
a `Periodic` boundary condition.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `matrices` | no | the operators the run assembles | operators to write: `stiffness`, `mass`, `external force` |
| `steps` | no | `[0]` | stops at which to write: `all`, a stop index, a list of stop indices, or `{every: n}` |

```yaml
matrix output:
  matrices: [stiffness, mass, external force]
  steps: [0]
```

Files, prefixed with the simulation name and indexed by the stop:

| File | Contents |
|---|---|
| `<name>-stiffness-NNNN.mm` | stiffness matrix K (tangent for nonlinear materials) plus the Robin stiffness α W, Matrix Market |
| `<name>-mass.mm` | consistent mass matrix M, written once |
| `<name>-lumped-mass.mm` | lumped mass as a diagonal matrix, written once (explicit or matrix-free dynamics) |
| `<name>-external-force-NNNN.csv` | body force plus boundary (Neumann) force, one node per line |
| `<name>-free-dofs-NNNN.csv` | one line per DOF: `true` free, `false` Dirichlet-constrained |

Matrices include every DOF, constrained or not. DOF `3(n-1)+c` is component
`c` (1 = x, 2 = y, 3 = z) of node `n` in Exodus order. The stiffness requires an
implicit integrator (quasi-static or Newmark) with the Hessian minimizer, and
the mass a dynamic integrator (Newmark or central difference). Without
`matrices`, the stiffness is written when the run assembles it, the mass when
the run is dynamic, and the external force always; an operator listed in
`matrices` that the run does not assemble aborts at setup. A quasi-static run
assembles the stiffness in its first solve, so at stop 0 the stiffness is not
written unless `initial equilibrium` is set; the log says so. With HHT-α the
Newton solve uses (1 − α) K; the file holds the unscaled K, and the external
solver applies α. Initial conditions and reference solutions come from the
regular CSV output. Example:
`examples/ahead/single/clamped/dynamic-linear-matrix-output/`.

## Canonical examples

- Basic Exodus output: any file under `examples/single/`
- CSV snapshot output: cases with `CSV output interval` under
  `examples/overlap/` and `examples/ahead/`
