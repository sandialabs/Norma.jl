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
These field names are produced by the output writer and are not input keys.

## Canonical examples

- Basic Exodus output: any file under `examples/single/`
- CSV snapshot output: cases with `CSV output interval` under
  `examples/overlap/` and `examples/ahead/`
