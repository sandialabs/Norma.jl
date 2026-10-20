# Examples

The `examples/` directory holds hundreds of runnable input files spanning every
capability. To run one, activate the project and call `Norma.run` on its YAML
file:

```bash
cd /path/to/Norma.jl/examples/contact/implicit-dynamic/2-bars
julia
]
activate .
using Norma
Norma.run("bars.yaml")
```

## How the examples are organized

| Directory | Contents |
|---|---|
| `examples/single/` | Single-domain runs: `static-solid`, `implicit-dynamic-solid`, `explicit-dynamic-solid` |
| `examples/materials/` | Material-model cases: J2 plasticity, anisotropic hyperelastic damage |
| `examples/element-types/` | Element-library cases (`tet4`, `tet10`, …) |
| `examples/crystal-plasticity/` | Polycrystals with a lattice orientation per element: LCM's RubiksCube and MicroStructPy Voronoi polycrystals, with their mesh generators, on rollers |
| `examples/overlap/` | Overlapping Schwarz coupling (same and different time steps) |
| `examples/nonoverlap/` | Non-overlapping Schwarz coupling: Dirichlet–Neumann (including the matching exchange), Robin–Robin, and matrix inclusion |
| `examples/contact/` | Contact via Schwarz (static, implicit, explicit) |
| `examples/adaptive-time-stepping/` | Adaptive time stepping and mesh swapping |
| `examples/ems/` | Mesh smoothing toward size and metric fields, and smoothing alternated with topological operations |
| `examples/ahead/` | Application-motivated cases, including the reduced-order-model workflow |

Each example is a self-contained set of a YAML input file and its Exodus mesh.
The input, read with the [Input File Reference](reference/index.md), shows how
a feature is configured.

## Choosing a starting point

- New to Norma: start from `examples/single/static-solid` and
  `examples/single/implicit-dynamic-solid`.
- Explicit dynamics: `examples/single/explicit-dynamic-solid`.
- Schwarz coupling: `examples/overlap/` and `examples/nonoverlap/`; see also the
  theory note in `docs/notes/schwarz-coupling`.
