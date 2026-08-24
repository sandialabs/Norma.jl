# Input file reference

Norma is driven by a single YAML input file, passed on the command line:

```bash
bin/norma input.yaml
```

This reference documents every section and key of that file. It covers the
full-order solid-mechanics capabilities; the reduced-order (Operator Inference)
model types are documented with the `norma-opinf` companion package and are out
of scope here.

## Two kinds of input file

Every input file declares its `type`:

- **`type: single`** — a single-domain simulation. Contains the mesh, model,
  time integrator, solver, and conditions. This is the common case and the basis
  for everything else.
- **`type: multi`** — a multidomain (Schwarz-coupled) simulation. A short
  controller file that lists subdomain files and governs the coupling. Each
  subdomain is itself a `type: single` file. See
  [Multidomain and Schwarz coupling](multidomain.md).

## Structure of a single-domain file

A single-domain input has the following top-level sections:

| Section | Required | Purpose |
|---|---|---|
| `type` | yes | `single` or `multi` |
| `input mesh file`, `output mesh file` | yes | [Mesh and output](mesh-and-io.md) |
| `model` | yes | physics and material assignment — [Model](model.md), [Materials](materials.md) |
| `time integrator` | yes | time-stepping scheme — [Time integrators](time-integrators.md) |
| `solver` | yes | solution algorithm — [Solvers](solvers.md) |
| `boundary conditions` | no | [Boundary conditions](boundary-conditions.md) |
| `initial conditions` | no | [Initial conditions](initial-conditions.md) |
| `Exodus output interval`, `CSV output interval`, `CSV write sidesets` | no | output cadence — [Mesh and output](mesh-and-io.md) |
| `restart`, `swaps` | no | restart and mesh swapping — [Multidomain and Schwarz coupling](multidomain.md) |

## A complete minimal example

An explicit dynamic simulation of a cantilever with a fixed end and an initial
displacement profile:

```yaml
type: single
input mesh file: cantilever.g
output mesh file: cantilever.e
model:
  type: solid mechanics
  material:
    blocks:
      cantilever: steel
    steel:
      model: linear elastic
      elastic modulus: 6.895e+09
      Poisson's ratio: 0.25
      density: 2768.0
time integrator:
  type: central difference
  initial time: 0.0
  final time: 3.0e-6
  CFL: 0.2
  γ: 0.5
initial conditions:
  displacement:
    - node set: nsall
      component: y
      function: "0.393701 * x * x"
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsx-
      component: y
      function: "0.0"
    - node set: nsx-
      component: z
      function: "0.0"
solver:
  type: explicit solver
  step: explicit
```

## Reference pages

- [Mesh and output](mesh-and-io.md)
- [Model](model.md)
- [Materials](materials.md)
- [Time integrators](time-integrators.md)
- [Solvers](solvers.md)
- [Boundary conditions](boundary-conditions.md)
- [Initial conditions](initial-conditions.md)
- [Function expressions](functions.md)
- [Multidomain and Schwarz coupling](multidomain.md)

## Conventions

- Keys are literal and case-sensitive, including Unicode keys (`β`, `γ`) and
  punctuation (`Poisson's ratio`, `Lamé's first constant`).
- Values that vary in space or time are [function expressions](functions.md) in
  `t, x, y, z`.
- Every reference page links to a canonical file under `examples/` that
  exercises the feature.

## Keyword validation

Every loaded input file is checked against the set of keys Norma actually
reads. A key that no parser recognizes — a misspelling, a stale key from an
older input format, or a key placed in the wrong section — produces a warning
with a suggestion:

```text
[WARNING] Input file 'cuboid-1.yaml': unknown key "adjoint paring" in
          Schwarz impedance nonoverlap boundary condition 1.
          Did you mean "adjoint pairing"?
```

Unknown keys warn rather than abort, and the run proceeds with the unknown
key ignored, exactly as before. Missing required keys and invalid values
remain hard errors reported by the section that reads them.
