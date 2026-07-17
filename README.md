[![CI](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml/badge.svg)](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/sandialabs/Norma.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sandialabs/Norma.jl)
[![License: BSD 3-Clause](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE.md)
[![Julia version](https://img.shields.io/badge/Julia-%E2%89%A5%201.10-blueviolet)](https://julialang.org/downloads/)

# Norma.jl

**Norma** is a Julia prototype for testing algorithms and ideas for
  coupling and multiphysics, primarily in solid mechanics and heat
  conduction.

![Norma Contact Simulation](https://github.com/sandialabs/Norma.jl/blob/main/doc/bars.gif)
*Impact simulation of two bars using different time integrators and
mesh types: the left bar uses hexahedral elements with an implicit
time integrator; the right bar uses tetrahedral elements with an
explicit time integrator. Each subdomain advances independently with
its own time step. (~100,000 elements)*

![Norma Torsion Simulation](https://github.com/sandialabs/Norma.jl/blob/main/doc/torsion.gif)
*Dynamic torsion of a solid cylinder undergoing large deformations. (~160,000 elements)*

![Norma Sphere Simulation](https://github.com/sandialabs/Norma.jl/blob/main/doc/sphere.gif)
*Large-deformation vibration of a soft rubber ball. The animation
plays in real time (100 frames per second) to match simulation
time. (~300,000 elements)*

---

## Quick Start

```bash
# Self-activating CLI wrapper (recommended)
bin/norma input.yaml

# Multi-threaded
bin/norma input.yaml --threads 8

# Or directly with julia
julia --project=. src/Norma.jl input.yaml
```

Interactive:
```julia
using Pkg; Pkg.activate(".")
using Norma
Norma.run("input.yaml")
```

---

## Documentation

Full documentation lives in the [`docs/`](docs/) directory and is built with
[Documenter.jl](https://documenter.juliadocs.org/). Build and view it locally
with:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
# then open docs/build/index.html
```

The documentation covers:

- [Installation](docs/src/installation.md)
- [Running Norma](docs/src/running.md)
- [Testing](docs/src/testing.md)
- [Examples](docs/src/examples.md)
- [Profiling and debugging](docs/src/profiling-and-debugging.md)
- [Troubleshooting](docs/src/troubleshooting.md)
- **[Input File Reference](docs/src/reference/index.md)** — every section and key
  of the YAML input file

Runnable input files for every capability are in the [`examples/`](examples/)
directory.

---

## Features

A prototyping framework for solid mechanics, multiphysics, and coupling
algorithms, with an emphasis on extensibility and experimentation:

- **Material models:** linear elastic, Saint-Venant–Kirchhoff, neohookean,
  Seth-Hill, and finite-deformation J2 plasticity.
- **Time integration:** quasi-static (with adaptive stepping), Newmark
  (implicit dynamic), and central difference (explicit dynamic).
- **Solvers:** Newton (Hessian minimizer) with optional line search,
  matrix-free steepest descent, and an explicit lumped-mass solver.
- **Elements:** BAR2; TRI3, TRI6, QUAD4; TETRA4, TETRA10, HEX8; with Gauss,
  Dunavant, and Keast quadrature.
- **Multidomain coupling:** overlapping and non-overlapping Schwarz, impedance
  (absorbing) coupling, subcycled multi-time-step coupling, and Schwarz contact.
- **Boundary conditions:** Dirichlet, Neumann traction and follower pressure,
  Robin, and inclined-surface (roller) support.
- **Adaptive mesh swapping**, **nodal field recovery**, and **Exodus/CSV output**.
- **Reduced-order models:** Operator Inference and RBF kernel ROMs, documented
  with the [`norma-opinf`](https://github.com/sandialabs/norma-opinf) companion
  package.

See the [documentation](docs/) for the full list and details.

---

## License

Norma is licensed under the BSD 3-Clause License. See [LICENSE.md](LICENSE.md)
for details.
