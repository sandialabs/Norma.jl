[![CI](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml/badge.svg)](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml)
[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sandialabs.github.io/Norma.jl/dev/)
[![codecov](https://codecov.io/gh/sandialabs/Norma.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sandialabs/Norma.jl)
[![License: BSD 3-Clause](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE.md)
[![Julia version](https://img.shields.io/badge/Julia-%E2%89%A5%201.10-blueviolet)](https://julialang.org/downloads/)

# Norma.jl

**Norma** is a Julia prototype for testing algorithms and ideas for
  coupling and multiphysics, primarily in solid mechanics and heat
  conduction.

**Documentation:** [https://sandialabs.github.io/Norma.jl/dev/](https://sandialabs.github.io/Norma.jl/dev/)

![Norma Contact Simulation](https://github.com/sandialabs/Norma.jl/blob/main/docs/assets/bars.gif)
*Impact simulation of two bars using different time integrators and
mesh types: the left bar uses hexahedral elements with an implicit
time integrator; the right bar uses tetrahedral elements with an
explicit time integrator. Each subdomain advances independently with
its own time step. (~100,000 elements)*

![Norma Torsion Simulation](https://github.com/sandialabs/Norma.jl/blob/main/docs/assets/torsion.gif)
*Dynamic torsion of a solid cylinder undergoing large deformations. (~160,000 elements)*

![Norma Sphere Simulation](https://github.com/sandialabs/Norma.jl/blob/main/docs/assets/sphere.gif)
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

## Formatting

The Julia sources are formatted with JuliaFormatter in the blue style with a
line width of 120, as set in `.JuliaFormatter.toml`. The formatter version is
pinned in `format/Manifest.toml`, and CI rejects a change whose files are not
formatted with that version. Format the sources before committing with:

```bash
julia --project=format -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=format format/format.jl
```

The commit that formatted the whole tree is listed in `.git-blame-ignore-revs`;
`git config blame.ignoreRevsFile .git-blame-ignore-revs` makes `git blame`
skip it.

---

## Documentation

The documentation is published at [https://sandialabs.github.io/Norma.jl/dev/](https://sandialabs.github.io/Norma.jl/dev/). It is rebuilt
from the latest commit on `main` and covers:

- [Installation](https://sandialabs.github.io/Norma.jl/dev/installation/)
- [Running Norma](https://sandialabs.github.io/Norma.jl/dev/running/)
- [Testing](https://sandialabs.github.io/Norma.jl/dev/testing/)
- [Examples](https://sandialabs.github.io/Norma.jl/dev/examples/)
- [Profiling and debugging](https://sandialabs.github.io/Norma.jl/dev/profiling-and-debugging/)
- [Troubleshooting](https://sandialabs.github.io/Norma.jl/dev/troubleshooting/)
- **[Input File Reference](https://sandialabs.github.io/Norma.jl/dev/reference/)**: every section and key
  of the YAML input file

The sources are in the [`docs/`](docs/) directory and are built with
[Documenter.jl](https://documenter.juliadocs.org/). To build and view them
locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
# then open docs/build/index.html
```

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
- **Multidomain coupling:** overlapping and non-overlapping Schwarz, the
  matching Dirichlet–Neumann exchange for dynamics, subcycled
  multi-time-step coupling, and Schwarz contact.
- **Boundary conditions:** Dirichlet, Neumann traction and follower pressure,
  Robin, and inclined-surface (roller) support.
- **Adaptive mesh swapping**, **nodal field recovery**, and **Exodus/CSV output**.
- **Mesh improvement:** energy-minimization smoothing of tetrahedral meshes
  toward a prescribed size or anisotropic metric field, alternated with edge
  swaps, collapses, and splits.
- **Reduced-order models:** Operator Inference and RBF kernel ROMs, documented
  with the [`norma-opinf`](https://github.com/sandialabs/norma-opinf) companion
  package.

See the [documentation](docs/) for the full list and details.

---

## Citing Norma

If you use Norma in published work, please cite it. The file
[CITATION.cff](CITATION.cff) lists the authors and the version; GitHub
renders it under "Cite this repository" in BibTeX and APA formats. Each
release is archived on Zenodo with its own DOI; cite the release used to
produce the results.

## License

Norma is licensed under the BSD 3-Clause License. See [LICENSE.md](LICENSE.md)
for details.
