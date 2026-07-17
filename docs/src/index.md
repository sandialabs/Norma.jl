# Norma

**Norma** is a Julia prototype for testing algorithms and ideas for coupling and
multiphysics, primarily in solid mechanics and heat conduction. It is a
prototyping framework built for extensibility and experimentation, supporting
implicit and explicit time integration, large deformations, contact, and domain
decomposition by Schwarz methods.

```@raw html
<img src="https://raw.githubusercontent.com/sandialabs/Norma.jl/main/docs/assets/bars.gif" alt="Contact simulation of two bars" style="max-width:100%;">
```
*Impact of two bars using different time integrators and mesh types: the left
bar uses hexahedral elements with an implicit integrator; the right bar uses
tetrahedral elements with an explicit integrator. Each subdomain advances with
its own time step.*

```@raw html
<img src="https://raw.githubusercontent.com/sandialabs/Norma.jl/main/docs/assets/torsion.gif" alt="Torsion simulation" style="max-width:100%;">
```
*Dynamic torsion of a solid cylinder undergoing large deformations.*

```@raw html
<img src="https://raw.githubusercontent.com/sandialabs/Norma.jl/main/docs/assets/sphere.gif" alt="Vibrating sphere simulation" style="max-width:100%;">
```
*Large-deformation vibration of a soft rubber ball, playing in real time.*

## Quick start

```bash
# self-activating CLI wrapper
bin/norma input.yaml

# multi-threaded
bin/norma input.yaml --threads 8

# or directly with Julia
julia --project=. src/Norma.jl input.yaml
```

Interactive:

```julia
using Pkg; Pkg.activate(".")
using Norma
Norma.run("input.yaml")
```

See [Installation](installation.md) to set up the environment and
[Running Norma](running.md) for the full set of ways to launch a simulation.

## Documentation

- [Installation](installation.md)
- [Running Norma](running.md)
- [Testing](testing.md)
- [Examples](examples.md)
- [Profiling and debugging](profiling-and-debugging.md)
- [Troubleshooting](troubleshooting.md)
- **[Input File Reference](reference/index.md)** — every section and key of the YAML input file

## Features

Capabilities currently implemented:

**Material models** (infinitesimal- and finite-strain kinematics)

- Linear elastic
- Saint-Venant–Kirchhoff
- Neohookean
- Seth-Hill (generalized hyperelastic with tunable exponents)
- J2 plasticity (finite deformation, multiplicative split, radial return, linear isotropic hardening)

**Time integration**

- Quasi-static (with adaptive stepping)
- Newmark (implicit dynamic)
- Central difference (explicit dynamic)

**Solvers**

- Newton's method (Hessian minimizer) with optional backtracking line search
- Matrix-free steepest descent
- Explicit (lumped-mass) solver

**Element library**

- 1D: BAR2
- 2D: TRI3, TRI6, QUAD4
- 3D: TETRA4, TETRA10, HEX8
- Gauss, Dunavant, and Keast quadrature rules

**Multidomain coupling — Schwarz alternating method**

- Overlapping and non-overlapping (Dirichlet–Neumann) Schwarz
- Robin–Robin and impedance-matching (absorbing) Schwarz, including the
  adjoint-paired non-overlapping impedance coupling (see
  `docs/notes/schwarz-coupling`)
- Multi-time-step (subcycled) coupling with time interpolation of the partner state
- Frictionless and tied contact via Schwarz
- Fixed and adaptive (Aitken) relaxation with interface prediction

**Boundary conditions**

- Dirichlet (time-dependent displacement/velocity/acceleration)
- Neumann traction and follower pressure
- Robin (mixed displacement–traction)
- Inclined-surface (roller) support

**Adaptive mesh swapping**

- Triggered by time, stress-recovery error, overlap displacement error, or the elastic-to-plastic transition

**Nodal recovery and field transfer**

- L2-projection recovery (lumped and consistent) of stress, von Mises stress, deformation gradient, and internal variables
- Mesh-to-mesh field transfer for multidomain and remeshing workflows

**Input/output**

- YAML input configuration
- Exodus II output (nodal, element, and global variables)
- CSV output of nodal and side-set quantities
- Overlap-corrected (Arlequin-blended) energy output for multidomain runs

**Reduced-order models** — Operator Inference and RBF kernel ROMs are available
and documented with the [`norma-opinf`](https://github.com/sandialabs/norma-opinf)
companion package; they are outside the scope of this guide.

## License

Norma is licensed under the BSD 3-Clause License. See `LICENSE.md` in the
repository for details.
