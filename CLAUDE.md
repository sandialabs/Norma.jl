# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is Norma.jl?

A Julia prototype for testing algorithms for coupling and multiphysics, primarily in solid mechanics and heat conduction. Supports:
- Full Order Models (FOM): nonlinear FEM, implicit/explicit time integrators, contact, large deformations
- Reduced Order Models (ROM): Operator Inference (OpInf) via Python companion package `norma-opinf`
- Domain decomposition: Schwarz methods (overlapping, non-overlapping, contact-based)

## Development Setup

```julia
# In Julia REPL, enter package manager with ]
pkg> activate .
pkg> registry update
pkg> update
pkg> instantiate
```

## Running Norma

```bash
# Command-line
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml

# Multi-threaded
JULIA_NUM_THREADS=4 julia --project=. src/Norma.jl input.yaml

# Debug logging
JULIA_DEBUG=Norma julia --project=. src/Norma.jl input.yaml
```

```julia
# Interactive
using Pkg; Pkg.activate(".")
using Norma
cd("examples/contact/implicit-dynamic/2-bars")
Norma.run("bars.yaml")
# After code changes, reload: using Norma
```

## Testing

```bash
# Full test suite
julia --project=. test/runtests.jl

# Quick subset only
julia --project=. test/runtests.jl --quick

# Run tests by index
julia --project=. test/runtests.jl 1 3 5

# Filter by name (case-insensitive)
julia --project=. test/runtests.jl --filter cube
julia --project=. test/runtests.jl --quick --filter dynamic

# List all tests
julia --project=. test/runtests.jl --list

# With neural-network OpInf (requires norma-opinf Python package)
julia --project=. test/runtests.jl --with-nnopinf
```

## Code Architecture

### Entry Points
- `src/Norma.jl` — module entry; defines `Norma.run(input_file)` which creates and evolves a simulation
- Input files are YAML; examples live in `examples/`
- Meshes use Exodus format (`.g`, `.exo`)

### Paired `_types.jl` / `.jl` Convention
Every major component splits into:
- `*_types.jl` — struct/type definitions
- `*.jl` — methods operating on those types

### Core Components

| File(s) | Role |
|---|---|
| `simulation_types.jl` / `simulation.jl` | `SingleDomainSimulation` and `MultiDomainSimulation`; top-level `create_simulation` and `evolve` |
| `model_types.jl` / `model.jl` | `SolidMechanics` struct (mesh, DOFs, forces, stiffness/mass matrices, kinematics) |
| `constitutive_types.jl` / `constitutive.jl` | Material hierarchy: `Material → Solid → {Elastic, Inelastic}`; stress computation |
| `solver_types.jl` / `solver.jl` | `HessianMinimizer` (Newton), `ExplicitSolver`, `SteepestDescent`; line search |
| `time_integrator_types.jl` / `time_integrator.jl` | Implicit (Newmark) and explicit time stepping schemes |
| `ics_bcs_types.jl` / `ics_bcs.jl` | Initial and boundary conditions |
| `interpolation_types.jl` / `interpolation.jl` | FEM shape functions and quadrature |
| `io.jl` | Exodus mesh I/O, CSV output for ROM snapshots, mesh smoothing |
| `minitensor.jl` | Tensor algebra utilities (used in constitutive models and tests) |
| `utils.jl` | `norma_log` (colored logging), `norma_abort` (error with exceptions), FP exception traps |

### OpInf ROM Submodule (`src/opinf/`)
Mirrors the FOM structure with `opinf_*` prefixed files. ROM workflow:
1. Run FOM with `CSV output interval` and `CSV write sidesets` enabled to collect snapshots
2. Run `norma-opinf` (Python) to build operator `.npz` file
3. Run Norma with ROM input (e.g., `model: type: linear opinf rom; model-file: opinf-operator.npz`)

### Multi-Domain Simulation
`MultiDomainSimulation` holds multiple sub-simulations plus a `SolidMultiDomainTimeController` for Schwarz iteration management. Each subdomain can use a different time integrator and time step.

### Key Data Flow
```
Norma.run(yaml) → create_simulation → SingleDomain or MultiDomain
                                            ↓
                               evolve(sim) loop:
                               - advance time controller
                               - solver.solve → model residual/stiffness
                               - time_integrator.update fields
                               - io.write output (Exodus / CSV)
```
