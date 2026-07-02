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

## What is Norma.jl?

- A prototyping framework for multiphysics and coupling algorithms
- Focused on solid mechanics and heat conduction
- Designed for high extensibility and experimentation
- Supports implicit and explicit time integrators

---

## **Table of Contents**
1. [Features](#features)
2. [Installation](#installation)
3. [Running the Code](#running-the-code)
4. [Testing](#testing)
   - [Selective Test Execution](#selective-test-execution)
   - [Filtering by Name](#filtering-by-name)
5. [Examples](#examples)
6. [Profiling](#profiling)
7. [Debugging](#debugging)
8. [Troubleshooting](#troubleshooting)
9. [License](#license)
10. [Documentation](#documentation)

---

## **Features**

A prototyping framework for solid mechanics, multiphysics, and coupling
algorithms, with an emphasis on extensibility and experimentation. The
capabilities currently implemented are:

**Material models** (infinitesimal- and finite-strain kinematics)
- Linear elastic
- Saint-Venant–Kirchhoff
- Neohookean
- Seth-Hill (generalized hyperelastic with tunable exponents)
- J2 plasticity (finite deformation, multiplicative split, radial return, linear isotropic hardening)

**Time integration**
- Quasi-static (with adaptive stepping)
- Newmark (implicit dynamic)
- Central difference (explicit dynamic, CFL-controlled)

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
- Single-domain and multi-domain simulations
- Overlapping and non-overlapping (Dirichlet–Neumann) Schwarz
- Robin–Robin and impedance-matching (absorbing) Schwarz
- Frictionless/tied contact via Schwarz
- Fixed and adaptive Aitken relaxation — Aitken recursive (Irons–Tuck) and Aitken secant — with interface predictor acceleration

**Boundary conditions**
- Dirichlet (time-dependent displacement/velocity/acceleration)
- Neumann traction and pressure
- Robin (mixed displacement–traction)
- Inclined-surface support via rotation matrices

**Reduced-order models**
- Operator Inference (OpInf): linear, quadratic, cubic, and neural-network ROMs (the neural-network variant needs the optional PyCall + PyTorch backend; the rest are pure Julia)
- RBF kernel ROMs (Gaussian, inverse quadratic, inverse multiquadric, thin-plate-spline, Matérn variants)

**Adaptive mesh swapping**
- Triggered by time, stress-recovery error, overlap L2 displacement error, or elastic-to-plastic transition

**Nodal recovery and field transfer**
- L2-projection recovery (lumped and consistent mass) of stress, von Mises stress, deformation gradient, and internal variables
- Mesh-to-mesh field transfer for multidomain and remeshing workflows

**Input/output**
- YAML input configuration
- Exodus II output (nodal, element, and global variables)
- CSV output of nodal and side-set quantities

---

## **Installation**

### Clone the Repository
```bash
cd /path/to
git clone git@github.com:sandialabs/Norma.jl.git
cd Norma.jl
julia
```

### Set Up the Environment
Within the Julia package manager (enter by pressing `]` in the Julia REPL):
```julia
pkg> activate .
pkg> registry update
pkg> update
pkg> instantiate
```
Press `Backspace` or `Delete` to exit the package manager.

This installs **no Python dependencies**. Norma's full-order solver, Schwarz
coupling, the linear/quadratic/cubic OpInf ROMs, and the RBF kernel ROMs are all
pure Julia. Only the **neural-network OpInf ROM** needs Python (PyTorch via
PyCall); it is an optional [package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).
To enable it, add PyCall to the environment and make sure its Python has
`torch` installed:
```julia
pkg> add PyCall
```
Julia then loads the backend automatically. Without it, requesting a
`neural network opinf rom` aborts with a message telling you to install PyCall.

---

### If Interested in Reduced Order Model (ROM) Capabilities in Norma: Clone and Build Norma-OpInf Repository
```bash
cd /path/to
git clone git@gitlab-ex.sandia.gov:ejparis/norma-opinf.git
cd /path/to/norma-opinf
pip install -e .
```
For this, Python3 is required.  For more details and troubleshooting, please see the [norma-opinf repository](https://github.com/sandialabs/norma-opinf).

## **Running the Code**

### Running the Code in Full Order Model (FOM) Mode (Default)

To run the main program, assuming Julia is in your executable path:
```bash
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

To run Norma interactively from a Julia session:
```julia
cd /path/to/Norma.jl
julia
using Pkg
Pkg.activate(".")
using Norma
```
Then, navigate to your desired example folder and run the simulation. For example:
```julia
cd("examples/contact/implicit-dynamic/2-bars")
Norma.run("bars.yaml")
```

**Note**: If you make changes to the Norma code, you need to reload the Norma module (`using Norma`) for those changes to take effect.

### Running the Code with Operator Inference (OpInf) Reduced Order Models (ROMs)

Running Norma with OpInf ROMs is a process consisting of three steps.  More details can be found in the README file found [here](https://github.com/sandialabs/Norma.jl/blob/main/examples/ahead/overlap/cuboid/dynamic-opinf-fom/README.md).

#### Step 1: Run Norma in FOM Mode to Generate Training Data for ROM

Run the main program in FOM mode, assuming Julia is in your executable path:
```bash
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```
after modifying ```input.yaml``` to enable ```CSV output``` and ```CSV write sidesets```, e.g.,
```
CSV output interval: {time_step}
CSV write sidesets: true
```
An example input file can be found [here](https://github.com/sandialabs/Norma.jl/blob/main/examples/ahead/single/cuboid/dynamic-opinf-fom/cuboid.yaml).

#### Step 2: Run Norma-OpInf to Build ROM from Snapshot Data Generated in Step 1

Please see the [Norma-OpInf repo Wiki Page](https://github.com/sandialabs/norma-opinf) for details.

#### Step 3: Run Norma in ROM Mode to Obtain ROM Solution

Run the main program in ROM mode, assuming Julia is in your executable path:
```bash
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input_rom.yaml
```
after modifying ```input_rom.yaml``` to utilize a ROM model type and read in the ```.npz``` file produced in Step 2, e.g.,
```
model:
  type: linear opinf rom
  model-file: opinf-operator.npz
```
An example ROM input file can be found [here](https://github.com/sandialabs/Norma.jl/blob/main/examples/ahead/single/cuboid/dynamic-opinf-rom/cuboid.yaml).

### Running with Multiple Threads

```bash
bin/norma input.yaml --threads 8
```

Or via the environment variable:
```bash
JULIA_NUM_THREADS=8 julia --project=. src/Norma.jl input.yaml
```

Inside Julia:
```julia
using Pkg
Pkg.activate("/path/to/Norma.jl")
using Norma
Norma.run("input.yaml")
```

---

## **Testing**

To run the test suite using the Julia REPL:
```julia
using Pkg
Pkg.test()
```

Or from the command line:
```bash
cd /path/to/Norma.jl/test
julia --project=@/path/to/Norma.jl ./runtests.jl
```

By default, all tests are run.

---

### Selective Test Execution

You can control which tests are run using command-line arguments.

#### Run a Quick Subset (Fast Tests Only)
```bash
julia --project=@/path/to/Norma.jl ./runtests.jl --quick
```
Use this when you want to verify functionality quickly without running the full suite.

#### Run Specific Tests by Index
```bash
julia --project=@/path/to/Norma.jl ./runtests.jl 1 3 5
```

#### List All Available Tests
```bash
julia --project=@/path/to/Norma.jl ./runtests.jl --list
```

---

### Filtering by Name

To run tests whose filenames contain a given string (case-insensitive):

```bash
julia --project=@/path/to/Norma.jl ./runtests.jl --filter cube
```

You can combine filters with specific indices or `--quick`:

```bash
julia --project=@/path/to/Norma.jl ./runtests.jl --quick --filter dynamic
julia --project=@/path/to/Norma.jl ./runtests.jl 2 4 --filter static
```

---

## **Examples**

To run the `examples/contact/implicit-dynamic/2-bars` example:
```bash
cd /path/to/Norma.jl/examples/contact/implicit-dynamic/2-bars
julia
]
activate .
using Norma
Norma.run("bars.yaml")
```

---

## **Profiling**

To identify performance bottlenecks in Norma.jl, you can use Julia's built-in `Profile` module and visualize results with [`ProfileCanvas.jl`](https://github.com/SciML/ProfileCanvas.jl):

### Step 1: Enable Profiling
```julia
using Profile
include("/path/to/Norma.jl/src/Norma.jl")
cd("/path/to/Norma.jl/examples/contact/implicit-dynamic/2-bars")
@profile Norma.run("bars.yaml")
```

### Step 2: Visualize with ProfileCanvas
```julia
using Pkg; Pkg.add("ProfileCanvas")
using ProfileCanvas
ProfileCanvas.canvas()
```

### Command-Line Workflow
```bash
julia --project=@/path/to/Norma.jl -e 'using Profile; using Norma; cd("examples/contact/implicit-dynamic/2-bars"); @profile Norma.run("bars.yaml")' -E 'using ProfileCanvas; ProfileCanvas.canvas()'
```

---

## **Debugging**

To enable debug-level logging in Norma.jl, use the `JULIA_DEBUG` environment variable:

```bash
JULIA_DEBUG=Norma julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

To add debug messages in code:
```julia
@debug "Starting simulation with input file: input.yaml"
```

To disable debug printing:
```bash
unset JULIA_DEBUG
```

Or suppress it at launch:
```bash
JULIA_DEBUG= julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

---

## **Troubleshooting**

If you are on Sandia's network and run into SSL/TLS certificate errors when
installing or fetching packages, see [README-sandia.md](README-sandia.md) for a
complete setup guide.

---

## **License**

Norma.jl is licensed under the BSD 3-Clause License. See [LICENSE.md](LICENSE.md) for details.

## **Documentation**

While we do not have documentation of the Norma.jl code, we have documented a number of application-motivated test cases available in the Norma.jl repository in the [following slides](https://docs.google.com/presentation/d/1UI1dOn5RU4L_nNXgEQoaSKIHXUv97TPn/edit?usp=drive_link&ouid=118047221648309868055&rtpof=true&sd=true).  The input and mesh files for these test cases can be found in the ```Norma.jl/examples/ahead``` directory.
