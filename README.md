[![CI](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml/badge.svg)](https://github.com/sandialabs/Norma.jl/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/sandialabs/Norma.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sandialabs/Norma.jl)
[![License: BSD 3-Clause](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE.md)
[![Julia version](https://img.shields.io/badge/Julia-%E2%89%A5%201.10-blueviolet)](https://julialang.org/downloads/)

# Norma.jl

**Norma** is a Julia prototype for testing algorithms and ideas for coupling and multiphysics, primarily in solid mechanics and heat conduction.

![Norma Contact Simulation 1](https://github.com/sandialabs/Norma.jl/blob/main/doc/bars.gif)
*Simulation of the impact of two bars: one using hexahedral elements with an implicit time integrator, and the other using tetrahedral elements with an explicit time integrator, each with different time steps.*

![Norma Torsion Simulation](https://github.com/sandialabs/Norma.jl/blob/main/doc/torsion.gif)
*Dynamic simulation of torsion with large deformations.*

---

## Quick Start

```bash
julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

Or run it interactively:
```julia
using Pkg; Pkg.activate("/path/to/Norma.jl")
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
5. [Examples](#examples)
6. [Profiling](#profiling)
7. [Debugging](#debugging)
8. [Troubleshooting](#troubleshooting)
9. [License](#license)

---

## **Features**
- Prototyping of coupling and multiphysics algorithms.
- Applications in solid mechanics and heat conduction.
- Designed for extensibility and experimentation.

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

---

## **Running the Code**

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

### Running with Multiple Threads

To run Norma.jl using multiple threads, set the `JULIA_NUM_THREADS` environment variable before launching Julia. For example, to use 4 threads:

```bash
JULIA_NUM_THREADS=4 julia --project=@/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl input.yaml
```

Or for interactive usage:
```bash
JULIA_NUM_THREADS=4 julia
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

### Selective Test Execution

To run only a subset of tests, pass one or more test indices (shown at startup) as arguments:

```bash
julia --project=@/path/to/Norma.jl ./runtests.jl 1 3 5
```

This will run tests 1, 3, and 5 only. The list of available tests and their corresponding indices is printed automatically at the start of the run.

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

### SSL Certificate Issues
```bash
cd ~/.julia/registries
git clone https://github.com/JuliaRegistries/General.git
export JULIA_SSL_CA_ROOTS_PATH=/etc/ssl/certs/ca-bundle.crt
```
Then retry installation.

---

## **License**

Norma.jl is licensed under the BSD 3-Clause License. See [LICENSE.md](LICENSE.md) for details.

