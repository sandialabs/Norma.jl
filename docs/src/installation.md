# Installation

## Requirements

- Julia 1.10 or newer.

Norma's full-order solver, Schwarz coupling, and kernel/Operator-Inference
reduced-order models are pure Julia and require no Python. Only the optional
neural-network Operator-Inference model needs Python (see below).

## Clone the repository

```bash
cd /path/to
git clone git@github.com:sandialabs/Norma.jl.git
cd Norma.jl
julia
```

## Set up the environment

In the Julia package manager (press `]` in the REPL):

```julia
pkg> activate .
pkg> registry update
pkg> update
pkg> instantiate
```

Press `Backspace` or `Delete` to leave the package manager.

## Optional: neural-network reduced-order models

Only the neural-network Operator-Inference model needs Python (PyTorch via
PyCall); it is an optional package extension. To enable it, add PyCall to the
environment and make sure its Python has `torch` installed:

```julia
pkg> add PyCall
```

Julia then loads the backend automatically. Without it, requesting a neural-network
model aborts with a message telling you to install PyCall. The reduced-order
capabilities as a whole are documented with the
[`norma-opinf`](https://github.com/sandialabs/norma-opinf) companion package and
are outside the scope of this guide.

## Sandia network

If you are on Sandia's network and hit SSL/TLS certificate errors when
installing or fetching packages, see `README-sandia.md` in the repository for a
complete setup guide, and [Troubleshooting](troubleshooting.md).
