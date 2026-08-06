# Testing

Run the test suite from the Julia REPL:

```julia
using Pkg
Pkg.test()
```

or from the command line, in the `test/` directory (test files use relative
paths such as `../examples/`):

```bash
cd /path/to/Norma.jl/test
julia --project=@/path/to/Norma.jl ./runtests.jl
```

By default every test runs except the neural-network Operator-Inference tests,
which depend on an optional Python stack and are enabled with a flag (see
[Neural-network reduced-order tests](#neural-network-reduced-order-tests)
below). The suite reports which tests it skipped for this reason.

## Selective execution

Run a quick subset (fast tests only, excluding the slow coupling cases):

```bash
julia --project=.. runtests.jl --quick
```

Run specific tests by index:

```bash
julia --project=.. runtests.jl 1 3 5
```

List all available tests:

```bash
julia --project=.. runtests.jl --list
```

## Filtering by name

Run tests whose filenames contain a string (case-insensitive):

```bash
julia --project=.. runtests.jl --filter cube
```

Combine filters with indices or `--quick`:

```bash
julia --project=.. runtests.jl --quick --filter dynamic
julia --project=.. runtests.jl 2 4 --filter static
```

## Neural-network reduced-order tests

The neural-network Operator-Inference tests require the optional `norma-opinf`
Python package, along with PyCall built against a Python that has `torch`
installed. They load PyCall at include time to reach the `NormaPyTorchExt`
package extension, so they would error rather than fail without that stack.
For that reason they are excluded from the default suite and enabled with a
flag:

```bash
julia --project=.. runtests.jl --with-nnopinf
```

This flag runs the entire suite, the neural-network tests included. Continuous
integration provisions the Python stack and always uses it, so these tests are
covered there on every pull request. Selecting a test by index also bypasses
the flag, so `runtests.jl 52` runs the neural-network test on its own.
