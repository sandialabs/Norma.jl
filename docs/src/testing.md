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

By default all tests run.

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
Python package and are enabled with a flag:

```bash
julia --project=.. runtests.jl --with-nnopinf
```
