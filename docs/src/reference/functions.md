# Function expressions

Many input values — boundary-condition data, initial conditions, Neumann
tractions and pressures, level-set surfaces, and mesh-smoothing size fields —
are given as `function` strings rather than constants. Each string is a
symbolic expression evaluated at each node.

## Available symbols

Exactly four symbols are available in every `function` expression:

| Symbol | Meaning |
|---|---|
| `t` | current simulation time |
| `x` | reference (undeformed) x-coordinate of the node |
| `y` | reference y-coordinate of the node |
| `z` | reference z-coordinate of the node |

The spatial coordinates are always the **reference** configuration
coordinates, not the deformed ones.

## Syntax

Expressions are parsed and compiled with the Symbolics package, so any Julia
numeric syntax that Symbolics understands is valid: arithmetic (`+ - * / ^`),
elementary functions (`sin`, `cos`, `exp`, `log`, `sqrt`, `abs`, `tanh`, …),
the constant `pi`, and `ifelse(condition, a, b)` for piecewise definitions.
Multiple statements separated by semicolons are allowed, with the last
expression supplying the value:

```yaml
function: "a = 0.001; s = 0.02; a * exp(-z*z / s/s / 2)"
```

Constants are a single number written as a string:

```yaml
function: "0.0"
```

A linear ramp in time:

```yaml
function: "0.01 * t"
```

A smooth pressure pulse:

```yaml
function: "-937.5e3 * 1500 * (0.5 - 0.5 * cos(pi * t / 1500.0))"
```

## Time derivatives (Dirichlet and initial conditions)

For a `Dirichlet` boundary condition and for a `displacement` initial
condition, the `function` gives the **displacement**. The corresponding
velocity and acceleration are obtained automatically by symbolic
differentiation in `t` — you do not supply them separately. A `velocity`
initial condition sets the velocity directly from its `function`.

Because differentiation is symbolic, a displacement expression that depends on
`t` produces the physically consistent velocity and acceleration; a
time-independent displacement produces zero velocity and acceleration.
