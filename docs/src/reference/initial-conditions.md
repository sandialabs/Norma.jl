# Initial conditions

The `initial conditions` block sets the initial displacement and velocity
fields. It groups conditions by type, each mapping to a **list** of entries. The
two types are `displacement` and `velocity`:

```yaml
initial conditions:
  displacement:
    - node set: nsall
      component: y
      function: "0.393701 * x * x"
  velocity:
    - node set: nsall
      component: z
      function: "10.0"
```

Each entry uses the same fields:

| Key | Required | Meaning |
|---|---|---|
| `node set` | yes | node set to initialize |
| `component` | yes | `x`, `y`, or `z` |
| `function` | yes | [function expression](functions.md) in `t, x, y, z` |

For a `displacement` entry, the `function` gives the initial displacement; the
initial velocity is taken from its time derivative when that is nonzero (see
[Function expressions](functions.md)). For a `velocity` entry, the `function`
sets the velocity directly.

The block is optional. It cannot be combined with a `restart` (the restart
snapshot already carries the initial state). Conflicting velocity assignments to
the same node and component are rejected.

## Canonical examples

- Initial displacement profile: `examples/single/static-solid`,
  `examples/single/implicit-dynamic-solid`
