# Multidomain and Schwarz coupling

A multidomain simulation couples two or more single-domain models with the
Schwarz alternating method. It has two parts: a **top-level controller file**
(`type: multi`) that lists the subdomains and governs the coupling, and one
**subdomain input file** per domain — each an ordinary single-domain input with
Schwarz coupling boundary conditions added.

```yaml
type: multi
domains: ["cuboid-1.yaml", "cuboid-2.yaml"]
initial time: 0.0
final time: 1.0
time step: 1.0e-3
minimum iterations: 1
maximum iterations: 16
relative tolerance: 1.0e-12
absolute tolerance: 1.0e-08
```

## Controller file (`type: multi`)

### Domains and time window

| Key | Required | Meaning |
|---|---|---|
| `domains` | yes | list of subdomain input files, each a full single-domain input |
| `initial time` | yes | coupled simulation start time |
| `final time` | yes | coupled simulation end time |
| `time step` | yes | controller (Schwarz coupling) time step |

The controller supplies the time window to every subdomain, overriding any
`initial time`/`final time` in the subdomain files. Each subdomain keeps its own
time integrator and its own `time step`, which is clamped to at most the
controller step; a subdomain with a smaller step subcycles within each
controller step.

The controller file does **not** contain `model`, `solver`, or `time integrator`
blocks — those live in the subdomain files.

### Schwarz iteration control

Each controller step performs Schwarz iterations until the interface converges.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `minimum iterations` | yes | — | minimum Schwarz iterations per step |
| `maximum iterations` | yes | — | maximum Schwarz iterations per step |
| `absolute tolerance` | yes | — | absolute interface convergence tolerance |
| `relative tolerance` | yes | — | relative interface convergence tolerance |

### Relaxation and acceleration

| Key | Required | Default | Meaning |
|---|---|---|---|
| `relaxation` | no | fixed | `aitken recursive` (Irons–Tuck) or `aitken secant` adaptive relaxation; omit for a fixed factor |
| `relaxation parameter` | no | `1.0` | relaxation factor θ; the constant factor under fixed relaxation, and under either Aitken method the factor used wherever an adaptive one is not yet available |
| `aitken N0 parameter` | no | `1` | Schwarz iteration, counting from zero, at which the adaptive factor takes over from `relaxation parameter` (Aitken methods only) |
| `interface predictor` | no | `false` | extrapolate the interface state at the start of each step |
| `naive stabilized` | no | `false` | naive interface stabilization |

Relaxation is not tied to a particular transmission condition: it applies to
whatever datum a coupling boundary condition transmits. For `Schwarz overlap`
and `Schwarz DN nonoverlap` the relaxed quantity is the interface displacement;
for `Schwarz impedance nonoverlap`, and therefore for its `Schwarz RR
nonoverlap` alias, it is the impedance right-hand side. Both Aitken forms work
with all of them, and `relaxation parameter` and `aitken N0 parameter` mean the
same thing in each. On the impedance and Robin conditions the acceleration is
substantial: on the cantilever benchmark either Aitken form converges in about
a tenth of the sweeps that a fixed factor needs.

`relaxation parameter` is not ignored when `relaxation` names an Aitken method.
It is the factor applied for Schwarz iterations below `aitken N0 parameter`,
and the fallback whenever an adaptive factor cannot be formed, which includes
the first iterations of every step before two iterates exist to compare.

Aitken acceleration is applied only to stops with a single substep, that is
when the controller `time step` equals the relaxed subdomain's time step. A
windowed stop couples all of its time slots in one sweep, where the adaptive
factors were measured to diverge or to lose to a fixed factor, so such stops
use `relaxation parameter` throughout. This is automatic and needs no input.

If a relaxation factor ever becomes small enough to leave an interface iterate
unchanged, the sweep carries no information: every subdomain re-solves against
the data it already had and returns the solution it already had. The
displacement-based convergence test cannot distinguish that from convergence,
so such a sweep is refused as evidence of convergence and the run logs
`Relaxation factor near zero froze an interface iterate`. Seeing that message
repeatedly means the coupling is not advancing; check `relaxation parameter`.

### Output

| Key | Required | Default | Meaning |
|---|---|---|---|
| `Exodus output interval` | no | the controller `time step` | applied uniformly to all subdomains |
| `CSV output interval` | no | `0.0` (disabled) | applied uniformly to all subdomains |
| `blended energy output` | no | `false` | write an Arlequin-blended kinetic/stored/total-energy CSV each stop, removing the double count in overlapping regions |

## Schwarz coupling boundary conditions (subdomain files)

Inside each subdomain's `boundary conditions` block, coupling to a partner
subdomain is expressed with one of the Schwarz condition types below. All share
these keys:

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source` | yes | — | name of the coupled subdomain (its input-file basename) |
| `side set` | yes | — | this subdomain's coupling surface |
| `source side set` | for non-overlapping | `""` | the partner's coupling surface |
| `source block` | for overlapping | `""` | the partner element block searched for the overlap |
| `search tolerance` | no | `1.0e-6` | geometric search tolerance for locating partner points |

### `Schwarz overlap`

Dirichlet-overlap coupling: the subdomain reads its partner's displacement in
the shared overlap region.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block covering the overlap |
| `weak` | no | `false` | weak (integrated) rather than pointwise transfer |
| `compute overlap L2 relative error` | no | `""` | report the overlap L2 relative error of `disp`, `velo`, or `acce` (drives the overlap-error mesh-swap criterion) |

### `Schwarz DN nonoverlap`

Non-overlapping Dirichlet–Neumann coupling across a shared interface. The two
sides must take opposite roles.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `default BC type` | no | `Dirichlet` | this side's role: `Dirichlet` or `Neumann` (the two sides must be opposite) |
| `swap BC types` | no | `false` | swap the Dirichlet/Neumann roles between Schwarz iterations |

### `Schwarz impedance nonoverlap`

Non-overlapping impedance (absorbing) coupling — the default and recommended
non-overlapping method. The interface transmits the partner traction plus a
dashpot term, making the interface energy exchange dissipative. See
`docs/notes/schwarz-coupling` for the theory.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `robin parameter` | no | `0.0` | Robin coefficient α; must be identical on both sides; affects convergence rate only, not the converged solution |
| `impedance scale` | no | `1.0` | scalar scaling of the dashpot impedance (`0.0` disables it); must be ≥ 0 |
| `adjoint pairing` | no | `true` | use the adjoint-paired shared cross-mass transfer (recommended); `false` restores the legacy per-side transfer |

`Schwarz RR nonoverlap` is a deprecated alias for this condition.

### `Schwarz impedance overlap`

Overlapping impedance coupling with recovered or consistent partner tractions.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block covering the overlap |
| `robin parameter` | no | `0.0` | Robin coefficient α |
| `impedance scale` | no | `1.0` | dashpot impedance scaling (scalar, or a P/S-split schedule) |
| `partner traction` | no | `auto` | `auto`, `consistent traction`, or `recovered stress` |
| `transfer` | no | `variational` | partner-field transfer: `pointwise` or `variational` |
| `transfer quadrature subdivisions` | no | `1` | quadrature refinement for variational transfer (integer ≥ 1) |
| `representable dashpot` | no | `false` | restrict the dashpot to the representable subspace |
| `content aware absorption` | no | `false` | content-aware absorption variant |

Requesting `Schwarz impedance overlap` forces consistent nodal stress recovery
on for the coupled model.

### `Schwarz contact`

Frictionless or tied contact enforced through Schwarz coupling.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner contact surface |
| `friction type` | yes | — | `frictionless` or `tied` |
| `swap BC types` | no | `false` | swap the interface roles between iterations |

## Mesh swapping (`swaps`)

Both single- and multi-domain runs can replace a mesh mid-run when a criterion
fires. `swaps` is a list; each entry names a `replacement` file and a
`criterion` (in a multidomain run it also names the `subsim` to replace). Swaps
cannot be combined with `restart`.

| Key | Required | Meaning |
|---|---|---|
| `subsim` | multidomain | subdomain to replace |
| `replacement` | yes | replacement input file |
| `criterion` | yes | trigger (see below) |

Criterion `type` values: `time` (with `t_swap`); `stress recovery` (with
`tolerance`, default `1.0e-2`, and `direction` `refine`/`coarsen`);
`elastic to plastic transition` (with `tolerance`, default `0.05`);
`overlap l2 relative error` (with `tolerance`, default `1.0e-6`, and `direction`
`refine`/`coarsen`).

## Restart (`restart`)

A run can resume from an Exodus snapshot by adding a `restart` block. Its one
key is `index`, the snapshot time-step index (negative values count back from
the end, so `-1` is the last snapshot). Restart is incompatible with `swaps`,
`initial conditions`, `j2 plasticity`, `Schwarz contact`, and mesh smoothing.

```yaml
restart:
  index: -1
```

## Canonical examples

- Overlapping Schwarz: `examples/overlap/`
- Non-overlapping impedance (same and subcycled steps):
  `examples/nonoverlap/dynamic-same-step/`,
  `examples/nonoverlap/dynamic-different-steps/`
- Contact: `examples/contact/`
- Adaptive mesh swapping: `examples/adaptive-time-stepping/`,
  `examples/ahead/`

See `docs/notes/schwarz-coupling` for the theory and stability analysis of the
impedance coupling.
